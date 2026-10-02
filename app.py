"""
o8G-miRNA Retargeting Explorer  —  Streamlit app
=================================================
Browse human miRNA seeds, enumerate 8-oxoguanine (o8G) seed-oxidation states,
inspect the predicted mRNA target list for each state, and compare target
lists by pathway enrichment (volcano + heatmap).

Run:  streamlit run app.py
Data: o8g_targets.db, o8g_states.parquet, genesets/  (same folder)
Engine modules: o8g_engine.py, o8g_scanner.py, o8g_enrich.py, o8g_plots.py, o8g_db.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Prefer ViennaRNA CLIs from the mirna_viewer conda env when present
_VRNA_BIN = "/opt/homebrew/anaconda3/envs/mirna_viewer/bin"
if os.path.isdir(_VRNA_BIN):
    os.environ["PATH"] = _VRNA_BIN + os.pathsep + os.environ.get("PATH", "")
import numpy as np
import pandas as pd
import streamlit as st

from o8g_engine import extract_seed, enumerate_states, g_positions, SeedState
import o8g_db as _o8g_db
import o8g_precision as _o8g_precision
import importlib as _importlib
_o8g_precision = _importlib.reload(_o8g_precision)  # pick up new PrecisionMode members
_o8g_db = _importlib.reload(_o8g_db)  # pick up ConservationUnavailable across hot-reloads
from o8g_db import TargetDB, ConservationUnavailable
from o8g_enrich import enrich, enrich_within_pool, compare_states, available_libraries
from o8g_genes import ID_TYPES, GeneResolver
from o8g_precision import (
    PrecisionMode,
    PrecisionConfig,
    mode_value,
    is_high_stringency,
)
from o8g_thermo import METRIC_CAPTION, vienna_available
from o8g_pubthermo import annotate_gene_mirna_hits, PROVENANCE_CAPTION as PUBTHERMO_CAPTION
import o8g_refsets as refsets
import o8g_plots as plots
import o8g_sections as sections
plots = _importlib.reload(plots)   # pick up edits to the plotting module on rerun
refsets = _importlib.reload(refsets)
sections = _importlib.reload(sections)
import o8g_oboe as _o8g_oboe
import o8g_oboe_model as _o8g_oboe_model

_o8g_oboe_model = _importlib.reload(_o8g_oboe_model)
_o8g_oboe = _importlib.reload(_o8g_oboe)  # OBOE RNABERT ranking API

st.set_page_config(page_title="o8G-miRNA Retargeting Explorer", layout="wide",
                   initial_sidebar_state="expanded")

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "o8g_targets.db")
RANK_LABEL = {1: "6mer", 2: "7mer-A1", 3: "7mer-m8", 4: "8mer"}


# Bump when TargetDB / GeneResolver constructor API changes so Streamlit
# drops stale @cache_resource instances across hot-reloads.
_CACHE_VER = 26


@st.cache_resource
def get_db(_ver: int = _CACHE_VER):
    return TargetDB(DB_PATH)

@st.cache_resource
def get_scanner(_ver: int = _CACHE_VER):
    """Full in-memory UTR k-mer index (~5 GB). The hosted UI must not call this.

    De novo reads o8g_targets.db. Antagomir collateral streams utr3_human.parquet.
    """
    from o8g_scanner import TargetScanner
    utr = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utr3_human.parquet")
    if not os.path.exists(utr):
        return None
    return TargetScanner.from_parquet(utr)

@st.cache_resource
def get_resolver(_ver: int = _CACHE_VER):
    try:
        return GeneResolver()
    except FileNotFoundError:
        return None

db = get_db()
# Heal stale cached TargetDB instances from before reverse-index / precision APIs.
if (
    not hasattr(db, "targets_filtered")
    or not hasattr(db, "retarget_partition")
    or not hasattr(db, "states_lost_on_oxidation")
    or not hasattr(db, "_targetscan_for")
    or not hasattr(db, "_anchor_symbols_for")
):
    get_db.clear()
    db = get_db()
resolver = get_resolver()


# ---------- helpers ----------
def _ox_positions_from_label(label: str) -> tuple[int, ...]:
    if not label or label == "none":
        return ()
    part = label.replace("o8G@", "")
    return tuple(int(x) for x in part.split(",") if x.strip().isdigit())


def _pathway_gene_peek(enrich_df: pd.DataFrame, *, key: str, title: str = "Pathway genes") -> None:
    """Compact gene-list peek for one selected enrichment term (avoids wide tables)."""
    if enrich_df is None or enrich_df.empty or "term" not in enrich_df.columns:
        return
    with st.expander(title, expanded=False):
        terms = enrich_df["term"].astype(str).tolist()
        pick = st.selectbox("Pathway", terms, key=key)
        row = enrich_df[enrich_df["term"].astype(str) == pick].iloc[0]
        genes = str(row.get("genes", "") or "")
        n = int(row.get("overlap", 0) or 0)
        st.caption(f"{n} overlapping genes · q={row.get('q_value', float('nan')):.2e}")
        if genes:
            # soft wrap as comma-separated text (not a huge dataframe)
            glist = [g for g in genes.split(";") if g]
            st.text(", ".join(glist[:80]) + (" …" if len(glist) > 80 else ""))
            st.download_button(
                "⬇ Genes in this pathway (CSV)",
                pd.DataFrame({"symbol": glist}).to_csv(index=False),
                file_name=f"pathway_genes_{pick[:40].replace(' ', '_')}.csv",
                mime="text/csv",
                key=f"{key}_dl",
            )
        else:
            st.caption("No gene list stored for this term.")


def enrich_binding(
    df: pd.DataFrame,
    *,
    label: str,
    sort: bool = True,
    with_thermo: bool = False,
    scanner=None,
) -> pd.DataFrame:
    """Optional RNAduplex / RNAup / TargetScan context++ columns; sort by site_rank."""
    out = df.copy()
    if with_thermo and scanner is not None and info.get("seq_dna"):
        try:
            from o8g_thermo import score_targets_thermo

            out = score_targets_thermo(
                out,
                scanner=scanner,
                mature_dna=info["seq_dna"],
                oxidized_positions=_ox_positions_from_label(label),
                mirna=mirna,
                is_unmodified=(label == "none"),
            )
        except Exception:
            pass
    if sort and not out.empty:
        cols, asc = [], []
        if "site_rank" in out.columns:
            cols.append("site_rank")
            asc.append(False)
        if "symbol" in out.columns:
            cols.append("symbol")
            asc.append(True)
        if cols:
            out = out.sort_values(cols, ascending=asc).reset_index(drop=True)
    return out


def seed_html(seed: str, oxidized: set) -> str:
    """Render seed with G colored; oxidized Gs boxed red (o8G), normal G blue."""
    out = []
    for i, b in enumerate(seed):
        pos = i + 2
        if b == "G" and pos in oxidized:
            out.append(f"<span style='background:#c44e52;color:white;padding:2px 5px;"
                       f"border-radius:3px;font-weight:700' title='o8G at position {pos}'>{b}<sub>o8</sub></span>")
        elif b == "G":
            out.append(f"<span style='background:#dbe4f0;color:#4c72b0;padding:2px 5px;"
                       f"border-radius:3px;font-weight:700' title='G at position {pos}'>{b}</span>")
        else:
            out.append(f"<span style='padding:2px 4px;color:#333'>{b}</span>")
    return ("<div style='font-family:monospace;font-size:26px;letter-spacing:2px'>"
            "<span style='color:#999;font-size:13px'>5'-</span>" + "".join(out) +
            "<span style='color:#999;font-size:13px'>-3'</span></div>"
            "<div style='font-family:monospace;font-size:10px;color:#999;letter-spacing:2px'>"
            "&nbsp;&nbsp;&nbsp;&nbsp;" + "".join(f"&nbsp;{p}&nbsp;&nbsp;" for p in range(2,9)) + "</div>")


# ---------- sidebar: miRNA selection ----------
st.sidebar.title("🧬 o8G-miRNA Explorer")
st.sidebar.caption("Predicting miRNA target retargeting by 8-oxoguanine seed oxidation")

# Apply a jump requested from the Gene → miRNA/oxomiR tab (previous run).
_jump = st.session_state.pop("_jump_mirna", None)
if _jump:
    st.session_state["mirna_search"] = ""  # clear filter so the target is in the dropdown
    st.session_state["_pending_mirna_select"] = _jump
    _ox = st.session_state.pop("_jump_ox_positions", None)
    if _ox is not None:
        st.session_state["single_state_ox_positions"] = list(_ox)

mir_df = db.mirnas()
query = st.sidebar.text_input(
    "Search miRNA (optional filter)",
    value="",
    placeholder="e.g. miR-1, miR-124, let-7 — leave blank for all",
    help="Type part of a miRNA name to narrow the list. "
         "Leave blank to browse all miRNAs in the dropdown below.",
    key="mirna_search",
)
matches = mir_df[mir_df["mirna"].str.contains(query, case=False, na=False)] if query else mir_df
if len(matches) == 0:
    st.sidebar.warning("No miRNA matches that search.")
    st.stop()
# preselect hsa-miR-1-3p when present so the app opens on a worked example
opts = matches["mirna"].tolist()
default_idx = 0
for _pref in ("hsa-miR-1-3p", "hsa-miR-1-2-3p", "hsa-miR-1"):
    if _pref in opts:
        default_idx = opts.index(_pref)
        break

_pending = st.session_state.pop("_pending_mirna_select", None)
if _pending is not None:
    if _pending in opts:
        st.session_state["mirna_select"] = _pending
        st.sidebar.success(f"Loaded **{_pending}** from Gene → miRNA/oxomiR")
    else:
        st.sidebar.warning(f"`{_pending}` is not in the current miRNA list.")

# Keep selectbox value valid when the search filter shrinks the options
if st.session_state.get("mirna_select") not in opts:
    st.session_state["mirna_select"] = opts[default_idx]

mirna = st.sidebar.selectbox(
    f"Select miRNA ({len(matches)} of {len(mir_df)} shown)",
    opts,
    key="mirna_select",
)
info = db.mirna_info(mirna)
seed = info["seed"]
gpos = g_positions(seed)

# ----- sidebar gene → open Gene / oxomiR reverse lookup -----
@st.cache_data(show_spinner=False)
def _sidebar_gene_symbols(_ver=_CACHE_VER):
    return sorted(str(s) for s in db.symbols if s)


st.sidebar.markdown("---")
st.sidebar.markdown("### Select gene")
_gene_type = st.sidebar.text_input(
    "Gene symbol / Ensembl / Entrez",
    value="",
    placeholder="e.g. HDAC4 or ENSG…",
    key="sidebar_gene_type",
    help="Type a gene, then pick a match (or an exact hit). Opens Gene → miRNA/oxomiR.",
)
_sym_list = _sidebar_gene_symbols()
_gene_choice = None
_q = (_gene_type or "").strip()
if _q:
    _qu = _q.upper()
    _exact = [s for s in _sym_list if s.upper() == _qu]
    if _exact:
        _gene_choice = _exact[0]
        st.sidebar.caption(f"Exact match: `{_gene_choice}`")
    else:
        _pref = [s for s in _sym_list if s.upper().startswith(_qu)][:60]
        _sub = [s for s in _sym_list if _qu in s.upper() and s not in _pref][:40]
        _opts = _pref + _sub
        if _opts:
            _gene_choice = st.sidebar.selectbox(
                f"Matches ({len(_opts)} shown)",
                options=["—"] + _opts,
                index=0,
                key="sidebar_gene_match",
            )
            if _gene_choice == "—":
                _gene_choice = None
        else:
            # allow Ensembl / Entrez / novel symbols through to the Gene tab resolver
            st.sidebar.caption("No UTR symbol match — will resolve on Gene tab.")
            _gene_choice = _q

if _gene_choice:
    _norm = str(_gene_choice).strip()
    _prev = st.session_state.get("_sidebar_gene_loaded")
    if _prev != _norm.upper():
        st.session_state["_sidebar_gene_loaded"] = _norm.upper()
        st.session_state["gene_tab_query"] = _norm
        st.session_state["main_section"] = "Gene → miRNA/oxomiR"
        st.sidebar.success(f"Loading **{_norm}** in Gene → miRNA/oxomiR")
        st.rerun()

st.sidebar.markdown("---")
st.sidebar.markdown(f"**Accession:** `{info['accession']}`")
st.sidebar.markdown(f"**Mature (DNA):** `{info['seq_dna']}`")
st.sidebar.markdown(f"**Seed (pos 2–8):** `{seed}`")
st.sidebar.markdown(f"**Guanines in seed:** {len(gpos)} at position(s) {gpos or '—'}")
st.sidebar.markdown(f"**Oxidation states:** {2**len(gpos)}")

# OBOE oxo-G ranking of which seed Gs to oxidize
with st.sidebar.expander("OBOE oxo-G prior", expanded=False):
    try:
        _ok, _msg = _o8g_oboe.model_status()
        st.caption(
            f"{'Local OBOE RNABERT' if _ok else 'Fallback GC prior'} — {_msg}. "
            "Ranks seed Gs / ox states by P(o8G); not a causal LoF call."
        )
        _g_tbl = _o8g_oboe.seed_g_table(info["seq_dna"])
        _seed_only = _g_tbl[_g_tbl["in_seed"]] if len(_g_tbl) else _g_tbl
        if len(_seed_only):
            _cols = [c for c in ["mature_pos", "window", "oboe_prior", "source"] if c in _seed_only.columns]
            st.dataframe(
                _seed_only[_cols],
                hide_index=True,
                height=min(160, 28 * (len(_seed_only) + 1)),
            )
            _ranked = _o8g_oboe.rank_oxidation_states(info["seq_dna"])
            if len(_ranked):
                st.caption("Top ox states by mean OBOE prior")
                st.dataframe(
                    _ranked.head(5)[["ox_label", "mean_oboe_prior", "n_ox"]],
                    hide_index=True,
                )
        else:
            st.caption("No guanines on this mature.")
    except Exception as _e:
        st.caption(f"Prior unavailable: {_e}")

lib = st.sidebar.selectbox("Pathway library", available_libraries(), index=0)

st.sidebar.markdown("### Prediction mode")
_MODE_OPTIONS = [m.value for m in PrecisionMode]
mode_name = st.sidebar.radio(
    "Target-list filter",
    _MODE_OPTIONS,
    index=0,  # Sequence-based (low stringency) — discovery default
    key=f"prediction_mode_v{_CACHE_VER}",
    help=(
        "Sequence-based (low stringency): 7mer-m8+8mer (explorer DB). "
        "Sequence-based (high stringency): 8mer only. "
        "TargetScan: explorer ∩ catalog WT predictions (WT baseline anchor). "
        "TargetScan de novo: precomputed TargetScanS 8mer/7mer-m8 lists on WT and "
        "oxidized seeds (o8G→T WC encoding; not the web catalog). "
        "Consensus: explorer ∩ TargetScan conserved WT families."
    ),
)
precision_cfg = PrecisionConfig.from_mode(str(mode_name))
effective_mode = str(mode_name)
st.sidebar.caption(
    "Paper / claims: **Sequence-based (high stringency)** or **Consensus**. "
    "**TargetScan** = catalog WT anchor for lost/gained (ox display lists stay rank≥3). "
    "**TargetScan de novo** = precomputed TargetScanS strong sites on oxidized seeds. "
    "Catalog TargetScan and de novo can show the **same oxidized gene count** — "
    "that is expected when both use the same rank≥3 site types on the ox arm. "
    "Gained/lost are always computed after filtering both states."
)


@st.cache_data(show_spinner=False)
def _conservation_status(_seed, _mirna, _ver=_CACHE_VER):
    """Return "" if Consensus is usable, else the reason it is not."""
    try:
        db._conserved_for(_seed, _mirna)
        return ""
    except ConservationUnavailable as e:
        return str(e)


@st.cache_data(show_spinner=False)
def _targetscan_status(_mirna, _ver=_CACHE_VER):
    """Return "" if TargetScan prediction mode is usable."""
    try:
        db._targetscan_for(_mirna)
        return ""
    except ConservationUnavailable as e:
        return str(e)


@st.cache_data(show_spinner=False)
def _utr_universe(_ver=_CACHE_VER):
    """Gene symbols with a 3'UTR in our index — correct ORA population."""
    return set(str(s).upper() for s in db.symbols)


_consensus_problem = ""
if mode_name == "Consensus":
    _consensus_problem = _conservation_status(seed, mirna)
    if _consensus_problem:
        precision_cfg = PrecisionConfig.from_mode(PrecisionMode.SENSITIVE)
        effective_mode = f"{PrecisionMode.SENSITIVE.value} (Consensus unavailable)"
        st.sidebar.error(
            "Consensus unavailable — showing sequence-based (low stringency). "
            + _consensus_problem
        )
elif mode_name == "TargetScan":
    _ts_problem = _targetscan_status(mirna)
    if _ts_problem:
        precision_cfg = PrecisionConfig.from_mode(PrecisionMode.SENSITIVE)
        effective_mode = f"{PrecisionMode.SENSITIVE.value} (TargetScan unavailable)"
        st.sidebar.error(
            "TargetScan mode unavailable — showing sequence-based (low stringency). "
            + _ts_problem
        )

universe = _utr_universe()

if mode_name == "TargetScan de novo":
    st.sidebar.caption(
        "De novo uses the precomputed TargetScanS 8mer/7mer-m8 scan "
        "(o8G encoded as T). Same lists as the explorer database."
    )


def filtered_targets(label, **kw):
    """Single entry point so no call site bypasses the precision ladder."""
    global precision_cfg, effective_mode
    try:
        return db.targets_filtered(seed, label, precision_cfg, mirna=mirna, **kw)
    except ConservationUnavailable as e:
        # Defensive: never crash a view if TS/Consensus/de novo anchor is empty mid-render
        precision_cfg = PrecisionConfig.from_mode(PrecisionMode.SENSITIVE)
        effective_mode = f"{PrecisionMode.SENSITIVE.value} (anchor unavailable)"
        st.sidebar.error(
            "Prediction-mode anchor unavailable — showing sequence-based (low stringency). "
            + str(e)
        )
        return db.targets_filtered(seed, label, precision_cfg, mirna=mirna, **kw)


# Do NOT build TargetScanner (~40s) on every page load — only when thermo / live
# UTR scans are explicitly requested. Streamlit re-runs the whole script on each
# widget change, so eager loading freezes every tab.
scanner = None
# Legacy slider kept as an additional floor within the mode (rarely needed)
min_rank = st.sidebar.select_slider(
    "Extra min site quality floor",
    options=[1, 2, 3, 4],
    value=3,
    format_func=lambda r: {1: "6mer", 2: "7mer-A1", 3: "7mer-m8", 4: "8mer"}[r],
)
st.sidebar.caption("Applied on top of the prediction mode (usually leave at 7mer-m8).")


# ---------- main ----------
st.title(f"{mirna}")
if _consensus_problem:
    st.error(
        "Consensus unavailable — showing sequence-based (low stringency). "
        + _consensus_problem
    )
if len(gpos) == 0:
    st.info("This miRNA's seed contains no guanine, so 8-oxoG oxidation cannot alter its "
            "seed pairing — only one target list exists. Choose a G-containing seed to explore retargeting.")

st.markdown("#### Seed sequence")
st.caption("Guanines are oxidation-eligible. Blue = normal G (pairs C); red = o8G (pairs A).")

states = enumerate_states(seed)
state_labels = [s.label for s in states]

# Conditional sections (not st.tabs): Streamlit executes *all* tab bodies on every
# rerun, so expensive Targets / External-DB work was freezing Gene → miRNA.
# Drop retired views from session (Statistics / Loss of function removed)
if st.session_state.get("main_section") in ("Statistics", "Loss of function"):
    st.session_state["main_section"] = "All states"

_SECTION = st.radio(
    "View",
    [
        "Single state — targets",
        "Compare two states",
        "All states",
        "Overlap (Venn/UpSet)",
        "Transcription factors",
        "Antagomir design",
        "RNA-seq / DEG upload",
        "Gene → miRNA/oxomiR",
        "External DB comparison",
    ],
    horizontal=True,
    key="main_section",
)

_SECTION_DISPATCH = {
    "Overlap (Venn/UpSet)": "render_overlap",
    "Transcription factors": "render_tf",
    "Antagomir design": "render_antagomir",
    "RNA-seq / DEG upload": "render_deg_upload",
}

# ===== TAB 1: single state target list =====
if _SECTION == "Single state — targets":
    c1, c2 = st.columns([3, 2])
    with c1:
        # Gene-tab jump may prefill ox positions; drop any that are not Gs on this seed.
        if "single_state_ox_positions" in st.session_state:
            st.session_state["single_state_ox_multiselect"] = [
                int(p)
                for p in st.session_state.pop("single_state_ox_positions")
                if int(p) in gpos
            ]
        _cur_ox = st.session_state.get("single_state_ox_multiselect", [])
        if any(p not in gpos for p in _cur_ox):
            st.session_state["single_state_ox_multiselect"] = [p for p in _cur_ox if p in gpos]
        ox_choice = st.multiselect(
            "Oxidize which guanine position(s)?",
            gpos,
            key="single_state_ox_multiselect",
            help="Select seed G positions to model as 8-oxoG (o8G). None = fully normal seed.",
        )
        cur = SeedState(seed, tuple(sorted(ox_choice)))
        st.markdown(seed_html(seed, set(ox_choice)), unsafe_allow_html=True)
        st.markdown(f"**State:** `{cur.label}`  •  target motifs (5'→3' on mRNA):")
        st.code(f"8mer   : {cur.motifs['8mer']}\n7mer-m8: {cur.motifs['7mer-m8']}\n"
                f"7mer-A1: {cur.motifs['7mer-A1']}\n6mer   : {cur.motifs['6mer']}")
    with c2:
        meta = db.states_for_seed(seed)
        row = meta[meta["label"] == cur.label]
        if len(row):
            r = row.iloc[0]
            st.metric("Strong-site target genes (7mer-m8 + 8mer)", int(r["n_strong"]))
            st.caption(f"Site inventory: {int(r['n_8mer'])} 8mer · {int(r['n_7mer_m8'])} 7mer-m8 · "
                       f"{int(r['n_7mer_A1'])} 7mer-A1 · {int(r['n_6mer'])} 6mer")

    tdf = filtered_targets(cur.label, scanner=None, mature_dna=info["seq_dna"])
    if "site_rank" in tdf.columns:
        tdf = tdf[tdf["site_rank"] >= min_rank]
    score_targets_thermo = st.checkbox(
        "Add ViennaRNA / TargetScan energetics (slow)",
        value=False,
        help="Joins TargetScan context++ for unmodified miRNAs. "
             "RNAduplex / RNAup need the in-memory UTR index and stay blank here. "
             "Off by default — Streamlit re-runs the script on every click.",
        key="targets_tab_thermo",
    )
    if score_targets_thermo:
        with st.spinner("Joining TargetScan context++…"):
            tdf = enrich_binding(
                tdf, label=cur.label, sort=True, with_thermo=True, scanner=None
            )
        scanner = None
    else:
        tdf = enrich_binding(
            tdf, label=cur.label, sort=True, with_thermo=False, scanner=None
        )
    st.markdown(
        f"**{len(tdf)} predicted target genes** · mode `{effective_mode}` "
        f"(extra floor ≥ {RANK_LABEL[min_rank]})"
    )
    st.caption("Ranked by site type (8mer > 7mer-m8 > 7mer-A1 > 6mer), then symbol.")
    if score_targets_thermo:
        st.caption(METRIC_CAPTION)
        if scanner is None:
            st.caption(
                "RNAduplex / RNAup need per-site coordinates from the in-memory UTR index, "
                "so those columns stay blank. context++ still fills for unmodified miRNAs "
                "from the TargetScan table."
            )
        elif not vienna_available():
            st.caption("ViennaRNA Python package not installed — duplex/RNAup columns unavailable.")
    show_cols = [c for c in [
        "symbol", "gene_id",
        "dG_RNAduplex", "dG_RNAup", "contextpp_TargetScan", "context_analog",
        "site_type", "site_rank",
        "n_8mer", "n_7mer_m8", "score", "context_score", "is_conserved",
    ] if c in tdf.columns]
    st.dataframe(tdf[show_cols] if show_cols else tdf, width='stretch', height=320)
    _mode_fn = effective_mode.replace(" ", "_").replace("(", "").replace(")", "")
    st.download_button("⬇ Download target list (CSV)", tdf.to_csv(index=False),
                       file_name=f"{mirna}_{cur.label}_{_mode_fn}_targets.csv", mime="text/csv")

    if len(tdf) >= 5:
        with st.spinner("Running pathway enrichment…"):
            e = enrich(tdf["symbol"].tolist(), lib, background=universe, top=20)
        st.markdown(f"**Top {lib} pathways** (hypergeometric; background = 3′UTR universe, N={len(universe):,})")
        st.dataframe(
            e[["term", "overlap", "term_size", "odds_ratio", "p_value", "q_value"]],
            width="stretch",
            height=280,
        )
        _pathway_gene_peek(e, key="single_pw_genes", title="Genes in a selected pathway")

# ===== TAB 2: compare two states → volcano =====
elif _SECTION == "Compare two states":
    if len(gpos) == 0:
        st.info("No oxidation states to compare for a G-free seed.")
    else:
        cc1, cc2 = st.columns(2)
        sA = cc1.selectbox("State A", state_labels, index=0)
        sB = cc2.selectbox("State B", state_labels, index=min(3, len(state_labels)-1))
        stA = next(s for s in states if s.label == sA)
        stB = next(s for s in states if s.label == sB)
        cc1.markdown(seed_html(seed, set(stA.oxidized_positions)), unsafe_allow_html=True)
        cc2.markdown(seed_html(seed, set(stB.oxidized_positions)), unsafe_allow_html=True)

        tA = filtered_targets(sA, scanner=None, mature_dna=info["seq_dna"])
        tB = filtered_targets(sB, scanner=None, mature_dna=info["seq_dna"])
        if "site_rank" in tA.columns:
            tA = tA[tA["site_rank"] >= max(min_rank, 3 if not is_high_stringency(precision_cfg) else 4)]
            tB = tB[tB["site_rank"] >= max(min_rank, 3 if not is_high_stringency(precision_cfg) else 4)]
        gA = tA["symbol"].tolist()
        gB = tB["symbol"].tolist()
        sa, sb = set(gA), set(gB)
        m1,m2,m3 = st.columns(3)
        m1.metric(f"Targets in {sA}", len(sa))
        m2.metric(f"Targets in {sB}", len(sb))
        m3.metric("Shared", len(sa & sb), delta=f"{len(sa^sb)} state-specific")

        try:
            fig_v = plots.venn2_plotly(sa, sb, label_a=sA, label_b=sB)
            st.plotly_chart(fig_v, width="stretch")
        except Exception as e:
            st.caption(f"Venn skipped ({type(e).__name__}: {e})")

        # ---- gene-level differential (parallel to the pathway volcano) ----
        st.markdown("#### Differential target genes")
        st.caption(
            f"**Lost** = targeted in `{sA}` but not `{sB}` · "
            f"**Gained** = new targets in `{sB}` · **Shared** = targeted in both. "
            f"Prediction mode `{effective_mode}` applied to each state before the set difference."
        )
        # merge on symbol, carrying each state's best site type
        mA = tA.drop_duplicates("symbol").set_index("symbol")
        mB = tB.drop_duplicates("symbol").set_index("symbol")
        merged = pd.DataFrame(index=sorted(sa | sb))
        merged[f"site_{sA}"] = mA["site_type"].reindex(merged.index)
        merged[f"site_{sB}"] = mB["site_type"].reindex(merged.index)
        merged["gene_id"] = mA["gene_id"].reindex(merged.index).fillna(
            mB["gene_id"].reindex(merged.index))
        merged["category"] = np.where(
            merged.index.isin(sa & sb), "Shared",
            np.where(merged.index.isin(sa - sb), "Lost", "Gained"))
        merged = merged.reset_index().rename(columns={"index": "symbol"})

        lost   = merged[merged["category"] == "Lost"][["symbol", "gene_id", f"site_{sA}"]]
        gained = merged[merged["category"] == "Gained"][["symbol", "gene_id", f"site_{sB}"]]
        shared = merged[merged["category"] == "Shared"][["symbol", "gene_id", f"site_{sA}", f"site_{sB}"]]
        # optional thermo columns when annotate path is used elsewhere; compare stays lean
        tA_ann = enrich_binding(tA, label=sA, sort=False).drop_duplicates("symbol")
        tB_ann = enrich_binding(tB, label=sB, sort=False).drop_duplicates("symbol")
        for col in ("dG_RNAduplex", "dG_RNAup", "contextpp_TargetScan"):
            if col in tA_ann.columns:
                lost[col] = lost["symbol"].map(tA_ann.set_index("symbol")[col])
            if col in tB_ann.columns:
                gained[col] = gained["symbol"].map(tB_ann.set_index("symbol")[col])
                shared[col] = shared["symbol"].map(tB_ann.set_index("symbol")[col]).fillna(
                    shared["symbol"].map(tA_ann.set_index("symbol")[col])
                    if col in tA_ann.columns
                    else np.nan
                )
        # rank by site type string order via optional site_rank maps
        def _sort_diff(frame: pd.DataFrame, site_col: str) -> pd.DataFrame:
            if site_col in frame.columns:
                order = {"8mer": 4, "7mer-m8": 3, "7mer-A1": 2, "6mer": 1}
                return frame.assign(
                    _rk=frame[site_col].map(order)
                ).sort_values(["_rk", "symbol"], ascending=[False, True], na_position="last").drop(
                    columns="_rk"
                )
            return frame.sort_values("symbol")

        lost = _sort_diff(lost, f"site_{sA}")
        gained = _sort_diff(gained, f"site_{sB}")
        shared = _sort_diff(shared, f"site_{sB}")
        st.caption(METRIC_CAPTION)

        gc1, gc2, gc3 = st.columns(3)
        with gc1:
            st.markdown(f"**🔴 Lost on {sB}** ({len(lost)})")
            st.dataframe(lost, width='stretch', height=300, hide_index=True)
        with gc2:
            st.markdown(f"**🔵 Gained on {sB}** ({len(gained)})")
            st.dataframe(gained, width='stretch', height=300, hide_index=True)
        with gc3:
            st.markdown(f"**⚪ Shared** ({len(shared)})")
            st.dataframe(shared, width='stretch', height=300, hide_index=True)

        st.download_button(
            "⬇ Download differential gene table (CSV)",
            merged[["symbol", "gene_id", f"site_{sA}", f"site_{sB}", "category"]].to_csv(index=False),
            file_name=f"{mirna}_{sA}_vs_{sB}_differential_genes.csv", mime="text/csv")

        st.markdown("#### Differential pathway enrichment")
        st.info(
            "Genome-wide (UTR-universe) enrichment answers *what this miRNA targets*. "
            "The **within-pool** control uses background = union(state A, state B) and "
            "answers *what oxidation changed*. Never read the genome panel alone."
        )
        chart = st.radio("Display", ["Diverging bar", "Volcano"], horizontal=True,
            help="Diverging bar names each pathway and shows direction directly — clearest "
                 "when few terms are significant. Volcano is better with many terms.")
        if len(sa) >= 5 and len(sb) >= 5:
            pool = sa | sb
            with st.spinner("Computing differential pathway enrichment…"):
                cmp = compare_states(
                    gA, gB, library=lib, label_A=sA, label_B=sB,
                    background=universe, min_overlap=3,
                )
                cmp["log2_or"] = (-cmp["log2_or_ratio"]).clip(-8, 8)   # x>0 → toward A
                cmp["q_value"] = cmp[[f"q_{sA}", f"q_{sB}"]].min(axis=1)
                cmp["neglog10_q"] = cmp["neglog10_q_best"]
                lost_genes = list(sa - sb)
                gained_genes = list(sb - sa)
                within_lost = enrich_within_pool(lost_genes, pool, library=lib, min_overlap=2, top=30)
                within_gained = enrich_within_pool(gained_genes, pool, library=lib, min_overlap=2, top=30)
            n_sig_genome = int((cmp["q_value"] < 0.05).sum()) if len(cmp) else 0
            n_sig_pool = int(
                ((within_lost["q_value"] < 0.05).sum() if len(within_lost) else 0)
                + ((within_gained["q_value"] < 0.05).sum() if len(within_gained) else 0)
            )
            c_g, c_p = st.columns(2)
            c_g.metric(f"Significant terms vs UTR universe (q<0.05)", n_sig_genome)
            c_p.metric(f"Significant terms within target pool (q<0.05)", n_sig_pool)
            if chart == "Diverging bar":
                fig = plots.diverging_bar_plotly(cmp, label_A=sA, label_B=sB, x_col="log2_or")
                st.plotly_chart(fig, width='stretch')
                st.caption(f"Each bar is a pathway; length = −log10 q. Right (red) = enriched in "
                           f"**{sA}**; left (blue) = enriched in **{sB}**. Dashed lines = q<0.05. "
                           f"Background = 3′UTR universe (N={len(universe):,}).")
            else:
                fig = plots.volcano_plotly(cmp, label_A=sA, label_B=sB, x_col="log2_or")
                st.plotly_chart(fig, width='stretch')
                st.caption(f"Each point is a pathway. Right (red) = enriched in **{sA}**; "
                           f"left (blue) = enriched in **{sB}**. Dashed line = q<0.05. "
                           f"Background = 3′UTR universe (N={len(universe):,}).")
            st.dataframe(cmp.sort_values("neglog10_q_best", ascending=False)
                         [["term","library",f"overlap_{sA}",f"overlap_{sB}",
                           f"odds_{sA}",f"odds_{sB}","q_value"]].head(40),
                         width='stretch', height=300)
            st.markdown("##### Within-pool control (lost / gained vs union of both states)")
            wp1, wp2 = st.columns(2)
            with wp1:
                st.caption(f"Lost genes vs pool (n_pool={len(pool):,})")
                st.dataframe(
                    within_lost[["term", "overlap", "odds_ratio", "q_value"]] if len(within_lost) else within_lost,
                    width="stretch", height=220,
                )
                _pathway_gene_peek(
                    within_lost, key="cmp_lost_pw_genes", title="Genes in a lost pathway"
                )
            with wp2:
                st.caption(f"Gained genes vs pool (n_pool={len(pool):,})")
                st.dataframe(
                    within_gained[["term", "overlap", "odds_ratio", "q_value"]] if len(within_gained) else within_gained,
                    width="stretch", height=220,
                )
                _pathway_gene_peek(
                    within_gained, key="cmp_gain_pw_genes", title="Genes in a gained pathway"
                )
        else:
            st.warning("Need ≥5 strong-site targets in each state for enrichment.")

# ===== TAB 3: all states heatmap =====
elif _SECTION == "All states":
    if len(gpos) == 0:
        st.info("Only one state exists for a G-free seed.")
    else:
        st.markdown(
            f"#### All oxidation states · `{mirna}` · mode `{effective_mode}`"
        )
        # LoF / gained-lost summary (moved from former Loss of function tab)
        def _strong_set_all(label: str) -> set[str]:
            df = filtered_targets(label, scanner=None, mature_dna=info["seq_dna"])
            if "site_rank" in df.columns:
                floor = 4 if is_high_stringency(precision_cfg) else max(min_rank, 3)
                df = df[df["site_rank"] >= floor]
            return set(df["symbol"].astype(str))

        _lof_ctx = sections.SectionContext(
            db=db,
            mirna=mirna,
            info=info,
            seed=seed,
            gpos=gpos,
            state_labels=state_labels,
            library=lib,
            precision_mode=effective_mode,
            strong_set=_strong_set_all,
            universe=universe,
            matched_background=universe,
            external_refs=refsets,
            precision_cfg=precision_cfg,
            scanner=None,
        )
        from o8g_lof import render_state_summary as _render_lof_summary

        _render_lof_summary(_lof_ctx)

        st.markdown("---")
        st.markdown(
            f"Pathway enrichment (−log10 q) across all **{len(states)}** seed-oxidation "
            f"states, library **{lib}** · background = 3′UTR universe (N={len(universe):,})."
        )
        try:
            with st.expander("OBOE ranking of oxidation states", expanded=False):
                _ok, _msg = _o8g_oboe.model_status()
                st.caption(
                    f"Prioritize G→o8G combinations by local OBOE site probability "
                    f"({'RNABERT' if _ok else 'GC fallback'}: {_msg}). "
                    "Consequence tables below remain the LoF layer."
                )
                _rt = _o8g_oboe.rank_oxidation_states(info["seq_dna"])
                st.dataframe(_rt, hide_index=True, height=260, use_container_width=True)
                if st.checkbox("Probe live OBOE server (best-effort)", value=False, key="oboe_remote"):
                    with st.spinner("Contacting rnamd.org OBOE…"):
                        rem = _o8g_oboe.try_remote_oboe(info["seq_dna"])
                    if rem is None or rem.get("status") == "error":
                        st.warning(
                            "Live OBOE unavailable (server-side model down). "
                            "Using the local OBOE RNABERT / prior table above."
                        )
                    else:
                        st.json(rem)
        except Exception as e:
            st.caption(f"OBOE prior skipped: {e}")
        view = st.radio("Display", ["Heatmap", "Dot plot"], horizontal=True,
            help="Heatmap = enrichment intensity per state. Dot plot adds gene-overlap as dot "
                 "size, so you also see how many genes support each term in each state.")
        top_terms = st.slider("Pathways to show", 8, 40, 22)
        def _enrich_utr(genes, library=lib, **kw):
            return enrich(genes, library=library, background=universe, **kw)
        sg = {}
        for s in states:
            df = filtered_targets(s.label, scanner=None, mature_dna=info["seq_dna"])
            if "site_rank" in df.columns:
                df = df[df["site_rank"] >= max(min_rank, 3)]
            sg[s.label] = df["symbol"].tolist()
        states_order = [s.label for s in states]
        if view == "Heatmap":
            with st.spinner("Enriching all states…"):
                mat = plots.enrichment_matrix(sg, _enrich_utr, library=lib, top_terms=top_terms)
            if mat.shape[0] >= 2:
                fig = plots.heatmap_plotly(mat)
                st.plotly_chart(fig, width='stretch')
                st.download_button("⬇ Download enrichment matrix (CSV)", mat.to_csv(),
                                   file_name=f"{mirna}_state_enrichment_{lib}.csv", mime="text/csv")
            else:
                st.warning("Not enough enriched pathways to build a heatmap at this setting.")
        else:
            with st.spinner("Enriching all states…"):
                long = plots.state_dotplot_data(sg, _enrich_utr, library=lib, top_terms=top_terms)
            if long["term"].nunique() >= 2:
                fig = plots.dotplot_plotly(long, states_order=states_order)
                st.plotly_chart(fig, width='stretch')
                st.caption("Dot size = number of overlapping genes; color = −log10 q. "
                           "Compare a row across columns to see how oxidation shifts each pathway.")
                st.download_button("⬇ Download enrichment table (CSV)", long.to_csv(index=False),
                                   file_name=f"{mirna}_state_dotplot_{lib}.csv", mime="text/csv")
            else:
                st.warning("Not enough enriched pathways to build a dot plot at this setting.")

elif _SECTION in _SECTION_DISPATCH:
    def _strong_set(label: str) -> set[str]:
        df = filtered_targets(label, scanner=None, mature_dna=info["seq_dna"])
        if "site_rank" in df.columns:
            floor = 4 if is_high_stringency(precision_cfg) else max(min_rank, 3)
            df = df[df["site_rank"] >= floor]
        return set(df["symbol"].astype(str))

    ctx = sections.SectionContext(
        db=db,
        mirna=mirna,
        info=info,
        seed=seed,
        gpos=gpos,
        state_labels=state_labels,
        library=lib,
        precision_mode=effective_mode,
        strong_set=_strong_set,
        universe=universe,
        matched_background=universe,
        external_refs=refsets,
        precision_cfg=precision_cfg,
        scanner=None,
    )
    getattr(sections, _SECTION_DISPATCH[_SECTION])(ctx)

# ===== TAB 4: reverse gene → miRNA/oxomiR =====
elif _SECTION == "Gene → miRNA/oxomiR":
    st.markdown("#### Find miRNA / oxomiR states that target a gene")
    st.caption(
        "Offline reverse lookup against strong sites (7mer-m8 + 8mer) in the precomputed "
        "database. Gene IDs resolve via a local NCBI + HGNC map (no live API)."
    )
    if resolver is None:
        st.error(
            "Gene alias map missing. Build it with:\n\n"
            "`python scripts/fetch_gene_aliases.py`\n\n"
            "Then rebuild the reverse index if needed:\n\n"
            "`python scripts/build_reverse_index.py`"
        )
    elif getattr(db, "reverse_path", None) is None:
        st.error(
            "Reverse index missing. Build it with:\n\n"
            "`python scripts/build_reverse_index.py`"
        )
    else:
        gc1, gc2 = st.columns([3, 1])
        if "gene_tab_query" not in st.session_state:
            st.session_state["gene_tab_query"] = "HDAC4"
        gene_q = gc1.text_input(
            "Gene",
            placeholder="HDAC4 · ENSG00000068024 · 9759",
            help="Symbol, Ensembl gene ID, or Entrez Gene ID. Prefill from the sidebar Select gene.",
            key="gene_tab_query",
        )
        id_type = gc2.selectbox("ID type", ID_TYPES, index=0)

        if gene_q.strip():
            res = resolver.resolve(gene_q, id_type)
            if res.error:
                st.warning(res.error)
            else:
                if res.ambiguous:
                    st.info(
                        f"Ambiguous {res.detected_type} match ({len(res.hits)} genes). "
                        "Showing the first; refine the query or pick Ensembl/Entrez."
                    )
                    pick_labels = [
                        f"{h.symbol} · {h.ensembl}"
                        + (f" · Entrez {h.entrez}" if h.entrez else "")
                        for h in res.hits
                    ]
                    pick = st.selectbox("Choose gene", pick_labels, index=0)
                    hit = res.hits[pick_labels.index(pick)]
                else:
                    hit = res.hits[0]

                meta_bits = [
                    f"**{hit.symbol}**",
                    f"Ensembl `{hit.ensembl}`",
                ]
                if hit.entrez is not None:
                    meta_bits.append(f"Entrez `{hit.entrez}`")
                meta_bits.append(f"matched as {res.detected_type if id_type == 'Auto' else id_type}")
                st.markdown(" · ".join(meta_bits))
                built = resolver.meta.get("built_at", "?")
                st.caption(
                    f"ID map built {built} · sources: {resolver.meta.get('sources', 'local')}"
                )

                with st.spinner("Querying reverse index…"):
                    hits = db.states_targeting_gene(hit.gene_idx)

                if hits.empty:
                    st.info("No strong-site miRNA/oxomiR states target this gene in the database.")
                else:
                    filt = st.radio(
                        "Show states",
                        [
                            "All",
                            "Unmodified only",
                            "Oxidized only",
                            "Gained on oxidation only",
                            "Lost on oxidation only",
                            "Also in unmodified",
                        ],
                        horizontal=True,
                        help=(
                            "Lost on oxidation = gene is a strong-site target of the unmodified "
                            "seed but not of that oxidized state (inferred; not stored as a hit)."
                        ),
                    )
                    if filt == "Lost on oxidation only":
                        with st.spinner("Inferring losses vs oxidized states…"):
                            view = db.states_lost_on_oxidation(hit.gene_idx)
                        if view.empty:
                            st.info(
                                "No losses inferred: every oxidized state of seeds that target "
                                "this gene when unmodified still targets it (or the gene has no "
                                "unmodified strong-site hits)."
                            )
                    else:
                        view = hits
                        if filt == "Unmodified only":
                            view = hits[hits["state_label"] == "none"]
                        elif filt == "Oxidized only":
                            view = hits[hits["state_label"] != "none"]
                        elif filt == "Gained on oxidation only":
                            view = hits[hits["vs_unmodified"] == "gained on oxidation"]
                        elif filt == "Also in unmodified":
                            view = hits[hits["vs_unmodified"] == "also in unmodified"]

                    if not view.empty:
                        n_seeds = view["seed"].nunique()
                        n_ox = int((view["state_label"] != "none").sum())
                        m1, m2, m3 = st.columns(3)
                        m1.metric("Matching states", len(view))
                        m2.metric("Unique seeds", int(n_seeds))
                        m3.metric(
                            "Lost ox-state rows" if filt.startswith("Lost") else "Oxidized-state rows",
                            n_ox,
                        )

                        score_thermo = st.checkbox(
                            "Add published energetics (RNAduplex / RNAup / TargetScan context++)",
                            value=False,
                            help="Same methods as the interaction viewer: Bartel site type is already "
                                 "in site_type; ΔG from ViennaRNA; context++ from TargetScan 8 tables. "
                                 "Off by default — scoring runs RNAduplex/RNAup on the first N rows. "
                                 "Not available for Lost rows (gene is not a target of that ox state).",
                            key="gene_tab_pubthermo_v2",
                            disabled=filt.startswith("Lost"),
                        )
                        max_thermo = st.slider(
                            "Max rows to score with RNAduplex/RNAup",
                            min_value=20,
                            max_value=min(500, max(20, len(view))),
                            value=min(50, max(20, len(view))),
                            help="Full gene hits can be thousands of rows; energetics are computed "
                                 "for the first N rows of the filtered table (results are cached).",
                            key="gene_tab_thermo_n_v2",
                            disabled=not score_thermo or filt.startswith("Lost"),
                        )

                        show = view[
                            [
                                "mirna",
                                "seed",
                                "state_label",
                                "oxidized_positions",
                                "site_type",
                                "site_rank",
                                "motif_7mer_m8",
                                "motif_8mer",
                                "vs_unmodified",
                                "all_mirnas",
                            ]
                        ].rename(
                            columns={
                                "state_label": "oxidation_state",
                                "vs_unmodified": "vs_unmodified_seed",
                            }
                        )

                        if score_thermo and not filt.startswith("Lost"):
                            with st.spinner(
                                f"Scoring up to {max_thermo} rows (RNAduplex/RNAup) + "
                                "one-pass TargetScan context++…"
                            ):
                                scored = annotate_gene_mirna_hits(
                                    view,
                                    gene_symbol=hit.symbol,
                                    gene_idx=hit.gene_idx,
                                    scanner=None,
                                    db=db,
                                    max_thermo_rows=max_thermo,
                                )
                            for col in (
                                "dG_hybrid",
                                "dG_open",
                                "ddG",
                                "contextpp_TargetScan",
                            ):
                                if col in scored.columns:
                                    show[col] = scored[col].values
                            # Prefer a readable column order
                            lead = [
                                "mirna",
                                "oxidation_state",
                                "site_type",
                                "dG_hybrid",
                                "dG_open",
                                "ddG",
                                "contextpp_TargetScan",
                                "seed",
                                "oxidized_positions",
                                "site_rank",
                                "motif_7mer_m8",
                                "motif_8mer",
                                "vs_unmodified_seed",
                                "all_mirnas",
                            ]
                            show = show[[c for c in lead if c in show.columns]]
                            st.caption(PUBTHERMO_CAPTION)
                            _utr_pq = os.path.join(
                                os.path.dirname(os.path.abspath(__file__)), "utr3_human.parquet"
                            )
                            if not os.path.exists(_utr_pq):
                                st.caption(
                                    "UTR parquet unavailable — motif windows missing; "
                                    "ΔG columns may be empty."
                                )

                        st.dataframe(show, width="stretch", height=420, hide_index=True)
                        st.caption(
                            "**vs_unmodified_seed:** *unmodified* = targeted by the normal seed; "
                            "*also in unmodified* = still targeted when oxidized; "
                            "*gained on oxidation* = only appears under an o8G state for that seed; "
                            "*lost on oxidation* = targeted by unmodified but not by that o8G state "
                            "(site_type/rank are from the unmodified hit that was lost)."
                        )

                        st.markdown("##### Open a hit in Explorer")
                        st.caption(
                            "Loads the miRNA into the sidebar (and optionally its oxidation state "
                            "into **Single state — targets**). Then switch tabs to browse / compare / "
                            "enrich as usual."
                        )
                        jump_opts = [
                            f"{r.mirna}  ·  {r.oxidation_state}"
                            for r in show.itertuples(index=False)
                        ]
                        # Prefer unique labels while preserving order
                        seen = set()
                        jump_opts_u = []
                        for lab in jump_opts:
                            if lab not in seen:
                                seen.add(lab)
                                jump_opts_u.append(lab)
                        if jump_opts_u:
                            jump_lab = st.selectbox(
                                "miRNA · oxidation state",
                                jump_opts_u,
                                key="gene_tab_jump_pick",
                            )
                            j_mirna, j_state = [x.strip() for x in jump_lab.split("·", 1)]
                            # Parse oxidized positions from the matching row
                            j_row = show[
                                (show["mirna"] == j_mirna)
                                & (show["oxidation_state"] == j_state)
                            ].iloc[0]
                            raw_pos = str(j_row.get("oxidized_positions", "") or "")
                            digits = [
                                int(x) for x in raw_pos.split(",") if x.strip().isdigit()
                            ]
                            if j_state == "none":
                                digits = []
                            cja, cjb = st.columns([1, 2])
                            with cja:
                                if st.button(
                                    "Load into Explorer",
                                    type="primary",
                                    key="gene_tab_jump_btn",
                                    help="Set sidebar miRNA (+ Single-state oxidation positions).",
                                ):
                                    st.session_state["_jump_mirna"] = j_mirna
                                    st.session_state["_jump_ox_positions"] = digits
                                    st.rerun()
                            with cjb:
                                st.caption(
                                    f"Will set sidebar → **{j_mirna}**"
                                    + (
                                        f", Single state ox → positions {digits}"
                                        if digits
                                        else ", Single state ox → none"
                                    )
                                )

                        st.download_button(
                            "⬇ Download miRNA/oxomiR hits (CSV)",
                            show.to_csv(index=False),
                            file_name=f"{hit.symbol}_targeting_oxomirs.csv",
                            mime="text/csv",
                        )

# ===== TAB 5: external DB comparison + UpSet =====
elif _SECTION == "External DB comparison":
    st.markdown("#### Multi-database target comparison")
    st.caption(
        "Compare **any Explorer oxidation state** (unmodified or o8G) to open external "
        "resources. External DBs only catalog **unmodified** miRNAs — useful for asking "
        "which oxomiR targets are novel vs already known wild-type targets. "
        "**Explorer** always follows the sidebar prediction mode "
        f"(currently `{effective_mode}`). "
        "Optional **TargetScan de novo** is the precomputed TargetScanS strong-site "
        "list for that same oxidation state (independent of the sidebar mode). "
        "Predicted catalogs: TargetScan 8 / miRDB≥80 / DIANA-microT≥0.7 / miRmap. "
        "Experimental: ENCORI CLIP; miRTarBase if a local file is present. "
        "Note: **catalog TargetScan** only lists unmodified miRNAs; on an oxidized "
        "Explorer state, catalog ∩ explorer can look like **TargetScan de novo** "
        "because both ox arms use the same rank≥3 site types — contrast catalog WT "
        "TargetScan vs de novo for a true algorithm difference."
    )
    db_state = st.selectbox(
        "Explorer oxidation state",
        state_labels,
        index=0,
        key="ext_db_explorer_state",
        help="External catalogs stay on the unmodified mature; Explorer and optional "
             "TargetScan de novo follow this oxidation state.",
    )
    st.markdown(
        seed_html(
            seed,
            set(next(s for s in states if s.label == db_state).oxidized_positions),
        ),
        unsafe_allow_html=True,
    )
    include_unmod_explorer = st.checkbox(
        "Also include Explorer unmodified (`none`) as a separate set",
        value=(db_state != "none"),
        help="When comparing an oxidized state, keep the wild-type Explorer list in the UpSet "
             "so gained/lost vs external DBs is visible.",
    )
    include_ts_denovo = st.checkbox(
        "Also include TargetScan de novo for this oxidation state",
        value=False,
        key="ext_db_ts_denovo",
        help="Precomputed TargetScanS 8mer/7mer-m8 list for the selected state "
             "(o8G→T WC encoding). Contrasts catalog TargetScan (WT only) with the "
             "ox-aware strong-site list, even when the sidebar mode is sequence-based.",
    )
    avail = refsets.available_tools()
    default_tools = [t for t, ok in avail.items() if ok and t != "miRTarBase"]
    if avail.get("miRTarBase"):
        default_tools.append("miRTarBase")
    tool_opts = ["Explorer"] + list(refsets.LOADERS.keys())
    chosen = st.multiselect(
        "Databases to compare",
        tool_opts,
        default=["Explorer"] + default_tools[:4],
        help=(
            f"Explorer uses sidebar prediction mode `{effective_mode}` on the oxidation "
            "state selected above. Catalog TargetScan is always unmodified WT."
        ),
    )
    n_est = len(chosen) + (
        1 if include_unmod_explorer and db_state != "none" and "Explorer" in chosen else 0
    ) + (1 if include_ts_denovo else 0)
    # Streamlit raises when min_value == max_value; two selected sets → fixed min_dbs=2.
    if n_est > 2:
        min_dbs = st.slider(
            "Master list: gene must appear in at least this many selected DBs",
            min_value=2,
            max_value=n_est,
            value=2,
            help="Raise to the number of sets for a strict intersection (present in every selected DB).",
        )
    else:
        min_dbs = 2
        st.caption("Master list requires agreement of both selected databases.")

    if len(chosen) < 2 and not (
        (include_unmod_explorer and "Explorer" in chosen) or include_ts_denovo
    ):
        st.info(
            "Select at least two databases (or Explorer + unmodified / TargetScan de novo)."
        )
    else:
        with st.spinner(f"Loading reference sets for {mirna} [{db_state}]…"):
            sets: dict[str, set] = {}
            ours = pd.DataFrame()
            ours_unmod = pd.DataFrame()

            def _load_explorer(label: str) -> pd.DataFrame:
                df = filtered_targets(label, scanner=None, mature_dna=info["seq_dna"])
                if "site_rank" in df.columns:
                    df = df[df["site_rank"] >= min_rank]
                return df

            if "Explorer" in chosen:
                ours = _load_explorer(db_state)
                label_key = f"Explorer ({db_state})"
                sets[label_key] = set(ours["symbol"])
            if include_unmod_explorer and db_state != "none":
                ours_unmod = _load_explorer("none")
                sets["Explorer (none)"] = set(ours_unmod["symbol"])
            if include_ts_denovo:
                try:
                    _dn = db.targets(seed, db_state)
                    if "site_rank" in _dn.columns and len(_dn):
                        _dn = _dn[_dn["site_rank"] >= max(min_rank, 3)]
                    sets[f"TargetScan de novo ({db_state})"] = (
                        set(_dn["symbol"].astype(str)) if len(_dn) else set()
                    )
                except Exception as _e:
                    st.warning(f"TargetScan de novo set unavailable: {_e}")
            ext = refsets.load_selected(mirna, [t for t in chosen if t != "Explorer"])
            sets.update(ext)

        n_sets = len(sets)
        if n_sets < 2:
            st.info("Need at least two non-empty sets to compare.")
        else:
            # Set sizes — one row per set. Do not pack metrics into 4 columns with
            # i%4 (Streamlit overwrites cards and makes fixed external DBs look like
            # they changed when only Explorer's label/order shifts).
            size_rows = []
            for name, s in sets.items():
                base = name.split(" (")[0]
                meta = refsets.TOOL_META.get(base, {})
                if name.startswith("Explorer (") and name != "Explorer (none)":
                    ox_note = f"yes — Explorer under precision `{effective_mode}`"
                elif name == "Explorer (none)":
                    ox_note = "fixed — Explorer unmodified (`none`)"
                elif name.startswith("TargetScan de novo"):
                    ox_note = "yes — live TargetScanS on this ox state (not the catalog)"
                else:
                    ox_note = "no — unmodified miRNA catalog (same for every o8G state)"
                size_rows.append(
                    {
                        "set": name,
                        "n_genes": len(s),
                        "changes_with_oxidation": ox_note,
                        "threshold": meta.get("threshold", ""),
                    }
                )
            st.markdown("#### Set sizes")
            st.caption(
                "Only **Explorer (o8G…)** and **TargetScan de novo** change with oxidation. "
                "Catalog TargetScan / miRDB / DIANA / miRmap / ENCORI stay fixed for this miRNA. "
                "UpSet bars and the master list change because intersections with Explorer change."
            )
            st.dataframe(pd.DataFrame(size_rows), hide_index=True, width="stretch")

            empty_local = [
                t for t in chosen
                if t not in ("Explorer", "ENCORI") and not avail.get(t, False)
            ]
            if empty_local:
                st.warning(
                    "Local files missing for: **" + ", ".join(empty_local) + "**. "
                    "Place downloads under `paper/data/` (see `o8g_refsets.py` docstring)."
                )

            mat = refsets.membership_matrix(sets)
            if mat.empty:
                st.warning("No genes loaded from the selected databases.")
            else:
                # clamp min_dbs to available sets
                min_use = min(min_dbs, n_sets)
                core = refsets.consensus_intersection(sets, min_dbs=min_use)
                st.markdown(
                    f"#### Master comparison list "
                    f"({len(core)} genes in ≥{min_use} of {n_sets} sets) · Explorer state `{db_state}`"
                )
                master = mat[mat["n_dbs"] >= min_use].copy()
                # optional thermo / site columns from the selected Explorer state
                be_src = ours if not ours.empty else ours_unmod
                if not be_src.empty:
                    be_lab = db_state if not ours.empty else "none"
                    be_ann = enrich_binding(be_src, label=be_lab, sort=False).drop_duplicates(
                        "symbol"
                    )
                    if "site_type" in be_ann.columns:
                        master["site_type"] = master["symbol"].map(
                            be_ann.set_index("symbol")["site_type"]
                        )
                    for col in ("dG_RNAduplex", "dG_RNAup", "contextpp_TargetScan"):
                        if col in be_ann.columns:
                            master[col] = master["symbol"].map(be_ann.set_index("symbol")[col])
                    sort_cols = ["n_dbs"]
                    asc = [False]
                    if "site_type" in master.columns:
                        order = {"8mer": 4, "7mer-m8": 3, "7mer-A1": 2, "6mer": 1}
                        master["_rk"] = master["site_type"].map(order)
                        sort_cols.append("_rk")
                        asc.append(False)
                    sort_cols.append("symbol")
                    asc.append(True)
                    master = master.sort_values(
                        sort_cols, ascending=asc, na_position="last"
                    )
                    if "_rk" in master.columns:
                        master = master.drop(columns="_rk")
                else:
                    master = master.sort_values(["n_dbs", "symbol"], ascending=[False, True])

                # flag genes gained on oxidation vs unmodified Explorer (when both present)
                if "Explorer (none)" in sets and any(k.startswith("Explorer (") and k != "Explorer (none)" for k in sets):
                    ox_key = next(k for k in sets if k.startswith("Explorer (") and k != "Explorer (none)")
                    master["vs_unmodified_explorer"] = np.where(
                        master["symbol"].isin(sets[ox_key] - sets["Explorer (none)"]),
                        "gained on oxidation",
                        np.where(
                            master["symbol"].isin(sets["Explorer (none)"] - sets[ox_key]),
                            "lost on oxidation",
                            np.where(
                                master["symbol"].isin(sets[ox_key] & sets["Explorer (none)"]),
                                "shared",
                                "external only",
                            ),
                        ),
                    )

                st.dataframe(master, width="stretch", height=360, hide_index=True)
                st.download_button(
                    "⬇ Download master comparison list (CSV)",
                    master.to_csv(index=False),
                    file_name=f"{mirna}_{db_state}_multidb_master_min{min_use}.csv",
                    mime="text/csv",
                )

                st.markdown("#### UpSet plot (exclusive intersections)")
                fig_u = plots.upset_plotly(sets)
                st.plotly_chart(fig_u, width="stretch")

                with st.expander("Database versions & thresholds"):
                    meta_rows = []
                    for name in sets:
                        base = name.split(" (")[0]
                        m = refsets.TOOL_META.get(base, {})
                        meta_rows.append(
                            {
                                "set": name,
                                "kind": m.get("kind", "predicted"),
                                "version": m.get("version", ""),
                                "threshold": m.get("threshold", ""),
                                "citation": m.get("citation", ""),
                                "n_genes": len(sets[name]),
                            }
                        )
                    st.dataframe(pd.DataFrame(meta_rows), hide_index=True, width="stretch")

                # ---- Additive: prioritize Explorer-lost genes with external WT support ----
                ox_keys = [k for k in sets if k.startswith("Explorer (") and k != "Explorer (none)"]
                if db_state != "none" and "Explorer (none)" in sets and ox_keys:
                    ox_key = ox_keys[0]
                    lost_syms = sets["Explorer (none)"] - sets[ox_key]
                    ext_only = {k: v for k, v in sets.items() if not k.startswith("Explorer")}
                    if lost_syms and ext_only:
                        st.markdown("#### Predicted lost targets (external corroboration)")
                        st.caption(
                            "Same comparison as above, filtered to genes **lost on oxidation** "
                            "(in Explorer unmodified, absent after the selected o8G state). "
                            "`n_external` counts how many selected WT catalogs list each gene. "
                            "Explorer is strongest at loss-of-binding; use this list for follow-up."
                        )
                        min_ext = st.slider(
                            "Lost list: require support from at least this many external DBs",
                            min_value=0,
                            max_value=max(1, len(ext_only)),
                            value=1,
                            key="ext_lost_min_support",
                        )
                        unmod_src = ours_unmod if not ours_unmod.empty else ours
                        unmod_ann = (
                            enrich_binding(unmod_src, label="none", sort=False)
                            .drop_duplicates("symbol")
                            .set_index("symbol")
                        )
                        support = refsets.external_support_for_genes(lost_syms, ext_only)
                        lost_tbl = support[support["n_external"] >= min_ext].copy()
                        if "gene_id" in unmod_ann.columns:
                            lost_tbl["gene_id"] = lost_tbl["symbol"].map(unmod_ann["gene_id"])
                        if "site_type" in unmod_ann.columns:
                            lost_tbl["site_unmodified"] = lost_tbl["symbol"].map(
                                unmod_ann["site_type"]
                            )
                        for col in ("dG_RNAduplex", "dG_RNAup", "contextpp_TargetScan"):
                            if col in unmod_ann.columns:
                                lost_tbl[col] = lost_tbl["symbol"].map(unmod_ann[col])
                        lost_tbl = lost_tbl.sort_values(
                            ["n_external", "symbol"],
                            ascending=[False, True],
                            na_position="last",
                        )
                        lead = [
                            "symbol", "gene_id", "site_unmodified",
                            "n_external",
                            "dG_RNAduplex", "dG_RNAup", "contextpp_TargetScan",
                        ]
                        flags = [c for c in ext_only if c in lost_tbl.columns]
                        lost_tbl = lost_tbl[[c for c in lead if c in lost_tbl.columns] + flags]
                        c1, c2, c3 = st.columns(3)
                        c1.metric("Explorer lost", len(lost_syms))
                        c2.metric("Lost ∩ ≥1 external", int((support["n_external"] >= 1).sum()))
                        c3.metric(f"Shown (≥{min_ext})", len(lost_tbl))
                        st.dataframe(lost_tbl, width="stretch", height=320, hide_index=True)
                        st.download_button(
                            "⬇ Download lost targets + external support (CSV)",
                            lost_tbl.to_csv(index=False),
                            file_name=f"{mirna}_{db_state}_lost_external_min{min_ext}.csv",
                            mime="text/csv",
                            key="dl_lost_ext",
                        )
                    elif lost_syms and not ext_only:
                        st.info(
                            "Select at least one external database above to annotate "
                            f"the {len(lost_syms)} Explorer-lost genes."
                        )

st.markdown("---")
st.caption("Prediction: seed positions 2–8; unoxidized G pairs C, 8-oxoG (o8G) pairs A "
           "(Hoogsteen). Targets = TargetScan-style 6mer/7mer/8mer sites in human 3′UTRs "
           "(longest per gene, Ensembl). Enrichment = hypergeometric ORA, BH-corrected. "
           "Large seed-match target lists yield broad enrichment — interpret pathway shifts, "
           "not absolute significance, and validate candidates experimentally. "
           "Gene ID resolution uses a locally cached NCBI/HGNC map (no runtime API). "
           "External DB comparison uses local TargetScan/miRDB/DIANA/miRmap files when present "
           "plus the ENCORI open API; optional **TargetScan de novo** adds the precomputed "
           "ox-aware TargetScanS strong-site set. Explorer follows the sidebar prediction mode. An optional lost-gene "
           "list ranks Explorer losses by external WT support. Target tables are ranked by "
           "site type (8mer > 7mer-m8 > …).")
