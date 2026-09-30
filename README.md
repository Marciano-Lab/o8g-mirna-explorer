# o8G-miRNA Retargeting Explorer

Predict how **8-oxoguanine (o8G) oxidation of a microRNA seed** rewires its mRNA
target list, compare oxidation states by pathway enrichment, and cross-check
predictions against external **unmodified**-miRNA catalogs.

**Live demo:** [oxomir.samnti.com](https://oxomir.samnti.com)
([o8g.samnti.com](https://o8g.samnti.com))

**Repository:** [Marciano-Lab/o8g-mirna-explorer](https://github.com/Marciano-Lab/o8g-mirna-explorer).
Clone and pull from this repository. It is the project folder in the
[Marciano-Lab](https://github.com/Marciano-Lab) organization.

> **Biology.** A seed guanine normally pairs **C**; oxidized to 8-oxoG it pairs
> **A** (Hoogsteen). Each of the *k* guanines in a seed can independently be
> normal or o8G, giving **2^k** target repertoires per miRNA.
>
> **What is unique here.** Only this explorer **rewrites seed motifs under
> oxidation** and rescans 3′UTRs (or serves precomputed states). TargetScan /
> miRDB / DIANA / miRmap / ENCORI / miRTarBase are **lookup catalogs for
> unmodified matures** — they do not accept o8G sequences. Optional ViennaRNA
> energetics use a **G→U pairing-code proxy** on oxidized matures (no native
> o8G parameters); site gain/loss uses biological **G·C → o8G·A**, not G→U.

## Quick start

```bash
git clone https://github.com/Marciano-Lab/o8g-mirna-explorer.git
cd o8g-mirna-explorer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./start.sh
# open http://localhost:8501
```

### Required data (not in git)

Large SQLite / parquet assets are gitignored. Place them next to the code
(or rebuild — see below). Views that need a missing file report which path
is required.

| Asset | Needed for |
|---|---|
| `o8g_targets.db` | All target browsing |
| `genesets/` | Pathway + TF enrichment (Enrichr GMTs) |
| `gene_aliases.parquet` + `o8g_reverse.db` | Gene → miRNA/oxomiR (+ sidebar gene jump) |
| `mirdb_ref.parquet` | Optional miRDB attrition column in **All states**; `python scripts/build_mirdb_ref.py` |
| `utr3_human.parquet` | **TargetScan de novo** precision mode; optional live UTR rescans / ViennaRNA windows |
| `paper/data/Predicted_Targets_Info.default_predictions.txt` | **TargetScan** precision mode (WT catalog anchor) |
| `paper/data/Conserved_Family_Info.txt` | **Consensus** precision mode |
| `paper/data/Conserved_Site_Context_Scores.txt` | TargetScan context++ columns |
| `models/oboe_rnabert/checkpoint/` | OBOE RNABERT oxo-G prior (optional; GC fallback without it) |
| ViennaRNA (Python `ViennaRNA`) | Opt-in energetics / antagomir ΔG |
| `pip install -r requirements-oboe.txt` | Optional torch stack for OBOE RNABERT inference |

Fetch extra TF motif libraries (if missing):

```bash
python scripts/fetch_tf_genesets.py
```

## App views

Single Streamlit app (`app.py`). Navigation is a **radio switch** (not tabs) so
only the active section runs on each rerun.

| View | Description | Optional data |
|---|---|---|
| **Single state — targets** | One oxidation state; site-type ranking; optional energetics; pathway gene peek | `utr3_human.parquet`, ViennaRNA |
| **Compare two states** | Venn + lost/gained/shared + pathway volcano / within-pool control + gene lists | `genesets/` |
| **All states** | Per-state lost/**gained**/shared (+ optional miRDB attrition); pathway heatmap/dotplot | `genesets/`, `mirdb_ref.parquet` |
| **Overlap (Venn/UpSet)** | 2–3 states → Venn; 4+ → UpSet; region gene table | — |
| **Transcription factors** | Enriched TF regulons (TRANSFAC/JASPAR, JASPAR 2025, Genome Browser PWMs, ChEA/ENCODE/TRRUST) | TF GMTs; `scripts/fetch_tf_genesets.py` |
| **Antagomir design** | Oxo-selective oligo discrimination (G→U proxy); ranks single-G states; explains weak full-length folds | ViennaRNA optional |
| **RNA-seq / DEG upload** | Exploratory UP↔lost / DOWN↔gained concordance for a miRNA panel (CSV/TSV/XLSX; multi-sheet Excel picker) | session upload only |
| **Gene → miRNA/oxomiR** | Reverse lookup (symbol / Ensembl / Entrez); also reachable from sidebar **Select gene** | `gene_aliases.parquet`, `o8g_reverse.db` |
| **External DB comparison** | Explorer (sidebar prediction mode) vs WT catalogs (TargetScan / miRDB / DIANA / miRmap / ENCORI / miRTarBase); optional **TargetScan de novo** set for the selected ox state | local ref files; `utr3_human.parquet` for de novo |

Sidebar also has **Select gene** (typeahead → opens Gene → miRNA/oxomiR) and an
**OBOE oxo-G prior** — local fine-tuned RNABERT (Xia et al. Figshare data; held-out
test AUC ~0.98) ranks which seed Gs / ox states are sequence-plausible; GC motif
fallback if the checkpoint or torch stack is missing.

### Prediction modes

**Sequence-based (low stringency)** / **Sequence-based (high stringency)** /
**TargetScan** / **TargetScan de novo** / **Consensus** (sidebar; UI label
**Prediction mode**). Legacy labels `Sensitive` / `Stringent` still resolve in
code. All of these apply to Explorer target lists used across every view
(including External DB’s Explorer arm, DEG concordance, overlap, etc.).

| Mode | Rule |
|---|---|
| **Sequence-based (low stringency)** | Strong sites (7mer-m8 + 8mer) from the explorer DB. Discovery default. Formerly Sensitive. |
| **Sequence-based (high stringency)** | 8mer only. Formerly Stringent. |
| **TargetScan** | Rank ≥ 3 ∩ TargetScanHuman 8.0 *predicted* strong sites on the **unmodified** baseline (`Predicted_Targets_Info`). Catalog lookup — WT anchor for lost/gained; oxidized *display* lists stay rank≥3. let-7*-5p / miR-98-5p → family `let-7-5p/98-5p`. |
| **TargetScan de novo** | Live TargetScanS site finding on **both** WT and oxidized seeds (o8G encoded as **T** for WC search so sites require **A**). Uses `TargetScanner` / vendored `third_party/targetscan/targetscan_70.pl`; **not** the web catalog. Requires `utr3_human.parquet`. |
| **Consensus** | Rank ≥ 3 ∩ TargetScan *conserved* families on the unmodified baseline. |

**Catalog TargetScan vs TargetScan de novo.** On an **oxidized** state, gene
counts often **match** because catalog mode only intersects the *unmodified*
baseline for lost/gained; the oxidized arm is still rank≥3 — the same site types
de novo uses. That is **not a bug**. Prefer External DB: catalog WT TargetScan
vs de novo ox for a true algorithm contrast. Low-stringency sequence-based vs
de novo also usually match (same TargetScanS rules on the same 3′UTRs).

If TargetScan/Consensus data are missing **or** the mature has no TS entry
(e.g. some poorly annotated miRNAs), the app **refuses an empty baseline** and
falls back to sequence-based (low stringency) — it does not silently claim
“everything is gained.”

Filters always apply to each oxidation state **before** gained/lost.

## Binding metrics

| Column | Source | Notes |
|---|---|---|
| `site_type` / `site_rank` | Engine (8mer > 7mer-m8 > 7mer-A1 > 6mer) | Default UI ranking |
| `dG_RNAduplex` | ViennaRNA `duplexfold` | Opt-in |
| `dG_RNAup` | Duplex + seed opening (RNAup-style) | Opt-in; slower |
| `contextpp_TargetScan` | TargetScan 8 weighted context++ | Unmodified mature only |

**o8G folding caveat.** ViennaRNA has no native 8-oxoG parameters. For oxidized
states, oxidized seed Gs are temporarily replaced with **U** so Watson–Crick
U·A stands in for Hoogsteen o8G·A. That is a **pairing-code proxy for ranking**,
not a chemical force field. Site-type / gained–lost calls use **G·C → o8G·A**
and do **not** depend on the U proxy.

## Repository layout

| Path | Role |
|---|---|
| `app.py` | Streamlit explorer |
| `o8g_engine.py` | Seed extraction, 2^k state enumeration, motifs (G·C → o8G·A) |
| `o8g_scanner.py` | k-mer index over 3′UTRs |
| `o8g_db.py` | Read-only SQLite accessor (+ reverse gene query) |
| `o8g_precision.py` | Prediction modes: sequence-based low/high / TargetScan / de novo / Consensus |
| `o8g_ts_denovo.py` | Offline TargetScanS de novo oxomiR site finding |
| `third_party/targetscan/` | Vendored `targetscan_70.pl` (Bartel) |
| `o8g_sections.py` | Section context + dispatch |
| `o8g_venn.py`, `o8g_tf.py`, `o8g_lof.py`, `o8g_anti.py`, `o8g_energy.py` | Views / helpers |
| `o8g_deg.py`, `o8g_deg_score.py`, `o8g_deg_upload.py` | DEG upload parse + scoring + UI |
| `o8g_oboe.py` / `o8g_oboe_model.py` | OBOE RNABERT site prior + ranking |
| `models/oboe_rnabert/` | Fine-tuned checkpoint (from Figshare CSVs) |
| `third_party/oboe/` | Vendored Xia et al. Figshare code + data |
| `scripts/train_oboe_rnabert.py` | Retrain OBOE RNABERT |
| `requirements-oboe.txt` | Optional torch/transformers stack |
| `o8g_thermo.py` | ViennaRNA + TargetScan context++ |
| `o8g_pubthermo.py` | Duplex / RNAup / context++ for Gene → miRNA |
| `o8g_refsets.py` | External catalog loaders (TargetScan family resolution, etc.) |
| `o8g_enrich.py` | Hypergeometric ORA + BH (+ within-pool; pathway `genes` column) |
| `o8g_plots.py` | Volcano, heatmap, diverging bar, dot plot, Venn, UpSet |
| `o8g_genes.py` | Local gene ID resolver (NCBI/HGNC map) |
| `conservation.py` | TargetScan conserved-family helpers for Consensus |
| `precompute_db.py` | Rebuild `o8g_targets.db` |
| `start.sh` | Local Streamlit launcher |
| `scripts/fetch_tf_genesets.py` | Download Enrichr TF/motif GMTs |
| `scripts/build_mirdb_ref.py` | Build `mirdb_ref.parquet` |

## Rebuild

```bash
python scripts/fetch_mirbase.py
python scripts/fetch_utr3.py       # slow if not cached
python precompute_db.py
python scripts/validate_db.py
```

### Gene aliases + reverse index

```bash
python scripts/fetch_gene_aliases.py   # NCBI gene_info + gene2ensembl + HGNC
python scripts/build_reverse_index.py  # invert o8g_targets.db → o8g_reverse.db
```

Smoke check (HDAC4 targeted by unmodified miR-1 / seed `GGAATGT`, **absent** at
`o8G@7`):

```python
from o8g_genes import GeneResolver
from o8g_db import TargetDB

hit = GeneResolver().resolve("HDAC4").hits[0]
df = TargetDB("o8g_targets.db").states_targeting_gene(hit.gene_idx)
mir1 = df[df["seed"] == "GGAATGT"]
assert (mir1["state_label"] == "none").any()
assert not (mir1["state_label"] == "o8G@7").any()
```

## Programmatic use

```python
from o8g_db import TargetDB
from o8g_precision import PrecisionConfig
from o8g_scanner import TargetScanner

db = TargetDB("o8g_targets.db")
info = db.mirna_info("hsa-miR-1-3p")
parts = db.retarget_partition(
    info["seed"], "o8G@7", PrecisionConfig.from_mode("Sequence-based (low stringency)"), mirna=info["mirna"]
)
print({k: len(v) for k, v in parts.items()})  # unmod / oxid / lost / gained / shared

# TargetScan de novo (live UTR scan; needs utr3_human.parquet)
scanner = TargetScanner.from_parquet("utr3_human.parquet")
parts_dn = db.retarget_partition(
    info["seed"],
    "o8G@7",
    PrecisionConfig.from_mode("TargetScan de novo"),
    mirna=info["mirna"],
    scanner=scanner,
)

# reverse: which oxomiR states target HDAC4?
from o8g_genes import GeneResolver
g = GeneResolver().resolve("HDAC4").hits[0]
hits = db.states_targeting_gene(g.gene_idx)
```

## Caveats

- Seed-match prediction produces large target lists; treat pathway output as
  **direction of retargeting**, not validated targets.
- Reverse lookup uses the same strong-site (7mer-m8 + 8mer) precomputed DB —
  not live TargetScan/miRDB consensus.
- External databases catalog **unmodified** mature miRNAs only; only Explorer
  (and optional TargetScan de novo) change with o8G state.
- **TargetScan** (catalog) / **Consensus** need TS8 files and a mature that TS
  actually annotates; otherwise the UI falls back to sequence-based (low stringency).
- **TargetScan de novo** needs `utr3_human.parquet`; first UI selection builds
  the UTR index once (~30–40s). Optional Perl backend:
  `O8G_TS_DENOVO_BACKEND=perl`.
- Oxidized-state ViennaRNA ΔG uses a **G→U proxy**. Prefer site type and
  gold/null benchmarks for claims; use duplex numbers exploratorily.
- OBOE scores are **site-likelihood priors** for ranking ox states, not causal
  LoF calls (window-classifier AUC ≠ per-seed-G ranking accuracy).
- DEG upload is an **exploratory live-site** concordance module (session-only
  uploads); not a DESeq2 pipeline and not a primary Database Issue claim.

## License

See repository metadata / `LICENSE` if present.
