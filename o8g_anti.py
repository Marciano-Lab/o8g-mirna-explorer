"""Antagomir design view — state discrimination + mRNA collateral hybridization."""
from __future__ import annotations

import importlib

import pandas as pd
import streamlit as st

import o8g_energy as _o8g_energy

_o8g_energy = importlib.reload(_o8g_energy)
from o8g_energy import (  # noqa: E402  — reload-first for Streamlit hot-reload
    collateral_mrna_offtargets,
    design_antagomir,
    feasibility_label,
    position2_warning,
    rank_single_g_designs,
    utr_parquet_path,
)
from o8g_sections import SectionContext


def render(ctx: SectionContext) -> None:
    st.markdown(f"#### Antagomir design · `{ctx.precision_mode}`")
    st.warning(
        "Oxidized-state ViennaRNA ΔG uses a **G→U proxy** (no o8G:A Hoogsteen parameter). "
        "`dG_normal` is the reliable arm; `dG_oxo` / `ddG` / fold-preference are "
        "**directional proxies**."
    )
    st.info(
        "An antagomir **is not AGO-loaded and has no seed**. This view scores "
        "**state discrimination** (oxo-selective oligo vs mature) and separately "
        "lists **mRNA 3′UTR hybridization off-targets** (collateral).\n\n"
        "**Why WT-selective is omitted:** a full-length reverse-complement only "
        "changes the base(s) opposite oxidized G(s). With 1–3 of ~22 nt differing, "
        "|ΔΔG| stays small and fold-preference near 1× **by construction**, so the "
        "UI always designs an **oxo-selective** oligo (A opposite oxidized G)."
    )

    ox_labs = [s for s in ctx.state_labels if s != "none"] or list(ctx.state_labels)
    # Prefer interior single-G default when available
    default_i = 0
    for pref in ("o8G@7", "o8G@6", "o8G@5", "o8G@4", "o8G@3"):
        if pref in ox_labs:
            default_i = ox_labs.index(pref)
            break
    ox_lab = st.selectbox(
        "Design against oxidation state",
        ox_labs,
        index=default_i,
        key="anti_state",
    )
    positions: tuple[int, ...] = ()
    if ox_lab and ox_lab != "none":
        positions = tuple(
            int(x) for x in ox_lab.replace("o8G@", "").split(",") if x.strip().isdigit()
        )

    if position2_warning(positions):
        st.error(
            "Position-2-only oxidation is a **poor antagomir target** in literature-style "
            "screens (weak fold preference vs interior Gs). Prefer an interior o8G state."
        )

    design = design_antagomir(
        ctx.info.get("seq_dna", ""),
        positions,
        kind="oxo-selective",
    )
    st.caption("Chemistry: **oxo-selective** (A opposite o8G; C opposite other Gs).")
    st.code(design.sequence)
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("dG_normal", f"{design.dG_normal:.2f}" if design.dG_normal is not None else "—")
    m2.metric("dG_oxo (proxy)", f"{design.dG_oxo:.2f}" if design.dG_oxo is not None else "—")
    m3.metric("ΔΔG (ox−wt)", f"{design.ddG:.2f}" if design.ddG is not None else "—")
    m4.metric(
        "Fold → oxidized",
        f"{design.fold_preference:.1f}×" if design.fold_preference is not None else "—",
    )
    st.caption(design.note)

    level, msg = feasibility_label(design.fold_preference)
    if level == "good":
        st.success(msg)
    elif level == "modest":
        st.warning(msg)
    elif level == "weak":
        st.error(msg)
    else:
        st.info(msg)

    # Rank single-G states so users see which G is least-bad
    with st.expander("Rank single-G states by oxo-selective fold preference", expanded=True):
        ranked = rank_single_g_designs(ctx.info.get("seq_dna", ""))
        if not ranked:
            st.caption("No seed guanines on this mature.")
        else:
            tab = pd.DataFrame(
                [
                    {
                        "ox_label": f"o8G@{d.target_state}",
                        "fold→ox": d.fold_preference,
                        "ΔΔG": d.ddG,
                        "dG_wt": d.dG_normal,
                        "dG_ox": d.dG_oxo,
                        "oligo": d.sequence,
                    }
                    for d in ranked
                ]
            )
            st.dataframe(tab, hide_index=True, width="stretch")
            best = ranked[0]
            st.caption(
                f"Best single-G on this mature: **o8G@{best.target_state}** "
                f"({best.fold_preference:.1f}×). Values <10× are still weak absolute "
                "discriminators — use for ranking states, not as a Kd claim."
            )

    st.markdown("##### Collateral mRNA off-targets (3′UTR hybridization)")
    st.caption(
        "Genes whose 3′UTR contains a long contiguous stretch complementary to the "
        "antagomir (sequence the oligo can bind). This is **mRNA off-target risk**, "
        "not a scan of other miRNAs."
    )
    min_match = st.slider(
        "Minimum contiguous match (nt)",
        min_value=8,
        max_value=16,
        value=10,
        key="anti_collateral_k",
    )
    try:
        with st.spinner("Scanning 3′UTRs for oligo hybridization sites…"):
            hit_df = collateral_mrna_offtargets(
                design.sequence,
                ctx.scanner,
                min_match=int(min_match),
                top_n=80,
            )
        if hit_df.empty and ctx.scanner is None and utr_parquet_path() is None:
            st.caption(
                "UTR table unavailable (`utr3_human.parquet`). Cannot predict mRNA off-targets."
            )
            return
        st.dataframe(hit_df, hide_index=True, width="stretch")
        st.caption(
            f"{len(hit_df)} genes with ≥{min_match}-nt contiguous match "
            "(showing top 80 by match length)."
        )
        if len(hit_df):
            st.download_button(
                "⬇ Download mRNA collateral hits (CSV)",
                hit_df.to_csv(index=False),
                file_name=f"{ctx.mirna}_antagomir_mrna_offtargets.csv",
                mime="text/csv",
                key="anti_collateral_dl",
            )
    except Exception as e:
        st.caption(f"Collateral mRNA scan skipped ({type(e).__name__}: {e}).")
