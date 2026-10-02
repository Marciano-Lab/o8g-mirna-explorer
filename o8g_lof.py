"""Loss-of-function / state attrition helpers (miRDB-anchored + Explorer gained/lost)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from o8g_sections import SectionContext

ROOT = Path(__file__).resolve().parent


def _load_mirdb(mirna: str, score_min: float) -> set[str] | None:
    """Return miRDB symbols for mirna at score >= cutoff, or None if unavailable."""
    candidates = [
        ROOT / "mirdb_ref.parquet",
        ROOT / "paper" / "data" / "mirdb_ref.parquet",
        ROOT / "mirdb_custom_cache.db",
    ]
    for p in candidates:
        if not p.exists():
            continue
        if p.suffix == ".parquet":
            try:
                df = pd.read_parquet(p)
                cols = {c.lower(): c for c in df.columns}
                mir_col = cols.get("mirna") or cols.get("mir") or list(df.columns)[0]
                score_col = cols.get("score") or cols.get("mirdb_score")
                sym_col = cols.get("symbol") or cols.get("gene") or cols.get("gene_symbol")
                sub = df[df[mir_col].astype(str) == mirna]
                if score_col is not None:
                    sub = sub[sub[score_col] >= score_min]
                return set(sub[sym_col].astype(str).str.upper())
            except Exception:
                continue
        if p.suffix == ".db" or p.name.endswith(".db"):
            try:
                import sqlite3

                con = sqlite3.connect(p)
                try:
                    rows = con.execute(
                        "SELECT symbol FROM mirdb WHERE mirna=? AND score>=?",
                        (mirna, float(score_min)),
                    ).fetchall()
                except sqlite3.Error:
                    rows = []
                con.close()
                if rows:
                    return {str(r[0]).upper() for r in rows}
            except Exception:
                continue
    try:
        import o8g_refsets as refsets

        if refsets.available_tools().get("miRDB"):
            return {g.upper() for g in refsets.load_mirdb(mirna, score_min=float(score_min))}
    except Exception:
        pass
    return None


def state_summary_table(
    ctx: SectionContext, *, mirdb_score: float = 80.0, include_mirdb: bool = True
) -> pd.DataFrame:
    """Per-ox-state lost / gained / shared vs unmodified, plus optional miRDB attrition."""
    unmod = {g.upper() for g in ctx.strong_set("none")}
    baseline = _load_mirdb(ctx.mirna, float(mirdb_score)) if include_mirdb else None
    ox_labels = [s for s in ctx.state_labels if s != "none"]
    rows = []
    for lab in ox_labels:
        ox = {g.upper() for g in ctx.strong_set(lab)}
        lost = unmod - ox
        gained = ox - unmod
        shared = unmod & ox
        row = {
            "state": lab,
            "n_unmod": len(unmod),
            "n_oxidized": len(ox),
            "lost": len(lost),
            "gained": len(gained),
            "shared": len(shared),
            "retained_frac": (len(shared) / len(unmod)) if unmod else 0.0,
        }
        if baseline is not None:
            row["mirdb_baseline"] = len(baseline)
            row["mirdb_lost"] = len(baseline - ox)
            row["mirdb_retained"] = len(baseline & ox)
        rows.append(row)
    return pd.DataFrame(rows)


def render_state_summary(ctx: SectionContext) -> None:
    """Compact LoF / gained table for embedding in All states."""
    st.markdown("##### Target attrition across oxidation states")
    st.caption(
        "Lost / gained / shared are Explorer set-differences vs unmodified under the "
        "current precision mode. Optional **miRDB_lost** is attrition from a WT miRDB "
        "baseline (miRDB has no oxomiR catalog — that column is loss-only)."
    )
    want_mirdb = st.checkbox(
        "Add miRDB attrition column",
        value=False,
        key="allstates_mirdb",
        help="Reads the hosted miRDB catalog. Off by default so opening this tab "
             "does not load that database.",
    )
    score = 80.0
    if want_mirdb:
        score = st.slider(
            "miRDB score cutoff (for mirdb_lost column)",
            min_value=50,
            max_value=100,
            value=80,
            step=5,
            key="allstates_mirdb_score",
        )
    tab = state_summary_table(
        ctx, mirdb_score=float(score), include_mirdb=bool(want_mirdb)
    )
    if tab.empty:
        st.info("No oxidized states for this seed.")
        return
    show = tab.copy()
    st.dataframe(show, width="stretch", hide_index=True)
    chart_cols = [c for c in ("lost", "gained") if c in show.columns]
    if chart_cols:
        st.bar_chart(show.set_index("state")[chart_cols])
    st.caption(
        f"{ctx.mirna}: median lost={show['lost'].median():.0f}, "
        f"median gained={show['gained'].median():.0f} "
        f"(n_unmod={int(show['n_unmod'].iloc[0]) if len(show) else 0})."
    )


def render(ctx: SectionContext) -> None:
    """Standalone view kept for back-compat; prefer All states embedding."""
    st.markdown(f"#### Loss of function / state summary · `{ctx.precision_mode}`")
    render_state_summary(ctx)
