"""
o8g_thermo.py
=============
Three complementary binding metrics for explorer target tables.

1. **RNAduplex (ViennaRNA)** — ``dG_RNAduplex`` (kcal/mol)
   ``RNA.duplexfold(miRNA, local UTR window)``. More negative = stronger
   hybridization. Approximates RNAhybrid / RNAduplex seed-centered duplex energy.

2. **RNAup / IntaRNA-style** — ``dG_RNAup`` (kcal/mol)
   ``dG_RNAduplex + ED``, where ED is the opening energy of the seed region
   from unpaired probabilities (ViennaRNA partition function on the local
   window): ED = −RT Σ ln P_unpaired(i). Same energy decomposition used by
   RNAup / IntaRNA (hybridization + accessibility). More negative = stronger
   *and* more accessible.

3. **TargetScan context++** — ``contextpp_TargetScan``
   Official weighted context++ score from TargetScanHuman 8.0
   ``Conserved_Site_Context_Scores.txt`` (best / most negative per gene).
   Available for **unmodified** catalogs only; NaN under o8G states (TargetScan
   has no oxomiR motifs). Also emits ``context_analog`` from our lightweight
   scanner context++ *signs* when live-scanned.

o8G approximation for ViennaRNA (no native 8-oxoG parameter): oxidized seed
guanines are temporarily replaced with U so WC U:A stands in for Hoogsteen
o8G:A when folding. Documented as a proxy, not a chemical force-field.

The engine-internal ``binding_efficiency`` column in ``o8g_binding.py`` is
unchanged and remains the default sort key.
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "paper" / "data"
TS_CONTEXT = DATA / "Conserved_Site_Context_Scores.txt"

# Boltzmann (kcal/mol) at 37 °C
RT = 0.001987 * 310.15

# Local window around the 6mer start for duplex / accessibility
FLANK_5 = 12
FLANK_3 = 20
SEED_LEN = 7  # positions 2–8

_VIENNA_OK: bool | None = None


def vienna_available() -> bool:
    global _VIENNA_OK
    if _VIENNA_OK is None:
        try:
            import RNA  # noqa: F401

            _VIENNA_OK = True
        except Exception:
            _VIENNA_OK = False
    return bool(_VIENNA_OK)


def _to_rna(seq: str) -> str:
    return (seq or "").upper().replace("T", "U").replace("N", "N")


def mature_for_thermo(mature_dna: str, oxidized_positions: tuple[int, ...] | list[int] = ()) -> str:
    """Mature miRNA as RNA; oxidized Gs (seed pos 2–8) → U as o8G:A WC proxy."""
    rna = list(_to_rna(mature_dna))
    ox = set(int(p) for p in oxidized_positions)
    # mature DNA is 5'→3'; position 1 = index 0
    for pos in ox:
        i = pos - 1
        if 0 <= i < len(rna) and rna[i] == "G":
            rna[i] = "U"
    return "".join(rna)


def _local_window(utr: str, site_start: int) -> tuple[str, int]:
    """Return (window RNA, seed offset within window)."""
    lo = max(0, int(site_start) - FLANK_5)
    hi = min(len(utr), int(site_start) + SEED_LEN + FLANK_3)
    win = _to_rna(utr[lo:hi])
    seed_off = int(site_start) - lo
    return win, seed_off


def dg_rnaduplex(mir_rna: str, utr_window_rna: str) -> float:
    """RNAduplex-style hybridization free energy (kcal/mol)."""
    import RNA

    if not mir_rna or not utr_window_rna or "N" in mir_rna or "N" in utr_window_rna:
        return float("nan")
    try:
        return float(RNA.duplexfold(mir_rna, utr_window_rna).energy)
    except Exception:
        return float("nan")


def opening_energy(utr_window_rna: str, seed_off: int, seed_len: int = SEED_LEN) -> float:
    """ED ≈ −RT Σ ln P_unpaired over the seed (RNAup / IntaRNA accessibility term)."""
    import RNA

    n = len(utr_window_rna)
    if n < seed_len + 2 or seed_off < 0 or seed_off + seed_len > n:
        return float("nan")
    if "N" in utr_window_rna:
        return float("nan")
    try:
        fc = RNA.fold_compound(utr_window_rna)
        fc.pf()
        bpp = fc.bpp()  # 1-indexed; bpp[i][j] for i < j
    except Exception:
        return float("nan")

    ed = 0.0
    for i0 in range(seed_off, seed_off + seed_len):
        i = i0 + 1  # 1-based
        paired = 0.0
        for j in range(1, n + 1):
            if j == i:
                continue
            a, b = (i, j) if i < j else (j, i)
            try:
                paired += float(bpp[a][b])
            except Exception:
                pass
        pu = max(1e-12, 1.0 - paired)
        ed += -RT * math.log(pu)
    return float(ed)


def dg_rnaup(mir_rna: str, utr_window_rna: str, seed_off: int) -> float:
    """IntaRNA/RNAup-style total: hybridization + opening energy."""
    hyb = dg_rnaduplex(mir_rna, utr_window_rna)
    if hyb != hyb:  # NaN
        return float("nan")
    ed = opening_energy(utr_window_rna, seed_off)
    if ed != ed:
        return hyb
    return float(hyb + ed)


@lru_cache(maxsize=8)
def _targetscan_contextpp_table(mirna: str) -> pd.DataFrame:
    """Best (most negative) weighted context++ score per gene symbol for one miRNA."""
    from conservation import ensure_context_scores_txt

    ensure_context_scores_txt(TS_CONTEXT.parent)
    if not TS_CONTEXT.exists():
        return pd.DataFrame(columns=["symbol", "contextpp_TargetScan"])
    best: dict[str, float] = {}
    usecols = [
        "Gene Symbol",
        "miRNA",
        "context++ score",
        "weighted context++ score",
    ]
    for chunk in pd.read_csv(
        TS_CONTEXT, sep="\t", usecols=usecols, chunksize=250_000, low_memory=False
    ):
        sub = chunk[chunk["miRNA"] == mirna]
        if sub.empty:
            continue
        # prefer weighted context++; fall back to context++
        score = pd.to_numeric(sub["weighted context++ score"], errors="coerce")
        fallback = pd.to_numeric(sub["context++ score"], errors="coerce")
        score = score.fillna(fallback)
        for sym, sc in zip(sub["Gene Symbol"].astype(str), score):
            if sc != sc:
                continue
            if sym not in best or sc < best[sym]:
                best[sym] = float(sc)
    if not best:
        return pd.DataFrame(columns=["symbol", "contextpp_TargetScan"])
    return pd.DataFrame(
        {"symbol": list(best.keys()), "contextpp_TargetScan": list(best.values())}
    )


def targetscan_contextpp_map(mirna: str) -> pd.Series:
    tab = _targetscan_contextpp_table(mirna)
    if tab.empty:
        return pd.Series(dtype=float)
    return tab.set_index("symbol")["contextpp_TargetScan"]


def score_targets_thermo(
    df: pd.DataFrame,
    *,
    scanner,
    mature_dna: str,
    oxidized_positions: tuple[int, ...] | list[int] = (),
    mirna: str | None = None,
    is_unmodified: bool = True,
    rnaup_max_genes: int = 800,
) -> pd.DataFrame:
    """Attach ``dG_RNAduplex``, ``dG_RNAup``, ``contextpp_TargetScan`` columns.

    Requires ``site_start`` (6mer index in UTR) and a live ``scanner`` for Vienna
    metrics. TargetScan context++ is joined by symbol when ``mirna`` is set and
    ``is_unmodified``.

    RNAup opening energies use a partition function and are slower — for lists
    longer than ``rnaup_max_genes``, RNAup is computed for the top-N by site
    rank / score and left NaN elsewhere (duplex + TargetScan still fill all rows).
    """
    out = df.copy().reset_index(drop=True)
    n = len(out)
    duplex = np.full(n, np.nan)
    rnaup = np.full(n, np.nan)

    do_rnaup_idx: set[int] = set(range(n))
    if n > rnaup_max_genes:
        # Prefer strong sites / high multiplicity for the expensive accessibility term
        rank = out["site_rank"] if "site_rank" in out.columns else pd.Series(0, index=out.index)
        sc = out["score"] if "score" in out.columns else pd.Series(0.0, index=out.index)
        order = (
            pd.DataFrame({"rank": rank, "score": sc})
            .sort_values(["rank", "score"], ascending=[False, False])
            .head(rnaup_max_genes)
            .index
        )
        do_rnaup_idx = set(int(i) for i in order)

    if scanner is not None and vienna_available() and mature_dna and n and "site_start" in out.columns:
        mir = mature_for_thermo(mature_dna, oxidized_positions)
        sym_to_gidx = {}
        if "gene_idx" not in out.columns and hasattr(scanner, "symbols"):
            sym_to_gidx = {s: i for i, s in enumerate(scanner.symbols)}

        for i in range(n):
            try:
                if "gene_idx" in out.columns and pd.notna(out.at[i, "gene_idx"]):
                    gi = int(out.at[i, "gene_idx"])
                else:
                    gi = sym_to_gidx.get(str(out.at[i, "symbol"]), None)
                if gi is None:
                    continue
                if pd.isna(out.at[i, "site_start"]):
                    continue
                site_start = int(out.at[i, "site_start"])
                utr = scanner.utrs[gi]
                win, seed_off = _local_window(utr, site_start)
                duplex[i] = dg_rnaduplex(mir, win)
                if i in do_rnaup_idx:
                    rnaup[i] = dg_rnaup(mir, win, seed_off)
            except Exception:
                continue

    out["dG_RNAduplex"] = duplex
    out["dG_RNAup"] = rnaup

    if mirna and is_unmodified:
        ts = targetscan_contextpp_map(mirna)
        out["contextpp_TargetScan"] = (
            out["symbol"].map(ts) if "symbol" in out.columns else np.nan
        )
    else:
        out["contextpp_TargetScan"] = np.nan

    if "context_score" in out.columns:
        out["context_analog"] = out["context_score"]

    return out


METRIC_CAPTION = (
    "Optional energetics: **dG_RNAduplex** = ViennaRNA duplexfold (kcal/mol; more negative → stronger). "
    "**dG_RNAup** = duplex + seed opening energy (RNAup/IntaRNA-style; more negative → stronger & accessible). "
    "**contextpp_TargetScan** = official TargetScan 8 weighted context++ (unmodified only; more negative → stronger). "
    "o8G duplex uses G→U proxy at oxidized seed positions. "
    "Default ranking is by site type (8mer > 7mer-m8 > 7mer-A1 > 6mer)."
)
