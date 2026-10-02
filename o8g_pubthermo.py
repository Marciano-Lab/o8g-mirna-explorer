"""
o8g_pubthermo.py
================
Publication-safe energetics for explorer tables (Gene → miRNA and others).

Uses the same tools as ``mirna_target_viewer`` — no invented NN tables:

* ``dG_hybrid`` — RNAduplex CLI or ViennaRNA ``RNA.duplexfold``
  (Lorenz et al. 2011; Rehmsmeier et al. 2004; Turner 2004 params)
* ``dG_open`` / ``ddG`` — RNAup (Mückstein et al. Bioinformatics 2006)
* ``contextpp_TargetScan`` — TargetScanHuman 8 weighted context++ when the
  local Conserved_Site_Context_Scores table is present (Agarwal eLife 2015)

o8G note: oxidized seed Gs are replaced with U for folding only (WC proxy for
Hoogsteen o8G:A); documented in provenance, not a chemical force field.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "paper" / "data"
TS_CONTEXT = DATA / "Conserved_Site_Context_Scores.txt"

# Prefer the conda mirna_viewer toolchain when present
_CONDA_BIN = Path("/opt/homebrew/anaconda3/envs/mirna_viewer/bin")
if _CONDA_BIN.is_dir():
    os.environ["PATH"] = str(_CONDA_BIN) + os.pathsep + os.environ.get("PATH", "")

FLANK_5 = 12
FLANK_3 = 20


def _to_rna(seq: str) -> str:
    return (seq or "").upper().replace("T", "U")


def mature_for_thermo(mature_dna: str, oxidized_positions: str | tuple | list = ()) -> str:
    rna = list(_to_rna(mature_dna))
    ox: set[int] = set()
    if isinstance(oxidized_positions, str):
        ox = {int(x) for x in oxidized_positions.split(",") if x.strip().isdigit()}
    else:
        ox = {int(x) for x in oxidized_positions}
    for pos in ox:
        i = pos - 1
        if 0 <= i < len(rna) and rna[i] == "G":
            rna[i] = "U"
    return "".join(rna)


def _which(name: str) -> Optional[str]:
    return shutil.which(name)


@lru_cache(maxsize=8192)
def rnaduplex_dg(mir_rna: str, tgt_rna: str) -> Optional[float]:
    """ΔG_hybrid (kcal/mol) via RNAduplex or ViennaRNA duplexfold."""
    if not mir_rna or not tgt_rna or "N" in mir_rna or "N" in tgt_rna:
        return None
    path = _which("RNAduplex")
    if path:
        try:
            proc = subprocess.run(
                [path, "-s"],
                input=f"{mir_rna}\n{tgt_rna}\n",
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            raw = (proc.stdout or "") + (proc.stderr or "")
            m = re.search(r"\(\s*([-+]?\d+\.\d+)\s*\)", raw)
            if m:
                return float(m.group(1))
        except Exception:
            pass
    try:
        import RNA

        return float(RNA.duplexfold(mir_rna, tgt_rna).energy)
    except Exception:
        return None


@lru_cache(maxsize=4096)
def rnaup_open_ddg(mir_rna: str, tgt_rna: str) -> tuple[Optional[float], Optional[float]]:
    """Return (dG_open, ddG=dGtot) from RNAup, or (None, None)."""
    path = _which("RNAup")
    if not path or not mir_rna or not tgt_rna:
        return None, None
    try:
        proc = subprocess.run(
            [path, "-b", "-o"],
            input=f"{tgt_rna}\n{mir_rna}\n",
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
            cwd="/tmp",
        )
        raw = (proc.stdout or "") + "\n" + (proc.stderr or "")
        m = re.search(
            r"\(\s*(-?\d+\.\d+)\s*=\s*(-?\d+\.\d+)\s*\+\s*(-?\d+\.\d+)(?:\s*\+\s*(-?\d+\.\d+))?\s*\)",
            raw,
        )
        if not m:
            return None, None
        dG_tot = float(m.group(1))
        o1 = float(m.group(3))
        o2 = float(m.group(4)) if m.group(4) is not None else 0.0
        return o1 + o2, dG_tot
    except Exception:
        return None, None


@lru_cache(maxsize=512)
def _ts_contextpp_map(mirna: str) -> dict[str, float]:
    """miRNA → {gene_symbol: best weighted context++}. Prefer `_ts_contextpp_for_gene`."""
    from conservation import ensure_context_scores_txt

    ensure_context_scores_txt(TS_CONTEXT.parent)
    if not TS_CONTEXT.exists() or not mirna:
        return {}
    best: dict[str, float] = {}
    usecols = ["Gene Symbol", "miRNA", "context++ score", "weighted context++ score"]
    for chunk in pd.read_csv(
        TS_CONTEXT, sep="\t", usecols=usecols, chunksize=250_000, low_memory=False
    ):
        sub = chunk[chunk["miRNA"] == mirna]
        if sub.empty:
            continue
        score = pd.to_numeric(sub["weighted context++ score"], errors="coerce")
        score = score.fillna(pd.to_numeric(sub["context++ score"], errors="coerce"))
        for sym, sc in zip(sub["Gene Symbol"].astype(str), score):
            if sc != sc:
                continue
            if sym not in best or sc < best[sym]:
                best[sym] = float(sc)
    return best


@lru_cache(maxsize=256)
def _ts_contextpp_for_gene(gene_symbol: str, mirnas_key: tuple[str, ...]) -> dict[str, float]:
    """One pass over TargetScan tables: miRNA → context++ for a single gene.

    Scanning the 138 MB Conserved_Site_Context_Scores once per gene (not once
    per miRNA) keeps Gene→miRNA energetics interactive even with hundreds of
    unmodified hits.
    """
    from conservation import ensure_context_scores_txt

    ensure_context_scores_txt(TS_CONTEXT.parent)
    if not TS_CONTEXT.exists() or not gene_symbol or not mirnas_key:
        return {}
    want = set(mirnas_key)
    best: dict[str, float] = {}
    usecols = ["Gene Symbol", "miRNA", "context++ score", "weighted context++ score"]
    for chunk in pd.read_csv(
        TS_CONTEXT, sep="\t", usecols=usecols, chunksize=250_000, low_memory=False
    ):
        sub = chunk[
            (chunk["Gene Symbol"].astype(str) == gene_symbol)
            & (chunk["miRNA"].isin(want))
        ]
        if sub.empty:
            continue
        score = pd.to_numeric(sub["weighted context++ score"], errors="coerce")
        score = score.fillna(pd.to_numeric(sub["context++ score"], errors="coerce"))
        for mir, sc in zip(sub["miRNA"].astype(str), score):
            if sc != sc:
                continue
            if mir not in best or sc < best[mir]:
                best[mir] = float(sc)
    return best


@lru_cache(maxsize=512)
def _utr_for_gene(gene_symbol: str = "", gene_idx: int = -1) -> str:
    """Fetch one 3′UTR without building the full TargetScanner index."""
    parquet = ROOT / "utr3_human.parquet"
    if not parquet.exists():
        return ""
    try:
        import pyarrow.parquet as pq

        if gene_symbol:
            table = pq.read_table(
                parquet,
                filters=[("symbol", "==", gene_symbol)],
                columns=["utr3"],
            )
            if table.num_rows:
                return str(table.column("utr3")[0].as_py() or "")
        if gene_idx is not None and int(gene_idx) >= 0:
            # Row order matches TargetDB gene_idx
            df = pd.read_parquet(parquet, columns=["utr3"])
            idx = int(gene_idx)
            if 0 <= idx < len(df):
                return str(df.iloc[idx]["utr3"] or "")
    except Exception:
        return ""
    return ""


def _local_window(utr: str, site_start: int) -> str:
    lo = max(0, int(site_start) - FLANK_5)
    hi = min(len(utr), int(site_start) + 8 + FLANK_3)
    return _to_rna(utr[lo:hi])


def _site_start_in_utr(utr: str, motif_8mer: str, motif_7m8: str, site_rank: int) -> Optional[int]:
    """Locate best site motif on the UTR (DNA alphabet)."""
    utr_d = (utr or "").upper().replace("U", "T")
    if int(site_rank) >= 4 and motif_8mer:
        m = str(motif_8mer).upper().replace("U", "T")
        i = utr_d.find(m)
        if i >= 0:
            return i
    if motif_7m8:
        m = str(motif_7m8).upper().replace("U", "T")
        i = utr_d.find(m)
        if i >= 0:
            # 7mer-m8 includes m8 base; 6mer core starts at i+1 for geometry,
            # but for windowing the motif start is fine.
            return i
    return None


def annotate_gene_mirna_hits(
    hits: pd.DataFrame,
    *,
    gene_symbol: str,
    gene_idx: int,
    scanner=None,
    db,
    max_thermo_rows: int = 200,
    include_rnaup: bool = True,
) -> pd.DataFrame:
    """Add dG_hybrid, dG_open, ddG, contextpp_TargetScan columns.

    Thermodynamics are scored for at most ``max_thermo_rows`` (filtered table
    order). Remaining rows get null energetics; context++ is filled for all
    unmodified rows when TargetScan tables exist (one file pass per gene).
    """
    out = hits.copy().reset_index(drop=True)
    n = len(out)
    out["dG_hybrid"] = pd.NA
    out["dG_open"] = pd.NA
    out["ddG"] = pd.NA
    out["contextpp_TargetScan"] = pd.NA
    if n == 0:
        return out

    # TargetScan context++ (unmodified only) — single pass over the big TSV
    mirnas_none = (
        tuple(
            sorted(
                str(m)
                for m in out.loc[out["state_label"] == "none", "mirna"].dropna().unique()
            )
        )
        if "state_label" in out.columns
        else ()
    )
    ts_for_gene = _ts_contextpp_for_gene(str(gene_symbol), mirnas_none)
    for i in range(n):
        if str(out.at[i, "state_label"]) != "none":
            continue
        sc = ts_for_gene.get(str(out.at[i, "mirna"]))
        if sc is not None:
            out.at[i, "contextpp_TargetScan"] = sc

    utr = ""
    if scanner is not None:
        try:
            utr = scanner.utrs[int(gene_idx)]
        except Exception:
            utr = ""
    if not utr:
        utr = _utr_for_gene(str(gene_symbol), int(gene_idx))
    if not utr:
        return out

    mature_cache: dict[str, str] = {}

    def mature(mir: str) -> str:
        if mir not in mature_cache:
            try:
                mature_cache[mir] = str(db.mirna_info(mir)["seq_dna"])
            except Exception:
                mature_cache[mir] = ""
        return mature_cache[mir]

    limit = min(n, max(0, int(max_thermo_rows)))
    for i in range(limit):
        mir = str(out.at[i, "mirna"])
        mat = mature(mir)
        if not mat:
            continue
        ox = out.at[i, "oxidized_positions"] if "oxidized_positions" in out.columns else ""
        mir_rna = mature_for_thermo(
            mat, ox if str(out.at[i, "state_label"]) != "none" else ""
        )
        site_start = _site_start_in_utr(
            utr,
            str(out.at[i, "motif_8mer"] or ""),
            str(out.at[i, "motif_7mer_m8"] or ""),
            int(out.at[i, "site_rank"]),
        )
        if site_start is None:
            win = _to_rna(utr[: min(len(utr), 80)])
        else:
            win = _local_window(utr, site_start)
        if len(win) < 8:
            continue
        hyb = rnaduplex_dg(mir_rna, win)
        if hyb is not None:
            out.at[i, "dG_hybrid"] = round(hyb, 2)
        if include_rnaup:
            opn, tot = rnaup_open_ddg(mir_rna, win)
            if opn is not None:
                out.at[i, "dG_open"] = round(opn, 2)
            if tot is not None:
                out.at[i, "ddG"] = round(tot, 2)

    return out


PROVENANCE_CAPTION = (
    "**dG_hybrid** = RNAduplex / ViennaRNA duplexfold (Lorenz 2011; Turner 2004). "
    "**dG_open** / **ddG** = RNAup (Mückstein 2006; ddG = dGtot). "
    "**contextpp_TargetScan** = TargetScan 8 weighted context++ (unmodified only; "
    "Agarwal eLife 2015). Thermodynamics scored for the first *N* rows of the "
    "filtered table (cached). o8G folding uses G→U at oxidized seed positions."
)
