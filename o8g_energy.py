"""Antagomir / oligo energetics helpers.

ViennaRNA has no native o8G:A Hoogsteen parameters. Oxidized-state ΔG uses a
G→U proxy so Watson–Crick U·A stands in for o8G·A. ``dG`` against the
unmodified mature is trustworthy; oxidized numbers are directional with
placeholder chemistry.

Why many designs look "weak"
----------------------------
A full-length reverse-complement of the mature only differs at the oxidized
position(s). With 1–3 of ~22 nt changed, |ΔΔG| is small by construction and
fold-preference stays near 1×. Discrimination requires an **oxo-selective**
oligo (A opposite oxidized G, not C) and is still a proxy, not a wet-lab Kd.
"""
from __future__ import annotations

from dataclasses import dataclass

from o8g_engine import clean_seq, g_positions, extract_seed


@dataclass
class OligoDesign:
    sequence: str
    kind: str  # "oxo-selective" | "wt-selective"
    target_state: str
    dG_normal: float | None
    dG_oxo: float | None
    ddG: float | None  # dG_oxo - dG_normal; negative ⇒ prefers oxidized
    fold_preference: float | None  # >1 ⇒ prefers oxidized form
    note: str


def _to_rna(seq: str) -> str:
    return clean_seq(seq).replace("T", "U")


def _rc_rna(rna: str) -> list[str]:
    comp = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}
    return [comp.get(b, "N") for b in reversed(rna)]


def _oxo_selective_rc(mature_dna: str, oxidized_positions: tuple[int, ...]) -> str:
    """RC of mature with A (not C) opposite oxidized Gs — prefers o8G:A pairing."""
    rna = _to_rna(mature_dna)
    rc = _rc_rna(rna)
    ox = {int(p) for p in oxidized_positions}
    n = len(rna)
    for p in ox:
        i = p - 1
        if 0 <= i < n and rna[i] == "G":
            rc[n - 1 - i] = "A"
    return "".join(rc)


def _wt_selective_rc(mature_dna: str) -> str:
    return "".join(_rc_rna(_to_rna(mature_dna)))


def _duplex(oligo_rna: str, mature_rna: str) -> float | None:
    try:
        from o8g_thermo import dg_rnaduplex, vienna_available

        if not vienna_available():
            return None
        e = float(dg_rnaduplex(oligo_rna, mature_rna))
        return e if e == e else None  # NaN check
    except Exception:
        return None


def _mature_proxy(mature_dna: str, oxidized_positions: tuple[int, ...]) -> str:
    try:
        from o8g_thermo import mature_for_thermo

        return mature_for_thermo(mature_dna, oxidized_positions)
    except Exception:
        rna = list(_to_rna(mature_dna))
        for p in oxidized_positions:
            i = int(p) - 1
            if 0 <= i < len(rna) and rna[i] == "G":
                rna[i] = "U"
        return "".join(rna)


def _fold_from_ddg(ddG: float) -> float:
    """Boltzmann-ish fold favoring oxidized when ddG = dG_ox - dG_wt is negative."""
    # RT ≈ 0.616 kcal/mol at 37 °C; report as preference for oxidized form
    import math

    return float(math.exp(-ddG / 0.616))


def reverse_complement_dna(seq: str) -> str:
    s = clean_seq(seq)
    table = str.maketrans("ACGT", "TGCA")
    return s.translate(table)[::-1]


_COLLATERAL_COLS = ["symbol", "gene_id", "match_nt", "site_in_utr"]


def utr_parquet_path():
    """Longest-UTR table shipped with the app. None when the file is absent."""
    from pathlib import Path

    path = Path(__file__).resolve().parent / "utr3_human.parquet"
    return path if path.is_file() else None


def _normalize_utr(seq: str) -> str:
    return str(seq).upper().replace("U", "T")


def _longest_oligo_match(site: str, utr: str, min_match: int) -> tuple[int, str]:
    """Longest contiguous piece of ``site`` found in ``utr``, or (0, '')."""
    if len(site) < min_match or len(utr) < min_match:
        return 0, ""
    if not any(site[i : i + min_match] in utr for i in range(len(site) - min_match + 1)):
        return 0, ""
    for length in range(len(site), min_match - 1, -1):
        for start in range(len(site) - length + 1):
            window = site[start : start + length]
            if window in utr:
                return length, window
    return 0, ""


def _collateral_frame(hits: list[dict], top_n: int) -> "pd.DataFrame":
    import pandas as pd

    if not hits:
        return pd.DataFrame(columns=_COLLATERAL_COLS)
    out = pd.DataFrame(hits).sort_values(
        ["match_nt", "symbol"], ascending=[False, True]
    )
    return out.head(top_n).reset_index(drop=True)


def collateral_mrna_offtargets(
    oligo_dna: str,
    scanner=None,
    *,
    parquet_path=None,
    min_match: int = 10,
    top_n: int = 80,
) -> "pd.DataFrame":
    """Predict mRNA 3′UTR off-targets for an antagomir (hybridization, not AGO).

    The oligo binds RNA stretches complementary to itself. Those sites are the
    reverse-complement of the oligo (≈ mature for a full RC). We report genes
    whose 3′UTR contains a contiguous match of length ≥ ``min_match``.

    A live ``scanner`` is used when one is already built. Otherwise the search
    streams ``utr3_human.parquet`` one batch at a time and does not build the
    in-memory k-mer index.
    """
    oligo = clean_seq(oligo_dna)
    if len(oligo) < min_match:
        return pd.DataFrame(columns=_COLLATERAL_COLS)

    site = reverse_complement_dna(oligo)
    if scanner is not None and getattr(scanner, "utrs", None):
        hits: list[dict] = []
        symbols = scanner.symbols
        genes = scanner.genes
        for i, utr in enumerate(scanner.utrs):
            best, best_site = _longest_oligo_match(site, _normalize_utr(utr), min_match)
            if best >= min_match:
                hits.append(
                    {
                        "symbol": str(symbols[i]),
                        "gene_id": str(genes[i]),
                        "match_nt": int(best),
                        "site_in_utr": best_site,
                    }
                )
        return _collateral_frame(hits, top_n)

    path = parquet_path or utr_parquet_path()
    if path is None:
        return pd.DataFrame(columns=_COLLATERAL_COLS)
    return _collateral_from_parquet(site, path, min_match=min_match, top_n=top_n)


def _collateral_from_parquet(site: str, path, *, min_match: int, top_n: int) -> "pd.DataFrame":
    import pyarrow.parquet as pq

    hits: list[dict] = []
    parquet = pq.ParquetFile(path)
    for batch in parquet.iter_batches(
        batch_size=512, columns=["gene_id", "symbol", "utr3"]
    ):
        gene_ids = batch.column("gene_id").to_pylist()
        symbols = batch.column("symbol").to_pylist()
        utrs = batch.column("utr3").to_pylist()
        for gene_id, symbol, utr in zip(gene_ids, symbols, utrs):
            if not utr:
                continue
            best, best_site = _longest_oligo_match(site, _normalize_utr(utr), min_match)
            if best >= min_match:
                hits.append(
                    {
                        "symbol": str(symbol),
                        "gene_id": str(gene_id),
                        "match_nt": int(best),
                        "site_in_utr": best_site,
                    }
                )
    return _collateral_frame(hits, top_n)


def design_antagomir(
    mature_dna: str,
    oxidized_positions: tuple[int, ...],
    *,
    kind: str = "oxo-selective",
) -> OligoDesign:
    """Score an oligo for state discrimination (not AGO seed targeting)."""
    ox = tuple(sorted(int(p) for p in oxidized_positions))
    if kind == "wt-selective":
        oligo = _wt_selective_rc(mature_dna)
        kind_l = "wt-selective"
        note = (
            "WT-selective RC (C opposite every G). Prefers unmodified mature; "
            "oxidized Gs become mismatches under the G→U proxy."
        )
    else:
        oligo = _oxo_selective_rc(mature_dna, ox)
        kind_l = "oxo-selective"
        note = (
            "Oxo-selective RC: A opposite oxidized G (o8G:A), C opposite other Gs. "
            "dG_oxo uses G→U proxy (no native o8G parameters in ViennaRNA)."
        )

    wt = _to_rna(mature_dna)
    ox_m = _mature_proxy(mature_dna, ox)
    dG_n = _duplex(oligo, wt)
    dG_o = _duplex(oligo, ox_m)

    # Fallback compositional estimate when ViennaRNA missing
    if dG_n is None or dG_o is None:
        def _gc_score(a: str, b: str) -> float:
            s = 0.0
            for x, y in zip(a, reversed(b) if len(b) == len(a) else b):
                # crude: reward WC pairs
                pair = {x, y}
                if pair in ({"G", "C"}, {"C", "G"}):
                    s -= 1.0
                elif pair in ({"A", "U"}, {"U", "A"}, {"A", "T"}, {"T", "A"}):
                    s -= 0.5
            return s

        # oligo is same length as mature; align 5'–3' both
        dG_n = _gc_score(oligo, wt)
        dG_o = _gc_score(oligo, ox_m)
        note += " ViennaRNA unavailable — compositional placeholder ΔG only."

    ddG = (dG_o - dG_n) if dG_n is not None and dG_o is not None else None
    fold = _fold_from_ddg(ddG) if ddG is not None else None

    if fold is not None and fold < 10 and kind_l == "oxo-selective":
        n_ox = max(len(ox), 1)
        note += (
            f" Weak discrimination expected: only {n_ox} of {len(wt)} nt encode the "
            f"oxidation difference (|ΔΔG|={abs(ddG):.2f} kcal/mol). "
            "Interior single-G states usually beat position-2-only; short seed-centered "
            "chemistries (not scored here) are typically needed for large fold preferences."
        )

    return OligoDesign(
        sequence=oligo.replace("U", "T"),
        kind=kind_l,
        target_state=",".join(str(p) for p in ox) or "none",
        dG_normal=dG_n,
        dG_oxo=dG_o,
        ddG=ddG,
        fold_preference=fold,
        note=note,
    )


def rank_single_g_designs(mature_dna: str) -> list[OligoDesign]:
    """Rank oxo-selective designs for every single seed-G oxidation."""
    seed = extract_seed(mature_dna)
    out = []
    for p in g_positions(seed):
        out.append(design_antagomir(mature_dna, (p,), kind="oxo-selective"))
    out.sort(key=lambda d: (d.fold_preference or 0.0), reverse=True)
    return out


def position2_warning(oxidized_positions: tuple[int, ...]) -> bool:
    """True when the only oxidized seed position is position 2 (poor antagomir)."""
    return tuple(oxidized_positions) == (2,)


def feasibility_label(fold: float | None) -> tuple[str, str]:
    """Return (level, message) level in {good, modest, weak, na}."""
    if fold is None:
        return "na", "No ΔG available."
    if fold >= 50:
        return "good", f"Promising discrimination ({fold:.1f}× for oxidized form)."
    if fold >= 10:
        return "modest", f"Modest discrimination ({fold:.1f}×) — directional only."
    return (
        "weak",
        f"Weak discrimination ({fold:.1f}×). Full-length oligos only differ at oxidized "
        "position(s); small |ΔΔG| is expected, not a ViennaRNA bug. Prefer interior Gs "
        "and treat numbers as proxy ranking.",
    )
