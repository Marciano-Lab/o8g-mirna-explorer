"""
o8g_enrich.py
=============
Self-contained over-representation (hypergeometric) pathway enrichment on
predicted miRNA target gene lists, using locally cached Enrichr .gmt libraries.
No network access is required at query time.

Enrichment model
-----------------
For a query gene set Q against a pathway gene set P, within a background of N
genes, the number of overlaps k follows a hypergeometric distribution.  The
enrichment p-value is the upper-tail P(X >= k).  We report:
    overlap, query_size, term_size, expected, odds_ratio (fold enrichment),
    p_value, and Benjamini-Hochberg adjusted q_value per library.

Background
----------
`background` may be:
  - None  → library-union population (legacy; inflated for this app — prefer
            the 3'UTR symbol universe from o8g_targets.db)
  - int   → fixed population size N (legacy; term sizes unrestricted)
  - iterable of symbols → restrict BOTH the population and the query to
            set(background) ∩ (union of library gene sets). Term sizes are
            likewise intersected with that pool. This is the correct path for
            UTR-indexed target lists.

Differential enrichments must also report a within-pool control
(`enrich_within_pool`) whose background is union(unmodified, oxidized) —
genome-wide ORA alone answers "what does this miRNA target", not "what did
oxidation change".
"""
from __future__ import annotations
import glob, os, re
from collections.abc import Iterable
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

GENESET_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "genesets")

LIBRARY_FILES = {
    "GO_BP": "GO_Biological_Process_2023.gmt",
    "GO_MF": "GO_Molecular_Function_2023.gmt",
    "GO_CC": "GO_Cellular_Component_2023.gmt",
    "KEGG": "KEGG_2021_Human.gmt",
    "Reactome": "Reactome_2022.gmt",
    "Hallmark": "MSigDB_Hallmark_2020.gmt",
    # TF motif / ChIP libraries (optional — views degrade if files absent)
    # Fetch: python scripts/fetch_tf_genesets.py
    "TRANSFAC_JASPAR": "TRANSFAC_and_JASPAR_PWMs.gmt",
    "JASPAR_2025": "JASPAR_PWM_Human_2025.gmt",
    "GenomeBrowser_PWM": "Genome_Browser_PWMs.gmt",  # UCSC/ENCODE ≈ Ensembl-adjacent PWMs
    "ENCODE_ChEA_Consensus": "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt",
    "ChEA_TF": "ChEA_2022.gmt",
    "ENCODE_TF": "ENCODE_TF_ChIP-seq_2015.gmt",
    "TRRUST_TF": "TRRUST_Transcription_Factors_2019.gmt",
    "TF_Perturb": "TF_Perturbations_Followed_by_Expression.gmt",
}

# Motif-first defaults, then ChIP / curated TF–target networks
TF_LIBRARIES = (
    "TRANSFAC_JASPAR",
    "JASPAR_2025",
    "GenomeBrowser_PWM",
    "ENCODE_ChEA_Consensus",
    "ChEA_TF",
    "ENCODE_TF",
    "TRRUST_TF",
    "TF_Perturb",
)

TF_MOTIF_LIBRARIES = ("TRANSFAC_JASPAR", "JASPAR_2025", "GenomeBrowser_PWM")


@lru_cache(maxsize=None)
def load_library(name: str) -> tuple:
    """Return (terms tuple, list-of-frozenset genes) for a cached library."""
    path = os.path.join(GENESET_DIR, LIBRARY_FILES[name])
    terms, sets = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = frozenset(g.split(",")[0].strip().upper() for g in parts[2:] if g.strip())
            if genes:
                terms.append(term); sets.append(genes)
    return tuple(terms), sets


def available_libraries() -> list[str]:
    # Primary pathway libs first; TF libs available when present
    primary = ["GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark"]
    out = [k for k in primary if k in LIBRARY_FILES
           and os.path.exists(os.path.join(GENESET_DIR, LIBRARY_FILES[k]))]
    for k in TF_LIBRARIES:
        if k in LIBRARY_FILES and os.path.exists(os.path.join(GENESET_DIR, LIBRARY_FILES[k])):
            out.append(k)
    return out


def available_tf_libraries() -> list[str]:
    return [
        k
        for k in TF_LIBRARIES
        if k in LIBRARY_FILES
        and os.path.exists(os.path.join(GENESET_DIR, LIBRARY_FILES[k]))
    ]


def _library_universe(sets) -> set[str]:
    allg: set[str] = set()
    for s in sets:
        allg |= s
    return allg


def enrich(query_genes, library: str = "GO_BP", background: int | Iterable | None = None,
           min_overlap: int = 2, top: int | None = None, *,
           include_genes: bool = True) -> pd.DataFrame:
    """Hypergeometric over-representation of `query_genes` in one library."""
    terms, sets = load_library(library)
    lib_u = _library_universe(sets)

    if background is None:
        pool = lib_u
        N = len(pool)
        restrict_terms = False
    elif isinstance(background, (int, np.integer)):
        pool = lib_u  # query still restricted to lib genes; N forced
        N = int(background)
        restrict_terms = False
    else:
        bg = {str(g).upper() for g in background}
        pool = bg & lib_u
        N = len(pool)
        restrict_terms = True

    Q = frozenset(g.upper() for g in query_genes) & pool
    q = len(Q)
    rows = []
    for term, P0 in zip(terms, sets):
        P = (P0 & pool) if restrict_terms else P0
        k = len(Q & P)
        if k < min_overlap:
            continue
        M = len(P)
        if M == 0 or N == 0 or q == 0:
            continue
        pval = hypergeom.sf(k - 1, N, M, q)
        expected = q * M / N if N else np.nan
        odds = (k / expected) if expected else np.nan
        genes = ";".join(sorted(Q & P)) if include_genes else ""
        rows.append((term, k, q, M, expected, odds, pval, N, genes))
    df = pd.DataFrame(
        rows,
        columns=[
            "term",
            "overlap",
            "query_size",
            "term_size",
            "expected",
            "odds_ratio",
            "p_value",
            "background_size",
            "genes",
        ],
    )
    if len(df):
        df = df.sort_values("p_value").reset_index(drop=True)
        m = len(df)
        ranks = np.arange(1, m + 1)
        q_raw = df["p_value"].to_numpy() * m / ranks
        df["q_value"] = np.minimum.accumulate(q_raw[::-1])[::-1].clip(0, 1)
        df["neglog10_q"] = -np.log10(df["q_value"].clip(lower=1e-300))
        df["library"] = library
    else:
        for c in ["q_value", "neglog10_q"]:
            df[c] = pd.Series(dtype=float)
        df["library"] = library
    if top:
        df = df.head(top)
    return df


def enrich_within_pool(query_genes, pool_genes, library: str = "GO_BP",
                       min_overlap: int = 2, top: int | None = None) -> pd.DataFrame:
    """ORA against an explicit gene pool (e.g. union of unmod + oxidized targets).

    This is the correct control for differential (lost/gained) enrichment:
    genome-wide ORA alone answers 'what does this miRNA target', not
    'what did oxidation change'.
    """
    return enrich(query_genes, library=library, background=pool_genes,
                  min_overlap=min_overlap, top=top)


def compare_states(genes_A, genes_B, library: str = "GO_BP",
                   label_A: str = "A", label_B: str = "B",
                   background: int | Iterable | None = None, min_overlap: int = 2) -> pd.DataFrame:
    """Merge enrichment of two gene lists for differential (volcano/heatmap) views.

    Volcano axes:
        x = log2( odds_ratio_B / odds_ratio_A )   (enrichment shift B vs A)
        y = -log10( min(q_A, q_B) )               (best significance of the term)
    """
    a = enrich(genes_A, library, background, min_overlap).set_index("term")
    b = enrich(genes_B, library, background, min_overlap).set_index("term")
    terms = a.index.union(b.index)
    def col(df, c):
        return df[c].reindex(terms)
    out = pd.DataFrame({
        "term": terms,
        f"overlap_{label_A}": col(a, "overlap").fillna(0).astype(int),
        f"overlap_{label_B}": col(b, "overlap").fillna(0).astype(int),
        f"odds_{label_A}": col(a, "odds_ratio"),
        f"odds_{label_B}": col(b, "odds_ratio"),
        f"q_{label_A}": col(a, "q_value").fillna(1.0),
        f"q_{label_B}": col(b, "q_value").fillna(1.0),
    })
    oa = out[f"odds_{label_A}"].fillna(0.0).to_numpy()
    ob = out[f"odds_{label_B}"].fillna(0.0).to_numpy()
    eps = 1e-3
    out["log2_or_ratio"] = np.log2((ob + eps) / (oa + eps))
    out["neglog10_q_best"] = -np.log10(np.minimum(out[f"q_{label_A}"], out[f"q_{label_B}"]).clip(lower=1e-300))
    out["library"] = library
    out["direction"] = np.where(out["log2_or_ratio"] > 0, f"up in {label_B}", f"up in {label_A}")
    return out.sort_values("neglog10_q_best", ascending=False).reset_index(drop=True)


def differential_pathways(genes_A, genes_B, library: str = "GO_BP",
                          label_A: str = "A", label_B: str = "B",
                          min_overlap: int = 3) -> pd.DataFrame:
    """Fisher's-exact differential pathway test between two target lists.

    For each pathway P, build the 2x2 table over the genes that are targets in
    exactly one state (the symmetric difference):

                       in P        not in P
        only-in-A       a             b
        only-in-B       c             d

    Fisher's exact tests whether P is over-represented among the A-specific
    (lost) vs B-specific (gained) targets.  This directly tests the *shift*
    of the target repertoire, so it is the natural volcano statistic:

        x = log2 odds-ratio   (>0 -> pathway skews toward A / lost)
        y = -log10 q          (BH-adjusted)

    Returns tidy DataFrame sorted by significance.
    """
    from scipy.stats import fisher_exact
    A = frozenset(g.upper() for g in genes_A)
    B = frozenset(g.upper() for g in genes_B)
    onlyA = A - B
    onlyB = B - A
    terms, sets = load_library(library)
    nA, nB = len(onlyA), len(onlyB)
    rows = []
    for term, P in zip(terms, sets):
        a = len(onlyA & P); c = len(onlyB & P)
        if a + c < min_overlap:
            continue
        b = nA - a; d = nB - c
        orr, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append((term, a, c, orr, p))
    df = pd.DataFrame(rows, columns=["term", f"n_{label_A}", f"n_{label_B}",
                                     "odds_ratio", "p_value"])
    if len(df):
        df = df.sort_values("p_value").reset_index(drop=True)
        m = len(df); ranks = np.arange(1, m + 1)
        q_raw = df["p_value"].to_numpy() * m / ranks
        df["q_value"] = np.minimum.accumulate(q_raw[::-1])[::-1].clip(0, 1)
        df["neglog10_q"] = -np.log10(df["q_value"].clip(lower=1e-300))
        # log2 OR, guarding zeros
        orr = df["odds_ratio"].to_numpy()
        orr = np.where(np.isfinite(orr) & (orr > 0), orr, np.nan)
        # cap infinities/zeros for display
        with np.errstate(divide="ignore"):
            l2 = np.log2(orr)
        # replace inf odds (c==0) with a large finite value, 0 odds (a==0) with negative
        l2 = np.where(np.isinf(df["odds_ratio"]) | (df["odds_ratio"] > 1e6), np.log2(orr[np.isfinite(orr)].max() if np.isfinite(orr).any() else 8)+2, l2)
        df["log2_or"] = pd.Series(l2).fillna(0.0)
        df["direction"] = np.where(df["odds_ratio"] >= 1, f"toward {label_A}", f"toward {label_B}")
        df["library"] = library
    return df


if __name__ == "__main__":
    libs = available_libraries()
    print("libraries:", libs)
    demo = ["HDAC4", "TWF1", "GJA1", "KCNJ2", "BDNF", "CDK6", "CCND2", "MET",
            "PTEN", "VEGFA", "IGF1", "FOXP1", "SRF", "MYOCD", "CACNA1C"]
    e = enrich(demo, "GO_BP", top=5)
    print(e[["term", "overlap", "term_size", "odds_ratio", "q_value"]].to_string(index=False))
