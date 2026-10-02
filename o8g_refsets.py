"""
o8g_refsets.py
==============
External miRNA–target reference sets for in-app multi-DB comparison.

Sources (downloadable / open API — cite thresholds before changing)
-------------------------------------------------------------------
Predicted:
  TargetScanHuman 8.0  — local Predicted_Targets_Info (8mer/7mer-m8)
  miRDB v6.0           — local gzip; score ≥ 80 (Chen & Wang 2020)
  DIANA-microT-CDS     — local flat file; interaction_score ≥ 0.7
  miRmap 202203        — local zst (bench cache) or on-demand filter; pct ≥ 80

Experimentally supported:
  ENCORI / starBase    — open REST API (CLIP ≥ 1 experiment)
  miRTarBase           — local hsa_MTI if present under paper/data/mirtarbase/;
                         else ENCORI is the live experimental fallback
                         (CUHK bulk download often 404)

Paths default to paper/data/ (local lab files; gitignored). The hosted app ships
paper/data/hosted_refsets.sqlite instead of the multi-gigabyte downloads.
Per-miRNA results are cached under paper/data/cache/refsets/.
"""
from __future__ import annotations

import gzip
import io
import json
import sqlite3
import urllib.parse
import urllib.request
from functools import lru_cache
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "paper" / "data"
CACHE = DATA / "cache" / "refsets"
HOSTED = DATA / "hosted_refsets.sqlite"
HOSTED_GZ = DATA / "hosted_refsets.sqlite.gz"


def _ensure_hosted_sqlite() -> None:
    """Expand the gzipped extract when a checkout has only the shipped archive."""
    if HOSTED.exists() or not HOSTED_GZ.exists():
        return
    HOSTED.parent.mkdir(parents=True, exist_ok=True)
    tmp = HOSTED.with_suffix(".sqlite.partial")
    with gzip.open(HOSTED_GZ, "rb") as src, open(tmp, "wb") as out:
        while True:
            chunk = src.read(8 * 1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(HOSTED)

MIRDB_SCORE_MIN = 80.0
DIANA_SCORE_MIN = 0.7
MIRMAP_PERCENTILE_MIN = 80.0
ENCORI_CLIP_MIN = 1

# Known TargetScan families for common miRNAs (extend via resolve_targetscan_family)
TS_FAMILIES = {
    "hsa-miR-1-3p": "miR-1-3p/206",
    "hsa-miR-206": "miR-1-3p/206",
    "hsa-miR-124-3p": "miR-124-3p.1",
    "hsa-let-7a-5p": "let-7-5p/98-5p",
    "hsa-let-7b-5p": "let-7-5p/98-5p",
    "hsa-let-7c-5p": "let-7-5p/98-5p",
    "hsa-let-7d-5p": "let-7-5p/98-5p",
    "hsa-let-7e-5p": "let-7-5p/98-5p",
    "hsa-let-7f-5p": "let-7-5p/98-5p",
    "hsa-let-7g-5p": "let-7-5p/98-5p",
    "hsa-let-7i-5p": "let-7-5p/98-5p",
    "hsa-miR-98-5p": "let-7-5p/98-5p",
    "hsa-miR-122-5p": "miR-122-5p",
    "hsa-miR-21-5p": "miR-21-5p",
    "hsa-miR-155-5p": "miR-155-5p",
    "hsa-miR-17-5p": "miR-17-5p/20-5p/93-5p/106-5p",
    "hsa-miR-34a-5p": "miR-34-5p/449-5p",
}

# Alternate family strings that appear in Predicted_Targets_Info across species rows
_TS_FAMILY_ALIASES: dict[str, set[str]] = {
    "let-7-5p/98-5p": {
        "let-7-5p/98-5p",
        "let-7/98",
        "let-7-5p/miR-98",
        "let-7-5p/miR-98-5p",
        "let-7-5p/miR-98/6134",
    },
    "miR-1-3p/206": {"miR-1-3p/206", "miR-1-3p/206/6132", "miR-1/206"},
}


def resolve_targetscan_family(mirna: str) -> str | None:
    """Map a mature miRNA name to its TargetScan family string when known."""
    if mirna in TS_FAMILIES:
        return TS_FAMILIES[mirna]
    stem = mirna.replace("hsa-", "")
    # All let-7*-5p + miR-98-5p share one conserved/predicted family in TS8
    if stem.lower() == "mir-98-5p":
        return "let-7-5p/98-5p"
    if stem.startswith("let-7") and stem.endswith("-5p"):
        return "let-7-5p/98-5p"
    if stem in ("miR-1-3p", "miR-206"):
        return "miR-1-3p/206"
    return None


def _family_column_matches(fam_col: str, wanted: str | None, mirna: str) -> bool:
    if not fam_col:
        return False
    if wanted:
        aliases = _TS_FAMILY_ALIASES.get(wanted, {wanted})
        if fam_col in aliases or fam_col == wanted:
            return True
        # substring only when wanted is a clean token (avoid "7" matching everything)
        if wanted in fam_col:
            return True
        # Known family: do not also accept a looser stem hit. "miR-124-3p" is a
        # substring of both miR-124-3p.1 and miR-124-3p.2, and only .1 is this miRNA.
        return False
    stem = mirna.replace("hsa-", "")
    return bool(stem and stem in fam_col)

MIRNA_TO_MIMAT = {
    "hsa-let-7a-5p": "MIMAT0000062",
    "hsa-miR-1-3p": "MIMAT0000416",
    "hsa-miR-122-5p": "MIMAT0000421",
    "hsa-miR-124-3p": "MIMAT0000422",
}

TOOL_META = {
    "Explorer": {
        "kind": "predicted",
        "version": "o8G seed engine (this app)",
        "threshold": "prediction mode (sequence-based low/high stringency / Consensus)",
        "citation": "local UTR seed scan",
    },
    "TargetScan": {
        "kind": "predicted",
        "version": "TargetScanHuman 8.0",
        "threshold": "site_type in {8mer, 7mer-m8}",
        "citation": "Agarwal eLife 2015; McGeary Science 2019",
    },
    "miRDB": {
        "kind": "predicted",
        "version": "miRDB v6.0",
        "threshold": f"score >= {MIRDB_SCORE_MIN:g}",
        "citation": "Chen & Wang NAR 2020",
    },
    "DIANA-microT": {
        "kind": "predicted",
        "version": "microT-CDS (miRBase 22.1 flat)",
        "threshold": f"interaction_score >= {DIANA_SCORE_MIN:g}",
        "citation": "Paraskevopoulou NAR 2013",
    },
    "miRmap": {
        "kind": "predicted",
        "version": "miRmap 1.2.0 / mirmap_202203",
        "threshold": f"within-miRNA percentile >= {MIRMAP_PERCENTILE_MIN:g}",
        "citation": "Vejnar NAR 2013",
    },
    "ENCORI": {
        "kind": "experimental",
        "version": "ENCORI/starBase API hg38",
        "threshold": f"clipExpNum >= {ENCORI_CLIP_MIN}",
        "citation": "Li et al. NAR 2014; ENCORI web API",
    },
    "miRTarBase": {
        "kind": "experimental",
        "version": "miRTarBase local hsa_MTI (strong evidence)",
        "threshold": "Luciferase / Western / qPCR (Functional MTI)",
        "citation": "Huang et al. NAR 2020/2025",
    },
}


def _cache_path(tool: str, mirna: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    safe = mirna.replace("/", "_")
    return CACHE / f"{tool}_{safe}.json"


def _load_cache(tool: str, mirna: str) -> set[str] | None:
    p = _cache_path(tool, mirna)
    if not p.exists():
        return None
    try:
        return set(json.loads(p.read_text()))
    except Exception:
        return None


def _save_cache(tool: str, mirna: str, genes: set[str]) -> None:
    _cache_path(tool, mirna).write_text(json.dumps(sorted(genes)))


def ensembl_to_symbol() -> dict[str, str]:
    utr = ROOT / "utr3_human.parquet"
    if not utr.exists():
        return {}
    df = pd.read_parquet(utr, columns=["gene_id", "symbol"])
    return dict(zip(df["gene_id"].astype(str), df["symbol"].astype(str)))


def _hosted_connect() -> sqlite3.Connection:
    from o8g_db import connect_sqlite

    return connect_sqlite(HOSTED, check_same_thread=True)


@lru_cache(maxsize=1)
def _hosted_tools() -> frozenset[str]:
    """Tools present in the shippable extract. Empty if that file is absent."""
    _ensure_hosted_sqlite()
    if not HOSTED.exists():
        return frozenset()
    con = _hosted_connect()
    try:
        rows = con.execute("SELECT tool FROM tools").fetchall()
    except sqlite3.OperationalError:
        rows = []
    finally:
        con.close()
    return frozenset(r[0] for r in rows)


@lru_cache(maxsize=256)
def _hosted_symbols(tool: str, mirna: str) -> set[str] | None:
    """Gene symbols from the hosted extract, or None if that tool was not shipped."""
    if tool not in _hosted_tools():
        return None
    con = _hosted_connect()
    try:
        rows = con.execute(
            "SELECT symbol FROM targets WHERE tool=? AND mirna=?",
            [tool, mirna],
        ).fetchall()
    finally:
        con.close()
    return frozenset(r[0] for r in rows)


def available_tools() -> dict[str, bool]:
    """Which external tools have data/API available right now."""
    hosted = _hosted_tools()
    return {
        "TargetScan": "TargetScan" in hosted
        or (DATA / "Predicted_Targets_Info.default_predictions.txt").exists(),
        "miRDB": "miRDB" in hosted
        or (
            (DATA / "miRDB_v6.0_prediction_result.txt.gz").exists()
            and (DATA / "refseq_to_symbol.tsv").exists()
        ),
        "DIANA-microT": "DIANA-microT" in hosted
        or (DATA / "interactions_human.microT.mirbase.txt.gz").exists()
        or (DATA / "diana_microt_4mirs.tsv").exists(),
        "miRmap": "miRmap" in hosted
        or (DATA / "mirmap_4mirs.parquet").exists()
        or (DATA / "mirmap_202203_homsap_targets_1to1.csv.zst").exists(),
        "ENCORI": True,  # live API
        "miRTarBase": "miRTarBase" in hosted
        or (
            any(
                p.stat().st_size > 1000
                for p in (DATA / "mirtarbase").glob("hsa_MTI*")
            )
            if (DATA / "mirtarbase").exists()
            else False
        ),
    }


def load_targetscan(mirna: str) -> set[str]:
    cached = _load_cache("TargetScan", mirna)
    # Never trust an empty cache entry (legacy poison from failed name match)
    if cached:
        return cached
    if cached is not None and len(cached) == 0:
        try:
            _cache_path("TargetScan", mirna).unlink(missing_ok=True)
        except Exception:
            pass

    hosted = _hosted_symbols("TargetScan", mirna)
    if hosted is not None:
        if hosted:
            _save_cache("TargetScan", mirna, hosted)
        return hosted

    path = DATA / "Predicted_Targets_Info.default_predictions.txt"
    if not path.exists():
        return set()
    fam = resolve_targetscan_family(mirna)
    out: set[str] = set()
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or p[4] != "9606":
                continue
            if p[9] not in ("8mer", "7mer-m8"):
                continue
            if _family_column_matches(p[0], fam, mirna):
                out.add(p[2])
    # Only cache non-empty hits (empty = unknown / missing, retry next time)
    if out:
        _save_cache("TargetScan", mirna, out)
    return out


def load_mirdb(mirna: str, score_min: float = MIRDB_SCORE_MIN) -> set[str]:
    cached = _load_cache("miRDB", mirna)
    if cached is not None:
        return cached
    hosted = _hosted_symbols("miRDB", mirna)
    if hosted is not None:
        _save_cache("miRDB", mirna, hosted)
        return hosted
    map_path = DATA / "refseq_to_symbol.tsv"
    pred = DATA / "miRDB_v6.0_prediction_result.txt.gz"
    if not map_path.exists() or not pred.exists():
        return set()
    mp = dict(
        pd.read_csv(map_path, sep="\t")[["refseq", "symbol"]].itertuples(index=False)
    )
    out: set[str] = set()
    with gzip.open(pred, "rt") as f:
        for line in f:
            mir, ref, score = line.rstrip().split("\t")
            if mir == mirna and float(score) >= score_min:
                sym = mp.get(ref.split(".")[0])
                if sym:
                    out.add(sym)
    _save_cache("miRDB", mirna, out)
    return out


def load_diana(mirna: str, score_min: float = DIANA_SCORE_MIN) -> set[str]:
    cached = _load_cache("DIANA-microT", mirna)
    if cached is not None:
        return cached
    hosted = _hosted_symbols("DIANA-microT", mirna)
    if hosted is not None:
        _save_cache("DIANA-microT", mirna, hosted)
        return hosted
    cache4 = DATA / "diana_microt_4mirs.tsv"
    raw = DATA / "interactions_human.microT.mirbase.txt.gz"
    e2s = ensembl_to_symbol()
    out: set[str] = set()

    def _consume(df: pd.DataFrame) -> None:
        cols = {c.lower(): c for c in df.columns}
        mir_c = cols.get("mirna") or df.columns[0]
        gene_c = cols.get("ensembl_gene_id") or df.columns[1]
        score_c = cols.get("interaction_score") or df.columns[2]
        scores = pd.to_numeric(df[score_c], errors="coerce")
        sub = df[(df[mir_c].astype(str) == mirna) & (scores >= score_min)]
        for gid in sub[gene_c]:
            sym = e2s.get(str(gid).split(".")[0]) or e2s.get(str(gid))
            if sym:
                out.add(sym)

    if cache4.exists():
        _consume(pd.read_csv(cache4, sep="\t"))
    if not out and raw.exists():
        rows = []
        with gzip.open(raw, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p and p[0] == mirna:
                    rows.append(p)
        if rows:
            _consume(pd.DataFrame(rows, columns=header))
    _save_cache("DIANA-microT", mirna, out)
    return out


def load_mirmap(mirna: str, pct_min: float = MIRMAP_PERCENTILE_MIN) -> set[str]:
    cached = _load_cache("miRmap", mirna)
    if cached is not None:
        return cached
    hosted = _hosted_symbols("miRmap", mirna)
    if hosted is not None:
        _save_cache("miRmap", mirna, hosted)
        return hosted
    path = DATA / "mirmap_4mirs.parquet"
    out: set[str] = set()
    if path.exists():
        df = pd.read_parquet(path)
        # expected columns from benchmark cache
        mir_c = "mirna" if "mirna" in df.columns else None
        if mir_c and "symbol" in df.columns:
            sub = df[df[mir_c] == mirna]
            if "percentile" in sub.columns:
                sub = sub[sub["percentile"] >= pct_min]
            out = set(sub["symbol"].astype(str))
    _save_cache("miRmap", mirna, out)
    return out


def load_encori(mirna: str, clip_min: int = ENCORI_CLIP_MIN) -> set[str]:
    """ENCORI/starBase open API — CLIP-supported targets."""
    cached = _load_cache("ENCORI", mirna)
    if cached is not None:
        return cached
    qs = urllib.parse.urlencode(
        {
            "assembly": "hg38",
            "geneType": "mRNA",
            "miRNA": mirna,
            "clipExpNum": str(clip_min),
            "degraExpNum": "0",
            "pancancerNum": "0",
            "programNum": "1",
            "program": "None",
            "target": "all",
            "cellType": "all",
        }
    )
    url = f"https://rnasysu.com/encori/api/miRNATarget/?{qs}"
    out: set[str] = set()
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "o8g-mirna-explorer/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        # skip comment lines; find header
        lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
        if len(lines) < 2:
            _save_cache("ENCORI", mirna, out)
            return out
        header = lines[0].split("\t")
        # geneName or geneName often present
        try:
            gi = header.index("geneName")
        except ValueError:
            gi = 1 if len(header) > 1 else 0
        for ln in lines[1:]:
            p = ln.split("\t")
            if len(p) > gi and p[gi]:
                out.add(p[gi])
    except Exception:
        pass
    _save_cache("ENCORI", mirna, out)
    return out


def load_mirtarbase(mirna: str) -> set[str]:
    cached = _load_cache("miRTarBase", mirna)
    if cached is not None:
        return cached
    hosted = _hosted_symbols("miRTarBase", mirna)
    if hosted is not None:
        _save_cache("miRTarBase", mirna, hosted)
        return hosted
    out: set[str] = set()
    mdir = DATA / "mirtarbase"
    files = list(mdir.glob("hsa_MTI*")) if mdir.exists() else []
    # Prefer TSV/CSV (fast) over xlsx; skip tiny 404 stubs
    candidates = [p for p in files if p.stat().st_size > 1000]
    candidates.sort(
        key=lambda p: (
            0 if p.suffix.lower() in {".tsv", ".txt", ".csv"} else 1,
            -p.stat().st_size,
        )
    )
    path = candidates[0] if candidates else None
    if path is None:
        return out
    if path.suffix.lower() in {".xlsx", ".xls"}:
        raw = pd.read_excel(path)
    else:
        raw = pd.read_csv(path, sep=None, engine="python")
    colmap = {c.lower().replace(" ", "_"): c for c in raw.columns}

    def col(*names):
        for n in names:
            if n in colmap:
                return colmap[n]
            if n in raw.columns:
                return n
        return None

    mir_c = col("mirna", "mirna")
    gene_c = col("target_gene", "gene", "target_symbol")
    exp_c = col("experiments", "experiment")
    if mir_c is None or gene_c is None:
        return out
    strong = ("luciferase", "reporter", "western", "qpcr", "qrt-pcr", "immunoblot")
    sub = raw[raw[mir_c].astype(str) == mirna]
    for _, r in sub.iterrows():
        exp = str(r[exp_c]).lower() if exp_c else ""
        if any(s in exp for s in strong) or not exp_c:
            out.add(str(r[gene_c]).upper() if str(r[gene_c]).islower() else str(r[gene_c]))
            # keep original case from file
            out.add(str(r[gene_c]))
    _save_cache("miRTarBase", mirna, out)
    return out


LOADERS = {
    "TargetScan": load_targetscan,
    "miRDB": load_mirdb,
    "DIANA-microT": load_diana,
    "miRmap": load_mirmap,
    "ENCORI": load_encori,
    "miRTarBase": load_mirtarbase,
}


def load_selected(mirna: str, tools: list[str]) -> dict[str, set[str]]:
    """Load gene-symbol sets for each requested tool."""
    avail = available_tools()
    out: dict[str, set[str]] = {}
    for t in tools:
        if t == "Explorer":
            continue
        if t not in LOADERS:
            continue
        if t != "ENCORI" and not avail.get(t, False):
            out[t] = set()
            continue
        out[t] = LOADERS[t](mirna)
    return out


def membership_matrix(sets: dict[str, set[str]]) -> pd.DataFrame:
    """Binary gene × database membership table."""
    all_genes = sorted(set().union(*sets.values()) if sets else [])
    rows = []
    for g in all_genes:
        row = {"symbol": g}
        for name, s in sets.items():
            row[name] = g in s
        row["n_dbs"] = sum(1 for name in sets if row[name])
        rows.append(row)
    return pd.DataFrame(rows)


def consensus_intersection(sets: dict[str, set[str]], min_dbs: int | None = None) -> set[str]:
    """Genes present in every set (or in at least min_dbs sets)."""
    if not sets:
        return set()
    if min_dbs is None or min_dbs >= len(sets):
        it = iter(sets.values())
        core = set(next(it))
        for s in it:
            core &= s
        return core
    mat = membership_matrix(sets)
    return set(mat.loc[mat["n_dbs"] >= min_dbs, "symbol"])


def external_support_for_genes(
    genes: set[str],
    external: dict[str, set[str]],
) -> pd.DataFrame:
    """Binary membership of *genes* in each external DB (+ ``n_external`` count).

    Used for loss-of-binding prioritization: pass Explorer-lost symbols and the
    unmodified catalogs (TargetScan / miRDB / …). Genes with high ``n_external``
    are known WT targets that Explorer predicts vanish after o8G.
    """
    if not genes:
        return pd.DataFrame(columns=["symbol", "n_external", *external.keys()])
    rows = []
    for g in sorted(genes):
        row: dict = {"symbol": g}
        for name, s in external.items():
            row[name] = g in s
        row["n_external"] = int(sum(1 for name in external if row[name]))
        rows.append(row)
    return pd.DataFrame(rows)


def lost_corroborated(
    lost: set[str],
    external: dict[str, set[str]],
    *,
    min_external: int = 1,
) -> set[str]:
    """Lost genes supported by at least ``min_external`` external catalogs."""
    if not lost or not external:
        return set() if min_external > 0 else set(lost)
    mat = external_support_for_genes(lost, external)
    return set(mat.loc[mat["n_external"] >= min_external, "symbol"])
