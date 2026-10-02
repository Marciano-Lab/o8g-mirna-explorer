#!/usr/bin/env python3
"""Build paper/data/hosted_refsets.sqlite for the Render deploy.

The raw downloads (DIANA ~307 MB, miRmap ~1.8 GB, plus TargetScan / miRDB /
miRTarBase) are gitignored and are not in the hosted image. This extract keeps
the same gene-symbol sets the in-app comparison uses, one row per
(tool, miRNA, symbol), small enough to ship.

Thresholds match o8g_refsets.py:
  TargetScan   site_type in {8mer, 7mer-m8}, human
  miRDB        score >= 80
  DIANA-microT interaction_score >= 0.7, Ensembl gene → symbol via utr3
  miRmap       best (most negative) mirmap_score per symbol on the indexed
               3'UTR transcript, within-miRNA percentile >= 80
  miRTarBase   strong-evidence experiments (luciferase / western / qPCR / …)
"""
from __future__ import annotations

import gzip
import io
import sqlite3
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import zstandard

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "data"
OUT = DATA / "hosted_refsets.sqlite"

import sys

sys.path.insert(0, str(ROOT))
from o8g_refsets import (  # noqa: E402
    DIANA_SCORE_MIN,
    MIRDB_SCORE_MIN,
    MIRMAP_PERCENTILE_MIN,
    _family_column_matches,
    resolve_targetscan_family,
)


def _connect() -> sqlite3.Connection:
    if OUT.exists():
        OUT.unlink()
    con = sqlite3.connect(str(OUT))
    con.execute("PRAGMA journal_mode=DELETE")
    con.execute("PRAGMA synchronous=OFF")
    con.executescript(
        """
        CREATE TABLE targets (
            tool TEXT NOT NULL,
            mirna TEXT NOT NULL,
            symbol TEXT NOT NULL
        );
        CREATE TABLE tools (
            tool TEXT PRIMARY KEY,
            n INTEGER NOT NULL
        );
        """
    )
    return con


def _flush(con: sqlite3.Connection, tool: str, rows: list[tuple[str, str]]) -> int:
    con.executemany(
        "INSERT INTO targets(tool, mirna, symbol) VALUES (?,?,?)",
        ((tool, m, s) for m, s in rows),
    )
    con.execute(
        "INSERT INTO tools(tool, n) VALUES (?,?)",
        (tool, len(rows)),
    )
    con.commit()
    print(f"  {tool}: {len(rows):,} rows", flush=True)
    return len(rows)


def build_targetscan(con: sqlite3.Connection, mirnas: list[str]) -> None:
    path = DATA / "Predicted_Targets_Info.default_predictions.txt"
    families: dict[str, set[str]] = defaultdict(set)
    t0 = time.time()
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or p[4] != "9606":
                continue
            if p[9] not in ("8mer", "7mer-m8"):
                continue
            families[p[0]].add(p[2])
    print(f"  TargetScan families {len(families):,} in {time.time()-t0:.0f}s", flush=True)
    rows: list[tuple[str, str]] = []
    fam_items = list(families.items())
    for mirna in mirnas:
        wanted = resolve_targetscan_family(mirna)
        genes: set[str] = set()
        for fam_col, syms in fam_items:
            if _family_column_matches(fam_col, wanted, mirna):
                genes |= syms
        rows.extend((mirna, g) for g in genes)
    _flush(con, "TargetScan", rows)


def build_mirdb(con: sqlite3.Connection) -> None:
    mp = dict(
        pd.read_csv(DATA / "refseq_to_symbol.tsv", sep="\t")[["refseq", "symbol"]].itertuples(
            index=False
        )
    )
    rows: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    with gzip.open(DATA / "miRDB_v6.0_prediction_result.txt.gz", "rt") as fh:
        for line in fh:
            mir, ref, score = line.rstrip("\n").split("\t")
            if float(score) < MIRDB_SCORE_MIN:
                continue
            sym = mp.get(ref.split(".")[0])
            if not sym:
                continue
            key = (mir, sym)
            if key in seen:
                continue
            seen.add(key)
            rows.append(key)
    _flush(con, "miRDB", rows)


def build_diana(con: sqlite3.Connection) -> None:
    utr = pd.read_parquet(ROOT / "utr3_human.parquet", columns=["gene_id", "symbol"])
    e2s = dict(zip(utr["gene_id"].astype(str), utr["symbol"].astype(str)))
    rows: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    path = DATA / "interactions_human.microT.mirbase.txt.gz"
    n = 0
    t0 = time.time()
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            n += 1
            mir, gid, score = line.rstrip("\n").split("\t")
            if float(score) < DIANA_SCORE_MIN:
                continue
            sym = e2s.get(gid.split(".")[0])
            if not sym:
                continue
            key = (mir, sym)
            if key in seen:
                continue
            seen.add(key)
            rows.append(key)
            if n % 5_000_000 == 0:
                print(f"  DIANA scanned {n:,} lines, kept {len(rows):,}", flush=True)
    print(f"  DIANA scanned {n:,} lines in {time.time()-t0:.0f}s", flush=True)
    _flush(con, "DIANA-microT", rows)


def build_mirmap(con: sqlite3.Connection, accession_to_mirnas: dict[str, list[str]]) -> None:
    utr = pd.read_parquet(
        ROOT / "utr3_human.parquet", columns=["transcript_id", "symbol"]
    )
    tx_to_symbol = dict(
        zip(utr["transcript_id"].astype(str), utr["symbol"].astype(str))
    )
    path = DATA / "mirmap_202203_homsap_targets_1to1.csv.zst"
    rows: list[tuple[str, str]] = []
    t0 = time.time()
    n_mir = 0

    def emit(mimat: str, best: dict[str, float]) -> None:
        nonlocal n_mir
        names = accession_to_mirnas.get(mimat) or []
        if not names or not best:
            return
        scores = pd.Series(list(best.values()), dtype="float64")
        symbols = list(best.keys())
        pct = scores.rank(method="average", ascending=False, pct=True).to_numpy() * 100.0
        kept = [sym for sym, p in zip(symbols, pct) if p >= MIRMAP_PERCENTILE_MIN]
        for mirna in names:
            rows.extend((mirna, sym) for sym in kept)
        n_mir += 1
        if n_mir % 250 == 0:
            print(
                f"  miRmap {n_mir} miRNAs, kept {len(rows):,}, {time.time()-t0:.0f}s",
                flush=True,
            )

    dctx = zstandard.ZstdDecompressor()
    with open(path, "rb") as fh, dctx.stream_reader(fh) as reader:
        text = io.TextIOWrapper(reader, encoding="utf-8")
        header = text.readline().rstrip("\n").split(",")
        mi = header.index("mirna_id")
        ti = header.index("transcript_stable_id")
        si = header.index("mirmap_score")
        current: str | None = None
        best: dict[str, float] = {}
        for line in text:
            p = line.rstrip("\n").split(",")
            mimat = p[mi]
            if current is None:
                current = mimat
            elif mimat != current:
                emit(current, best)
                current = mimat
                best = {}
            sym = tx_to_symbol.get(p[ti])
            if not sym:
                continue
            score = float(p[si])
            prev = best.get(sym)
            if prev is None or score < prev:
                best[sym] = score
        if current is not None:
            emit(current, best)
    print(f"  miRmap finished {n_mir} miRNAs in {time.time()-t0:.0f}s", flush=True)
    _flush(con, "miRmap", rows)


def build_mirtarbase(con: sqlite3.Connection) -> None:
    path = DATA / "mirtarbase" / "hsa_MTI.tsv"
    raw = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    colmap = {c.lower().replace(" ", "_"): c for c in raw.columns}
    mir_c = colmap["mirna"]
    gene_c = colmap["target_gene"]
    exp_c = colmap.get("experiments")
    strong = ("luciferase", "reporter", "western", "qpcr", "qrt-pcr", "immunoblot")
    rows: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for mir, gene, exp in zip(raw[mir_c], raw[gene_c], raw[exp_c] if exp_c else [""] * len(raw)):
        if not isinstance(mir, str) or not isinstance(gene, str):
            continue
        exp_l = str(exp).lower() if isinstance(exp, str) else ""
        if exp_c and not any(s in exp_l for s in strong):
            continue
        symbols = {gene}
        if gene.islower():
            symbols.add(gene.upper())
        for sym in symbols:
            key = (mir, sym)
            if key not in seen:
                seen.add(key)
                rows.append(key)
    _flush(con, "miRTarBase", rows)


def main() -> None:
    t0 = time.time()
    src = sqlite3.connect(str(ROOT / "o8g_targets.db"))
    mirnas = [r[0] for r in src.execute("SELECT mirna FROM mirnas ORDER BY mirna")]
    accession: dict[str, list[str]] = defaultdict(list)
    for mir, acc in src.execute(
        "SELECT mirna, accession FROM mirnas WHERE accession IS NOT NULL AND accession != ''"
    ):
        if mir not in accession[acc]:
            accession[acc].append(mir)
    src.close()
    con = _connect()
    print("TargetScan…", flush=True)
    build_targetscan(con, mirnas)
    print("miRDB…", flush=True)
    build_mirdb(con)
    print("DIANA-microT…", flush=True)
    build_diana(con)
    print("miRTarBase…", flush=True)
    build_mirtarbase(con)
    print("miRmap…", flush=True)
    build_mirmap(con, accession)
    print("indexing…", flush=True)
    con.execute("CREATE INDEX idx_targets_tool_mirna ON targets(tool, mirna)")
    con.commit()
    con.close()
    size = OUT.stat().st_size
    print(f"wrote {OUT} ({size/1e6:.1f} MB) in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
