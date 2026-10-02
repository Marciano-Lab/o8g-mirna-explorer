#!/usr/bin/env python3
"""Invert o8g_targets.db state→genes into gene→(seed,state,site) reverse index.

Writes o8g_reverse.db with:
  gene_targets(gene_idx, seed, label, site_rank, oxidized_positions)
  seed_mirnas(seed, mirna)   -- all mature names sharing a seed
"""
from __future__ import annotations

import argparse
import sqlite3
import time
import zlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def build(db_path: Path, out_path: Path):
    src = sqlite3.connect(str(db_path))
    if out_path.exists():
        out_path.unlink()
    dst = sqlite3.connect(str(out_path))
    dst.execute("PRAGMA journal_mode=WAL")
    dst.executescript(
        """
        CREATE TABLE gene_targets(
            gene_idx INTEGER NOT NULL,
            seed TEXT NOT NULL,
            label TEXT NOT NULL,
            site_rank INTEGER NOT NULL,
            oxidized_positions TEXT NOT NULL
        );
        CREATE TABLE seed_mirnas(
            seed TEXT NOT NULL,
            mirna TEXT NOT NULL,
            PRIMARY KEY(seed, mirna)
        );
        """
    )

    mir = src.execute("SELECT seed, mirna FROM mirnas").fetchall()
    dst.executemany("INSERT OR IGNORE INTO seed_mirnas(seed, mirna) VALUES (?,?)", mir)
    dst.commit()
    print(f"seed_mirnas: {len(mir)} rows", flush=True)

    cur = src.execute(
        "SELECT seed, label, oxidized_positions, gene_blob, rank_blob FROM states"
    )
    batch = []
    n_states = 0
    n_pairs = 0
    t0 = time.time()
    while True:
        rows = cur.fetchmany(200)
        if not rows:
            break
        for seed, label, ox, gblob, rblob in rows:
            gidx = np.frombuffer(zlib.decompress(gblob), dtype=np.int32)
            rnk = np.frombuffer(zlib.decompress(rblob), dtype=np.int8)
            ox = ox or ""
            for gi, rk in zip(gidx.tolist(), rnk.tolist()):
                batch.append((int(gi), seed, label, int(rk), ox))
            n_states += 1
            n_pairs += len(gidx)
        if batch:
            dst.executemany(
                "INSERT INTO gene_targets(gene_idx, seed, label, site_rank, oxidized_positions) "
                "VALUES (?,?,?,?,?)",
                batch,
            )
            dst.commit()
            batch.clear()
        if n_states % 1000 == 0:
            print(
                f"  {n_states} states, {n_pairs:,} pairs, {time.time()-t0:.0f}s",
                flush=True,
            )

    print("creating indexes…", flush=True)
    dst.execute("CREATE INDEX idx_gt_gene ON gene_targets(gene_idx)")
    dst.execute("CREATE INDEX idx_gt_seed ON gene_targets(seed)")
    dst.execute("CREATE INDEX idx_sm_seed ON seed_mirnas(seed)")
    dst.commit()
    src.close()
    dst.close()
    print(f"done → {out_path}: {n_states} states, {n_pairs:,} gene-state pairs, {time.time()-t0:.0f}s")


def build_compact(db_path: Path, out_path: Path):
    """Gene-keyed zlib payloads instead of one SQLite row per gene–state pair.

    The row-oriented reverse index is ~3 GB (about 51e6 strong-site pairs).
    This form keeps the same pairs in a file small enough to ship with the app.
    Payload records are little-endian uint16 state_id + uint8 site_rank.
    """
    import struct
    import zlib

    src = sqlite3.connect(str(db_path))
    if out_path.exists():
        out_path.unlink()
    dst = sqlite3.connect(str(out_path))
    dst.execute("PRAGMA journal_mode=DELETE")
    dst.executescript(
        """
        CREATE TABLE rev_states(
            state_id INTEGER PRIMARY KEY,
            seed TEXT NOT NULL,
            label TEXT NOT NULL,
            oxidized_positions TEXT NOT NULL
        );
        CREATE TABLE gene_rev(
            gene_idx INTEGER PRIMARY KEY,
            payload BLOB NOT NULL
        );
        CREATE TABLE seed_mirnas(
            seed TEXT NOT NULL,
            mirna TEXT NOT NULL,
            PRIMARY KEY(seed, mirna)
        );
        CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
        """
    )
    mir = src.execute("SELECT seed, mirna FROM mirnas").fetchall()
    dst.executemany(
        "INSERT OR IGNORE INTO seed_mirnas(seed, mirna) VALUES (?,?)", mir
    )
    dst.commit()
    print(f"seed_mirnas: {len(mir)} rows", flush=True)

    pack = struct.Struct("<HB").pack
    buckets: dict[int, bytearray] = {}
    n_states = 0
    n_pairs = 0
    t0 = time.time()
    cur = src.execute(
        "SELECT state_id, seed, label, oxidized_positions, gene_blob, rank_blob "
        "FROM states ORDER BY state_id"
    )
    state_rows = []
    for pack_id, (state_id, seed, label, ox, gblob, rblob) in enumerate(cur):
        if int(state_id) > 65535:
            raise SystemExit(f"state_id {state_id} does not fit in uint16")
        state_rows.append((pack_id, seed, label or "", ox or ""))
        gidx = np.frombuffer(zlib.decompress(gblob), dtype=np.int32)
        rnk = np.frombuffer(zlib.decompress(rblob), dtype=np.int8)
        for gi, rk in zip(gidx.tolist(), rnk.tolist()):
            buf = buckets.get(gi)
            if buf is None:
                buf = bytearray()
                buckets[gi] = buf
            buf += pack(pack_id, int(rk))
        n_states += 1
        n_pairs += len(gidx)
        if n_states % 2000 == 0:
            print(
                f"  {n_states} states, {n_pairs:,} pairs, {time.time()-t0:.0f}s",
                flush=True,
            )
    dst.executemany(
        "INSERT INTO rev_states(state_id, seed, label, oxidized_positions) VALUES (?,?,?,?)",
        state_rows,
    )
    print(f"compressing {len(buckets):,} gene payloads…", flush=True)
    payload_rows = [
        (int(gi), zlib.compress(bytes(buf), level=6)) for gi, buf in buckets.items()
    ]
    dst.executemany(
        "INSERT INTO gene_rev(gene_idx, payload) VALUES (?,?)", payload_rows
    )
    dst.executemany(
        "INSERT INTO meta(key, value) VALUES (?,?)",
        [
            ("format", "compact-v1"),
            ("n_states", str(n_states)),
            ("n_pairs", str(n_pairs)),
            ("n_genes", str(len(buckets))),
        ],
    )
    dst.commit()
    src.close()
    dst.close()
    print(
        f"done → {out_path}: {n_states} states, {n_pairs:,} pairs, "
        f"{len(buckets):,} genes, {time.time()-t0:.0f}s"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=ROOT / "o8g_targets.db")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument(
        "--compact",
        action="store_true",
        help="Write the shippable gene-keyed index (o8g_reverse.compact.db)",
    )
    args = ap.parse_args()
    if args.compact:
        out = args.out or (ROOT / "o8g_reverse.compact.db")
        build_compact(args.db, out)
    else:
        out = args.out or (ROOT / "o8g_reverse.db")
        build(args.db, out)


if __name__ == "__main__":
    main()
