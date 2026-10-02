"""
o8g_db.py
=========
Read-only data layer over the precomputed o8G target database (o8g_targets.db).

Schema (v1 → v2 precision migration)
------------------------------------
states columns (v1): gene_blob, rank_blob  (rank>=3 strong sites)
optional (v2): n8_blob, n7m8_blob, cons_blob  (parallel int8 arrays)

Backward compatible: missing optional blobs → columns filled with defaults /
on-the-fly enrichment when a TargetScanner is supplied.
"""
from __future__ import annotations

import sqlite3
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import o8g_precision as _o8g_precision
from o8g_precision import PrecisionConfig, PrecisionMode, partition_after_filter


def connect_sqlite(path, *, check_same_thread: bool = False) -> sqlite3.Connection:
    """Open SQLite without mmap.

    The hosted catalog is ~300 MB. On Linux, SQLite's default mmap counts
    toward the container memory limit and kills the 512 MB Render instance
    as soon as External DB comparison opens that file.
    """
    con = sqlite3.connect(str(path), check_same_thread=check_same_thread)
    con.execute("PRAGMA mmap_size=0")
    con.execute("PRAGMA cache_size=-512")
    return con


def _live_precision():
    """Use the already-imported precision module.

    Reloading it on every rerun duplicated module state on the small host.
    """
    return _o8g_precision


def _apply_precision_filter(*args, **kwargs):
    """Call through the live module so Streamlit reloads always take effect."""
    return _live_precision().apply_precision_filter(*args, **kwargs)

RANK_SITE = {4: "8mer", 3: "7mer-m8", 2: "7mer-A1", 1: "6mer"}

_ROOT = Path(__file__).resolve().parent
_FULL_REVERSE = _ROOT / "o8g_reverse.db"
_COMPACT_REVERSE = _ROOT / "o8g_reverse.compact.db"
SCHEMA_VERSION_KEY = "schema_version"


def _assemble_compact_reverse() -> None:
    """Join shipped parts when the compact index itself is not in the checkout.

    GitHub rejects the ~128 MB index, so the repo stores ``*.part*`` slices and
    the image concatenates them. A local checkout that already has the file
    is left alone.
    """
    if _COMPACT_REVERSE.exists():
        return
    parts = sorted(_ROOT.glob("o8g_reverse.compact.db.part*"))
    if not parts:
        return
    tmp = _COMPACT_REVERSE.with_suffix(".db.partial")
    with open(tmp, "wb") as out:
        for part in parts:
            with open(part, "rb") as fh:
                while True:
                    chunk = fh.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
    tmp.replace(_COMPACT_REVERSE)


def _default_reverse_path() -> Path:
    """Full row-oriented index when present; otherwise the shippable compact index.

    Hosted deploys only have ``o8g_reverse.compact.db`` (assembled from parts).
    A local checkout that still has the 3 GB ``o8g_reverse.db`` keeps using it.
    """
    if _FULL_REVERSE.exists():
        return _FULL_REVERSE
    _assemble_compact_reverse()
    return _COMPACT_REVERSE


class ConservationUnavailable(RuntimeError):
    """Raised when TargetScan conservation data cannot be loaded.

    Never substitute an empty conserved set: Consensus intersects the
    unmodified baseline with it, so an empty set silently collapses the
    baseline to zero and reports every oxidized target as gained.
    """


def _primary_mirna(names: list[str]) -> str:
    if not names:
        return ""
    threep = [n for n in names if n.endswith("-3p")]
    pool = threep or names
    for pref in ("hsa-miR-1-3p", "hsa-miR-1", "hsa-miR-124-3p"):
        if pref in pool:
            return pref
    return sorted(pool)[0]


class TargetDB:
    def __init__(self, path: str = "o8g_targets.db", reverse_path: str | Path | None = None):
        self.path = path
        self._con = connect_sqlite(path)
        g = pd.read_sql("SELECT gene_idx, gene_id, symbol FROM genes", self._con)
        self.symbols = g.sort_values("gene_idx")["symbol"].to_numpy()
        self.gene_ids = g.sort_values("gene_idx")["gene_id"].to_numpy()
        self._gene_table = g.set_index("gene_idx")
        self._reverse_candidate = Path(reverse_path) if reverse_path else _default_reverse_path()
        self._rev_compact: bool | None = None
        self._rev_seed: np.ndarray | None = None
        self._rev_label: np.ndarray | None = None
        self._rev_ox: np.ndarray | None = None
        self._rev: sqlite3.Connection | None = None
        self._state_cols = {
            r[1] for r in self._con.execute("PRAGMA table_info(states)").fetchall()
        }
        self.schema_version = self._read_schema_version()

    def _read_schema_version(self) -> int:
        tables = {
            r[0]
            for r in self._con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        if "meta" not in tables:
            return 1
        row = self._con.execute(
            "SELECT value FROM meta WHERE key=?", [SCHEMA_VERSION_KEY]
        ).fetchone()
        return int(row[0]) if row else 1

    @property
    def reverse_path(self) -> Path | None:
        return self._reverse_candidate if self._reverse_candidate.exists() else None

    def _rev_con(self) -> sqlite3.Connection:
        if self.reverse_path is None:
            raise FileNotFoundError(
                "o8g_reverse.db not found. Run: python scripts/build_reverse_index.py"
            )
        if self._rev is None:
            self._rev = connect_sqlite(self.reverse_path)
        return self._rev

    def mirnas(self) -> pd.DataFrame:
        return pd.read_sql(
            "SELECT mirna, accession, seq_dna, seed, n_G FROM mirnas ORDER BY mirna",
            self._con,
        )

    def search_mirnas(self, query: str) -> pd.DataFrame:
        q = f"%{query}%"
        return pd.read_sql(
            "SELECT mirna, accession, seq_dna, seed, n_G FROM mirnas "
            "WHERE mirna LIKE ? ORDER BY mirna",
            self._con,
            params=[q],
        )

    def mirna_info(self, mirna: str) -> dict | None:
        r = self._con.execute(
            "SELECT mirna, accession, seq_dna, seed, n_G FROM mirnas WHERE mirna=?",
            [mirna],
        ).fetchone()
        if not r:
            return None
        return dict(zip(["mirna", "accession", "seq_dna", "seed", "n_G"], r))

    def states_for_seed(self, seed: str) -> pd.DataFrame:
        return pd.read_sql(
            "SELECT state_id, label, oxidized_positions, motif_6mer, motif_7mer_m8, motif_8mer, "
            "n_6mer, n_7mer_A1, n_7mer_m8, n_8mer, n_strong FROM states WHERE seed=? ORDER BY state_id",
            self._con,
            params=[seed],
        )

    def _decode_state_blobs(self, seed: str, label: str) -> pd.DataFrame:
        cols = ["gene_blob", "rank_blob"]
        optional = []
        for c in ("n8_blob", "n7m8_blob", "cons_blob"):
            if c in self._state_cols:
                optional.append(c)
                cols.append(c)
        sql = f"SELECT {', '.join(cols)} FROM states WHERE seed=? AND label=?"
        r = self._con.execute(sql, [seed, label]).fetchone()
        empty = pd.DataFrame(
            columns=[
                "gene_idx",
                "gene_id",
                "symbol",
                "site_type",
                "site_rank",
                "n_8mer",
                "n_7mer_m8",
                "n_sites",
                "score",
                "is_conserved",
            ]
        )
        if not r:
            return empty
        gidx = np.frombuffer(zlib.decompress(r[0]), dtype=np.int32)
        rnk = np.frombuffer(zlib.decompress(r[1]), dtype=np.int8)
        n = len(gidx)
        n8 = np.zeros(n, dtype=np.int8)
        n7 = np.zeros(n, dtype=np.int8)
        cons = np.zeros(n, dtype=np.int8)
        # Map optional blobs by name
        colmap = {c: r[i] for i, c in enumerate(cols)}
        if colmap.get("n8_blob"):
            n8 = np.frombuffer(zlib.decompress(colmap["n8_blob"]), dtype=np.int8)
        if colmap.get("n7m8_blob"):
            n7 = np.frombuffer(zlib.decompress(colmap["n7m8_blob"]), dtype=np.int8)
        if colmap.get("cons_blob"):
            cons = np.frombuffer(zlib.decompress(colmap["cons_blob"]), dtype=np.int8)
        # If multiplicity blobs missing, approximate from best rank:
        # 8mer → n8=1; 7mer-m8 → n7=1 (underestimates multi-site genes).
        if "n8_blob" not in colmap or colmap["n8_blob"] is None:
            n8 = (rnk == 4).astype(np.int8)
            n7 = (rnk == 3).astype(np.int8)
        score = n8.astype(np.float64) * 1.0 + n7.astype(np.float64) * 0.7
        n_sites = n8.astype(np.int32) + n7.astype(np.int32)
        return pd.DataFrame(
            {
                "gene_idx": gidx,
                "gene_id": self.gene_ids[gidx],
                "symbol": self.symbols[gidx],
                "site_rank": rnk,
                "site_type": [RANK_SITE[int(x)] for x in rnk],
                "n_8mer": n8.astype(int),
                "n_7mer_m8": n7.astype(int),
                "n_sites": n_sites,
                "score": score,
                "is_conserved": cons.astype(bool),
            }
        )

    def targets(self, seed: str, label: str) -> pd.DataFrame:
        """Strong-site target genes for one seed-state (backward-compatible API)."""
        df = self._decode_state_blobs(seed, label)
        return (
            df.drop(columns=["gene_idx"], errors="ignore")
            .sort_values(["site_rank", "symbol"], ascending=[False, True])
            .reset_index(drop=True)
        )

    def targets_enriched(
        self,
        seed: str,
        label: str,
        *,
        scanner=None,
        mature_dna: str | None = None,
        conserved_symbols: set[str] | None = None,
    ) -> pd.DataFrame:
        """Targets with score / multiplicity / optional live context + conservation."""
        from o8g_engine import SeedState

        df = self._decode_state_blobs(seed, label)
        if conserved_symbols is not None:
            df["is_conserved"] = df["symbol"].isin(conserved_symbols)
        if scanner is not None:
            # Live multiplicity + context from UTR index (authoritative when present)
            ox = []
            if label != "none":
                # parse o8G@2,7 → (2,7)
                part = label.replace("o8G@", "")
                ox = [int(x) for x in part.split(",") if x]
            state = SeedState(seed, tuple(ox))
            live = scanner.scan_state_context(state, mature_dna=mature_dna, min_rank=3)
            if not live.empty:
                keep = [
                    "symbol",
                    "gene_idx",
                    "n_8mer",
                    "n_7mer_m8",
                    "n_sites",
                    "score",
                    "context_score",
                    "site_rank",
                    "site_type",
                    "site_start",
                    "gene_id",
                ]
                live = live[[c for c in keep if c in live.columns]]
                df = live.merge(
                    df[["symbol", "is_conserved"]],
                    on="symbol",
                    how="left",
                )
                df["is_conserved"] = df["is_conserved"].fillna(False)
        return df.sort_values(["site_rank", "symbol"], ascending=[False, True]).reset_index(
            drop=True
        )

    def _conserved_index(self):
        from conservation import get_conserved_index

        return get_conserved_index()

    def _conserved_for(self, seed, mirna):
        if not mirna:
            raise ConservationUnavailable(
                "Consensus mode needs the miRNA name to look up TargetScan "
                "conserved families; none was supplied."
            )
        try:
            from conservation import build_seed_family_map

            fam = build_seed_family_map(self.path).get(seed)
            syms = self._conserved_index().conserved_symbols_for_mirna(mirna, fam)
        except FileNotFoundError as e:
            raise ConservationUnavailable(
                "Consensus mode requires paper/data/Conserved_Family_Info.txt "
                "(TargetScanHuman 8.0), which is not installed. Use sequence-based "
                "(low/high stringency), or add the TargetScan release files."
            ) from e
        except ConservationUnavailable:
            raise
        except Exception as e:
            raise ConservationUnavailable(
                f"Consensus conservation lookup failed ({type(e).__name__}): {e}"
            ) from e
        if not syms:
            raise ConservationUnavailable(
                f"TargetScan returned no conserved families for {mirna}; "
                "refusing to run Consensus on an empty baseline."
            )
        return syms

    def _targetscan_for(self, mirna: str | None) -> set[str]:
        if not mirna:
            raise ConservationUnavailable(
                "TargetScan mode needs the miRNA name to look up Predicted_Targets_Info."
            )
        try:
            import o8g_refsets as refsets

            if not refsets.available_tools().get("TargetScan"):
                raise ConservationUnavailable(
                    "TargetScan mode requires paper/data/Predicted_Targets_Info."
                    "default_predictions.txt (TargetScanHuman 8.0 predictions)."
                )
            syms = {str(s).upper() for s in refsets.load_targetscan(mirna)}
        except ConservationUnavailable:
            raise
        except Exception as e:
            raise ConservationUnavailable(
                f"TargetScan prediction lookup failed ({type(e).__name__}): {e}"
            ) from e
        if not syms:
            raise ConservationUnavailable(
                f"TargetScan returned no predicted strong sites for {mirna}; "
                "refusing to run TargetScan mode on an empty baseline."
            )
        return syms

    def _anchor_symbols_for(self, cfg: PrecisionConfig, seed, mirna: str | None):
        _op = _live_precision()
        mode_s = _op.mode_value(cfg)
        if mode_s == "Consensus":
            return self._conserved_for(seed, mirna)
        if mode_s == "TargetScan":
            return self._targetscan_for(mirna)
        return None

    def targets_filtered(
        self,
        seed: str,
        label: str,
        cfg: PrecisionConfig | PrecisionMode | str,
        *,
        scanner=None,
        mature_dna: str | None = None,
        conserved_symbols: set[str] | None = None,
        mirna: str | None = None,
    ) -> pd.DataFrame:
        # Always re-hydrate via the live o8g_precision module (Streamlit reload-safe)
        _op = _live_precision()
        cfg = _op.PrecisionConfig.from_mode(cfg)
        mode_s = _op.mode_value(cfg)
        if mode_s == "TargetScan de novo":
            # Same TargetScanS 8mer/7mer-m8 calls already stored by precompute.
            # Reading them avoids rebuilding the ~5 GB in-memory UTR index.
            df = self.targets(seed, label)
            if "site_rank" in df.columns and len(df):
                df = df[df["site_rank"] >= 3]
            df = df.copy()
            df["source"] = "targetscan_denovo_precomputed"
            return df.reset_index(drop=True)
        if conserved_symbols is None:
            conserved_symbols = self._anchor_symbols_for(cfg, seed, mirna)
        df = self.targets_enriched(
            seed,
            label,
            scanner=scanner,
            mature_dna=mature_dna,
            conserved_symbols=conserved_symbols,
        )
        return _apply_precision_filter(
            df,
            cfg,
            conserved_symbols=conserved_symbols,
            is_unmodified_state=(label == "none"),
        )

    def retarget_partition(
        self,
        seed: str,
        ox_label: str,
        cfg: PrecisionConfig | PrecisionMode | str,
        **kwargs,
    ) -> dict[str, set[str]]:
        """Partition unmod vs oxidized after filtering both sides."""
        _op = _live_precision()
        cfg = _op.PrecisionConfig.from_mode(cfg)
        mode_s = _op.mode_value(cfg)
        if mode_s == "TargetScan de novo":
            su = {str(s).upper() for s in self.target_symbols(seed, "none", min_rank=3)}
            so = {str(s).upper() for s in self.target_symbols(seed, ox_label, min_rank=3)}
            return {
                "unmod": su,
                "oxid": so,
                "shared": su & so,
                "lost": su - so,
                "gained": so - su,
            }
        mirna = kwargs.get("mirna")
        conserved_symbols = kwargs.get("conserved_symbols")
        if conserved_symbols is None:
            conserved_symbols = self._anchor_symbols_for(cfg, seed, mirna)
        unmod = self.targets_enriched(
            seed,
            "none",
            scanner=kwargs.get("scanner"),
            mature_dna=kwargs.get("mature_dna"),
            conserved_symbols=conserved_symbols,
        )
        oxid = self.targets_enriched(
            seed,
            ox_label,
            scanner=kwargs.get("scanner"),
            mature_dna=kwargs.get("mature_dna"),
            conserved_symbols=conserved_symbols,
        )
        return _live_precision().partition_after_filter(
            unmod, oxid, cfg, conserved_symbols=conserved_symbols
        )

    def target_symbols(self, seed: str, label: str, min_rank: int = 3) -> list[str]:
        df = self.targets(seed, label)
        return df.loc[df["site_rank"] >= min_rank, "symbol"].tolist()

    def gene_info(self, gene_idx: int) -> dict | None:
        if gene_idx not in self._gene_table.index:
            return None
        row = self._gene_table.loc[gene_idx]
        return {
            "gene_idx": int(gene_idx),
            "gene_id": str(row["gene_id"]),
            "symbol": str(row["symbol"]),
        }

    def _load_compact_lut(self) -> None:
        if self._rev_seed is not None:
            return
        lut = pd.read_sql(
            "SELECT state_id, seed, label, oxidized_positions FROM rev_states ORDER BY state_id",
            self._rev_con(),
        )
        ids = lut["state_id"].to_numpy()
        if len(ids) == 0 or not np.array_equal(ids, np.arange(len(ids))):
            raise ValueError("compact reverse rev_states.state_id is not 0..n-1")
        self._rev_seed = lut["seed"].to_numpy()
        self._rev_label = lut["label"].to_numpy()
        self._rev_ox = lut["oxidized_positions"].fillna("").to_numpy()

    def _reverse_is_compact(self) -> bool:
        if self._rev_compact is None:
            rev = self._rev_con()
            self._rev_compact = (
                rev.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='gene_rev'"
                ).fetchone()
                is not None
            )
        return self._rev_compact

    def _reverse_gene_hits(self, gene_idx: int) -> pd.DataFrame:
        """Strong-site rows for one gene: seed, state_label, site_rank, oxidized_positions."""
        cols = ["seed", "state_label", "site_rank", "oxidized_positions"]
        rev = self._rev_con()
        if self._reverse_is_compact():
            row = rev.execute(
                "SELECT payload FROM gene_rev WHERE gene_idx=?", [int(gene_idx)]
            ).fetchone()
            if not row or not row[0]:
                return pd.DataFrame(columns=cols)
            raw = zlib.decompress(row[0])
            b = np.frombuffer(raw, dtype=np.uint8)
            if len(b) % 3 != 0:
                raise ValueError("compact reverse payload length is not a multiple of 3")
            state_ids = b[0::3].astype(np.uint32) | (b[1::3].astype(np.uint32) << 8)
            ranks = b[2::3].astype(int)
            self._load_compact_lut()
            return pd.DataFrame(
                {
                    "seed": self._rev_seed[state_ids],
                    "state_label": self._rev_label[state_ids],
                    "site_rank": ranks,
                    "oxidized_positions": self._rev_ox[state_ids],
                }
            )
        return pd.read_sql(
            "SELECT seed, label AS state_label, site_rank, oxidized_positions "
            "FROM gene_targets WHERE gene_idx=?",
            rev,
            params=[int(gene_idx)],
        )

    def states_targeting_gene(self, gene_idx: int) -> pd.DataFrame:
        rev = self._rev_con()
        gt = self._reverse_gene_hits(gene_idx)
        empty_cols = [
            "mirna",
            "all_mirnas",
            "seed",
            "state_label",
            "oxidized_positions",
            "site_rank",
            "site_type",
            "motif_7mer_m8",
            "motif_8mer",
            "vs_unmodified",
        ]
        if gt.empty:
            return pd.DataFrame(columns=empty_cols)

        seeds = gt["seed"].unique().tolist()
        ph = ",".join("?" * len(seeds))
        sm = pd.read_sql(
            f"SELECT seed, mirna FROM seed_mirnas WHERE seed IN ({ph})",
            rev,
            params=seeds,
        )
        mir_map: dict[str, list[str]] = {}
        for seed, mirna in sm.itertuples(index=False):
            mir_map.setdefault(seed, []).append(mirna)

        motifs = pd.read_sql(
            f"SELECT seed, label, motif_7mer_m8, motif_8mer FROM states WHERE seed IN ({ph})",
            self._con,
            params=seeds,
        ).rename(columns={"label": "state_label"})

        none_seeds = set(gt.loc[gt["state_label"] == "none", "seed"])

        def vs_unmod(seed: str, label: str) -> str:
            if label == "none":
                return "unmodified"
            if seed in none_seeds:
                return "also in unmodified"
            return "gained on oxidation"

        gt = gt.merge(motifs, on=["seed", "state_label"], how="left")
        gt["mirna"] = gt["seed"].map(lambda s: _primary_mirna(mir_map.get(s, [])))
        gt["all_mirnas"] = gt["seed"].map(lambda s: ";".join(sorted(mir_map.get(s, []))))
        gt["site_type"] = gt["site_rank"].map(RANK_SITE)
        gt["vs_unmodified"] = [
            vs_unmod(s, lab) for s, lab in zip(gt["seed"], gt["state_label"])
        ]
        gt["oxidized_positions"] = gt["oxidized_positions"].fillna("")
        return (
            gt[empty_cols]
            .sort_values(["mirna", "state_label", "site_rank"], ascending=[True, True, False])
            .reset_index(drop=True)
        )

    def states_lost_on_oxidation(self, gene_idx: int) -> pd.DataFrame:
        """Unmodified strong-site hits that are absent under an oxidized state of the same seed.

        Reverse index only stores states that *target* the gene, so loss is inferred:
        seed has ``none`` → gene, and an oxidized label for that seed is not in
        ``gene_targets`` for this gene. Each row is one (seed, oxidized state) loss.
        ``site_rank`` / ``site_type`` are from the unmodified hit that was lost.
        """
        empty_cols = [
            "mirna",
            "all_mirnas",
            "seed",
            "state_label",
            "oxidized_positions",
            "site_rank",
            "site_type",
            "motif_7mer_m8",
            "motif_8mer",
            "vs_unmodified",
        ]
        rev = self._rev_con()
        gt = self._reverse_gene_hits(gene_idx)
        if gt.empty:
            return pd.DataFrame(columns=empty_cols)

        none_hits = gt[gt["state_label"] == "none"].copy()
        if none_hits.empty:
            return pd.DataFrame(columns=empty_cols)

        none_seeds = none_hits["seed"].unique().tolist()
        ph = ",".join("?" * len(none_seeds))

        # All precomputed states for seeds that target this gene when unmodified
        all_states = pd.read_sql(
            f"SELECT seed, label AS state_label, motif_7mer_m8, motif_8mer "
            f"FROM states WHERE seed IN ({ph})",
            self._con,
            params=none_seeds,
        )
        present = set(zip(gt["seed"].astype(str), gt["state_label"].astype(str)))
        ox_states = all_states[all_states["state_label"] != "none"].copy()
        ox_states["_key"] = list(zip(ox_states["seed"].astype(str), ox_states["state_label"].astype(str)))
        ox_states = ox_states[~ox_states["_key"].isin(present)].drop(columns=["_key"])
        if ox_states.empty:
            return pd.DataFrame(columns=empty_cols)

        # Attach unmodified site quality (what was lost)
        none_site = none_hits.drop_duplicates("seed").set_index("seed")["site_rank"]
        ox_states["site_rank"] = ox_states["seed"].map(none_site)

        sm = pd.read_sql(
            f"SELECT seed, mirna FROM seed_mirnas WHERE seed IN ({ph})",
            rev,
            params=none_seeds,
        )
        mir_map: dict[str, list[str]] = {}
        for seed, mirna in sm.itertuples(index=False):
            mir_map.setdefault(seed, []).append(mirna)

        def _ox_pos(label: str) -> str:
            if not label or label == "none":
                return ""
            return label.replace("o8G@", "")

        ox_states["mirna"] = ox_states["seed"].map(lambda s: _primary_mirna(mir_map.get(s, [])))
        ox_states["all_mirnas"] = ox_states["seed"].map(
            lambda s: ";".join(sorted(mir_map.get(s, [])))
        )
        ox_states["site_type"] = ox_states["site_rank"].map(RANK_SITE)
        ox_states["vs_unmodified"] = "lost on oxidation"
        ox_states["oxidized_positions"] = ox_states["state_label"].map(_ox_pos)
        return (
            ox_states[empty_cols]
            .sort_values(["mirna", "state_label", "site_rank"], ascending=[True, True, False])
            .reset_index(drop=True)
        )

    def close(self):
        self._con.close()
        if self._rev is not None:
            self._rev.close()
            self._rev = None
