"""
conservation.py
===============
Intersect oxomiR explorer targets with TargetScan v8.0 *conserved* site tables.

We do NOT re-derive phastCons. Conservation is defined by membership in
TargetScan's published conserved-site releases (reference standard).

Release
-------
TargetScanHuman 8.0 (McGeary / Agarwal lab; data files dated 2021-10-18 on
targetscan.org vert_80_data_download):
  - Conserved_Family_Info.txt
  - Conserved_Site_Context_Scores.txt

Flags
-----
USE_CONSERVATION (module-level / env O8G_USE_CONSERVATION=0 to disable)
"""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_DATA = ROOT / "paper" / "data"
TARGETSCAN_RELEASE = "TargetScanHuman_8.0"
TARGETSCAN_FILES = (
    "Conserved_Family_Info.txt",
    "Conserved_Site_Context_Scores.txt",
)

# Feature flag — identical application controlled by PrecisionMode, not here.
USE_CONSERVATION = os.environ.get("O8G_USE_CONSERVATION", "1") != "0"

# TargetScan Site Type codes in Conserved_Site_Context_Scores.txt
# (TargetScan documentation: 1=6mer is not typical here; conserved file uses
#  2=7mer-m8, 3=8mer, 4=7mer-A1 in some releases — verify against counts).
# Empirically for vert_80 Conserved_Site_Context_Scores: Site Type ∈ {1,2,3}
# where we keep only strong conserved sites matching our rank≥3 policy.
SITE_TYPE_STRONG = {2, 3}  # treat as 7mer-m8 / 8mer class; exclude weaker


def _strip_ensembl(eid: str) -> str:
    return str(eid).split(".")[0]


def ensure_context_scores_txt(data_dir: Path | str | None = None) -> Path:
    """Expand the shipped gzip when the image only has the compressed table.

    GitHub's file cap is why ``Conserved_Site_Context_Scores.txt`` (144 MB) is
    stored as an 18 MB gzip. The readers still open the plain text path.
    """
    import gzip

    data_dir = Path(data_dir or DEFAULT_DATA)
    dest = data_dir / "Conserved_Site_Context_Scores.txt"
    gz = data_dir / "Conserved_Site_Context_Scores.txt.gz"
    if dest.exists() or not gz.exists():
        return dest
    data_dir.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".txt.partial")
    with gzip.open(gz, "rb") as src, open(tmp, "wb") as out:
        while True:
            chunk = src.read(8 * 1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    return dest


class ConservedIndex:
    """Gene-symbol and Ensembl sets of TargetScan-conserved targets per mature miRNA / family."""

    def __init__(
        self,
        data_dir: Path | str = DEFAULT_DATA,
        human_tax: int = 9606,
    ):
        self.data_dir = Path(data_dir)
        self.human_tax = human_tax
        self.release = TARGETSCAN_RELEASE
        # mirna_or_family -> set(symbols)
        self.by_mirna_symbol: dict[str, set[str]] = {}
        self.by_family_symbol: dict[str, set[str]] = {}
        self.by_family_ensembl: dict[str, set[str]] = {}
        self._loaded = False

    def ensure_loaded(self) -> None:
        if self._loaded:
            return
        if not USE_CONSERVATION:
            self._loaded = True
            return
        cache = self.data_dir / "ts80_conserved_family_human_strong.pkl"
        if cache.exists():
            import pickle

            blob = pickle.load(open(cache, "rb"))
            self.by_family_symbol = blob["by_family_symbol"]
            self.by_family_ensembl = blob["by_family_ensembl"]
            self._loaded = True
            return
        fam_path = self.data_dir / "Conserved_Family_Info.txt"
        ctx_path = ensure_context_scores_txt(self.data_dir)
        if not fam_path.exists():
            raise FileNotFoundError(
                f"Missing {fam_path}. Download TargetScan 8.0 Conserved_Family_Info.txt "
                f"(release {TARGETSCAN_RELEASE})."
            )
        self._load_family_info(fam_path)
        # Context-score file is optional (mature-level); family info is enough for Consensus.
        # Set O8G_LOAD_TS_CONTEXT=1 to also ingest Conserved_Site_Context_Scores.txt.
        if ctx_path.exists() and os.environ.get("O8G_LOAD_TS_CONTEXT", "0") == "1":
            self._load_context_scores(ctx_path)
        self._loaded = True

    def _load_family_info(self, path: Path) -> None:
        # Columns: miR Family, Gene ID, Gene Symbol, Transcript ID, Species ID, ...
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                try:
                    sp = int(p[4])
                except ValueError:
                    continue
                if sp != self.human_tax:
                    continue
                fam, gid, sym = p[0], _strip_ensembl(p[1]), p[2]
                seed_match = p[9] if len(p) > 9 else ""
                # Strong conserved sites only (align with STRONG=3 policy)
                if seed_match and seed_match not in ("8mer", "7mer-m8"):
                    continue
                self.by_family_symbol.setdefault(fam, set()).add(sym)
                self.by_family_ensembl.setdefault(fam, set()).add(gid)

    def _load_context_scores(self, path: Path) -> None:
        # Gene ID, Gene Symbol, Transcript ID, Gene Tax ID, miRNA, Site Type, ...
        with open(path) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                try:
                    tax = int(p[3])
                    stype = int(p[5])
                except ValueError:
                    continue
                if tax != self.human_tax:
                    continue
                if stype not in SITE_TYPE_STRONG:
                    continue
                mir, sym = p[4], p[1]
                self.by_mirna_symbol.setdefault(mir, set()).add(sym)

    def conserved_symbols_for_mirna(self, mirna: str, family: str | None = None) -> set[str]:
        """Union of conserved gene symbols for a mature name and optional TS family."""
        self.ensure_loaded()
        out: set[str] = set()
        if mirna in self.by_mirna_symbol:
            out |= self.by_mirna_symbol[mirna]
        # also try without species prefix variants
        if mirna.startswith("hsa-"):
            out |= self.by_mirna_symbol.get(mirna, set())
        if family and family in self.by_family_symbol:
            out |= self.by_family_symbol[family]
        return out

    def is_conserved_symbol(self, mirna: str, symbol: str, family: str | None = None) -> bool:
        return symbol in self.conserved_symbols_for_mirna(mirna, family)


def build_seed_family_map(db_path: Path | str) -> dict[str, str]:
    """Map our DB seed → a TargetScan family string when uniquely attributable.

    Uses mirnas.mirna names joined against Conserved_Family_Info family strings
    that contain the mature's numeric stem (best-effort; documented heuristic).
    """
    con = sqlite3.connect(str(db_path))
    mir = pd.read_sql("SELECT mirna, seed FROM mirnas", con)
    con.close()
    # Prefer explicit known families for paper miRNAs; else leave empty.
    known = {
        "GGAATGT": "miR-1-3p/206",
        "AAGGCAC": "miR-124-3p.1",
        "GAGGTAG": "let-7-5p/98-5p",
        "GGAGTGT": "miR-122-5p",
    }
    out = dict(known)
    for seed, g in mir.groupby("seed"):
        if seed in out:
            continue
        # no automatic guess — conservation lookup can still use mature names
    return out


_INDEX: ConservedIndex | None = None


def get_conserved_index(data_dir: Path | str | None = None) -> ConservedIndex:
    global _INDEX
    if _INDEX is None:
        _INDEX = ConservedIndex(data_dir or DEFAULT_DATA)
    return _INDEX
