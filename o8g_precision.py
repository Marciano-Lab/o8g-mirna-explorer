"""
o8g_precision.py
================
Precision ladder for oxomiR target lists.

CORE CONSTRAINT
---------------
Every filter is applied independently to the unmodified and oxidized target
tables; gained / lost / shared are set differences *after* filtering.
Never filter the delta itself.

Modes (PredictionMode / PrecisionMode)
--------------------------------------
UI label: **Prediction mode**.

Sequence-based (low stringency)  — formerly Sensitive; rank >= 3 (7mer-m8 + 8mer).
Sequence-based (high stringency) — formerly Stringent; rank == 4 (8mer only).
TargetScan — rank >= 3 AND TargetScan *predicted* strong sites on the *unmodified*
             baseline (Predicted_Targets_Info). Catalog WT anchor for lost/gained;
             oxidized-state *display* lists remain rank≥3 (same site types as de novo).
TargetScan de novo — precomputed TargetScanS 8mer/7mer-m8 lists on WT and oxidized seeds (o8G→T WC). Same scan as the explorer database; no live UTR index.
Consensus — rank >= 3 AND TargetScan-conserved on the unmodified baseline.

Deprecated aliases
----------------
``Sensitive`` / ``Stringent`` / ``HighConf`` still resolve via ``from_mode``.

Conservation note
-----------------
TargetScan conserved tables annotate *canonical* (unoxidized) seed matches.
Applying is_conserved==True identically to oxidized motifs would zero oxomiR
gains (those complementary words are absent from TargetScan). Consensus
therefore applies conservation to the unmodified baseline set and the same
*rank* filter to both states.

Feature flags (env or constructors)
-----------------------------------
O8G_USE_CONSERVATION=1|0
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum

import pandas as pd


class PrecisionMode(Enum):
    """Prediction-mode ladder (enum name kept for API stability).

    Not a ``str`` Enum — Streamlit reloads break ``isinstance(member, str)`` paths.
    """

    SENSITIVE = "Sequence-based (low stringency)"
    STRINGENT = "Sequence-based (high stringency)"
    TARGETSCAN = "TargetScan"
    TS_DENOVO = "TargetScan de novo"
    CONSENSUS = "Consensus"

    def __str__(self) -> str:
        return self.value


# Legacy UI / CSV labels → current mode value strings
_MODE_ALIASES = {
    "HighConf": "Sequence-based (high stringency)",
    "highconf": "Sequence-based (high stringency)",
    "Sensitive": "Sequence-based (low stringency)",
    "Stringent": "Sequence-based (high stringency)",
    "Sequence-based (low stringency)": "Sequence-based (low stringency)",
    "Sequence-based (high stringency)": "Sequence-based (high stringency)",
}


def mode_value(mode: object) -> str:
    """Extract the canonical mode label string from enum / config / str.

    Uses ``isinstance(..., Enum)`` (works across Streamlit-reloaded enum classes)
    rather than ``isinstance(..., PrecisionMode)`` (fails across reloads).
    """
    if mode is None:
        return PrecisionMode.SENSITIVE.value
    # Any Enum member (including stale PrecisionMode from a prior import)
    if isinstance(mode, Enum):
        return str(getattr(mode, "value", mode))
    if isinstance(mode, str):
        return _MODE_ALIASES.get(mode, mode)
    if type(mode).__name__ == "PrecisionConfig" and hasattr(mode, "mode"):
        return mode_value(getattr(mode, "mode"))
    if hasattr(mode, "mode") and not hasattr(mode, "value"):
        return mode_value(getattr(mode, "mode"))
    if hasattr(mode, "value"):
        return str(getattr(mode, "value"))
    # Last resort: string form often looks like "PrecisionMode.TS_DENOVO"
    text = str(mode)
    if "TargetScan de novo" in text:
        return "TargetScan de novo"
    for m in PrecisionMode:
        if m.value in text or m.name in text:
            return m.value
    raise TypeError(f"Cannot read precision mode value from {type(mode)!r}: {mode!r}")


def is_high_stringency(mode: object) -> bool:
    """True for Sequence-based (high stringency) / legacy Stringent."""
    return mode_value(mode) == PrecisionMode.STRINGENT.value


def is_low_stringency(mode: object) -> bool:
    return mode_value(mode) == PrecisionMode.SENSITIVE.value


def _coerce_mode(mode: PrecisionMode | str | object) -> PrecisionMode:
    """Normalize to *this* module's PrecisionMode by label string only.

    Never call ``PrecisionMode(x)`` — that fails across Streamlit hot-reloads when
    a stale enum class is missing newer members (e.g. TS_DENOVO).
    """
    val = mode_value(mode)
    for m in PrecisionMode:
        if m.value == val or m.name == val:
            return m
    known = ", ".join(m.value for m in PrecisionMode)
    raise ValueError(f"Unknown precision mode: {val!r}. Known: {known}")


@dataclass(frozen=True)
class PrecisionConfig:
    mode: PrecisionMode = PrecisionMode.SENSITIVE
    use_conservation: bool = True

    @staticmethod
    def from_mode(mode: PrecisionMode | str | "PrecisionConfig", **kwargs) -> "PrecisionConfig":
        # Streamlit hot-reload can produce a *different* PrecisionConfig class object;
        # always rebuild onto *this* module's classes.
        if type(mode).__name__ == "PrecisionConfig" and hasattr(mode, "mode"):
            return PrecisionConfig(
                mode=_coerce_mode(getattr(mode, "mode")),
                use_conservation=kwargs.get(
                    "use_conservation",
                    bool(getattr(mode, "use_conservation", True)),
                ),
            )
        return PrecisionConfig(mode=_coerce_mode(mode), **kwargs)

    @property
    def mode_label(self) -> str:
        return mode_value(self.mode)

    @staticmethod
    def ui_default() -> "PrecisionConfig":
        return PrecisionConfig(mode=PrecisionMode.SENSITIVE)

    @staticmethod
    def paper_default() -> "PrecisionConfig":
        return PrecisionConfig(mode=PrecisionMode.STRINGENT)

    @staticmethod
    def claim_default() -> "PrecisionConfig":
        return PrecisionConfig(mode=PrecisionMode.CONSENSUS)


def env_flag(name: str, default: bool = True) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    return v not in ("0", "false", "False", "no")


def apply_precision_filter(
    df: pd.DataFrame,
    cfg: PrecisionConfig,
    *,
    conserved_symbols: set[str] | None = None,
    is_unmodified_state: bool = True,
) -> pd.DataFrame:
    """Filter one state's target table. Does not touch deltas.

    Expected columns (subset OK): symbol, site_rank, is_conserved.
    """
    if df is None or df.empty:
        return df.iloc[0:0].copy() if df is not None else pd.DataFrame()

    cfg = PrecisionConfig.from_mode(cfg)
    out = df
    use_cons = cfg.use_conservation and env_flag("O8G_USE_CONSERVATION", True)

    mode_s = cfg.mode_label
    if mode_s == PrecisionMode.SENSITIVE.value:
        out = out[out["site_rank"] >= 3]
    elif mode_s == PrecisionMode.STRINGENT.value:
        out = out[out["site_rank"] >= 4]
    elif mode_s == PrecisionMode.TS_DENOVO.value:
        if "site_rank" in out.columns:
            out = out[out["site_rank"] >= 3]
    elif mode_s in (PrecisionMode.CONSENSUS.value, PrecisionMode.TARGETSCAN.value):
        out = out[out["site_rank"] >= 3]
        if use_cons and is_unmodified_state:
            if mode_s == PrecisionMode.CONSENSUS.value and "is_conserved" in out.columns:
                out = out[out["is_conserved"].astype(bool)]
            elif conserved_symbols is not None:
                syms = {str(s).upper() for s in conserved_symbols}
                out = out[out["symbol"].astype(str).str.upper().isin(syms)]
    else:
        raise ValueError(f"Unknown precision mode: {mode_s!r}")

    return out.reset_index(drop=True)


def partition_after_filter(
    unmod: pd.DataFrame,
    oxid: pd.DataFrame,
    cfg: PrecisionConfig,
    *,
    conserved_symbols: set[str] | None = None,
) -> dict[str, set[str]]:
    """Apply filters to each state; deltas = setdiff after filtering.

    Consensus special-case (documented): conservation shrinks the *unmodified
    baseline* for claims, but **gained** uses the same rank gate on both states
    without conservation on the oxidized side *and* subtracts the full
    rank-filtered unmodified set — so gains remain true o8G retargeting events,
    not artifacts of dropping non-conserved unmodified targets.
    """
    cfg = PrecisionConfig.from_mode(cfg)
    mode_s = cfg.mode_label
    # Consensus / TargetScan catalog: special gain/loss anchoring
    # TargetScan de novo + sequence-based modes: plain setdiff after filters
    if mode_s not in (PrecisionMode.CONSENSUS.value, PrecisionMode.TARGETSCAN.value):
        u = apply_precision_filter(
            unmod, cfg, conserved_symbols=conserved_symbols, is_unmodified_state=True
        )
        o = apply_precision_filter(
            oxid, cfg, conserved_symbols=conserved_symbols, is_unmodified_state=False
        )
        su = set(u["symbol"].astype(str).str.upper()) if len(u) else set()
        so = set(o["symbol"].astype(str).str.upper()) if len(o) else set()
        return {
            "unmod": su,
            "oxid": so,
            "shared": su & so,
            "lost": su - so,
            "gained": so - su,
        }

    # Consensus / TargetScan: anchor shrinks unmodified baseline; gains use rank-only
    u_rank = apply_precision_filter(
        unmod,
        PrecisionConfig(mode=PrecisionMode.SENSITIVE),
        is_unmodified_state=True,
    )
    o_rank = apply_precision_filter(
        oxid,
        PrecisionConfig(mode=PrecisionMode.SENSITIVE),
        is_unmodified_state=False,
    )
    su_rank = set(u_rank["symbol"].astype(str).str.upper()) if len(u_rank) else set()
    so_rank = set(o_rank["symbol"].astype(str).str.upper()) if len(o_rank) else set()
    u_anch = apply_precision_filter(
        unmod, cfg, conserved_symbols=conserved_symbols, is_unmodified_state=True
    )
    su_anch = set(u_anch["symbol"].astype(str).str.upper()) if len(u_anch) else set()
    return {
        "unmod": su_anch,
        "oxid": so_rank,
        "shared": su_anch & so_rank,
        "lost": su_anch - so_rank,
        "gained": so_rank - su_rank,
    }


def assert_retargeting_signal(
    parts: dict[str, set[str]],
    *,
    min_gained: int = 50,
    min_lost: int = 50,
    min_shared: int = 50,
    label: str = "miR-1 o8G@7",
) -> None:
    """Regression: oxidation still yields substantial partitions after filtering."""
    g, l, s = len(parts["gained"]), len(parts["lost"]), len(parts["shared"])
    if g < min_gained or l < min_lost or s < min_shared:
        raise AssertionError(
            f"{label} retargeting collapsed under filter: "
            f"gained={g}, lost={l}, shared={s} "
            f"(need ≥{min_gained}/{min_lost}/{min_shared})"
        )
