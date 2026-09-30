# Fix-up prompt for `Marciano-Lab/o8g-mirna-explorer`

Paste everything below the line into a coding agent working in a clone of
`https://github.com/Marciano-Lab/o8g-mirna-explorer` (branch `main`). It is ordered
so the two correctness bugs land first; the feature ports are optional and
independent of each other.

Three of these are live defects in the deployed copy at `oxomir.samnti.com`:
Task 1 and Task 2 are crashes or wrong output, and Task 6 inflates every
enrichment number the site reports. Fix those three even if you skip everything
else. The feature ports (Task 4) are optional.

---

You are working in a clone of the **o8G-miRNA Retargeting Explorer** — a
Streamlit app that predicts how 8-oxoguanine oxidation of microRNA seed
guanines rewires human mRNA target lists. Seed = mature positions 2–8. An
unoxidized seed G pairs C; an 8-oxoG pairs A (Hoogsteen). A seed with *k*
guanines therefore has 2^k target repertoires. Targets are precomputed into
`o8g_targets.db`; `o8g_precision.py` implements a Sensitive / Stringent /
Consensus filter ladder whose governing rule is **never filter the delta
itself** — filter each state's target list, then partition into
lost/gained/shared.

Do the following. After each task, run the regression check at the bottom.

## Task 1 (BUG, ship first) — Consensus mode silently reports garbage

**Symptom.** With `paper/data/` absent — which is the case in every fresh
clone, since `paper/` is gitignored — selecting precision mode **Consensus**
produces `unmod=0, oxid=3909, lost=0, gained=3909, shared=0` for
hsa-miR-1-3p. Every oxidized target is reported as newly gained and the
unmodified baseline vanishes. The sidebar simultaneously tells the user
"Paper / claims: use Stringent or Consensus."

**Cause.** `conservation.ConservedIndex.conserved_symbols_for_mirna` raises
`FileNotFoundError` when `paper/data/Conserved_Family_Info.txt` is missing.
Both `TargetDB.targets_filtered` and `TargetDB.retarget_partition` in
`o8g_db.py` wrap that call in `except Exception: conserved_symbols = set()`.
Consensus intersects the unmodified baseline with the conserved set, so an
empty set zeroes the baseline. The failure is invisible.

**Fix, in `o8g_db.py`:**

1. Define, near `RANK_SITE`:

   ```python
   class ConservationUnavailable(RuntimeError):
       """Raised when TargetScan conservation data cannot be loaded.

       Never substitute an empty conserved set: Consensus intersects the
       unmodified baseline with it, so an empty set silently collapses the
       baseline to zero and reports every oxidized target as gained.
       """
   ```

2. Add one helper used by both call sites:

   ```python
   def _conserved_for(self, seed, mirna):
       if not mirna:
           raise ConservationUnavailable(
               "Consensus mode needs the miRNA name to look up TargetScan "
               "conserved families; none was supplied.")
       try:
           syms = self._conserved_index().conserved_symbols_for_mirna(mirna)
       except FileNotFoundError as e:
           raise ConservationUnavailable(
               "Consensus mode requires paper/data/Conserved_Family_Info.txt "
               "(TargetScanHuman 8.0), which is not installed. Use Sensitive "
               "or Stringent, or add the TargetScan release files.") from e
       except Exception as e:
           raise ConservationUnavailable(
               f"Consensus conservation lookup failed ({type(e).__name__}): {e}") from e
       if not syms:
           raise ConservationUnavailable(
               f"TargetScan returned no conserved families for {mirna}; "
               "refusing to run Consensus on an empty baseline.")
       return syms
   ```

3. Replace **both** `except Exception: conserved_symbols = set()` blocks (in
   `targets_filtered` and in `retarget_partition`) with a call to
   `self._conserved_for(seed, mirna)`. In `retarget_partition` the mode check
   must not be gated on `and mirna` — a missing miRNA name currently skips the
   lookup silently; it should now raise.

**Fix, in `app.py`:** decide up front, not inside the render path.

```python
from o8g_db import TargetDB, ConservationUnavailable

@st.cache_data(show_spinner=False)
def _conservation_status(_seed, _mirna, _ver=_CACHE_VER):
    """Return "" if Consensus is usable, else the reason it is not."""
    try:
        db.retarget_partition(_seed, "none", "none",
                              PrecisionConfig.from_mode("Consensus"), mirna=_mirna)
        return ""
    except ConservationUnavailable as e:
        return str(e)

if mode_name == PrecisionMode.CONSENSUS.value:
    problem = _conservation_status(seed, mirna)
    if problem:
        precision_cfg = PrecisionConfig.from_mode("Sensitive")
        effective_mode = "Sensitive (Consensus unavailable)"
        st.sidebar.error("Consensus unavailable — showing Sensitive. " + problem)
```

and immediately after `st.title(...)`, repeat it as a main-pane `st.error`.

Two things matter here. The banner must be **persistent**, not shown once
behind a `st.session_state` dedup flag — a once-only banner disappears on the
next rerun and leaves Sensitive numbers sitting under a "Consensus" label.
And `effective_mode` (not `mode_name`) must be what captions, figure titles,
and CSV download headers report.

Route every target query through a single wrapper so no call site can bypass
the ladder:

```python
def filtered_targets(label, **kw):
    return db.targets_filtered(seed, label, precision_cfg, mirna=mirna, **kw)
```

Then rewrite all direct `db.targets_filtered(` call sites to use it. Verify
zero remain: `grep -n 'db\.targets_filtered' app.py` should match only the
wrapper body.

**Compare mode by value string, never by enum identity** — `mode_name ==
PrecisionMode.CONSENSUS.value`. Reloaded modules produce distinct enum objects
and identity comparison fails silently.

## Task 2 (BUG) — External DB comparison crashes with exactly two sets selected

In the External DB comparison view, the master-list slider is:

```python
min_dbs = st.slider("gene must appear in at least this many selected DBs",
                    min_value=2, max_value=max(2, n_sets), value=2)
```

When `n_sets == 2` — two reference files installed, or the user deselects down
to two — `min_value == max_value` and Streamlit raises. Guard it:

```python
if n_sets > 2:
    min_dbs = st.slider(..., min_value=2, max_value=n_sets, value=2)
else:
    min_dbs = 2
    st.caption("Master list requires agreement of both selected databases.")
```

## Task 3 — `scripts/benchmark_refsets.py` is a dangling symlink

It points at `../paper/scripts/benchmark_refsets.py`, and `paper/` is
gitignored, so the file is broken in every clone (and breaks `ls -L`, tarball
builds, and CI checkouts). Replace it with a small stub whose docstring says
where the real driver lives and that the reusable benchmark logic is importable
from `o8g_compare.py`. `sys.exit(2)` when run.

## Task 4 (feature) — five additional views

Add `o8g_sections.py` with a `SectionContext` dataclass and one `render_*(ctx)`
per view. The context carries `db, mirna, info, seed, gpos, state_labels,
library, precision_mode, strong_set, universe, matched_background,
external_refs`.

The design constraint that makes this safe: **sections never query the DB
directly.** Every target set arrives via `ctx.strong_set(state_label)`, which
routes through the precision ladder. Sections therefore inherit
Sensitive/Stringent/Consensus without knowing the filter rules, and no section
can compute a delta before filtering.

| View | Module | Content |
|---|---|---|
| Overlap (Venn/UpSet) | `o8g_venn.py` | Venn ≤3 sets, auto-switch to UpSet above; per-region gene table + CSV |
| Transcription factors | `o8g_tf.py` | TFs as direct targets (gained/lost/retained); optional one-hop regulon amplification behind a checkbox |
| Loss of function | `o8g_lof.py` | miRDB-anchored baseline with a score-cutoff slider; per-state attrition plot; vulnerability summary |
| Statistics | `o8g_stats.py` | Length-matched ORA or preranked GSEA, behind a "Run (resampling is slow)" checkbox |
| Antagomir design | `o8g_anti.py` + `o8g_energy.py` | Oligo design, discriminability, feasibility verdict, collateral miRNA scan |

Extend the router radio (keep `key="main_section"`; **keep the radio — do not
convert to `st.tabs`**, because Streamlit executes every tab body on every
rerun) to ten options and dispatch the five new labels through a dict.

Three caveats must appear in the UI, not only in docs:

- The oxidized-state ViennaRNA ΔG uses a **G→U proxy** (ViennaRNA has no o8G:A
  Hoogsteen parameter). `dG_normal` and collateral tables are trustworthy;
  `dG_oxo`/`ddG`/fold-preference are directionally right with placeholder
  magnitudes. `st.warning`.
- **Position-2-only oxidation states are poor antagomir targets**: median fold
  preference ≈8×, ~64% of designs below 10×, versus ≈150× better for interior
  positions. Dedicated `st.error` when such a state is selected.
- An antagomir **is not AGO-loaded and has no seed** — the view scores state
  discrimination and must never present an mRNA target list derived from the
  oligo.

Also add gene-set files `ChEA_2022.gmt`,
`ENCODE_TF_ChIP-seq_2015.gmt`, `TRRUST_Transcription_Factors_2019.gmt`,
`TF_Perturbations_Followed_by_Expression.gmt` to `genesets/`, and register them
in `o8g_enrich.py` with `TF_LIBRARIES = ("ChEA_TF","ENCODE_TF","TRRUST_TF","TF_Perturb")`.

`requirements.txt` additions: `requests>=2.31`, `lxml>=5.0`, `gseapy>=1.1`.

## Task 5 — repository hygiene

- `o8g_targetscan.db` is 102 MB and **exceeds GitHub's 100 MB hard limit**;
  `o8g_targets.db` at 87 MB is over the 50 MB warning. Do not commit either.
  Add them plus `o8g_confident.db`, `mirdb_ref.parquet`, and
  `mirdb_custom_cache.db` to `.gitignore`, and have the app degrade gracefully
  when each is absent (the views that need them should say which file is
  missing, not raise).
- Document in `README.md`: the ten-section router, that Consensus needs
  `paper/data/Conserved_Family_Info.txt` (TargetScanHuman 8.0) and refuses
  rather than falling back silently, and which optional data file each view
  requires.

## Task 6 (STATISTICS, affects every reported enrichment)

`o8g_enrich.enrich(query_genes, library, background=None, ...)` takes
`background` as an **integer population size**, and when it is `None` uses the
union of all genes in the library. Both behaviours are wrong for this
application, and the second is what the deployed app currently does.

Two distinct errors:

1. Genes with no 3'UTR in our index are in the library union but could never
   have been predicted as targets. Counting them in the population `N` inflates
   significance.
2. Query genes that no term in the library annotates can never be a success,
   but they are counted in the draw count `q`. This also inflates significance.

Fix: let `background` accept an iterable of symbols (keep the int and `None`
paths for compatibility), and restrict **both** the population and the query to
`set(background) & (union of library gene sets)`. Pool sizes for our
19,133-symbol 3'UTR universe come out as GO_BP 14,280, Reactome 10,230,
KEGG 7,607, Hallmark 4,353.

This is not cosmetic. For hsa-miR-1-3p `none` vs `o8G@7` the count of terms at
q < 0.05 changes by more than an order of magnitude in several layers.

**Also add a control that must be displayed next to every differential
enrichment.** Enriching a lost/gained set against a genome-wide background
answers "what does this miRNA target", not "what did oxidation change". The
correct background for the differential question is the miRNA's own target set
— `union(unmodified, oxidized)` for the layer in question.
`o8g_lof.enrich_within_baseline` already implements this shape; generalise it to
an explicit pool and call it alongside the genome-background result. Measured
on our side:

| miRNA | layer | sig. vs genome | sig. within pool |
|---|---|---|---|
| miR-1-3p | miRDB-anchored | 37 | 0 |
| miR-124-3p | miRDB-anchored | 130 | 0 |
| miR-1-3p | confident ≥2.0 | 37 | 0 |
| miR-124-3p | confident ≥2.0 | 49 | 0 |
| miR-1-3p | confident ≥3.0 | 19 | 1 |
| miR-124-3p | confident ≥3.0 | 33 | 7 |

A genome-background enrichment number must never be presented on its own.

**Do not make miRDB-anchored targets the default layer.** The anchored oxidized
set is `baseline ∩ oxidized_strong`, a subset of the baseline by construction,
so `gained` is identically 0 for every miRNA. miRDB catalogs only unmodified
matures, so no external anchor for gained targets exists. Anchoring belongs in
the Loss of function view only, with `gained = 0` labelled as structural rather
than measured. See `ANCHORING_FINDINGS.md` for the full argument and the
permutation null.

## Regression check

`git` and a display are not required; use Streamlit's headless harness. This
must print `fails=0`:

```python
import sys, os; sys.path.insert(0, os.getcwd())
from streamlit.testing.v1 import AppTest
SECTIONS = ["Single state — targets","Compare two states","All states",
            "Overlap (Venn/UpSet)","Transcription factors","Loss of function",
            "Statistics","Antagomir design","Gene → miRNA/oxomiR",
            "External DB comparison"]
fails = []
for mode in ["Sensitive","Stringent","Consensus"]:
    for sec in SECTIONS:
        at = AppTest.from_file("app.py", default_timeout=300).run()
        for r in at.radio:
            if r.options and mode in r.options:
                r.set_value(mode); break
        at.run()
        at.radio(key="main_section").set_value(sec); at.run()
        at.run()   # second rerun: the Consensus banner must still be there
        if at.exception:
            fails.append((mode, sec, str(at.exception[0].value)[:120]))
print("fails=%d" % len(fails)); [print(f) for f in fails]
```

Construct a fresh `AppTest` per combination — reusing one across two
`set_value` calls gives spurious failures.

Expected values for hsa-miR-1-3p (seed `GGAATGT`, G at seed positions 2, 3, 7;
8 states), comparing `none` against `o8G@7`:

| Mode | unmod | oxid | lost | gained | shared |
|---|---|---|---|---|---|
| Sensitive (rank ≥ 3) | 3317 | 3909 | 2251 | 2843 | 1066 |
| Stringent (8mer only) | 1355 | 1693 | 1146 | 1484 | 209 |

HDAC4 (a canonical miR-1 target) must be present at `none` and absent at
`o8G@7`. Consensus without the TargetScan files must raise
`ConservationUnavailable` from `o8g_db`, and the app must show the fallback
banner — it must **not** report `unmod=0`.

Clear `__pycache__` before each test run; stale `.pyc` files mask edits.
