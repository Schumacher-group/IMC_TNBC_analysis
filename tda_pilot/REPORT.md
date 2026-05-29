# Day 0 smoke test — Chromatic TDA on (tumour, CD8) — REPORT

**Date:** 2026-05-29
**Machine:** macOS (Darwin 24.6.0, arm64), 10 logical CPUs
**Decision in one line:** 🟢 **GREEN — full pilot is comfortably feasible.** A single
(tumour, CD8) chromatic six-pack on a median ROI takes **~0.9 s** (single-threaded);
the whole revision cohort for 5 pairs projects to **~1.1–2.2 h single-threaded, ~10–15 min on 10 cores.**

---

## Headline numbers

| Quantity | Value |
|---|---|
| Cohort | **808 ROIs, 63 patients** (revision cohort — see "Cohort" below) |
| Representative ROI | `Leap066_11` (pre-treatment **responder**, cohort-**median** pair size) |
| Points in ROI | 2008 (Tumour 1302, CD8 706) |
| **Per-pair compute (filtration + six-pack), 1 thread** | **0.874 s** |
| + persistent-statistics vectorisation | 0.137 s → **~1.0 s per pair end-to-end** |
| Peak process RSS (incl. chalc native) | **~338 MB** |
| Full pipeline wall-clock for the ROI (1 worker, 2 singletons + pair + I/O + startup) | 3.1 s |
| Raw six-pack diagrams on disk | 427 KB (`.h5`) + 234 KB (matrices `.pkl`) per ROI·pair |

### Projected full-cohort runtime
Formula: `runtime_per_pair_per_ROI × N_ROI × N_pairs`, with **N_ROI = 808** (revision
cohort) and `N_pairs = 5` (tumour–CD8 + four controls: tumour–CD4, tumour–B,
tumour–macrophage, tumour–CAF).

- **Single-threaded, median-ROI basis (~1.0 s/pair):** 808 × 5 × 1.0 s ≈ **67 min ≈ 1.1 h**
- **Conservative** (×2 for the right-hand size tail + any super-linear scaling): **≈ 2.2 h single-threaded**
- **With 10-core parallelism** (the pipeline parallelises across label-combinations/ROIs):
  **~7–15 min** for the whole cohort.

All comfortably under the 30 s/ROI "green" threshold (per-ROI per-pair is ~1 s).

> Caveat on the projection: timing is from **one median ROI** (the plan permits only one).
> Pair sizes across the cohort range 73 → 8958 points (median 2008, mean 2359, p90 4543).
> Runtime tracks point count; the median basis slightly under-projects the heavy tail,
> hence the conservative ×2 line above. A 2-point scaling check (median + p90 ROI) would
> tighten this but needs a second ROI run (out of Day-0 scope).

---

## Did the install need non-trivial fixes?

Yes — three, **all confined to `tda_pilot/`; chalc and M2S2_demo core were not modified.**

1. **Python 3.13 required (plan assumed 3.11/3.12).** The current `M2S2_demo`
   pins `requires-python >= 3.13` and `chalc~=14.0.0`; on this machine's Python 3.11
   the newest installable `chalc` is only 6.0.1. Resolved using the repo's own `uv`
   tooling (`uv sync`, which auto-provisioned **Python 3.13.5**). Installed: `chalc 14.0.0`,
   `gudhi 3.11.0`, `scikit-tda 1.1.1`, `pyarrow 20.0.0`.
2. **macOS multiprocessing (worth flagging to the author).** `utils/parallel.py`'s
   `PoolWithLogger` passes a **local closure** (`new_init`) as the worker-pool
   initialiser. This works under Linux's default `fork` start method but crashes under
   macOS/Python 3.13's default `spawn` (`Can't get local object ...new_init`). Worked
   around by forcing `multiprocessing.set_start_method("fork")` **in our launcher**.
   The author may want to make this portable (e.g. a module-level initialiser) for macOS users.
3. **`SixPack.save()` takes an h5py group, not a path string** (minor; our usage).
   Fixed by opening an `h5py.File` and passing the handle. Native `.h5` round-trip works.

## Did the M2S2 example/test script run cleanly?

The bundled examples (`colorectal_cancer_dataset_analysis`, `abm_dataset_analysis`)
require downloading the authors' datasets over the network (not our data). Instead, the
install was verified **end-to-end** by pushing a synthetic 2-colour point cloud through the
**exact** `make_cmdline_parser` + `generate_stats` machinery (chalc filtration → six-pack →
persistent stats → partitioned Parquet). It completed cleanly (exit 0, full six-pack with
non-empty kernel). This is the same code path the cohort run uses. See `install_test.py`.

## Does the output look biologically sensible?

Yes. On `Leap066_11` the six-pack is non-degenerate and shows clear loop structure:

```
        dim0   dim1
domain  2008   2089
codomn  2008   2390
kernel   320    450
cokern     0    730     (dim-0 cokernel is always empty by construction)
image   2008   1712
relativ    1    857
```

- **Degree-1 features are abundant** (e.g. `dom1=2089`, `rel1=857`) → real loop/void
  structure in the tumour+CD8 arrangement, not noise.
- **The inclusion kernel is non-trivial** (`ker0=320`, `ker1=450`) → there is CD8 topology
  that is "filled in" once tumour cells are added, i.e. genuine tumour–CD8 spatial coupling —
  exactly the signal that distinguishes infiltration from a rim/exclusion pattern.
- Stats feature vector: **456 values (6 diagrams × 2 dims × 38 stats), all finite.**

---

## Decisions used (confirmed with user)

- **Tumour aggregate** = `Cancer cell` + `B7H4 Cancer cell` → label `Tumour`.
- **CD8 (primary)** = `CD8 T cell` + `Memory CD8 T cell` → label `CD8` (NK/CD8 excluded).
  Two sensitivity definitions for the full pilot (CD8 alone; +NK/CD8) are just alternate
  `--labels-include` sets — no extra machinery.
- **No `qc_pass` filter** (use all cells).
- **Coordinates**: `centroid-1`→x, `centroid-0`→y, treated as micrometres (Hyperion = 1 px/µm).

### Cohort

This pilot uses the **revision cohort: N = 808 ROIs across 63 patients**, defined in
`cohort.py` and applied via `cohort.is_cohort_member` in both `adapter.py` and
`generate_stats.py`. It is the raw cell table minus:
- **patients `LEAP149` (4 ROIs) and `LEAP150` (12 ROIs)** — the originally submitted
  manuscript already excluded these (→ 813-ROI manuscript cohort); and
- **five collaborator-flagged ROIs** (→ 808): `Leap005_2_1`, `Leap005_2_2` (duplicate
  acquisitions of `Leap005_1`/`Leap005_2`), `Leap010_7`, `Leap094_7` (broken, retaken as
  `Leap010_8`/`Leap094_10`), and `Leap095_13` (non-tissue palette region).

## Cohort facts found (one full scan of the 8.2 GB cell table, 4.24 M cells, ~22 s)

- **Raw cell table: 829 ROIs.**
- **Manuscript cohort: 813 ROIs / 63 patients** (829 − LEAP149/LEAP150 = 829 − 16).
- **Revision cohort: 808 ROIs / 63 patients** (813 − 5 collaborator-flagged ROIs).
- **Cell counts are not a constraint:** every ROI (full table) has ≥20 of each species;
  807/829 have ≥50 of both tumour and CD8-primary. **0 ROIs** fall below 20 CD8 under any of
  the three CD8 definitions.
- Pair size (tumour+CD8) distribution: min 73, median 2008, mean 2359, p90 4543, max 8958.
- Per-ROI counts (with `cohort_member` + `Patient_ID`) saved to `tda_pilot/per_roi_counts.csv`
  (full 829-row table retained).

### Note on metadata coverage (the "unmatched" ROIs)
11 ROIs in the cell table carry no responder label after the metadata join. They all
belong to `LEAP084A` — a `LEAP_ID` present in the cell table (alongside `Leap084_b_*`
acquisitions) but **absent from `CleanCohort_Metadata.csv`**, i.e. a `LEAP084` / `LEAP084A`
naming mismatch, **not** the excluded LEAP149/LEAP150. These 11 ROIs *are* inside the 808
revision cohort but lack a `Response` value; the metadata key should be reconciled before
any response-stratified analysis. (Does not affect the Day-0 timing.)

## What `cell_table_bothneighbourhoods.csv` is
No notebook references it. Given the sibling files `Neighbourhoods_PreSub.csv` /
`Neighbourhoods_NonResSub.csv`, it is almost certainly `updated_cell_table` plus two extra
columns assigning each cell to a spatial-neighbourhood cluster (one per analysis subset).
Not needed for TDA.

---

## Storage note for the full pilot
Saving raw six-pack diagrams (per user request) costs ~0.4 MB (`.h5`) per ROI·pair →
≈ **1.7 GB** for 808 ROIs × 5 pairs. The persistent-statistics Parquet is tiny by comparison.

## Reproduce
```bash
cd tda_pilot
M2S2_demo/.venv/bin/python adapter.py --fov Leap066_11        # extract ROI -> data/Leap066_11.csv (+ overlay PNG)
M2S2_demo/.venv/bin/python run_smoke_test.py Leap066_11       # timing + raw diagrams + sanity
M2S2_demo/.venv/bin/python generate_stats.py --num-workers 1  # full pipeline -> output/stats (Parquet)
```
(Environment setup — `uv sync`, Python 3.13, the macOS fork note, and the `M2S2_demo`
clone command — is in `README.md`.)

## Files produced
- `data/Leap066_11.csv`, `data/Leap066_11_overlay.png` — ROI point cloud + sanity plot
- `output/diagrams/Leap066_11.h5` (native six-pack), `Leap066_11_matrices.pkl` (birth/death arrays)
- `output/stats/...` — partitioned persistent-statistics Parquet
- `per_roi_counts.csv` — per-ROI tumour/CD8 counts + responder status + `cohort_member` (whole cell table)
