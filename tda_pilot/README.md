# TDA pilot — chromatic persistence on (tumour, CD8)

Feasibility pilot for running M2S2 chromatic-TDA (`chalc` six-pack persistence)
on the TNBC IMC cohort, to distinguish CD8 infiltration from CD8 rim/exclusion
around tumour — a distinction squidpy neighbourhood enrichment cannot make.

See **`REPORT.md`** for the Day-0 smoke-test results and the go/no-go decision.

## What is tracked here

Small, hand-written, or reference files only. The heavy, regeneratable artefacts
(the upstream `M2S2_demo` clone, the Python venv, per-ROI point clouds, persistence
diagrams, and the Parquet stats) are git-ignored — see `.gitignore`. Everything
below can be regenerated from the cell table with the commands in this README.

| File | Purpose |
|---|---|
| `cohort.py` | Revision-cohort definition (N=808). Single source of truth for exclusions. |
| `adapter.py` | Extract one ROI's (tumour, CD8) point cloud → `data/<fov>.csv` (+ overlay PNG). |
| `run_smoke_test.py` | Time one chromatic six-pack; save raw diagrams; sanity-check output. |
| `generate_stats.py` | Run the M2S2 `generate_stats` pipeline → partitioned Parquet stats. |
| `roi_counts.py` | One pass over the 8 GB cell table → `per_roi_counts.csv`. |
| `install_test.py` | Synthetic end-to-end install check (no network). |
| `per_roi_counts.csv` | Per-ROI tumour/CD8 counts + responder status + `cohort_member`. |
| `data/Leap066_11_overlay.png` | Sanity plot of the smoke-test ROI. |

## Environment setup

Requires **Python 3.13** (the current `M2S2_demo` pins `requires-python >= 3.13`
and `chalc ~= 14.0.0`; older Pythons can only install `chalc <= 6`). We use
[`uv`](https://docs.astral.sh/uv/), which auto-provisions the right interpreter.

```bash
cd tda_pilot

# 1. Clone the upstream demo repo (NOT vendored here; it is a separate project).
git clone https://github.com/abhinavnatarajan/M2S2_demo

# 2. Create the env + install all deps (chalc 14, gudhi, scikit-tda, pyarrow, ...).
#    This provisions Python 3.13 automatically and creates M2S2_demo/.venv .
cd M2S2_demo && uv sync && cd ..
```

The scripts call the interpreter directly as `M2S2_demo/.venv/bin/python` and add
`M2S2_demo/` to `sys.path`, so no `activate` is needed.

### macOS note (required)

`M2S2_demo/utils/parallel.py` passes a local closure as the worker-pool
initialiser. That works under Linux's default `fork` start method but fails under
macOS / Python 3.13's default `spawn` (`Can't get local object ...new_init`). Our
launchers therefore force `fork` (`multiprocessing.set_start_method("fork")`) — this
is in our scripts only; **`M2S2_demo` and `chalc` are never modified.** If you adapt
other M2S2 entry points on macOS, do the same.

## Reproduce the smoke test

```bash
# from the repository root
P=tda_pilot/M2S2_demo/.venv/bin/python

$P tda_pilot/roi_counts.py                       # -> per_roi_counts.csv (≈22 s, scans 8 GB)
$P tda_pilot/adapter.py --fov Leap066_11         # -> data/Leap066_11.csv + overlay PNG
$P tda_pilot/run_smoke_test.py Leap066_11        # timing + raw diagrams + sanity report
$P tda_pilot/generate_stats.py --num-workers 1   # -> output/stats/ (partitioned Parquet)

$P tda_pilot/install_test.py                      # optional: synthetic end-to-end check
```

## Cohort

All pilot analyses use the **revision cohort (N = 808 ROIs, 63 patients)**, defined
in `cohort.py`: the raw cell table (829 ROIs) minus patients `LEAP149`/`LEAP150`
(→ 813, the originally submitted manuscript cohort) minus five collaborator-flagged
ROIs (duplicates / broken / palette → 808). `adapter.py` and `generate_stats.py`
both apply `cohort.is_cohort_member` so the filter is never duplicated ad hoc.

## Data inputs (not in this repo)

- `../CellTable_CleanCohort/updated_cell_table.csv` — per-cell table (8.2 GB).
- `../CellTable_CleanCohort/CleanCohort_Metadata.csv` — per-`LEAP_ID` response metadata.
