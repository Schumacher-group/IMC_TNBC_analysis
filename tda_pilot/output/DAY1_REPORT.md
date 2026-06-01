# Day 1 — Full chromatic TDA pilot: results

**Cohort:** 808 ROIs / 63 patients (revision cohort).
**Output:** clean per-(ROI, definition) feature table for Day-2 modelling. **No model fit** (per plan).

## Phase 2 — full run
- 7 pair definitions × 808 ROIs = **5,656 ROI-pair six-packs**, ROI-level parallelism (10 workers).
- **Total compute ≈ 36 min** wall: data build 4.5 min + run 15.9 min + validation <1 min + feature extraction ~15 min.
- **Zero failures.** Per-definition (completed / degenerate / failed), degenerate = a control species had <3 cells in that ROI:

| definition | completed | degenerate | failed | run time |
|---|---|---|---|---|
| CD8_primary | 808 | 0 | 0 | 189 s (incl. .h5 diagrams) |
| CD8_strict | 808 | 0 | 0 | 172 s |
| CD8_extended | 808 | 0 | 0 | 179 s |
| CD4 | 803 | 5 | 0 | 102 s |
| Bcell | 808 | 0 | 0 | 82 s |
| Macrophage | 806 | 2 | 0 | 102 s |
| Fibroblast | 807 | 1 | 0 | 127 s |

## Phase 3 — validation (see VALIDATION.md) — **ALLPASS**
- Completeness: completed+degenerate+failed = 808 for every definition; roi_ids unique.
- Finiteness: all 456 stat columns finite for every ok row.
- Non-degeneracy: kernel & image have >0 bars in dim 0 **and** 1 for **100%** of ok ROIs (≥90% required).
- Raw-diagram spot check (responder Leap066_11 vs non-responder Leap001_10): sensible degree-1 diagrams.

## Phase 4 — feature table
`output/per_roi_summary.parquet` — 5,656 rows × {roi_id, patient_id, response, pair_definition,
n_tumour, n_other, ker/im/cok degree-1 summaries}. 5,648 ok, 8 degenerate (NaN features), 0 failed.
Features computed exactly from degree-1 diagrams (5 µm lifetime threshold for n_features).

## Decision checkpoint — qualitative read (NOT a result)

**Focal vibes check** — `ker_dim1_total_persistence`, CD8_primary, responder vs non-responder
(`output/vibes_ker_dim1.png`): **flat.** Per-patient medians NR 1185 vs R 1179 (ratio 1.00);
per-ROI NR 1060 vs R 955 — heavy overlap, no visual separation.

**Wider univariate descriptive scan** (median R/NR ratio, CD8_primary; size-confounded, no model):

| feature | NR med | R med | R/NR |
|---|---|---|---|
| ker_dim1_total_persistence | 1060 | 955 | 0.90 |
| ker_dim1_n_features | 48 | 46 | 0.96 |
| im_dim1_total_persistence | 1843 | 1403 | **0.76** |
| im_dim1_n_features | 55 | 49 | 0.89 |
| cok_dim1_total_persistence | 862 | 835 | 0.97 |
| cok_dim1_n_features | 26 | 28 | 1.08 |
| ker/im/cok per-cell (size-normalised) | — | — | 0.94–1.10 |

**Read:** No single topological feature separates the groups univariately. The largest (still modest)
signal is in the **image** diagram — responders carry lower image total-persistence / fewer image loops.
The cokernel "encirclement/rim" signature is flat. This is a **weak/unencouraging** univariate result.

**Caveats (why this isn't a verdict):**
1. These are univariate medians on size-confounded features; the planned Day-2 test is a *multivariate
   hierarchical* model with partial pooling across 63 patients, which can surface joint/weak signal.
2. Signal could live in feature combinations or specific contrasts (e.g. CD8_primary vs CD8_strict/extended),
   not examined here.
3. But: a strong, obvious topological signal is **not** present in the simple read.

**Per plan, paused for go/no-go before Day-2 modelling.**
