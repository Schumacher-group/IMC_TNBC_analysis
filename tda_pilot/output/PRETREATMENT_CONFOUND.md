# Decisive confound: the responder signal is treatment-stage sampling, not topology

**Date:** 2026-05-30. Definition `CD8_primary`, feature `im_dim1_total_persistence`.
Cliff's δ = R vs NR (−ve = responders lower).

## Pre-treatment composition is flat
Per-ROI CD8:tumour ratio, **pre-treatment only** (`composition_ratio.png`):
- per-ROI: median NR 0.77 vs R 0.64, Cliff δ = −0.077, MWU p = 0.15
- **per-patient (honest unit): median NR 0.70 vs R 0.65, δ = +0.008, p = 0.97 — no difference.**

## The image-persistence signal is a pre/post artifact
| subset | nROI (NR/R) | raw δ | size-adj δ (ROI) | size-adj δ (patient) |
|---|---|--:|--:|--:|
| ALL (pre+post) | 797 (330/467) | −0.168 | −0.138 | −0.225 (p=0.13) |
| **PRE only** | 582 (160/422) | −0.061 | **+0.056** | +0.146 (p=0.34) |
| POST only | 177 (**170/7**) | −0.605 | +0.055 | +0.143 (p=0.91) |

**Mechanism:** post-treatment resections are ~96% non-responder (170 vs 7) — responders have
little/no residual tumour to image. Those residual-disease NR ROIs are tumour-imbalanced →
high image persistence → they inflate the NR average in the pooled cohort, manufacturing the
"responders lower" signal. On PRE-treatment biopsies the size-adjusted signal is null/slightly
reversed.

## Conclusion
- The −0.138 signal was **confounded by treatment stage**, not biology.
- Pre-treatment, neither **composition** (δ=+0.008, p=0.97) nor **size-adjusted topology**
  (δ=+0.146, p=0.34) distinguishes R from NR.
- **No-go** for TDA on the responder contrast. Honest statement: at baseline, tumour–CD8 spatial
  organisation does not distinguish responders from non-responders by composition or chromatic topology
  in this cohort.

Caveats: pre-treatment per-patient comparison is modestly powered (25 NR / 37 R patients) — rules out a
strong effect, not a tiny one. **Any pooled pre+post analysis (incl. cross-PCF) inherits the same
treatment-stage confound** and should be redone pre-treatment-only.

Script: `day1_composition_check.py` (+ inline pre/post decomposition).

---

## Addendum (2026-08-05) — corrections and scope of the above

**1. The subset rows do not sum to ALL.** 582 (PRE) + 177 (POST) = 759, not 797. The missing 38 ROIs
are `Sample_Type = mid-early` (23) and `mid-late` (15). **All 38 are responders.** They are included
in the pooled ALL row but silently dropped from the pre/post split. This does not change the
conclusion — it makes the pooled row's responder-weighting slightly worse than the split suggests —
but the table should not be read as a partition.

**2. "Null" here means "no large effect", quantified.** With 25 NR / 37 R patients, a two-sided
Mann–Whitney at α = 0.05 has ~80% power only at Cliff's |δ| ≈ **0.42**, and ~33% power at |δ| = 0.22.
The reported pre-treatment δ = +0.146 (p = 0.34) is entirely consistent with a real small-to-moderate
effect. State the floor alongside the no-go.

**3. The pre-treatment check covered one feature of one definition.** Only
`im_dim1_total_persistence` for `CD8_primary` was decomposed pre/post. Re-running the pre-treatment
contrast across all 8 pair definitions × 10 features × 4 analysis levels (320 comparisons) confirms
the conclusion — nothing reaches nominal significance, min p = 0.056 — so the no-go is better
supported than this note originally established. See `LENGTH_FEATURES.md`.

**4. But `im_dim1_total_persistence` was a poor choice of headline feature.** Median degree-1 image
bar lifetime in this cohort is ~0.33 µm; only ~3% of bars exceed the 5 µm threshold. Summing all bars
makes the statistic close to a rescaled bar count, i.e. it re-encodes the very cell-number and
class-imbalance confounds this note is about. The bar-*length* statistics (already computed, in
`output/stats/*.parquet`) are the scale-sensitive alternative and were never examined — see
`LENGTH_FEATURES.md` and, for the decisive test, `PERMUTATION_NULL.md`.

**5. The cross-PCF warning above does not apply — that pipeline already avoids the confound.**
The original note warned that "any pooled pre+post analysis (incl. cross-PCF) inherits the same
treatment-stage confound". Checked: `results_analysis/spatial_stats/export_for_spoox.py` builds
`conditions_pretreat.json` as responder vs non-responder **restricted to pre-treatment**, and
`conditions_nonresponder.json` as pre vs post **within non-responders only**. Both are correctly
stratified, and neither contains any ROI outside the 808 revision cohort. The warning stands for any
*future* pooled analysis, not for the existing SpOOx runs.

**6. Metadata bug: 11 ROIs lost their response label.** `Leap084a_1..11` parse to LEAP_ID `LEAP084A`,
which has no metadata row, so they were dropped from every response-stratified Day-1 analysis (this is
why ALL is 797 and not 808). They are the second tissue block of the LEAP084 core, whose other block is
named `Leap084_b_12..21` — the separator differs, nothing else. `conditions_pretreat.json` lists 433
pre-treatment responder ROIs against this table's 422, a difference of exactly these 11, confirming
they are pre-treatment responders.

**Applied 2026-08-05.** `roi_counts.py` carries the alias and `per_roi_counts.csv` has been
regenerated. The labelled cohort is now **808 ROIs / 63 patients**, and pre-treatment is
**593 ROIs — 433 Responder / 160 Non-Responder**, which matches `conditions_pretreat.json`
exactly (two independent paths through the metadata now agree). Patient counts are unchanged at
25 NR / 37 R, since LEAP084A joins existing patient 44. All Day-1 and Day-2 outputs in this
directory have been regenerated on the corrected cohort; effect sizes moved by ≤ 0.02.

**7. Neither ROI area nor cell density was ever controlled.** The Day-1 "size adjustment" regressed
on cell counts only. Convex-hull area spans 3.7e4–2.0e6 µm² (54×) across the cohort and pre-treatment
responder ROIs are significantly less dense (per-ROI Cliff δ = −0.28 on density, p < 1e-3). This cuts
both ways: geometry was uncontrolled when the signal was claimed *and* when it was dismissed.
`roi_geometry.py` now provides these covariates.
