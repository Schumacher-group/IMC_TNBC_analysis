# Specification sensitivity: directional statistic, pre-treatment

Endpoint: per-patient median of `ker1-avg_length`, CD8-directional inclusion, 593 ROIs / 25 NR + 37 R patients. Cliff's δ is R vs NR (+ve = responders higher).

## 1. Regression adjustment — answer depends on the covariate set

| adjustment | z_delta | z_p | reldev_delta | reldev_p |
| --- | --- | --- | --- | --- |
| none (pre-specified endpoint) | 0.380 | 0.012 | 0.356 | 0.019 |
| cd8_frac | 0.356 | 0.019 | 0.373 | 0.014 |
| log_cells | 0.362 | 0.017 | 0.408 | 0.007 |
| log_density | 0.323 | 0.033 | 0.358 | 0.018 |
| log_area | 0.250 | 0.099 | 0.360 | 0.017 |
| log_cells + log_area (minimal, non-redundant) | 0.226 | 0.136 | 0.185 | 0.223 |
| log_density + log_cells | 0.343 | 0.023 | 0.397 | 0.009 |
| log_cells + log_area + cd8_frac | 0.191 | 0.207 | 0.217 | 0.151 |
| all four (near-collinear) | 0.008 | 0.966 | -0.044 | 0.774 |

The four covariates are near-collinear (density = cells / area; condition number **448**), and together they explain only ~9% of the variance of z, yet the full set drives δ from +0.379 to +0.008. A 9%-of-variance covariate block extinguishing a moderate group difference is a signature of over-adjustment under collinearity, not of confounding removed.

## 2. Matching — assumption-free, and it does not collapse

| area+density band | n_ROI | z_delta | z_p | reldev_delta | reldev_p | n_NR_pat | n_R_pat |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all ROIs | 593 | 0.380 | 0.012 | 0.356 | 0.019 | 25 | 37 |
| 10–90 pct | 384 | 0.270 | 0.084 | 0.346 | 0.026 | 24 | 34 |
| 20–80 pct | 213 | 0.376 | 0.027 | 0.355 | 0.037 | 20 | 29 |
| IQR (25–75 pct) | 146 | 0.326 | 0.098 | 0.297 | 0.132 | 14 | 25 |

Restricting to a common band of area **and** density needs no model and cannot suffer collinearity. δ stays at **+0.270 … +0.379** in every band — it does not shrink toward zero, which is what a pure geometry artefact would do (compare `IMBALANCE_CONFOUND.md`, where the Day-1 signal fell from −0.138 to +0.035 as imbalanced ROIs were dropped). p-values widen because matching costs patients.

## Verdict

- Across reasonable specifications δ ranges **+0.008 … +0.379** (p from 0.009 to 0.97). The result is **not robust to analytic choice**, so it is neither established nor refuted.
- Matching, the more trustworthy arbiter, leaves the effect intact (δ ≈ +0.32) but underpowered (p ≈ 0.03–0.10 at n = 14–24 NR patients).
- **Therefore: indeterminate.** The honest statement is that this design cannot separate a moderate arrangement effect from ROI geometry. Resolving it needs more patients or ROI sampling matched on area and density by design, not more analysis of these data.
- **The direction is stable and matters more than the significance.** Responders score *higher* on the exclusion-like axis under every specification. The immune-exclusion hypothesis predicts the opposite (non-responders higher). So no specification supports 'reduced infiltration of cytotoxic T cells in non-responders'; the unresolved question is whether there is a weak effect running the other way.

Script: `day2_sensitivity.py`.
