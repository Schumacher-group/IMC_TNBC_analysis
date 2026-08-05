# Permutation-null test: CD8_primary, pre-treatment, inclusion = `other_only`

Inclusion: `SubChromaticInclusion(filt, [[CD8]])` — DIRECTIONAL: domain is the CD8-only complex, so degree-1 kernel bars are CD8 loops that tumour fills in, the signature an immune-exclusion hypothesis names.

Cohort: 593 ROIs / 25 NR + 37 R patients. B = 99 label permutations per ROI, cell positions held fixed.

Each ROI's z = (observed − null mean) / null sd. Composition, cell counts, density, hull area and dispersion are identical between an ROI and its own null, so they cannot bias the null MEAN. They do, however, set the null SD, so z is comparable within an ROI but not across ROIs of different geometry — see the cross-ROI robustness section, which is the decisive one. Cliff's δ is R vs NR (+ve = responders higher).

## Sanity checks

- `cod0-num_bars` is invariant under relabelling: max null sd = **0.000e+00** (must be 0). PASS
- All 50 z columns finite: **True**
- Resolving power: **45%** of all ROI × statistic z-values exceed |z| = 2, and the most deviant statistics sit enormously far from their nulls (median |z| 13 for `cok1-num_bars`, 13 for `ker0-entropy`; overall median |z| = 2.1). Tumour and CD8 are therefore *strongly* non-randomly arranged in both groups — the assay has power, and a null response contrast is informative rather than vacuous.

## Primary endpoint

`ker1-avg_length` z, per-patient median:

- Cliff's δ = **+0.379**
- patient-label permutation p (20,000 draws) = **0.0106**
- bootstrap 95% CI over patients = [+0.094, +0.643]
- leave-one-patient-out δ range = +0.356 … +0.435
- median z: NR +2.11, R +2.91
- size-free companion (deviation as a fraction of the null mean): δ = **+0.356**, p = 0.0186 — guards against |z| growing with ROI size, since the null sd shrinks as an ROI gains cells


## Cross-ROI geometry robustness (decisive)

The null holds composition and geometry fixed **within** each ROI, so the null mean is unbiased. But z = (obs − mean)/sd and the null sd shrinks as an ROI gains cells or density, so the arrangement-to-z mapping is itself geometry-dependent: z is comparable within an ROI, **not across** ROIs of different size and density. Responder and non-responder ROIs do differ in geometry, so a raw-z response contrast can be produced by that alone. Residualising on CD8 fraction, log(cells), log(density) and log(hull area):

| estimator | raw_delta | geom_adj_delta | geom_adj_p |
| --- | --- | --- | --- |
| z | 0.3795 | 0.0076 | 0.9657 |
| relative deviation | 0.3557 | -0.0443 | 0.7741 |

**Survives geometry adjustment: False.**

**Primary endpoint is NOT significant.** The decision rule stands: the pilot's no-go is confirmed on the strongest available evidence. Even after conditioning each ROI on its own composition- and geometry-matched null, baseline tumour-CD8 topology does not distinguish responders from non-responders. The bar-length signal in `LENGTH_FEATURES.md` is therefore attributable to tissue geometry (responder ROIs are less dense and more dispersed), not to spatial arrangement of T cells relative to tumour.

Pre-specified direction for exclusion was δ < 0 (non-responders higher); observed δ = +0.379, not significant either way.


## Balance sensitivity (primary endpoint)

The Day-1 signal collapsed as imbalanced ROIs were dropped (IMBALANCE_CONFOUND.md). Class imbalance is already held fixed within each ROI here, so this ladder is not testing the same thing — it asks whether the arrangement contrast is concentrated in ROIs where both species are well represented (where the segregation-vs-interspersion question is best posed).

| min_frac>= | n_ROI | n_NR_pat | n_R_pat | pat_delta | pat_p |
| --- | --- | --- | --- | --- | --- |
| 0 | 593 | 25 | 37 | 0.3795 | 0.012 |
| 0.1 | 554 | 24 | 37 | 0.3716 | 0.0152 |
| 0.15 | 521 | 23 | 37 | 0.4501 | 0.0037 |
| 0.2 | 463 | 22 | 37 | 0.4054 | 0.0099 |
| 0.25 | 388 | 20 | 34 | 0.3824 | 0.0204 |
| 0.3 | 307 | 20 | 33 | 0.3061 | 0.0652 |

Observed: delta rises monotonically from +0.379 (all ROIs) to +0.306 at minority fraction >= 0.30, reaching nominal p < 0.05 on the last 5 rung(s). Treat this as hypothesis-generating only — it is a post-hoc scan over six nested, highly correlated thresholds, patient counts fall as the threshold rises (25→20 NR), and the primary endpoint was pre-specified at the 0.00 rung. It is not evidence for an effect, but it is the one place worth a pre-registered look if this line is ever revisited.

## Secondary: all retained statistics (per-patient, BH-adjusted)

| statistic | roi_delta | pat_delta | pat_p | reldev_pat_delta | median_z_NR | median_z_R | pat_q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ker1-avg_length | 0.1685 | 0.3795 | 0.012 | 0.3557 | 2.114 | 2.914 | 0.3826 |
| ker0-med_length | 0.1857 | 0.3665 | 0.0153 | 0.3535 | 2.206 | 3.589 | 0.3826 |
| ker1-med_length | 0.1524 | 0.3254 | 0.0314 | 0.4141 | 0.2445 | 0.6397 | 0.4682 |
| ker1-p90_length | 0.1147 | 0.3146 | 0.0375 | 0.3384 | 1.552 | 1.932 | 0.4682 |
| ker0-avg_length | 0.1898 | 0.2778 | 0.0662 | 0.2843 | 5.567 | 7.954 | 0.6624 |
| ker0-entropy | -0.2072 | -0.2541 | 0.0932 | -0.2497 | -9.799 | -14.02 | 0.7764 |
| ker1-num_bars | -0.2012 | -0.2368 | 0.1178 | -0.2605 | -7.944 | -9.914 | 0.8413 |
| ker0-p90_length | 0.1475 | 0.2 | 0.1868 | 0.1762 | 5.81 | 7.368 | 0.9153 |
| ker0-num_bars | -0.1978 | -0.1827 | 0.2281 | -0.2649 | -9.835 | -12.91 | 0.9153 |
| ker1-entropy | -0.1287 | -0.1805 | 0.2336 | -0.1676 | -6.337 | -8.054 | 0.9153 |
| cod1-entropy | -0.0297 | -0.1751 | 0.248 | -0.1643 | 0 | -0.0841 | 0.9153 |
| cod1-med_length | -0.0504 | -0.1297 | 0.2544 | -0.1351 | 0 | 0 | 0.9153 |
| im1-entropy | 0.1464 | 0.1697 | 0.263 | 0.0097 | 4.334 | 5.36 | 0.9153 |
| cok0-entropy | 0.0788 | 0.1654 | 0.2755 | 0.2281 | 3.274 | 3.71 | 0.9153 |
| cod0-avg_length | -0.0493 | -0.1492 | 0.3023 | -0.1708 | 0 | 0 | 0.9153 |
| dom0-med_length | 0.0441 | 0.1438 | 0.3436 | 0.1762 | -3.051 | -3.053 | 0.9153 |
| cod1-avg_length | -0.0666 | -0.1243 | 0.3664 | -0.1222 | 0 | 0 | 0.9153 |
| dom1-avg_length | 0.0779 | 0.1308 | 0.3892 | 0.2086 | -0.6693 | -0.2207 | 0.9153 |
| im0-med_length | 0.1248 | 0.1265 | 0.4052 | 0.1438 | 0.311 | 1.647 | 0.9153 |
| im1-num_bars | 0.1158 | 0.1243 | 0.4134 | 0.0076 | 9.53 | 11.73 | 0.9153 |
| cok1-num_bars | -0.127 | -0.1135 | 0.4555 | -0.2 | -10.87 | -13.57 | 0.9153 |
| cod0-entropy | -0.0811 | -0.0995 | 0.5134 | -0.0778 | 0 | 0 | 0.9153 |
| im1-med_length | -0.0459 | -0.0962 | 0.5278 | -0.0249 | -0.3556 | -0.5029 | 0.9153 |
| dom0-p90_length | -0.0111 | -0.0919 | 0.5467 | -0.0011 | 1.275 | 1.121 | 0.9153 |
| cok1-avg_length | -0.076 | -0.0897 | 0.5563 | -0.1568 | 0.3347 | -0.6201 | 0.9153 |
| cok1-med_length | -0.0162 | -0.0897 | 0.5563 | -0.1243 | 0.1131 | 0.0206 | 0.9153 |
| cok0-med_length | -0.0468 | -0.0876 | 0.566 | -0.0162 | 1.76 | 1.336 | 0.9153 |
| cok0-p90_length | -0.09 | -0.0876 | 0.566 | -0.1719 | -2.302 | -3.395 | 0.9153 |
| cok0-avg_length | -0.0623 | -0.0811 | 0.5955 | -0.1568 | -2.036 | -3.162 | 0.9153 |
| im0-avg_length | 0.0623 | 0.0811 | 0.5955 | 0.0789 | 2.036 | 3.162 | 0.9153 |
| cod1-p90_length | 0.0301 | 0.0735 | 0.6106 | 0.1341 | 0 | 0 | 0.9153 |
| cod0-med_length | -0.0418 | -0.0692 | 0.6369 | -0.0724 | 0 | 0 | 0.9153 |
| dom1-p90_length | 0.0292 | 0.0703 | 0.6461 | 0.053 | -0.6403 | -0.5639 | 0.9153 |
| im1-avg_length | 0.0362 | 0.0681 | 0.6564 | 0.1027 | -0.7436 | -0.5991 | 0.9153 |
| cok1-entropy | -0.0346 | -0.0681 | 0.6564 | -0.0508 | -3.687 | -3.347 | 0.9153 |
| im0-entropy | 0.0007 | -0.0659 | 0.6668 | -0.0789 | -2.462 | -2.582 | 0.9153 |
| dom1-num_bars | -0.1118 | -0.0638 | 0.6773 | -0.1459 | -0.5947 | -0.7793 | 0.9153 |
| cod0-p90_length | 0.0158 | -0.0497 | 0.7197 | -0.0357 | 0 | 0 | 0.9228 |
| dom0-avg_length | -0.0114 | -0.0551 | 0.7198 | 0.0422 | -1.468 | -1.597 | 0.9228 |
| dom1-entropy | -0.095 | -0.0508 | 0.7414 | -0.0746 | -0.5362 | -0.9301 | 0.9267 |
| dom0-entropy | 0.038 | 0.0443 | 0.7741 | 0.0378 | -2.51 | -2.633 | 0.944 |
| im1-p90_length | 0.0079 | -0.0184 | 0.9086 | 0.0141 | -0.4024 | -0.4256 | 1 |
| dom1-med_length | -0.0911 | -0.0141 | 0.9314 | 0.0551 | -0.2979 | -0.5322 | 1 |
| im0-p90_length | 0.0301 | -0.0119 | 0.9428 | 0.0335 | 2.816 | 3.469 | 1 |
| cok1-p90_length | -0.0669 | 0.0032 | 0.9886 | -0.0724 | 0.3044 | -0.3665 | 1 |
| dom0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cod0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cod1-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| im0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cok0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |

## Limitations — what this does and does not rule out

- **Power.** 25 NR vs 37 R patients gives ~80% power only at Cliff's |δ| ≈ 0.42. The primary CI is [+0.094, +0.643]: a moderate responder-higher effect is not excluded. 'No difference' here means 'no large difference'.
- **The two estimators of the same quantity disagree in magnitude.** z gives δ = +0.379 (p = 0.011); the size-free relative deviation gives δ = +0.356 (p = 0.019). Both are non-significant, but the second is borderline, so this is a soft null rather than a clean one.
- **Scope.** One pair definition (`CD8_primary`), one inclusion (`KChromaticInclusion(filt, 1)`, which is colour-symmetric and therefore tumour-dominated), degree 0 and 1 only, one filtration (`delaunay_cech`). A directional `SubChromaticInclusion` on CD8 alone — the statistic an *exclusion* hypothesis actually calls for — has still not been run under this null.
- **B = 99** bounds the per-ROI null resolution; z is estimated from 99 draws, so per-ROI z carries noise of order 1/sqrt(2·99) ≈ 7% on the sd.

Figure: `perm_null_primary_other_only.png`. Scripts: `perm_null.py`, `day2_permutation_test.py`.
