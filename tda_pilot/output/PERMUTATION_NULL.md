# Permutation-null test: CD8_primary, pre-treatment, inclusion = `symmetric`

Inclusion: `KChromaticInclusion(filt, 1)` — colour-symmetric: domain is the monochromatic subcomplex of both colours, so the kernel is tumour-dominated at a ~2:1 ratio.

Cohort: 593 ROIs / 25 NR + 37 R patients. B = 99 label permutations per ROI, cell positions held fixed.

Each ROI's z = (observed − null mean) / null sd. Composition, cell counts, density, hull area and dispersion are identical between an ROI and its own null, so they cannot bias the null MEAN. They do, however, set the null SD, so z is comparable within an ROI but not across ROIs of different geometry — see the cross-ROI robustness section, which is the decisive one. Cliff's δ is R vs NR (+ve = responders higher).

## Sanity checks

- `cod0-num_bars` is invariant under relabelling: max null sd = **0.000e+00** (must be 0). PASS
- All 50 z columns finite: **True**
- Resolving power: **42%** of all ROI × statistic z-values exceed |z| = 2, and the most deviant statistics sit enormously far from their nulls (median |z| 18 for `ker0-entropy`, 17 for `cok1-num_bars`; overall median |z| = 1.5). Tumour and CD8 are therefore *strongly* non-randomly arranged in both groups — the assay has power, and a null response contrast is informative rather than vacuous.

## Primary endpoint

`ker1-avg_length` z, per-patient median:

- Cliff's δ = **+0.165**
- patient-label permutation p (20,000 draws) = **0.2733**
- bootstrap 95% CI over patients = [-0.150, +0.464]
- leave-one-patient-out δ range = +0.131 … +0.212
- median z: NR +3.25, R +4.24
- size-free companion (deviation as a fraction of the null mean): δ = **+0.254**, p = 0.0932 — guards against |z| growing with ROI size, since the null sd shrinks as an ROI gains cells


## Cross-ROI geometry robustness (decisive)

The null holds composition and geometry fixed **within** each ROI, so the null mean is unbiased. But z = (obs − mean)/sd and the null sd shrinks as an ROI gains cells or density, so the arrangement-to-z mapping is itself geometry-dependent: z is comparable within an ROI, **not across** ROIs of different size and density. Responder and non-responder ROIs do differ in geometry, so a raw-z response contrast can be produced by that alone. Residualising on CD8 fraction, log(cells), log(density) and log(hull area):

| estimator | raw_delta | geom_adj_delta | geom_adj_p |
| --- | --- | --- | --- |
| z | 0.1654 | 0.0486 | 0.7522 |
| relative deviation | 0.2541 | 0.0205 | 0.8972 |

**Survives geometry adjustment: False.**

**Primary endpoint is NOT significant.** The decision rule stands: the pilot's no-go is confirmed on the strongest available evidence. Even after conditioning each ROI on its own composition- and geometry-matched null, baseline tumour-CD8 topology does not distinguish responders from non-responders. The bar-length signal in `LENGTH_FEATURES.md` is therefore attributable to tissue geometry (responder ROIs are less dense and more dispersed), not to spatial arrangement of T cells relative to tumour.


## Balance sensitivity (primary endpoint)

The Day-1 signal collapsed as imbalanced ROIs were dropped (IMBALANCE_CONFOUND.md). Class imbalance is already held fixed within each ROI here, so this ladder is not testing the same thing — it asks whether the arrangement contrast is concentrated in ROIs where both species are well represented (where the segregation-vs-interspersion question is best posed).

| min_frac>= | n_ROI | n_NR_pat | n_R_pat | pat_delta | pat_p |
| --- | --- | --- | --- | --- | --- |
| 0 | 593 | 25 | 37 | 0.1654 | 0.2755 |
| 0.1 | 554 | 24 | 37 | 0.1689 | 0.2714 |
| 0.15 | 521 | 23 | 37 | 0.2009 | 0.1962 |
| 0.2 | 463 | 22 | 37 | 0.2113 | 0.1802 |
| 0.25 | 388 | 20 | 34 | 0.3206 | 0.052 |
| 0.3 | 307 | 20 | 33 | 0.3424 | 0.039 |

Observed: delta rises monotonically from +0.165 (all ROIs) to +0.342 at minority fraction >= 0.30, reaching nominal p < 0.05 on the last 1 rung(s). Treat this as hypothesis-generating only — it is a post-hoc scan over six nested, highly correlated thresholds, patient counts fall as the threshold rises (25→20 NR), and the primary endpoint was pre-specified at the 0.00 rung. It is not evidence for an effect, but it is the one place worth a pre-registered look if this line is ever revisited.

## Secondary: all retained statistics (per-patient, BH-adjusted)

| statistic | roi_delta | pat_delta | pat_p | reldev_pat_delta | median_z_NR | median_z_R | pat_q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| dom0-entropy | -0.2196 | -0.3081 | 0.0416 | -0.3081 | -3.226 | -5.01 | 0.5988 |
| ker0-entropy | -0.2097 | -0.2951 | 0.051 | -0.2086 | -14.73 | -19.51 | 0.5988 |
| ker1-med_length | 0.0945 | 0.2822 | 0.0621 | 0.3297 | -0.4113 | -0.089 | 0.5988 |
| dom0-avg_length | -0.2322 | -0.2778 | 0.0662 | -0.2497 | -8.531 | -11.76 | 0.5988 |
| dom0-med_length | -0.2405 | -0.2778 | 0.0662 | -0.2627 | -8.387 | -11.47 | 0.5988 |
| ker0-avg_length | 0.2031 | 0.267 | 0.0776 | 0.3232 | 6.623 | 9.164 | 0.5988 |
| ker0-num_bars | -0.2065 | -0.2497 | 0.0989 | -0.2476 | -12.42 | -15.91 | 0.5988 |
| dom1-med_length | -0.1833 | -0.2151 | 0.1554 | -0.2346 | -0.4851 | -1.032 | 0.5988 |
| cok1-entropy | -0.1452 | -0.2151 | 0.1554 | -0.2541 | -11.74 | -16.27 | 0.5988 |
| ker1-entropy | -0.1826 | -0.2086 | 0.1683 | -0.1914 | -8.148 | -11.95 | 0.5988 |
| cok1-med_length | 0.1619 | 0.2086 | 0.1683 | 0.1978 | 0.903 | 1.769 | 0.5988 |
| dom1-avg_length | -0.1676 | -0.2086 | 0.1683 | -0.2324 | -3.701 | -4.733 | 0.5988 |
| cok1-avg_length | 0.1344 | 0.1978 | 0.1916 | 0.2324 | 5.125 | 6.471 | 0.5988 |
| cok1-p90_length | 0.1061 | 0.1935 | 0.2016 | 0.2086 | 3.485 | 4.402 | 0.5988 |
| ker1-num_bars | -0.1719 | -0.1827 | 0.2281 | -0.2216 | -11.14 | -13.58 | 0.5988 |
| ker0-med_length | 0.0886 | 0.1805 | 0.2336 | 0.2519 | 1.311 | 2.498 | 0.5988 |
| ker0-p90_length | 0.1679 | 0.1805 | 0.2336 | 0.2843 | 6.006 | 7.911 | 0.5988 |
| im1-num_bars | 0.1732 | 0.1805 | 0.2336 | 0.187 | 12.15 | 16.31 | 0.5988 |
| cod0-avg_length | -0.0625 | -0.1676 | 0.2375 | -0.1632 | 0 | 0 | 0.5988 |
| cod1-med_length | -0.0504 | -0.1297 | 0.2544 | -0.1351 | 0 | 0 | 0.5988 |
| im1-entropy | 0.1456 | 0.1697 | 0.263 | 0.1632 | 7.926 | 10.4 | 0.5988 |
| cok1-num_bars | -0.1661 | -0.1676 | 0.2692 | -0.2497 | -13.85 | -18.28 | 0.5988 |
| ker1-avg_length | 0.1575 | 0.1654 | 0.2755 | 0.2541 | 3.248 | 4.245 | 0.5988 |
| dom1-num_bars | 0.0785 | 0.1611 | 0.2883 | 0.1762 | 0.8118 | 1.206 | 0.6006 |
| im1-avg_length | -0.1009 | -0.1503 | 0.3221 | -0.1676 | -2.775 | -3.47 | 0.6442 |
| im1-p90_length | -0.0951 | -0.1416 | 0.351 | -0.107 | -1.362 | -1.692 | 0.6535 |
| im0-avg_length | -0.0752 | -0.1286 | 0.3583 | -0.1459 | 0 | 0 | 0.6535 |
| im1-med_length | -0.0977 | -0.1373 | 0.366 | -0.0616 | -0.1901 | -0.3536 | 0.6535 |
| dom1-p90_length | -0.1494 | -0.1265 | 0.4052 | -0.2259 | -3.205 | -4.114 | 0.6987 |
| cod0-entropy | -0.0304 | -0.1157 | 0.4462 | -0.1178 | 0.1759 | 0 | 0.7437 |
| cod1-entropy | 0.0426 | 0.0822 | 0.5903 | 0.0822 | 0 | 0 | 0.9109 |
| cod1-p90_length | 0.0301 | 0.0735 | 0.6106 | 0.1341 | 0 | 0 | 0.9109 |
| im0-med_length | -0.0418 | -0.0692 | 0.6369 | -0.0724 | 0 | 0 | 0.9109 |
| cod0-med_length | -0.0418 | -0.0692 | 0.6369 | -0.0724 | 0 | 0 | 0.9109 |
| cod1-avg_length | 0.0167 | -0.0659 | 0.6376 | -0.0443 | 0 | 0 | 0.9109 |
| im0-p90_length | 0.0158 | -0.0497 | 0.7197 | -0.0357 | 0 | 0 | 0.9725 |
| cod0-p90_length | 0.0158 | -0.0497 | 0.7197 | -0.0357 | 0 | 0 | 0.9725 |
| ker1-p90_length | 0.0173 | 0.0249 | 0.8746 | 0.0984 | 2.316 | 2.471 | 1 |
| im0-entropy | -0.0175 | 0.0227 | 0.8858 | 0.0292 | 0.1849 | 0 | 1 |
| dom0-p90_length | -0.0048 | -0.0119 | 0.9428 | 0.0162 | -3.427 | -3.513 | 1 |
| dom1-entropy | 0.0012 | -0.0054 | 0.9771 | -0.0097 | -0.1545 | -0.0686 | 1 |
| cod0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cok0-num_bars | 0 | 0 | 1 | nan | 0 | 0 | 1 |
| im0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cok0-avg_length | 0 | 0 | 1 | nan | 0 | 0 | 1 |
| cok0-med_length | 0 | 0 | 1 | nan | 0 | 0 | 1 |
| cok0-p90_length | 0 | 0 | 1 | nan | 0 | 0 | 1 |
| dom0-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cod1-num_bars | 0 | 0 | 1 | 0 | 0 | 0 | 1 |
| cok0-entropy | 0 | 0 | 1 | nan | 0 | 0 | 1 |

## Limitations — what this does and does not rule out

- **Power.** 25 NR vs 37 R patients gives ~80% power only at Cliff's |δ| ≈ 0.42. The primary CI is [-0.150, +0.464]: a moderate responder-higher effect is not excluded. 'No difference' here means 'no large difference'.
- **The two estimators of the same quantity disagree in magnitude.** z gives δ = +0.165 (p = 0.273); the size-free relative deviation gives δ = +0.254 (p = 0.093). Both are non-significant, but the second is borderline, so this is a soft null rather than a clean one.
- **Scope.** One pair definition (`CD8_primary`), one inclusion (`KChromaticInclusion(filt, 1)`, which is colour-symmetric and therefore tumour-dominated), degree 0 and 1 only, one filtration (`delaunay_cech`). A directional `SubChromaticInclusion` on CD8 alone — the statistic an *exclusion* hypothesis actually calls for — has still not been run under this null.
- **B = 99** bounds the per-ROI null resolution; z is estimated from 99 draws, so per-ROI z carries noise of order 1/sqrt(2·99) ≈ 7% on the sd.

Figure: `perm_null_primary.png`. Scripts: `perm_null.py`, `day2_permutation_test.py`.
