# Permutation-null test: CD8_primary, pre-treatment

Cohort: 582 ROIs / 25 NR + 37 R patients. B = 99 label permutations per ROI, cell positions held fixed.

Each ROI's z = (observed − null mean) / null sd. Composition, cell counts, density, hull area and dispersion are identical between an ROI and its own null, so they cannot produce a difference in z. Cliff's δ is R vs NR (+ve = responders higher).

## Sanity checks

- `cod0-num_bars` is invariant under relabelling: max null sd = **0.000e+00** (must be 0). PASS
- All 50 z columns finite: **True**
- Resolving power: **42%** of all ROI × statistic z-values exceed |z| = 2, and the most deviant statistics sit enormously far from their nulls (median |z| 18 for `ker0-entropy`, 17 for `cok1-num_bars`; overall median |z| = 1.5). Tumour and CD8 are therefore *strongly* non-randomly arranged in both groups — the assay has power, and a null response contrast is informative rather than vacuous.

## Primary endpoint

`ker1-avg_length` z, per-patient median:

- Cliff's δ = **+0.174**
- patient-label permutation p (20,000 draws) = **0.2505**
- bootstrap 95% CI over patients = [-0.140, +0.472]
- leave-one-patient-out δ range = +0.140 … +0.221
- median z: NR +3.25, R +4.25
- size-free companion (deviation as a fraction of the null mean): δ = **+0.269**, p = 0.0752 — guards against |z| growing with ROI size, since the null sd shrinks as an ROI gains cells

**Primary endpoint is NOT significant.** The decision rule stands: the pilot's no-go is confirmed on the strongest available evidence. Even after conditioning each ROI on its own composition- and geometry-matched null, baseline tumour-CD8 topology does not distinguish responders from non-responders. The bar-length signal in `LENGTH_FEATURES.md` is therefore attributable to tissue geometry (responder ROIs are less dense and more dispersed), not to spatial arrangement of T cells relative to tumour.

## Balance sensitivity (primary endpoint)

The Day-1 signal collapsed as imbalanced ROIs were dropped (IMBALANCE_CONFOUND.md). Class imbalance is already held fixed within each ROI here, so this ladder is not testing the same thing — it asks whether the arrangement contrast is concentrated in ROIs where both species are well represented (where the segregation-vs-interspersion question is best posed).

| min_frac>= | n_ROI | n_NR_pat | n_R_pat | pat_delta | pat_p |
| --- | --- | --- | --- | --- | --- |
| 0 | 582 | 25 | 37 | 0.1741 | 0.251 |
| 0.1 | 543 | 24 | 37 | 0.1734 | 0.2587 |
| 0.15 | 510 | 23 | 37 | 0.2056 | 0.1859 |
| 0.2 | 453 | 22 | 37 | 0.2334 | 0.1385 |
| 0.25 | 381 | 20 | 34 | 0.3382 | 0.0403 |
| 0.3 | 302 | 20 | 33 | 0.3394 | 0.0408 |

Observed: delta rises monotonically from +0.174 (all ROIs) to +0.339 at minority fraction >= 0.30, reaching nominal p < 0.05 on the last 2 rung(s). Treat this as hypothesis-generating only — it is a post-hoc scan over six nested, highly correlated thresholds, patient counts fall as the threshold rises (25→20 NR), and the primary endpoint was pre-specified at the 0.00 rung. It is not evidence for an effect, but it is the one place worth a pre-registered look if this line is ever revisited.

## Secondary: all retained statistics (per-patient, BH-adjusted)

| statistic | roi_delta | pat_delta | pat_p | reldev_pat_delta | median_z_NR | median_z_R | pat_q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| dom0-entropy | -0.2153 | -0.3081 | 0.0416 | -0.3211 | -3.226 | -4.935 | 0.5353 |
| ker0-entropy | -0.2053 | -0.2951 | 0.051 | -0.2195 | -14.73 | -19.49 | 0.5353 |
| ker1-med_length | 0.0997 | 0.2865 | 0.0582 | 0.3341 | -0.4113 | -0.0861 | 0.5353 |
| dom0-avg_length | -0.2257 | -0.28 | 0.0641 | -0.2692 | -8.531 | -11.73 | 0.5353 |
| dom0-med_length | -0.2375 | -0.2778 | 0.0662 | -0.28 | -8.387 | -11.43 | 0.5353 |
| ker0-avg_length | 0.205 | 0.2714 | 0.0729 | 0.3405 | 6.623 | 9.162 | 0.5353 |
| ker0-num_bars | -0.2039 | -0.2519 | 0.096 | -0.2541 | -12.42 | -15.91 | 0.5353 |
| cok1-entropy | -0.1455 | -0.2303 | 0.1282 | -0.2735 | -11.74 | -16.28 | 0.5353 |
| ker1-entropy | -0.1833 | -0.2195 | 0.1472 | -0.213 | -8.148 | -12.03 | 0.5353 |
| dom1-med_length | -0.1802 | -0.2173 | 0.1513 | -0.2454 | -0.4851 | -1.034 | 0.5353 |
| dom1-avg_length | -0.1629 | -0.213 | 0.1596 | -0.2432 | -3.701 | -4.7 | 0.5353 |
| cok1-p90_length | 0.1048 | 0.2065 | 0.1728 | 0.2346 | 3.485 | 4.387 | 0.5353 |
| cok1-med_length | 0.1674 | 0.2022 | 0.182 | 0.2 | 0.903 | 1.802 | 0.5353 |
| im1-num_bars | 0.1744 | 0.1978 | 0.1916 | 0.1957 | 12.15 | 16.42 | 0.5353 |
| cok1-avg_length | 0.1385 | 0.1978 | 0.1916 | 0.2432 | 5.125 | 6.597 | 0.5353 |
| ker0-p90_length | 0.1717 | 0.1935 | 0.2016 | 0.293 | 6.006 | 8.119 | 0.5353 |
| cod0-avg_length | -0.062 | -0.1816 | 0.2032 | -0.1751 | 0 | 0 | 0.5353 |
| ker0-med_length | 0.0968 | 0.1914 | 0.2067 | 0.267 | 1.311 | 2.51 | 0.5353 |
| ker1-num_bars | -0.1701 | -0.1892 | 0.2119 | -0.2324 | -11.14 | -13.57 | 0.5353 |
| cok1-num_bars | -0.1671 | -0.187 | 0.2172 | -0.2649 | -13.85 | -18.45 | 0.5353 |
| im1-entropy | 0.1447 | 0.1784 | 0.2393 | 0.1849 | 7.926 | 10.47 | 0.5353 |
| ker1-avg_length | 0.1647 | 0.1741 | 0.251 | 0.2692 | 3.248 | 4.255 | 0.5353 |
| cod1-med_length | -0.0541 | -0.1297 | 0.2544 | -0.1351 | 0 | 0 | 0.5353 |
| dom1-num_bars | 0.0809 | 0.1719 | 0.2569 | 0.1892 | 0.8118 | 1.209 | 0.5353 |
| im0-avg_length | -0.0765 | -0.1427 | 0.3115 | -0.1578 | 0 | 0 | 0.6194 |
| im1-avg_length | -0.0965 | -0.1503 | 0.3221 | -0.1784 | -2.775 | -3.461 | 0.6194 |
| dom1-p90_length | -0.1467 | -0.1416 | 0.351 | -0.2497 | -3.205 | -4.116 | 0.64 |
| im1-p90_length | -0.088 | -0.1395 | 0.3584 | -0.1114 | -1.362 | -1.665 | 0.64 |
| cod0-entropy | -0.0351 | -0.133 | 0.3807 | -0.1373 | 0.1759 | 0 | 0.6487 |
| im1-med_length | -0.0961 | -0.1308 | 0.3892 | -0.0616 | -0.1901 | -0.3526 | 0.6487 |
| cod1-entropy | 0.0412 | 0.0822 | 0.5903 | 0.0822 | 0 | 0 | 0.9445 |
| cod1-p90_length | 0.0264 | 0.0735 | 0.6106 | 0.1341 | 0 | 0 | 0.9445 |
| cod1-avg_length | 0.0112 | -0.0659 | 0.6376 | -0.0443 | 0 | 0 | 0.9445 |
| im0-p90_length | 0.0137 | -0.0605 | 0.6611 | -0.0378 | 0 | 0 | 0.9445 |
| cod0-p90_length | 0.0137 | -0.0605 | 0.6611 | -0.0378 | 0 | 0 | 0.9445 |
| im0-med_length | -0.0397 | -0.0476 | 0.7477 | -0.0573 | 0 | 0 | 1 |
| cod0-med_length | -0.0397 | -0.0476 | 0.7477 | -0.0573 | 0 | 0 | 1 |
| ker1-p90_length | 0.0191 | 0.027 | 0.8633 | 0.1135 | 2.316 | 2.472 | 1 |
| dom0-p90_length | 0.0009 | -0.0162 | 0.92 | 0.0032 | -3.427 | -3.494 | 1 |
| im0-entropy | -0.0248 | 0.0097 | 0.9542 | 0.0097 | 0.1849 | -0.0658 | 1 |
| dom1-entropy | 0.0046 | -0.0054 | 0.9771 | -0.0097 | -0.1545 | -0.0511 | 1 |
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

- **Power.** 25 NR vs 37 R patients gives ~80% power only at Cliff's |δ| ≈ 0.42. The primary CI is [-0.140, +0.472]: a moderate responder-higher effect is not excluded. 'No difference' here means 'no large difference'.
- **The two estimators of the same quantity disagree in magnitude.** z gives δ = +0.174 (p = 0.251); the size-free relative deviation gives δ = +0.269 (p = 0.075). Both are non-significant, but the second is borderline, so this is a soft null rather than a clean one.
- **Scope.** One pair definition (`CD8_primary`), one inclusion (`KChromaticInclusion(filt, 1)`, which is colour-symmetric and therefore tumour-dominated), degree 0 and 1 only, one filtration (`delaunay_cech`). A directional `SubChromaticInclusion` on CD8 alone — the statistic an *exclusion* hypothesis actually calls for — has still not been run under this null.
- **B = 99** bounds the per-ROI null resolution; z is estimated from 99 draws, so per-ROI z carries noise of order 1/sqrt(2·99) ≈ 7% on the sd.

Figure: `perm_null_primary.png`. Scripts: `perm_null.py`, `day2_permutation_test.py`.
