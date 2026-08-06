# Chromatic vs Dowker on the B7H4 question

556 pre-treatment ROIs with both analyses. Paired differences are B7H4+ minus B7H4−, within ROI, at matched cell counts.

## 1. Where the methods agree

| signal | rho_between_methods | p |
| --- | --- | --- |
| count | 0.6788 | 0 |
| length | 0.056 | 0.1877 |

The **count** signal is strongly concordant between two structurally independent constructions — that is the robust result. The **length** signal is essentially uncorrelated between them, so at most one of the two can be a clean readout of CD8 accessibility.

## 2. Which signal tracks the clustering confound?

At matched counts B7H4+ cancer cells are more tightly packed (median NN 11.3 vs 14.2 µm, tighter in 95% of ROIs). Correlating each paired difference against the paired clustering difference:

| signal | method | statistic | rho_with_clustering | p |
| --- | --- | --- | --- | --- |
| count | chromatic | ker1-num_bars | 0.2131 | 0 |
| count | Dowker | dow1-num_bars | 0.3647 | 0 |
| length | chromatic | ker0-avg_length | 0.0537 | 0.2064 |
| length | Dowker | dow1-avg_length | 0.2332 | 0 |

## 3. Restricting to ROIs where clustering is matched

The 139 ROIs whose subtypes are most similarly packed (|Δ nearest-neighbour distance| below its 25th percentile, 1.34 µm):

| signal | method | statistic | median_all | p_all | median_matched | p_matched |
| --- | --- | --- | --- | --- | --- | --- |
| count | chromatic | ker1-num_bars | -41 | 0 | -18 | 0 |
| count | Dowker | dow1-num_bars | -5.7 | 0 | -1.8 | 0.0012 |
| length | chromatic | ker0-avg_length | 0.2874 | 0 | 0.2667 | 0.0305 |
| length | Dowker | dow1-avg_length | -0.4052 | 0 | -0.1479 | 0.2095 |

**Verdict: the **Dowker** length effect is the clustering-driven one — it correlates with the clustering difference and vanishes when clustering is matched, while the chromatic one does neither.**

## What can be said

- **Robust:** B7H4+ tumour carries fewer degree-1 topological features than B7H4− tumour in the same tissue, by two independent constructions that agree strongly (ρ ≈ 0.68). Part of this tracks the clustering difference.
- **Directly measured, not inferred:** B7H4+ cancer cells are more tightly clustered.
- **Provisional:** the longer-bar 'exclusion' signature appears in the chromatic construction, is not explained by clustering, and survives clustering-matching — but it is not reproduced by an independent construction, so it should be reported as method-dependent rather than established.
- **Not supported:** a clean claim that CD8 accessibility differs between the subtypes beyond the architectural difference. The two methods would need to agree for that, and they do not.

Script: `day2_method_comparison.py`.
