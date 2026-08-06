# Dowker persistent homology: independent check on the B7H4 finding

556 pre-treatment ROIs. Landmarks = cancer cells (cap 300), witnesses = CD8 cells, cutoff 100 µm, 5 matched subsamples averaged. A long degree-1 bar is a loop of cancer cells that no single CD8 cell is close to — a region of tumour CD8 does not reach.

**Why this is an independent check, not a re-summary:** the Dowker filtration uses only cross-distances between the two species. Within-species structure — the source of every confound in the chromatic analyses (bar count ~ cell number ρ = 0.83, bar length ~ composition ρ = −0.69, z coupled to density) — never enters it.

## Primary: paired B7H4+ vs B7H4−  (13/16 at BH q < 0.05)

| statistic | median_cancer | median_b7h4 | median_diff | frac_b7h4_higher | wilcoxon_p | rank_biserial | q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| dow0-entropy | 4.472 | 4.349 | -0.1145 | 0.1601 | 0 | -0.6799 | 0 |
| dow1-total_persistence | 347.4 | 282.8 | -48.83 | 0.196 | 0 | -0.6079 | 0 |
| dow1-num_bars | 50.1 | 43 | -5.7 | 0.196 | 0 | -0.6079 | 0 |
| dow0-num_bars | 120.8 | 111.1 | -8.8 | 0.1781 | 0 | -0.6439 | 0 |
| dow0-total_persistence | 2041 | 1809 | -166.5 | 0.2158 | 0 | -0.5683 | 0 |
| dow1-entropy | 3.478 | 3.33 | -0.1436 | 0.1906 | 0 | -0.6187 | 0 |
| dow1-med_length | 4.82 | 4.378 | -0.4464 | 0.3291 | 0 | -0.3417 | 0 |
| dow0-med_length | 13.98 | 13.08 | -0.5889 | 0.3579 | 0 | -0.2842 | 0 |
| dow1-avg_length | 7.004 | 6.619 | -0.4052 | 0.3507 | 0 | -0.2986 | 0 |
| dow1-p90_length | 16.22 | 15.3 | -0.8823 | 0.3957 | 0 | -0.2086 | 0 |
| dow0-p90_length | 35.18 | 36.04 | 1.408 | 0.6277 | 0 | 0.2554 | 0 |
| dow1-max_length | 31.55 | 29.66 | -2.011 | 0.3777 | 0 | -0.2446 | 0 |
| dow0-avg_length | 17.19 | 16.96 | -0.0993 | 0.4712 | 0.0058 | -0.0576 | 0.0071 |
| dow0-censored | 2 | 2 | 0 | 0.3597 | 0.1823 | -0.2806 | 0.2083 |
| dow0-max_length | 70.65 | 70.99 | 0.5393 | 0.5198 | 0.5364 | 0.0396 | 0.5722 |
| dow1-censored | 2 | 2 | 0 | 0.4191 | 0.6824 | -0.1619 | 0.6824 |

## Secondary: response contrast on whole-tumour Dowker  (0/16 at BH q < 0.05)

| statistic | cliff_R_vs_NR | p | q |
| --- | --- | --- | --- |
| dow1-entropy | -0.2562 | 0.0904 | 0.399 |
| dow0-p90_length | 0.2519 | 0.096 | 0.399 |
| dow0-max_length | 0.2432 | 0.108 | 0.399 |
| dow1-num_bars | -0.2422 | 0.1096 | 0.399 |
| dow0-avg_length | 0.2324 | 0.1247 | 0.399 |
| dow0-med_length | 0.1589 | 0.2948 | 0.5845 |
| dow1-total_persistence | -0.1438 | 0.3436 | 0.5845 |
| dow0-censored | 0.1341 | 0.3748 | 0.5845 |
| dow0-num_bars | -0.1286 | 0.3972 | 0.5845 |
| dow1-max_length | -0.12 | 0.43 | 0.5845 |
| dow1-censored | 0.1189 | 0.4334 | 0.5845 |
| dow0-entropy | -0.1178 | 0.4384 | 0.5845 |
| dow1-med_length | 0.0811 | 0.5955 | 0.6919 |
| dow1-avg_length | 0.0789 | 0.6054 | 0.6919 |
| dow0-total_persistence | 0.0076 | 0.9657 | 1 |
| dow1-p90_length | -0.0011 | 1 | 1 |

Figure: `dowker_b7h4.png`. Scripts: `dowker.py` (construction + validation), `day2_dowker.py`.
