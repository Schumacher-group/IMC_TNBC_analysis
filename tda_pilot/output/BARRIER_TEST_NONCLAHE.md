# Barrier test (nonclahe masks)

442 pre-treatment ROIs, 57 patients. For each CD8 cell, the fraction of the straight path to its nearest tumour cell that lies inside collagen; median over cells per ROI. Compared against 99 relabellings in which the same number of CD8 cells is drawn from the non-tumour population of that ROI, holding cell positions, tumour and collagen fixed.

Collagen is preferentially between CD8 and tumour overall: median excess over the within-ROI null = **-0.0253** (18% of ROIs positive).

| statistic | meaning | median_R | median_NR | patient_delta | patient_p | roi_delta | roi_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| obs | raw path fraction (confounded by collagen amount) | 0.1222 | 0.1391 | -0.0423 | 0.7976 | -0.0291 | 0.6523 |
| excess | excess over within-ROI null | -0.0317 | -0.0205 | -0.0688 | 0.6731 | -0.082 | 0.2042 |
| z | z vs within-ROI null | -4.742 | -3.508 | -0.045 | 0.7849 | -0.0199 | 0.7588 |

`obs` is reported only to show what it does WITHOUT the null -- it is a measure of how much collagen an ROI contains as much as of where that collagen sits. `excess` and `z` are the interpretable ones.

Figure: `barrier_test_nonclahe.png`. Script: `day2_barrier_test.py`.
