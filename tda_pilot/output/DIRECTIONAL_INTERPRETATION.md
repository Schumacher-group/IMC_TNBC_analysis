# Directional statistic: what does `ker1-avg_length` z actually measure?

Cohort reference: 593 pre-treatment ROIs, z median +2.75, IQR [+1.19, +5.03], range [-4.0, +188.3].

## Named pair

### Geometry / composition

| roi_id | n_all | hull_area | bbox_area | med_nn | density | response | sample_type | tumour | cd8_primary | cd8_frac |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Leap008_5 | 4969 | 9.438e+05 | 9.466e+05 | 10.26 | 0.005265 | Responder | pre | 1898 | 1155 | 0.3783 |
| Leap046_1 | 5253 | 1.025e+06 | 1.138e+06 | 9.73 | 0.005125 | Non-Responder | pre | 3215 | 611 | 0.1597 |

### z-scores (raw, geometry-adjusted, cohort percentile)

| inclusion | statistic | z[Leap008_5] | adj[Leap008_5] | pct[Leap008_5] | z[Leap046_1] | adj[Leap046_1] | pct[Leap046_1] | adj_gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| directional | ker1-avg_length | 2.56 | -1.05 | 47.22 | 10.03 | 7.68 | 94.27 | 8.73 |
| directional | ker0-avg_length | 5.14 | -4.91 | 25.30 | 17.17 | 6.73 | 86.17 | 11.63 |
| directional | ker1-med_length | -1.22 | -1.86 | 13.66 | 1.91 | 1.26 | 79.43 | 3.12 |
| directional | ker1-p90_length | 2.65 | 0.77 | 67.96 | 3.28 | 2.51 | 88.03 | 1.74 |
| directional | cok1-avg_length | -0.00 | 2.36 | 69.48 | -10.58 | -5.86 | 8.94 | -8.22 |
| directional | im1-avg_length | -1.02 | -0.48 | 44.35 | -0.41 | -0.54 | 43.00 | -0.06 |
| directional | ker1-num_bars | -7.33 | 4.42 | 79.43 | -17.60 | -5.42 | 13.49 | -9.85 |
| symmetric | ker1-avg_length | 2.50 | -3.73 | 17.03 | 17.91 | 10.22 | 96.29 | 13.95 |
| symmetric | ker0-avg_length | 4.93 | -7.58 | 18.55 | 24.69 | 10.58 | 88.70 | 18.17 |
| symmetric | ker1-med_length | -0.71 | -0.69 | 36.09 | 1.71 | 1.48 | 82.29 | 2.17 |
| symmetric | ker1-p90_length | 2.80 | -0.82 | 44.18 | 9.06 | 4.72 | 93.59 | 5.53 |
| symmetric | cok1-avg_length | 3.48 | -5.69 | 12.82 | 14.73 | 3.92 | 81.11 | 9.61 |
| symmetric | im1-avg_length | -2.85 | 1.52 | 77.07 | -4.75 | 0.38 | 55.99 | -1.14 |
| symmetric | ker1-num_bars | -8.40 | 8.42 | 91.57 | -30.43 | -11.33 | 6.24 | -19.74 |


Figure: `directional_named_pair.png`.

