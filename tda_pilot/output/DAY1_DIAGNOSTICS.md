# Day-1 descriptive diagnostics (no models)

Effect size = Cliff's delta on per-ROI values (+ve = responders higher; |δ|≈0.1 small, 0.3 medium, 0.5 large).

## (2) Separation across all 7 definitions (Cliff's δ, R vs NR)

| definition | cok_dim1_total_persistence | im_dim1_n_bars_gt5um | im_dim1_total_persistence | ker_dim1_n_bars_gt5um | ker_dim1_total_persistence |
| --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.003 | -0.104 | -0.155 | -0.007 | -0.018 |
| CD8_strict | -0.017 | -0.093 | -0.155 | -0.027 | -0.041 |
| CD8_extended | -0.006 | -0.111 | -0.156 | -0.015 | -0.028 |
| CD8_NK_only | -0.025 | -0.057 | -0.020 | -0.050 | -0.060 |
| CD4 | 0.037 | 0.008 | 0.017 | 0.022 | 0.020 |
| Bcell | -0.058 | -0.106 | -0.046 | -0.098 | -0.106 |
| Macrophage | -0.023 | -0.032 | -0.002 | -0.072 | -0.070 |
| Fibroblast | -0.027 | -0.006 | 0.016 | -0.063 | -0.048 |

## Per-definition x feature detail

| definition | feature | NR_med | R_med | ratio | cliffs_delta |
| --- | --- | --- | --- | --- | --- |
| CD8_primary | ker_dim1_total_persistence | 1060.360 | 970.514 | 0.915 | -0.018 |
| CD8_primary | im_dim1_total_persistence | 1842.864 | 1443.432 | 0.783 | -0.155 |
| CD8_primary | cok_dim1_total_persistence | 862.462 | 847.538 | 0.983 | 0.003 |
| CD8_primary | im_dim1_n_bars_gt5um | 55.000 | 49.500 | 0.900 | -0.104 |
| CD8_primary | ker_dim1_n_bars_gt5um | 48.000 | 47.000 | 0.979 | -0.007 |
| CD8_strict | ker_dim1_total_persistence | 1001.638 | 902.815 | 0.901 | -0.041 |
| CD8_strict | im_dim1_total_persistence | 1810.533 | 1405.295 | 0.776 | -0.155 |
| CD8_strict | cok_dim1_total_persistence | 807.642 | 780.812 | 0.967 | -0.017 |
| CD8_strict | im_dim1_n_bars_gt5um | 56.000 | 52.500 | 0.938 | -0.093 |
| CD8_strict | ker_dim1_n_bars_gt5um | 46.000 | 42.500 | 0.924 | -0.027 |
| CD8_extended | ker_dim1_total_persistence | 1202.717 | 1090.875 | 0.907 | -0.028 |
| CD8_extended | im_dim1_total_persistence | 2004.716 | 1588.466 | 0.792 | -0.156 |
| CD8_extended | cok_dim1_total_persistence | 968.334 | 920.647 | 0.951 | -0.006 |
| CD8_extended | im_dim1_n_bars_gt5um | 66.000 | 57.000 | 0.864 | -0.111 |
| CD8_extended | ker_dim1_n_bars_gt5um | 57.500 | 52.000 | 0.904 | -0.015 |
| CD8_NK_only | ker_dim1_total_persistence | 386.736 | 346.616 | 0.896 | -0.060 |
| CD8_NK_only | im_dim1_total_persistence | 1188.861 | 1152.314 | 0.969 | -0.020 |
| CD8_NK_only | cok_dim1_total_persistence | 313.548 | 314.898 | 1.004 | -0.025 |
| CD8_NK_only | im_dim1_n_bars_gt5um | 49.000 | 43.000 | 0.878 | -0.057 |
| CD8_NK_only | ker_dim1_n_bars_gt5um | 20.000 | 17.000 | 0.850 | -0.050 |
| CD4 | ker_dim1_total_persistence | 434.499 | 471.012 | 1.084 | 0.020 |
| CD4 | im_dim1_total_persistence | 1249.471 | 1286.491 | 1.030 | 0.017 |
| CD4 | cok_dim1_total_persistence | 391.199 | 434.040 | 1.110 | 0.037 |
| CD4 | im_dim1_n_bars_gt5um | 48.000 | 50.000 | 1.042 | 0.008 |
| CD4 | ker_dim1_n_bars_gt5um | 20.000 | 22.000 | 1.100 | 0.022 |
| Bcell | ker_dim1_total_persistence | 459.331 | 400.130 | 0.871 | -0.106 |
| Bcell | im_dim1_total_persistence | 1410.544 | 1272.933 | 0.902 | -0.046 |
| Bcell | cok_dim1_total_persistence | 383.376 | 367.382 | 0.958 | -0.058 |
| Bcell | im_dim1_n_bars_gt5um | 59.500 | 51.000 | 0.857 | -0.106 |
| Bcell | ker_dim1_n_bars_gt5um | 21.000 | 19.000 | 0.905 | -0.098 |
| Macrophage | ker_dim1_total_persistence | 649.953 | 548.823 | 0.844 | -0.070 |
| Macrophage | im_dim1_total_persistence | 1221.200 | 1301.611 | 1.066 | -0.002 |
| Macrophage | cok_dim1_total_persistence | 549.461 | 521.683 | 0.949 | -0.023 |
| Macrophage | im_dim1_n_bars_gt5um | 50.000 | 52.000 | 1.040 | -0.032 |
| Macrophage | ker_dim1_n_bars_gt5um | 31.000 | 27.000 | 0.871 | -0.072 |
| Fibroblast | ker_dim1_total_persistence | 626.648 | 564.138 | 0.900 | -0.048 |
| Fibroblast | im_dim1_total_persistence | 1510.465 | 1510.385 | 1.000 | 0.016 |
| Fibroblast | cok_dim1_total_persistence | 614.973 | 588.783 | 0.957 | -0.027 |
| Fibroblast | im_dim1_n_bars_gt5um | 61.500 | 58.000 | 0.943 | -0.006 |
| Fibroblast | ker_dim1_n_bars_gt5um | 29.000 | 27.000 | 0.931 | -0.063 |

## (3) Image total persistence: raw vs size-adjusted (CD8 definitions)

| definition | raw_cliffs_delta | size_adj_cliffs_delta | resid_med_NR | resid_med_R |
| --- | --- | --- | --- | --- |
| CD8_primary | -0.155 | -0.143 | -27.710 | -105.199 |
| CD8_strict | -0.155 | -0.124 | -18.003 | -85.084 |
| CD8_extended | -0.156 | -0.137 | -9.986 | -105.027 |
| CD8_NK_only | -0.020 | -0.121 | 4.639 | -78.860 |

Figures: `output/diag_separation.png`, `output/diag_image_residual.png`.

