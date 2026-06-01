# Day-1 descriptive diagnostics (no models)

Effect size = Cliff's delta on per-ROI values (+ve = responders higher; |δ|≈0.1 small, 0.3 medium, 0.5 large).

## (2) Separation across all 7 definitions (Cliff's δ, R vs NR)

| definition | cok_dim1_total_persistence | im_dim1_n_features | im_dim1_total_persistence | ker_dim1_n_features | ker_dim1_total_persistence |
| --- | --- | --- | --- | --- | --- |
| CD8_primary | -0.011 | -0.121 | -0.168 | -0.022 | -0.034 |
| CD8_strict | -0.032 | -0.108 | -0.167 | -0.044 | -0.059 |
| CD8_extended | -0.020 | -0.128 | -0.170 | -0.029 | -0.043 |
| CD8_NK_only | -0.035 | -0.068 | -0.026 | -0.063 | -0.073 |
| CD4 | 0.033 | -0.002 | 0.011 | 0.016 | 0.014 |
| Bcell | -0.069 | -0.117 | -0.052 | -0.111 | -0.120 |
| Macrophage | -0.043 | -0.038 | -0.005 | -0.093 | -0.091 |
| Fibroblast | -0.043 | -0.013 | 0.014 | -0.082 | -0.066 |

## Per-definition x feature detail

| definition | feature | NR_med | R_med | ratio | cliffs_delta |
| --- | --- | --- | --- | --- | --- |
| CD8_primary | ker_dim1_total_persistence | 1060.360 | 955.436 | 0.901 | -0.034 |
| CD8_primary | im_dim1_total_persistence | 1842.864 | 1403.422 | 0.762 | -0.168 |
| CD8_primary | cok_dim1_total_persistence | 862.462 | 834.693 | 0.968 | -0.011 |
| CD8_primary | im_dim1_n_features | 55.000 | 49.000 | 0.891 | -0.121 |
| CD8_primary | ker_dim1_n_features | 48.000 | 46.000 | 0.958 | -0.022 |
| CD8_strict | ker_dim1_total_persistence | 1001.638 | 883.962 | 0.883 | -0.059 |
| CD8_strict | im_dim1_total_persistence | 1810.533 | 1399.496 | 0.773 | -0.167 |
| CD8_strict | cok_dim1_total_persistence | 807.642 | 758.289 | 0.939 | -0.032 |
| CD8_strict | im_dim1_n_features | 56.000 | 51.000 | 0.911 | -0.108 |
| CD8_strict | ker_dim1_n_features | 46.000 | 42.000 | 0.913 | -0.044 |
| CD8_extended | ker_dim1_total_persistence | 1202.717 | 1068.077 | 0.888 | -0.043 |
| CD8_extended | im_dim1_total_persistence | 2004.716 | 1560.340 | 0.778 | -0.170 |
| CD8_extended | cok_dim1_total_persistence | 968.334 | 916.740 | 0.947 | -0.020 |
| CD8_extended | im_dim1_n_features | 66.000 | 57.000 | 0.864 | -0.128 |
| CD8_extended | ker_dim1_n_features | 57.500 | 51.000 | 0.887 | -0.029 |
| CD8_NK_only | ker_dim1_total_persistence | 386.736 | 335.394 | 0.867 | -0.073 |
| CD8_NK_only | im_dim1_total_persistence | 1188.861 | 1134.040 | 0.954 | -0.026 |
| CD8_NK_only | cok_dim1_total_persistence | 313.548 | 312.529 | 0.997 | -0.035 |
| CD8_NK_only | im_dim1_n_features | 49.000 | 42.000 | 0.857 | -0.068 |
| CD8_NK_only | ker_dim1_n_features | 20.000 | 16.000 | 0.800 | -0.063 |
| CD4 | ker_dim1_total_persistence | 434.499 | 462.848 | 1.065 | 0.014 |
| CD4 | im_dim1_total_persistence | 1249.471 | 1263.837 | 1.011 | 0.011 |
| CD4 | cok_dim1_total_persistence | 391.199 | 427.960 | 1.094 | 0.033 |
| CD4 | im_dim1_n_features | 48.000 | 49.000 | 1.021 | -0.002 |
| CD4 | ker_dim1_n_features | 20.000 | 22.000 | 1.100 | 0.016 |
| Bcell | ker_dim1_total_persistence | 459.331 | 396.471 | 0.863 | -0.120 |
| Bcell | im_dim1_total_persistence | 1410.544 | 1266.093 | 0.898 | -0.052 |
| Bcell | cok_dim1_total_persistence | 383.376 | 361.225 | 0.942 | -0.069 |
| Bcell | im_dim1_n_features | 59.500 | 50.000 | 0.840 | -0.117 |
| Bcell | ker_dim1_n_features | 21.000 | 18.000 | 0.857 | -0.111 |
| Macrophage | ker_dim1_total_persistence | 649.953 | 541.287 | 0.833 | -0.091 |
| Macrophage | im_dim1_total_persistence | 1221.200 | 1299.224 | 1.064 | -0.005 |
| Macrophage | cok_dim1_total_persistence | 549.461 | 516.980 | 0.941 | -0.043 |
| Macrophage | im_dim1_n_features | 50.000 | 52.000 | 1.040 | -0.038 |
| Macrophage | ker_dim1_n_features | 31.000 | 26.000 | 0.839 | -0.093 |
| Fibroblast | ker_dim1_total_persistence | 626.648 | 555.974 | 0.887 | -0.066 |
| Fibroblast | im_dim1_total_persistence | 1510.465 | 1510.991 | 1.000 | 0.014 |
| Fibroblast | cok_dim1_total_persistence | 614.973 | 583.466 | 0.949 | -0.043 |
| Fibroblast | im_dim1_n_features | 61.500 | 57.000 | 0.927 | -0.013 |
| Fibroblast | ker_dim1_n_features | 29.000 | 26.000 | 0.897 | -0.082 |

## (3) Image total persistence: raw vs size-adjusted (CD8 definitions)

| definition | raw_cliffs_delta | size_adj_cliffs_delta | resid_med_NR | resid_med_R |
| --- | --- | --- | --- | --- |
| CD8_primary | -0.168 | -0.138 | -30.421 | -102.054 |
| CD8_strict | -0.167 | -0.118 | -12.886 | -83.072 |
| CD8_extended | -0.170 | -0.133 | -8.572 | -104.754 |
| CD8_NK_only | -0.026 | -0.132 | 6.579 | -76.447 |

Figures: `output/diag_separation.png`, `output/diag_image_residual.png`.

