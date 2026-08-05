# Day-2 re-scan: bar-length statistics (pre-treatment, per-patient)

Cohort: pre-treatment ROIs only, n = 593 ROIs / 25 NR + 37 R patients.
Effect size = Cliff's delta (+ve = responders higher); `pat_*` = per-patient medians (primary unit), `roi_*` = per-ROI. `pat_q` = Benjamini-Hochberg across each variant's family. **Descriptive and post hoc** -- see the caveat at the end.

## Why bar length rather than total persistence

- Median degree-1 image bar lifetime: **0.33 um** (mean 0.94 um) -- well below one cell diameter.
- Median fraction of image degree-1 bars clearing the 5 um threshold used for `n_bars_gt5um`: **3.1%** (median 1621 bars per ROI, of which 51 clear it).
- `im_dim1_total_persistence` sums *all* finite bars, so the statistic the Day-1 go/no-go rested on is dominated by sub-cell-scale features. `_dim1_feats` in `extract_features.py` compounds this: `total_persistence` and `persistent_entropy` use every bar while `n_bars_gt5um` uses only bars > 5 um, so they summarise different bar populations.

## Degree-1 mixing diagrams (ker / im / cok), avg_length


### variant: `raw`

| definition | pat_delta_cok1 | pat_delta_im1 | pat_delta_ker1 | pat_p_cok1 | pat_p_im1 | pat_p_ker1 |
| --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.377 | 0.310 | 0.345 | 0.013 | 0.040 | 0.023 |
| CD8_strict | 0.436 | 0.302 | 0.375 | 0.004 | 0.046 | 0.013 |
| CD8_extended | 0.362 | 0.286 | 0.293 | 0.017 | 0.058 | 0.053 |
| CD8_NK_only | 0.103 | 0.062 | 0.051 | 0.500 | 0.688 | 0.741 |
| CD4 | 0.083 | 0.066 | 0.070 | 0.586 | 0.667 | 0.646 |
| Bcell | 0.178 | 0.053 | 0.198 | 0.239 | 0.731 | 0.192 |
| Macrophage | 0.198 | 0.118 | 0.209 | 0.192 | 0.438 | 0.168 |
| Fibroblast | -0.057 | 0.012 | -0.038 | 0.709 | 0.943 | 0.807 |

### variant: `/med_nn`

| definition | pat_delta_cok1 | pat_delta_im1 | pat_delta_ker1 | pat_p_cok1 | pat_p_im1 | pat_p_ker1 |
| --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.308 | 0.291 | 0.317 | 0.042 | 0.054 | 0.036 |
| CD8_strict | 0.392 | 0.280 | 0.358 | 0.009 | 0.064 | 0.018 |
| CD8_extended | 0.293 | 0.276 | 0.254 | 0.053 | 0.068 | 0.093 |
| CD8_NK_only | 0.053 | -0.012 | 0.003 | 0.731 | 0.943 | 0.989 |
| CD4 | 0.049 | 0.008 | 0.010 | 0.752 | 0.966 | 0.954 |
| Bcell | 0.161 | -0.001 | 0.092 | 0.288 | 1.000 | 0.547 |
| Macrophage | 0.178 | 0.049 | 0.200 | 0.239 | 0.752 | 0.187 |
| Fibroblast | -0.118 | 0.005 | -0.109 | 0.438 | 0.977 | 0.473 |

### variant: `geom_adj`

| definition | pat_delta_cok1 | pat_delta_im1 | pat_delta_ker1 | pat_p_cok1 | pat_p_im1 | pat_p_ker1 |
| --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.274 | 0.155 | 0.248 | 0.071 | 0.308 | 0.102 |
| CD8_strict | 0.356 | 0.109 | 0.306 | 0.019 | 0.473 | 0.043 |
| CD8_extended | 0.267 | 0.137 | 0.142 | 0.078 | 0.366 | 0.351 |
| CD8_NK_only | -0.077 | -0.222 | -0.157 | 0.615 | 0.143 | 0.302 |
| CD4 | -0.042 | -0.103 | -0.055 | 0.785 | 0.500 | 0.720 |
| Bcell | 0.085 | -0.364 | -0.085 | 0.576 | 0.016 | 0.576 |
| Macrophage | 0.312 | -0.181 | 0.170 | 0.039 | 0.234 | 0.263 |
| Fibroblast | -0.163 | -0.265 | -0.152 | 0.282 | 0.080 | 0.315 |

## Full scan (all diagrams, all dimensions)

| definition | diagram | variant | roi_delta | roi_p | pat_delta | pat_p | pat_q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bcell | cod0 | /med_nn | -0.023 | 0.666 | 0.077 | 0.616 | 0.415 |
| Bcell | cod1 | /med_nn | -0.035 | 0.511 | 0.034 | 0.830 | 0.489 |
| Bcell | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.250 |
| Bcell | cok1 | /med_nn | 0.055 | 0.300 | 0.161 | 0.288 | 1.000 |
| Bcell | dom0 | /med_nn | -0.037 | 0.492 | 0.072 | 0.636 | 0.530 |
| Bcell | dom1 | /med_nn | -0.054 | 0.311 | 0.005 | 0.977 | 0.340 |
| Bcell | im0 | /med_nn | -0.023 | 0.666 | 0.077 | 0.616 | 0.250 |
| Bcell | im1 | /med_nn | -0.070 | 0.193 | -0.001 | 1.000 | 1.000 |
| Bcell | ker0 | /med_nn | 0.036 | 0.496 | 0.077 | 0.616 | 0.250 |
| Bcell | ker1 | /med_nn | 0.003 | 0.960 | 0.092 | 0.547 | 1.000 |
| CD4 | cod0 | /med_nn | 0.001 | 0.992 | 0.049 | 0.752 | 1.000 |
| CD4 | cod1 | /med_nn | 0.001 | 0.981 | 0.040 | 0.796 | 1.000 |
| CD4 | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.709 |
| CD4 | cok1 | /med_nn | 0.059 | 0.266 | 0.049 | 0.752 | 1.000 |
| CD4 | dom0 | /med_nn | -0.042 | 0.428 | 0.046 | 0.763 | 1.000 |
| CD4 | dom1 | /med_nn | -0.063 | 0.239 | 0.036 | 0.818 | 1.000 |
| CD4 | im0 | /med_nn | 0.001 | 0.992 | 0.049 | 0.752 | 1.000 |
| CD4 | im1 | /med_nn | -0.018 | 0.738 | 0.008 | 0.966 | 1.000 |
| CD4 | ker0 | /med_nn | 0.012 | 0.820 | -0.085 | 0.576 | 1.000 |
| CD4 | ker1 | /med_nn | 0.051 | 0.340 | 0.010 | 0.954 | 1.000 |
| CD8_NK_only | cod0 | /med_nn | -0.014 | 0.792 | 0.072 | 0.636 | 1.000 |
| CD8_NK_only | cod1 | /med_nn | -0.008 | 0.875 | 0.040 | 0.796 | 0.504 |
| CD8_NK_only | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_NK_only | cok1 | /med_nn | 0.007 | 0.889 | 0.053 | 0.731 | 1.000 |
| CD8_NK_only | dom0 | /med_nn | -0.033 | 0.538 | 0.085 | 0.576 | 1.000 |
| CD8_NK_only | dom1 | /med_nn | -0.027 | 0.613 | 0.036 | 0.818 | 0.324 |
| CD8_NK_only | im0 | /med_nn | -0.014 | 0.792 | 0.072 | 0.636 | 0.273 |
| CD8_NK_only | im1 | /med_nn | -0.032 | 0.546 | -0.012 | 0.943 | 0.273 |
| CD8_NK_only | ker0 | /med_nn | -0.078 | 0.145 | -0.116 | 0.447 | 0.273 |
| CD8_NK_only | ker1 | /med_nn | -0.031 | 0.564 | 0.003 | 0.989 | 0.273 |
| CD8_extended | cod0 | /med_nn | 0.234 | 0.000 | 0.269 | 0.075 | 0.273 |
| CD8_extended | cod1 | /med_nn | 0.243 | 0.000 | 0.323 | 0.033 | 0.962 |
| CD8_extended | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.273 |
| CD8_extended | cok1 | /med_nn | 0.233 | 0.000 | 0.293 | 0.053 | 1.000 |
| CD8_extended | dom0 | /med_nn | 0.225 | 0.000 | 0.271 | 0.073 | 0.751 |
| CD8_extended | dom1 | /med_nn | 0.230 | 0.000 | 0.217 | 0.151 | 0.962 |
| CD8_extended | im0 | /med_nn | 0.234 | 0.000 | 0.269 | 0.075 | 0.962 |
| CD8_extended | im1 | /med_nn | 0.208 | 0.000 | 0.276 | 0.068 | 0.935 |
| CD8_extended | ker0 | /med_nn | 0.038 | 0.481 | 0.144 | 0.344 | 0.962 |
| CD8_extended | ker1 | /med_nn | 0.132 | 0.013 | 0.254 | 0.093 | 0.962 |
| CD8_primary | cod0 | /med_nn | 0.259 | 0.000 | 0.321 | 0.034 | 0.349 |
| CD8_primary | cod1 | /med_nn | 0.238 | 0.000 | 0.349 | 0.021 | 0.200 |
| CD8_primary | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.200 |
| CD8_primary | cok1 | /med_nn | 0.245 | 0.000 | 0.308 | 0.042 | 1.000 |
| CD8_primary | dom0 | /med_nn | 0.249 | 0.000 | 0.304 | 0.044 | 0.648 |
| CD8_primary | dom1 | /med_nn | 0.235 | 0.000 | 0.289 | 0.056 | 0.200 |
| CD8_primary | im0 | /med_nn | 0.259 | 0.000 | 0.321 | 0.034 | 0.276 |
| CD8_primary | im1 | /med_nn | 0.210 | 0.000 | 0.291 | 0.054 | 0.200 |
| CD8_primary | ker0 | /med_nn | 0.061 | 0.252 | 0.170 | 0.263 | 0.276 |
| CD8_primary | ker1 | /med_nn | 0.164 | 0.002 | 0.317 | 0.036 | 0.276 |
| CD8_strict | cod0 | /med_nn | 0.275 | 0.000 | 0.338 | 0.025 | 0.962 |
| CD8_strict | cod1 | /med_nn | 0.251 | 0.000 | 0.351 | 0.020 | 0.962 |
| CD8_strict | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.638 |
| CD8_strict | cok1 | /med_nn | 0.306 | 0.000 | 0.392 | 0.009 | 1.000 |
| CD8_strict | dom0 | /med_nn | 0.252 | 0.000 | 0.302 | 0.046 | 0.962 |
| CD8_strict | dom1 | /med_nn | 0.231 | 0.000 | 0.286 | 0.058 | 0.962 |
| CD8_strict | im0 | /med_nn | 0.275 | 0.000 | 0.338 | 0.025 | 0.962 |
| CD8_strict | im1 | /med_nn | 0.196 | 0.000 | 0.280 | 0.064 | 0.962 |
| CD8_strict | ker0 | /med_nn | 0.107 | 0.045 | 0.191 | 0.207 | 0.962 |
| CD8_strict | ker1 | /med_nn | 0.217 | 0.000 | 0.358 | 0.018 | 0.962 |
| Fibroblast | cod0 | /med_nn | -0.010 | 0.846 | 0.064 | 0.677 | 0.873 |
| Fibroblast | cod1 | /med_nn | -0.030 | 0.570 | 0.042 | 0.785 | 0.395 |
| Fibroblast | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.501 |
| Fibroblast | cok1 | /med_nn | -0.041 | 0.445 | -0.118 | 0.438 | 1.000 |
| Fibroblast | dom0 | /med_nn | -0.049 | 0.365 | 0.025 | 0.875 | 0.613 |
| Fibroblast | dom1 | /med_nn | -0.070 | 0.193 | -0.038 | 0.807 | 0.489 |
| Fibroblast | im0 | /med_nn | -0.010 | 0.846 | 0.064 | 0.677 | 1.000 |
| Fibroblast | im1 | /med_nn | -0.022 | 0.681 | 0.005 | 0.977 | 0.324 |
| Fibroblast | ker0 | /med_nn | -0.048 | 0.372 | -0.109 | 0.473 | 1.000 |
| Fibroblast | ker1 | /med_nn | -0.043 | 0.420 | -0.109 | 0.473 | 0.250 |
| Macrophage | cod0 | /med_nn | 0.007 | 0.899 | 0.131 | 0.389 | 0.425 |
| Macrophage | cod1 | /med_nn | 0.010 | 0.853 | 0.114 | 0.456 | 0.320 |
| Macrophage | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.938 |
| Macrophage | cok1 | /med_nn | 0.171 | 0.001 | 0.178 | 0.239 | 1.000 |
| Macrophage | dom0 | /med_nn | -0.045 | 0.405 | 0.070 | 0.646 | 0.504 |
| Macrophage | dom1 | /med_nn | -0.045 | 0.402 | 0.001 | 1.000 | 0.504 |
| Macrophage | im0 | /med_nn | 0.007 | 0.899 | 0.131 | 0.389 | 0.265 |
| Macrophage | im1 | /med_nn | -0.032 | 0.548 | 0.049 | 0.752 | 0.366 |
| Macrophage | ker0 | /med_nn | 0.112 | 0.037 | 0.042 | 0.785 | 0.265 |
| Macrophage | ker1 | /med_nn | 0.146 | 0.006 | 0.200 | 0.187 | 0.542 |
| Bcell | cod0 | geom_adj | -0.228 | 0.000 | -0.254 | 0.093 | 0.250 |
| Bcell | cod1 | geom_adj | -0.227 | 0.000 | -0.258 | 0.088 | 0.402 |
| Bcell | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.516 |
| Bcell | cok1 | geom_adj | 0.025 | 0.640 | 0.085 | 0.576 | 1.000 |
| Bcell | dom0 | geom_adj | -0.201 | 0.000 | -0.269 | 0.075 | 0.250 |
| Bcell | dom1 | geom_adj | -0.238 | 0.000 | -0.278 | 0.066 | 0.504 |
| Bcell | im0 | geom_adj | -0.228 | 0.000 | -0.254 | 0.093 | 0.265 |
| Bcell | im1 | geom_adj | -0.278 | 0.000 | -0.364 | 0.016 | 1.000 |
| Bcell | ker0 | geom_adj | 0.079 | 0.138 | 0.114 | 0.456 | 1.000 |
| Bcell | ker1 | geom_adj | -0.058 | 0.276 | -0.085 | 0.576 | 1.000 |
| CD4 | cod0 | geom_adj | -0.251 | 0.000 | -0.200 | 0.187 | 0.973 |
| CD4 | cod1 | geom_adj | -0.232 | 0.000 | -0.209 | 0.168 | 1.000 |
| CD4 | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD4 | cok1 | geom_adj | -0.059 | 0.272 | -0.042 | 0.785 | 1.000 |
| CD4 | dom0 | geom_adj | -0.253 | 0.000 | -0.196 | 0.197 | 0.973 |
| CD4 | dom1 | geom_adj | -0.253 | 0.000 | -0.189 | 0.212 | 1.000 |
| CD4 | im0 | geom_adj | -0.251 | 0.000 | -0.200 | 0.187 | 0.598 |
| CD4 | im1 | geom_adj | -0.187 | 0.001 | -0.103 | 0.500 | 1.000 |
| CD4 | ker0 | geom_adj | 0.078 | 0.143 | 0.152 | 0.315 | 1.000 |
| CD4 | ker1 | geom_adj | -0.075 | 0.161 | -0.055 | 0.720 | 1.000 |
| CD8_NK_only | cod0 | geom_adj | -0.264 | 0.000 | -0.308 | 0.042 | 1.000 |
| CD8_NK_only | cod1 | geom_adj | -0.242 | 0.000 | -0.200 | 0.187 | 0.273 |
| CD8_NK_only | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_NK_only | cok1 | geom_adj | -0.060 | 0.258 | -0.077 | 0.616 | 0.273 |
| CD8_NK_only | dom0 | geom_adj | -0.276 | 0.000 | -0.332 | 0.028 | 1.000 |
| CD8_NK_only | dom1 | geom_adj | -0.261 | 0.000 | -0.274 | 0.071 | 0.273 |
| CD8_NK_only | im0 | geom_adj | -0.264 | 0.000 | -0.308 | 0.042 | 1.000 |
| CD8_NK_only | im1 | geom_adj | -0.205 | 0.000 | -0.222 | 0.143 | 0.636 |
| CD8_NK_only | ker0 | geom_adj | -0.023 | 0.669 | -0.001 | 1.000 | 1.000 |
| CD8_NK_only | ker1 | geom_adj | -0.113 | 0.035 | -0.157 | 0.301 | 0.273 |
| CD8_extended | cod0 | geom_adj | 0.250 | 0.000 | 0.334 | 0.027 | 0.273 |
| CD8_extended | cod1 | geom_adj | 0.178 | 0.001 | 0.239 | 0.115 | 0.962 |
| CD8_extended | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.273 |
| CD8_extended | cok1 | geom_adj | 0.175 | 0.001 | 0.267 | 0.078 | 0.962 |
| CD8_extended | dom0 | geom_adj | 0.174 | 0.001 | 0.202 | 0.182 | 0.273 |
| CD8_extended | dom1 | geom_adj | 0.135 | 0.012 | 0.168 | 0.269 | 1.000 |
| CD8_extended | im0 | geom_adj | 0.250 | 0.000 | 0.334 | 0.027 | 0.273 |
| CD8_extended | im1 | geom_adj | 0.153 | 0.004 | 0.137 | 0.366 | 0.962 |
| CD8_extended | ker0 | geom_adj | -0.018 | 0.734 | 0.094 | 0.537 | 1.000 |
| CD8_extended | ker1 | geom_adj | 0.050 | 0.353 | 0.142 | 0.351 | 0.529 |
| CD8_primary | cod0 | geom_adj | 0.242 | 0.000 | 0.345 | 0.022 | 0.465 |
| CD8_primary | cod1 | geom_adj | 0.145 | 0.007 | 0.226 | 0.136 | 0.200 |
| CD8_primary | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.276 |
| CD8_primary | cok1 | geom_adj | 0.185 | 0.001 | 0.274 | 0.071 | 0.279 |
| CD8_primary | dom0 | geom_adj | 0.174 | 0.001 | 0.202 | 0.182 | 0.465 |
| CD8_primary | dom1 | geom_adj | 0.128 | 0.017 | 0.172 | 0.257 | 0.276 |
| CD8_primary | im0 | geom_adj | 0.242 | 0.000 | 0.345 | 0.022 | 0.276 |
| CD8_primary | im1 | geom_adj | 0.101 | 0.060 | 0.155 | 0.308 | 0.499 |
| CD8_primary | ker0 | geom_adj | -0.027 | 0.614 | 0.094 | 0.537 | 1.000 |
| CD8_primary | ker1 | geom_adj | 0.097 | 0.069 | 0.248 | 0.102 | 0.200 |
| CD8_strict | cod0 | geom_adj | 0.240 | 0.000 | 0.338 | 0.025 | 0.962 |
| CD8_strict | cod1 | geom_adj | 0.160 | 0.003 | 0.211 | 0.164 | 0.962 |
| CD8_strict | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_strict | cok1 | geom_adj | 0.248 | 0.000 | 0.356 | 0.019 | 0.962 |
| CD8_strict | dom0 | geom_adj | 0.167 | 0.002 | 0.200 | 0.187 | 0.962 |
| CD8_strict | dom1 | geom_adj | 0.130 | 0.015 | 0.148 | 0.329 | 0.962 |
| CD8_strict | im0 | geom_adj | 0.240 | 0.000 | 0.338 | 0.025 | 0.529 |
| CD8_strict | im1 | geom_adj | 0.080 | 0.134 | 0.109 | 0.473 | 0.962 |
| CD8_strict | ker0 | geom_adj | 0.006 | 0.915 | 0.142 | 0.351 | 1.000 |
| CD8_strict | ker1 | geom_adj | 0.150 | 0.005 | 0.306 | 0.043 | 0.962 |
| Fibroblast | cod0 | geom_adj | -0.139 | 0.009 | -0.021 | 0.897 | 1.000 |
| Fibroblast | cod1 | geom_adj | -0.215 | 0.000 | -0.174 | 0.251 | 0.756 |
| Fibroblast | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.504 |
| Fibroblast | cok1 | geom_adj | -0.089 | 0.097 | -0.163 | 0.282 | 0.320 |
| Fibroblast | dom0 | geom_adj | -0.160 | 0.003 | -0.057 | 0.709 | 1.000 |
| Fibroblast | dom1 | geom_adj | -0.210 | 0.000 | -0.155 | 0.308 | 0.479 |
| Fibroblast | im0 | geom_adj | -0.139 | 0.009 | -0.021 | 0.897 | 0.504 |
| Fibroblast | im1 | geom_adj | -0.235 | 0.000 | -0.265 | 0.080 | 0.651 |
| Fibroblast | ker0 | geom_adj | -0.088 | 0.099 | -0.122 | 0.422 | 1.000 |
| Fibroblast | ker1 | geom_adj | -0.104 | 0.052 | -0.152 | 0.315 | 0.743 |
| Macrophage | cod0 | geom_adj | -0.125 | 0.020 | -0.008 | 0.966 | 0.415 |
| Macrophage | cod1 | geom_adj | -0.172 | 0.001 | -0.081 | 0.596 | 0.415 |
| Macrophage | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.446 |
| Macrophage | cok1 | geom_adj | 0.196 | 0.000 | 0.312 | 0.039 | 0.415 |
| Macrophage | dom0 | geom_adj | -0.239 | 0.000 | -0.178 | 0.239 | 0.415 |
| Macrophage | dom1 | geom_adj | -0.228 | 0.000 | -0.230 | 0.128 | 0.409 |
| Macrophage | im0 | geom_adj | -0.125 | 0.020 | -0.008 | 0.966 | 0.873 |
| Macrophage | im1 | geom_adj | -0.241 | 0.000 | -0.180 | 0.234 | 0.716 |
| Macrophage | ker0 | geom_adj | 0.139 | 0.009 | 0.291 | 0.054 | 1.000 |
| Macrophage | ker1 | geom_adj | 0.124 | 0.020 | 0.170 | 0.263 | 0.320 |
| Bcell | cod0 | raw | 0.047 | 0.383 | 0.077 | 0.616 | 1.000 |
| Bcell | cod1 | raw | 0.031 | 0.566 | 0.055 | 0.720 | 0.320 |
| Bcell | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.664 |
| Bcell | cok1 | raw | 0.109 | 0.041 | 0.178 | 0.239 | 1.000 |
| Bcell | dom0 | raw | 0.030 | 0.575 | 0.059 | 0.698 | 0.415 |
| Bcell | dom1 | raw | 0.007 | 0.895 | 0.005 | 0.977 | 1.000 |
| Bcell | im0 | raw | 0.047 | 0.383 | 0.077 | 0.616 | 0.415 |
| Bcell | im1 | raw | 0.015 | 0.787 | 0.053 | 0.731 | 1.000 |
| Bcell | ker0 | raw | 0.089 | 0.094 | 0.114 | 0.456 | 0.716 |
| Bcell | ker1 | raw | 0.075 | 0.160 | 0.198 | 0.192 | 1.000 |
| CD4 | cod0 | raw | 0.050 | 0.354 | 0.096 | 0.528 | 1.000 |
| CD4 | cod1 | raw | 0.054 | 0.309 | 0.072 | 0.636 | 0.795 |
| CD4 | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD4 | cok1 | raw | 0.098 | 0.068 | 0.083 | 0.586 | 1.000 |
| CD4 | dom0 | raw | 0.012 | 0.828 | 0.064 | 0.677 | 1.000 |
| CD4 | dom1 | raw | -0.004 | 0.948 | 0.059 | 0.698 | 1.000 |
| CD4 | im0 | raw | 0.050 | 0.354 | 0.096 | 0.528 | 1.000 |
| CD4 | im1 | raw | 0.048 | 0.369 | 0.066 | 0.667 | 1.000 |
| CD4 | ker0 | raw | 0.036 | 0.499 | -0.088 | 0.566 | 1.000 |
| CD4 | ker1 | raw | 0.099 | 0.064 | 0.070 | 0.646 | 1.000 |
| CD8_NK_only | cod0 | raw | 0.046 | 0.387 | 0.103 | 0.500 | 1.000 |
| CD8_NK_only | cod1 | raw | 0.050 | 0.350 | 0.070 | 0.646 | 0.273 |
| CD8_NK_only | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_NK_only | cok1 | raw | 0.063 | 0.238 | 0.103 | 0.500 | 0.273 |
| CD8_NK_only | dom0 | raw | 0.027 | 0.619 | 0.090 | 0.556 | 1.000 |
| CD8_NK_only | dom1 | raw | 0.026 | 0.623 | 0.040 | 0.796 | 0.273 |
| CD8_NK_only | im0 | raw | 0.046 | 0.387 | 0.103 | 0.500 | 0.273 |
| CD8_NK_only | im1 | raw | 0.046 | 0.387 | 0.062 | 0.688 | 0.273 |
| CD8_NK_only | ker0 | raw | -0.027 | 0.615 | -0.040 | 0.796 | 0.916 |
| CD8_NK_only | ker1 | raw | 0.028 | 0.597 | 0.051 | 0.741 | 0.273 |
| CD8_extended | cod0 | raw | 0.306 | 0.000 | 0.224 | 0.139 | 1.000 |
| CD8_extended | cod1 | raw | 0.306 | 0.000 | 0.308 | 0.042 | 0.962 |
| CD8_extended | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.273 |
| CD8_extended | cok1 | raw | 0.311 | 0.000 | 0.362 | 0.017 | 0.851 |
| CD8_extended | dom0 | raw | 0.284 | 0.000 | 0.252 | 0.096 | 0.273 |
| CD8_extended | dom1 | raw | 0.278 | 0.000 | 0.312 | 0.039 | 0.851 |
| CD8_extended | im0 | raw | 0.306 | 0.000 | 0.224 | 0.139 | 1.000 |
| CD8_extended | im1 | raw | 0.275 | 0.000 | 0.286 | 0.058 | 0.962 |
| CD8_extended | ker0 | raw | 0.091 | 0.087 | 0.174 | 0.251 | 0.962 |
| CD8_extended | ker1 | raw | 0.220 | 0.000 | 0.293 | 0.053 | 0.499 |
| CD8_primary | cod0 | raw | 0.318 | 0.000 | 0.263 | 0.083 | 1.000 |
| CD8_primary | cod1 | raw | 0.305 | 0.000 | 0.345 | 0.022 | 0.200 |
| CD8_primary | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.276 |
| CD8_primary | cok1 | raw | 0.316 | 0.000 | 0.377 | 0.013 | 0.314 |
| CD8_primary | dom0 | raw | 0.299 | 0.000 | 0.278 | 0.066 | 0.276 |
| CD8_primary | dom1 | raw | 0.293 | 0.000 | 0.360 | 0.017 | 0.314 |
| CD8_primary | im0 | raw | 0.318 | 0.000 | 0.263 | 0.083 | 0.276 |
| CD8_primary | im1 | raw | 0.273 | 0.000 | 0.310 | 0.040 | 0.200 |
| CD8_primary | ker0 | raw | 0.113 | 0.035 | 0.213 | 0.160 | 0.499 |
| CD8_primary | ker1 | raw | 0.247 | 0.000 | 0.345 | 0.022 | 0.200 |
| CD8_strict | cod0 | raw | 0.328 | 0.000 | 0.282 | 0.062 | 1.000 |
| CD8_strict | cod1 | raw | 0.322 | 0.000 | 0.358 | 0.018 | 0.962 |
| CD8_strict | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.962 |
| CD8_strict | cok1 | raw | 0.373 | 0.000 | 0.436 | 0.004 | 0.962 |
| CD8_strict | dom0 | raw | 0.303 | 0.000 | 0.282 | 0.062 | 0.962 |
| CD8_strict | dom1 | raw | 0.294 | 0.000 | 0.351 | 0.020 | 0.962 |
| CD8_strict | im0 | raw | 0.328 | 0.000 | 0.282 | 0.062 | 0.962 |
| CD8_strict | im1 | raw | 0.262 | 0.000 | 0.302 | 0.046 | 0.962 |
| CD8_strict | ker0 | raw | 0.154 | 0.004 | 0.211 | 0.164 | 0.962 |
| CD8_strict | ker1 | raw | 0.291 | 0.000 | 0.375 | 0.013 | 0.962 |
| Fibroblast | cod0 | raw | 0.074 | 0.168 | 0.036 | 0.818 | 1.000 |
| Fibroblast | cod1 | raw | 0.059 | 0.267 | 0.044 | 0.774 | 0.265 |
| Fibroblast | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.320 |
| Fibroblast | cok1 | raw | 0.033 | 0.538 | -0.057 | 0.709 | 0.324 |
| Fibroblast | dom0 | raw | 0.041 | 0.443 | 0.021 | 0.897 | 0.489 |
| Fibroblast | dom1 | raw | 0.014 | 0.799 | 0.034 | 0.830 | 0.324 |
| Fibroblast | im0 | raw | 0.074 | 0.168 | 0.036 | 0.818 | 0.479 |
| Fibroblast | im1 | raw | 0.065 | 0.228 | 0.012 | 0.943 | 0.320 |
| Fibroblast | ker0 | raw | 0.001 | 0.986 | -0.105 | 0.491 | 0.311 |
| Fibroblast | ker1 | raw | 0.021 | 0.698 | -0.038 | 0.807 | 0.743 |
| Macrophage | cod0 | raw | 0.087 | 0.104 | 0.142 | 0.351 | 1.000 |
| Macrophage | cod1 | raw | 0.093 | 0.081 | 0.129 | 0.397 | 0.769 |
| Macrophage | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.690 |
| Macrophage | cok1 | raw | 0.241 | 0.000 | 0.198 | 0.192 | 0.250 |
| Macrophage | dom0 | raw | 0.034 | 0.524 | 0.070 | 0.646 | 0.415 |
| Macrophage | dom1 | raw | 0.031 | 0.566 | 0.038 | 0.807 | 0.250 |
| Macrophage | im0 | raw | 0.087 | 0.104 | 0.142 | 0.351 | 0.250 |
| Macrophage | im1 | raw | 0.050 | 0.346 | 0.118 | 0.438 | 0.489 |
| Macrophage | ker0 | raw | 0.148 | 0.005 | 0.062 | 0.688 | 1.000 |
| Macrophage | ker1 | raw | 0.214 | 0.000 | 0.209 | 0.168 | 0.530 |

## Read

- CD8 aggregate definitions (primary/strict/extended), `/med_nn`: median per-patient delta = **+0.293**, 4/9 tests at p < 0.05.
- Control species (CD4, Bcell, Macrophage, Fibroblast), same variant: median delta = **+0.029**, 0/12 at p < 0.05.
- After BH across each variant's family, 0 of 240 tests survive q < 0.05.

**Caveat that matters.** The effect attenuates under `geom_adj`, and the same CD8-specific pattern appears in `dom0`/`dom1` -- the *monochromatic* complex, which contains no cross-species mixing information at all. That is consistent with 'responder tissue is more dispersed' rather than 'T cells infiltrate'. Regression adjustment cannot separate these, because density, dispersion and arrangement are collinear here. Only a null that holds cell positions fixed and shuffles labels can: see `perm_null.py` / `output/PERMUTATION_NULL.md`.

Script: `day2_length_features.py` (reads output/stats/*.parquet, no recomputation).
