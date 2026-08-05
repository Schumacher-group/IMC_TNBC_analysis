# Day-2 re-scan: bar-length statistics (pre-treatment, per-patient)

Cohort: pre-treatment ROIs only, n = 582 ROIs / 25 NR + 37 R patients.
Effect size = Cliff's delta (+ve = responders higher); `pat_*` = per-patient medians (primary unit), `roi_*` = per-ROI. `pat_q` = Benjamini-Hochberg across each variant's family. **Descriptive and post hoc** -- see the caveat at the end.

## Why bar length rather than total persistence

- Median degree-1 image bar lifetime: **0.33 um** (mean 0.94 um) -- well below one cell diameter.
- Median fraction of image degree-1 bars clearing the 5 um threshold used for `n_bars_gt5um`: **3.1%** (median 1606 bars per ROI, of which 50 clear it).
- `im_dim1_total_persistence` sums *all* finite bars, so the statistic the Day-1 go/no-go rested on is dominated by sub-cell-scale features. `_dim1_feats` in `extract_features.py` compounds this: `total_persistence` and `persistent_entropy` use every bar while `n_bars_gt5um` uses only bars > 5 um, so they summarise different bar populations.

## Degree-1 mixing diagrams (ker / im / cok), avg_length


### variant: `raw`

| definition | pat_delta_cok1 | pat_delta_im1 | pat_delta_ker1 | pat_p_cok1 | pat_p_im1 | pat_p_ker1 |
| --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.395 | 0.310 | 0.366 | 0.009 | 0.040 | 0.015 |
| CD8_strict | 0.444 | 0.297 | 0.399 | 0.003 | 0.049 | 0.008 |
| CD8_extended | 0.377 | 0.286 | 0.308 | 0.013 | 0.058 | 0.042 |
| CD8_NK_only | 0.103 | 0.051 | 0.053 | 0.500 | 0.741 | 0.731 |
| CD4 | 0.083 | 0.055 | 0.068 | 0.586 | 0.720 | 0.656 |
| Bcell | 0.181 | 0.042 | 0.196 | 0.234 | 0.785 | 0.197 |
| Macrophage | 0.198 | 0.109 | 0.209 | 0.192 | 0.473 | 0.168 |
| Fibroblast | -0.068 | 0.005 | -0.042 | 0.656 | 0.977 | 0.785 |

### variant: `/med_nn`

| definition | pat_delta_cok1 | pat_delta_im1 | pat_delta_ker1 | pat_p_cok1 | pat_p_im1 | pat_p_ker1 |
| --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.334 | 0.291 | 0.336 | 0.027 | 0.054 | 0.026 |
| CD8_strict | 0.414 | 0.282 | 0.375 | 0.006 | 0.062 | 0.013 |
| CD8_extended | 0.319 | 0.278 | 0.265 | 0.035 | 0.066 | 0.080 |
| CD8_NK_only | 0.057 | -0.025 | 0.003 | 0.709 | 0.875 | 0.989 |
| CD4 | 0.049 | -0.010 | 0.008 | 0.752 | 0.954 | 0.966 |
| Bcell | 0.161 | -0.010 | 0.092 | 0.288 | 0.954 | 0.547 |
| Macrophage | 0.183 | 0.040 | 0.200 | 0.228 | 0.796 | 0.187 |
| Fibroblast | -0.135 | -0.001 | -0.109 | 0.374 | 1.000 | 0.473 |

### variant: `geom_adj`

| definition | pat_delta_cok1 | pat_delta_im1 | pat_delta_ker1 | pat_p_cok1 | pat_p_im1 | pat_p_ker1 |
| --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 0.271 | 0.133 | 0.245 | 0.073 | 0.381 | 0.105 |
| CD8_strict | 0.354 | 0.096 | 0.306 | 0.019 | 0.528 | 0.043 |
| CD8_extended | 0.261 | 0.135 | 0.144 | 0.085 | 0.374 | 0.344 |
| CD8_NK_only | -0.064 | -0.239 | -0.139 | 0.677 | 0.114 | 0.358 |
| CD4 | -0.042 | -0.107 | -0.057 | 0.785 | 0.482 | 0.709 |
| Bcell | 0.094 | -0.388 | -0.088 | 0.537 | 0.010 | 0.566 |
| Macrophage | 0.319 | -0.196 | 0.170 | 0.035 | 0.197 | 0.263 |
| Fibroblast | -0.163 | -0.280 | -0.150 | 0.282 | 0.064 | 0.322 |

## Full scan (all diagrams, all dimensions)

| definition | diagram | variant | roi_delta | roi_p | pat_delta | pat_p | pat_q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bcell | cod0 | /med_nn | -0.029 | 0.586 | 0.070 | 0.646 | 0.403 |
| Bcell | cod1 | /med_nn | -0.044 | 0.408 | 0.027 | 0.863 | 0.524 |
| Bcell | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.239 |
| Bcell | cok1 | /med_nn | 0.057 | 0.292 | 0.161 | 0.288 | 1.000 |
| Bcell | dom0 | /med_nn | -0.044 | 0.413 | 0.064 | 0.677 | 0.550 |
| Bcell | dom1 | /med_nn | -0.065 | 0.229 | -0.001 | 1.000 | 0.336 |
| Bcell | im0 | /med_nn | -0.029 | 0.586 | 0.070 | 0.646 | 0.239 |
| Bcell | im1 | /med_nn | -0.079 | 0.140 | -0.010 | 0.954 | 1.000 |
| Bcell | ker0 | /med_nn | 0.041 | 0.444 | 0.075 | 0.626 | 0.239 |
| Bcell | ker1 | /med_nn | 0.001 | 0.985 | 0.092 | 0.547 | 1.000 |
| CD4 | cod0 | /med_nn | -0.006 | 0.914 | 0.044 | 0.774 | 1.000 |
| CD4 | cod1 | /med_nn | -0.008 | 0.887 | 0.031 | 0.841 | 1.000 |
| CD4 | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.652 |
| CD4 | cok1 | /med_nn | 0.052 | 0.329 | 0.049 | 0.752 | 1.000 |
| CD4 | dom0 | /med_nn | -0.049 | 0.365 | 0.044 | 0.774 | 1.000 |
| CD4 | dom1 | /med_nn | -0.072 | 0.179 | 0.027 | 0.863 | 1.000 |
| CD4 | im0 | /med_nn | -0.006 | 0.914 | 0.044 | 0.774 | 1.000 |
| CD4 | im1 | /med_nn | -0.028 | 0.604 | -0.010 | 0.954 | 1.000 |
| CD4 | ker0 | /med_nn | 0.005 | 0.918 | -0.092 | 0.547 | 1.000 |
| CD4 | ker1 | /med_nn | 0.041 | 0.445 | 0.008 | 0.966 | 1.000 |
| CD8_NK_only | cod0 | /med_nn | -0.021 | 0.703 | 0.068 | 0.656 | 1.000 |
| CD8_NK_only | cod1 | /med_nn | -0.017 | 0.748 | 0.027 | 0.863 | 0.504 |
| CD8_NK_only | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_NK_only | cok1 | /med_nn | 0.007 | 0.900 | 0.057 | 0.709 | 1.000 |
| CD8_NK_only | dom0 | /med_nn | -0.040 | 0.461 | 0.083 | 0.586 | 1.000 |
| CD8_NK_only | dom1 | /med_nn | -0.037 | 0.493 | 0.025 | 0.875 | 0.278 |
| CD8_NK_only | im0 | /med_nn | -0.021 | 0.703 | 0.068 | 0.656 | 0.257 |
| CD8_NK_only | im1 | /med_nn | -0.043 | 0.425 | -0.025 | 0.875 | 0.220 |
| CD8_NK_only | ker0 | /med_nn | -0.073 | 0.171 | -0.114 | 0.456 | 0.257 |
| CD8_NK_only | ker1 | /med_nn | -0.033 | 0.543 | 0.003 | 0.989 | 0.257 |
| CD8_extended | cod0 | /med_nn | 0.250 | 0.000 | 0.276 | 0.068 | 0.257 |
| CD8_extended | cod1 | /med_nn | 0.252 | 0.000 | 0.334 | 0.027 | 1.000 |
| CD8_extended | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.220 |
| CD8_extended | cok1 | /med_nn | 0.250 | 0.000 | 0.319 | 0.035 | 1.000 |
| CD8_extended | dom0 | /med_nn | 0.235 | 0.000 | 0.274 | 0.071 | 0.628 |
| CD8_extended | dom1 | /med_nn | 0.232 | 0.000 | 0.217 | 0.151 | 0.997 |
| CD8_extended | im0 | /med_nn | 0.250 | 0.000 | 0.276 | 0.068 | 1.000 |
| CD8_extended | im1 | /med_nn | 0.210 | 0.000 | 0.278 | 0.066 | 0.992 |
| CD8_extended | ker0 | /med_nn | 0.048 | 0.371 | 0.150 | 0.322 | 1.000 |
| CD8_extended | ker1 | /med_nn | 0.145 | 0.007 | 0.265 | 0.080 | 0.997 |
| CD8_primary | cod0 | /med_nn | 0.275 | 0.000 | 0.328 | 0.030 | 0.349 |
| CD8_primary | cod1 | /med_nn | 0.247 | 0.000 | 0.364 | 0.016 | 0.200 |
| CD8_primary | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.200 |
| CD8_primary | cok1 | /med_nn | 0.262 | 0.000 | 0.334 | 0.027 | 1.000 |
| CD8_primary | dom0 | /med_nn | 0.259 | 0.000 | 0.304 | 0.044 | 0.603 |
| CD8_primary | dom1 | /med_nn | 0.236 | 0.000 | 0.286 | 0.058 | 0.200 |
| CD8_primary | im0 | /med_nn | 0.275 | 0.000 | 0.328 | 0.030 | 0.276 |
| CD8_primary | im1 | /med_nn | 0.211 | 0.000 | 0.291 | 0.054 | 0.200 |
| CD8_primary | ker0 | /med_nn | 0.073 | 0.171 | 0.189 | 0.212 | 0.276 |
| CD8_primary | ker1 | /med_nn | 0.176 | 0.001 | 0.336 | 0.026 | 0.256 |
| CD8_strict | cod0 | /med_nn | 0.289 | 0.000 | 0.345 | 0.022 | 0.997 |
| CD8_strict | cod1 | /med_nn | 0.260 | 0.000 | 0.366 | 0.015 | 0.997 |
| CD8_strict | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.603 |
| CD8_strict | cok1 | /med_nn | 0.323 | 0.000 | 0.414 | 0.006 | 1.000 |
| CD8_strict | dom0 | /med_nn | 0.260 | 0.000 | 0.302 | 0.046 | 0.997 |
| CD8_strict | dom1 | /med_nn | 0.231 | 0.000 | 0.284 | 0.060 | 0.997 |
| CD8_strict | im0 | /med_nn | 0.289 | 0.000 | 0.345 | 0.022 | 0.997 |
| CD8_strict | im1 | /med_nn | 0.197 | 0.000 | 0.282 | 0.062 | 0.997 |
| CD8_strict | ker0 | /med_nn | 0.122 | 0.023 | 0.211 | 0.164 | 0.997 |
| CD8_strict | ker1 | /med_nn | 0.231 | 0.000 | 0.375 | 0.013 | 0.997 |
| Fibroblast | cod0 | /med_nn | -0.018 | 0.736 | 0.059 | 0.698 | 0.808 |
| Fibroblast | cod1 | /med_nn | -0.041 | 0.442 | 0.031 | 0.841 | 0.292 |
| Fibroblast | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.524 |
| Fibroblast | cok1 | /med_nn | -0.044 | 0.416 | -0.135 | 0.374 | 1.000 |
| Fibroblast | dom0 | /med_nn | -0.059 | 0.275 | 0.023 | 0.886 | 0.625 |
| Fibroblast | dom1 | /med_nn | -0.084 | 0.118 | -0.051 | 0.741 | 0.501 |
| Fibroblast | im0 | /med_nn | -0.018 | 0.736 | 0.059 | 0.698 | 1.000 |
| Fibroblast | im1 | /med_nn | -0.031 | 0.568 | -0.001 | 1.000 | 0.285 |
| Fibroblast | ker0 | /med_nn | -0.044 | 0.410 | -0.109 | 0.473 | 1.000 |
| Fibroblast | ker1 | /med_nn | -0.046 | 0.395 | -0.109 | 0.473 | 0.239 |
| Macrophage | cod0 | /med_nn | 0.001 | 0.981 | 0.129 | 0.397 | 0.403 |
| Macrophage | cod1 | /med_nn | 0.001 | 0.988 | 0.114 | 0.456 | 0.285 |
| Macrophage | cok0 | /med_nn | 0.000 | 1.000 | 0.000 | 1.000 | 0.923 |
| Macrophage | cok1 | /med_nn | 0.166 | 0.002 | 0.183 | 0.228 | 1.000 |
| Macrophage | dom0 | /med_nn | -0.052 | 0.328 | 0.066 | 0.667 | 0.550 |
| Macrophage | dom1 | /med_nn | -0.057 | 0.288 | -0.012 | 0.943 | 0.562 |
| Macrophage | im0 | /med_nn | 0.001 | 0.981 | 0.129 | 0.397 | 0.239 |
| Macrophage | im1 | /med_nn | -0.040 | 0.458 | 0.040 | 0.796 | 0.349 |
| Macrophage | ker0 | /med_nn | 0.120 | 0.025 | 0.044 | 0.774 | 0.239 |
| Macrophage | ker1 | /med_nn | 0.141 | 0.009 | 0.200 | 0.187 | 0.564 |
| Bcell | cod0 | geom_adj | -0.235 | 0.000 | -0.261 | 0.085 | 0.239 |
| Bcell | cod1 | geom_adj | -0.238 | 0.000 | -0.280 | 0.064 | 0.385 |
| Bcell | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.563 |
| Bcell | cok1 | geom_adj | 0.032 | 0.555 | 0.094 | 0.537 | 1.000 |
| Bcell | dom0 | geom_adj | -0.207 | 0.000 | -0.265 | 0.080 | 0.239 |
| Bcell | dom1 | geom_adj | -0.251 | 0.000 | -0.295 | 0.051 | 0.565 |
| Bcell | im0 | geom_adj | -0.235 | 0.000 | -0.261 | 0.085 | 0.246 |
| Bcell | im1 | geom_adj | -0.289 | 0.000 | -0.388 | 0.010 | 1.000 |
| Bcell | ker0 | geom_adj | 0.086 | 0.110 | 0.111 | 0.464 | 1.000 |
| Bcell | ker1 | geom_adj | -0.055 | 0.307 | -0.088 | 0.566 | 0.963 |
| CD4 | cod0 | geom_adj | -0.250 | 0.000 | -0.209 | 0.168 | 0.963 |
| CD4 | cod1 | geom_adj | -0.240 | 0.000 | -0.224 | 0.139 | 1.000 |
| CD4 | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD4 | cok1 | geom_adj | -0.056 | 0.295 | -0.042 | 0.785 | 1.000 |
| CD4 | dom0 | geom_adj | -0.252 | 0.000 | -0.198 | 0.192 | 0.963 |
| CD4 | dom1 | geom_adj | -0.261 | 0.000 | -0.191 | 0.207 | 1.000 |
| CD4 | im0 | geom_adj | -0.250 | 0.000 | -0.209 | 0.168 | 0.575 |
| CD4 | im1 | geom_adj | -0.194 | 0.000 | -0.107 | 0.482 | 1.000 |
| CD4 | ker0 | geom_adj | 0.088 | 0.100 | 0.146 | 0.336 | 1.000 |
| CD4 | ker1 | geom_adj | -0.078 | 0.145 | -0.057 | 0.709 | 1.000 |
| CD8_NK_only | cod0 | geom_adj | -0.265 | 0.000 | -0.317 | 0.036 | 1.000 |
| CD8_NK_only | cod1 | geom_adj | -0.246 | 0.000 | -0.209 | 0.168 | 0.220 |
| CD8_NK_only | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_NK_only | cok1 | geom_adj | -0.049 | 0.363 | -0.064 | 0.677 | 0.257 |
| CD8_NK_only | dom0 | geom_adj | -0.279 | 0.000 | -0.347 | 0.022 | 1.000 |
| CD8_NK_only | dom1 | geom_adj | -0.267 | 0.000 | -0.286 | 0.058 | 0.257 |
| CD8_NK_only | im0 | geom_adj | -0.265 | 0.000 | -0.317 | 0.036 | 1.000 |
| CD8_NK_only | im1 | geom_adj | -0.213 | 0.000 | -0.239 | 0.115 | 0.525 |
| CD8_NK_only | ker0 | geom_adj | -0.013 | 0.813 | 0.001 | 1.000 | 1.000 |
| CD8_NK_only | ker1 | geom_adj | -0.107 | 0.045 | -0.140 | 0.358 | 0.220 |
| CD8_extended | cod0 | geom_adj | 0.250 | 0.000 | 0.334 | 0.027 | 0.220 |
| CD8_extended | cod1 | geom_adj | 0.181 | 0.001 | 0.237 | 0.118 | 1.000 |
| CD8_extended | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.257 |
| CD8_extended | cok1 | geom_adj | 0.181 | 0.001 | 0.261 | 0.085 | 0.997 |
| CD8_extended | dom0 | geom_adj | 0.170 | 0.002 | 0.204 | 0.177 | 0.220 |
| CD8_extended | dom1 | geom_adj | 0.133 | 0.013 | 0.157 | 0.301 | 1.000 |
| CD8_extended | im0 | geom_adj | 0.250 | 0.000 | 0.334 | 0.027 | 0.220 |
| CD8_extended | im1 | geom_adj | 0.150 | 0.005 | 0.135 | 0.374 | 0.997 |
| CD8_extended | ker0 | geom_adj | -0.015 | 0.779 | 0.090 | 0.556 | 1.000 |
| CD8_extended | ker1 | geom_adj | 0.060 | 0.264 | 0.144 | 0.344 | 0.542 |
| CD8_primary | cod0 | geom_adj | 0.242 | 0.000 | 0.332 | 0.028 | 0.453 |
| CD8_primary | cod1 | geom_adj | 0.148 | 0.006 | 0.226 | 0.136 | 0.200 |
| CD8_primary | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.256 |
| CD8_primary | cok1 | geom_adj | 0.191 | 0.000 | 0.271 | 0.073 | 0.279 |
| CD8_primary | dom0 | geom_adj | 0.170 | 0.002 | 0.196 | 0.197 | 0.453 |
| CD8_primary | dom1 | geom_adj | 0.127 | 0.018 | 0.159 | 0.295 | 0.276 |
| CD8_primary | im0 | geom_adj | 0.242 | 0.000 | 0.332 | 0.028 | 0.256 |
| CD8_primary | im1 | geom_adj | 0.098 | 0.068 | 0.133 | 0.381 | 0.453 |
| CD8_primary | ker0 | geom_adj | -0.028 | 0.609 | 0.085 | 0.576 | 1.000 |
| CD8_primary | ker1 | geom_adj | 0.105 | 0.050 | 0.245 | 0.105 | 0.200 |
| CD8_strict | cod0 | geom_adj | 0.240 | 0.000 | 0.336 | 0.026 | 0.997 |
| CD8_strict | cod1 | geom_adj | 0.163 | 0.002 | 0.209 | 0.168 | 0.997 |
| CD8_strict | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_strict | cok1 | geom_adj | 0.254 | 0.000 | 0.353 | 0.019 | 0.997 |
| CD8_strict | dom0 | geom_adj | 0.163 | 0.002 | 0.196 | 0.197 | 0.997 |
| CD8_strict | dom1 | geom_adj | 0.128 | 0.017 | 0.137 | 0.366 | 0.997 |
| CD8_strict | im0 | geom_adj | 0.240 | 0.000 | 0.336 | 0.026 | 0.542 |
| CD8_strict | im1 | geom_adj | 0.079 | 0.141 | 0.096 | 0.528 | 0.997 |
| CD8_strict | ker0 | geom_adj | 0.005 | 0.924 | 0.144 | 0.344 | 1.000 |
| CD8_strict | ker1 | geom_adj | 0.158 | 0.003 | 0.306 | 0.043 | 0.997 |
| Fibroblast | cod0 | geom_adj | -0.144 | 0.007 | -0.040 | 0.796 | 0.923 |
| Fibroblast | cod1 | geom_adj | -0.223 | 0.000 | -0.202 | 0.182 | 0.715 |
| Fibroblast | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.524 |
| Fibroblast | cok1 | geom_adj | -0.088 | 0.100 | -0.163 | 0.282 | 0.292 |
| Fibroblast | dom0 | geom_adj | -0.169 | 0.002 | -0.070 | 0.646 | 0.923 |
| Fibroblast | dom1 | geom_adj | -0.221 | 0.000 | -0.159 | 0.295 | 0.403 |
| Fibroblast | im0 | geom_adj | -0.144 | 0.007 | -0.040 | 0.796 | 0.548 |
| Fibroblast | im1 | geom_adj | -0.241 | 0.000 | -0.280 | 0.064 | 0.663 |
| Fibroblast | ker0 | geom_adj | -0.079 | 0.140 | -0.120 | 0.430 | 1.000 |
| Fibroblast | ker1 | geom_adj | -0.103 | 0.054 | -0.150 | 0.322 | 0.716 |
| Macrophage | cod0 | geom_adj | -0.134 | 0.013 | -0.010 | 0.954 | 0.403 |
| Macrophage | cod1 | geom_adj | -0.182 | 0.001 | -0.098 | 0.518 | 0.403 |
| Macrophage | cok0 | geom_adj | 0.000 | 1.000 | 0.000 | 1.000 | 0.413 |
| Macrophage | cok1 | geom_adj | 0.197 | 0.000 | 0.319 | 0.035 | 0.403 |
| Macrophage | dom0 | geom_adj | -0.252 | 0.000 | -0.180 | 0.234 | 0.403 |
| Macrophage | dom1 | geom_adj | -0.242 | 0.000 | -0.258 | 0.088 | 0.349 |
| Macrophage | im0 | geom_adj | -0.134 | 0.013 | -0.010 | 0.954 | 0.860 |
| Macrophage | im1 | geom_adj | -0.246 | 0.000 | -0.196 | 0.197 | 0.730 |
| Macrophage | ker0 | geom_adj | 0.150 | 0.005 | 0.312 | 0.039 | 1.000 |
| Macrophage | ker1 | geom_adj | 0.127 | 0.018 | 0.170 | 0.263 | 0.292 |
| Bcell | cod0 | raw | 0.041 | 0.449 | 0.072 | 0.636 | 1.000 |
| Bcell | cod1 | raw | 0.022 | 0.681 | 0.042 | 0.785 | 0.292 |
| Bcell | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.716 |
| Bcell | cok1 | raw | 0.109 | 0.041 | 0.180 | 0.234 | 1.000 |
| Bcell | dom0 | raw | 0.023 | 0.669 | 0.055 | 0.720 | 0.403 |
| Bcell | dom1 | raw | -0.003 | 0.956 | 0.003 | 0.989 | 1.000 |
| Bcell | im0 | raw | 0.041 | 0.449 | 0.072 | 0.636 | 0.403 |
| Bcell | im1 | raw | 0.006 | 0.915 | 0.042 | 0.785 | 1.000 |
| Bcell | ker0 | raw | 0.094 | 0.079 | 0.118 | 0.438 | 0.731 |
| Bcell | ker1 | raw | 0.073 | 0.175 | 0.196 | 0.197 | 1.000 |
| CD4 | cod0 | raw | 0.045 | 0.407 | 0.088 | 0.566 | 1.000 |
| CD4 | cod1 | raw | 0.047 | 0.379 | 0.064 | 0.677 | 0.795 |
| CD4 | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD4 | cok1 | raw | 0.091 | 0.089 | 0.083 | 0.586 | 1.000 |
| CD4 | dom0 | raw | 0.007 | 0.898 | 0.057 | 0.709 | 1.000 |
| CD4 | dom1 | raw | -0.011 | 0.837 | 0.049 | 0.752 | 1.000 |
| CD4 | im0 | raw | 0.045 | 0.407 | 0.088 | 0.566 | 1.000 |
| CD4 | im1 | raw | 0.040 | 0.456 | 0.055 | 0.720 | 1.000 |
| CD4 | ker0 | raw | 0.030 | 0.581 | -0.090 | 0.556 | 1.000 |
| CD4 | ker1 | raw | 0.090 | 0.092 | 0.068 | 0.656 | 1.000 |
| CD8_NK_only | cod0 | raw | 0.041 | 0.441 | 0.101 | 0.509 | 1.000 |
| CD8_NK_only | cod1 | raw | 0.042 | 0.430 | 0.059 | 0.698 | 0.233 |
| CD8_NK_only | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 1.000 |
| CD8_NK_only | cok1 | raw | 0.061 | 0.257 | 0.103 | 0.500 | 0.220 |
| CD8_NK_only | dom0 | raw | 0.021 | 0.697 | 0.083 | 0.586 | 1.000 |
| CD8_NK_only | dom1 | raw | 0.017 | 0.747 | 0.029 | 0.852 | 0.220 |
| CD8_NK_only | im0 | raw | 0.041 | 0.441 | 0.101 | 0.509 | 0.257 |
| CD8_NK_only | im1 | raw | 0.037 | 0.489 | 0.051 | 0.741 | 0.257 |
| CD8_NK_only | ker0 | raw | -0.023 | 0.663 | -0.044 | 0.774 | 0.859 |
| CD8_NK_only | ker1 | raw | 0.025 | 0.637 | 0.053 | 0.731 | 0.220 |
| CD8_extended | cod0 | raw | 0.314 | 0.000 | 0.224 | 0.139 | 1.000 |
| CD8_extended | cod1 | raw | 0.312 | 0.000 | 0.315 | 0.037 | 0.997 |
| CD8_extended | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.257 |
| CD8_extended | cok1 | raw | 0.325 | 0.000 | 0.377 | 0.013 | 0.925 |
| CD8_extended | dom0 | raw | 0.288 | 0.000 | 0.252 | 0.096 | 0.220 |
| CD8_extended | dom1 | raw | 0.278 | 0.000 | 0.310 | 0.040 | 0.925 |
| CD8_extended | im0 | raw | 0.314 | 0.000 | 0.224 | 0.139 | 1.000 |
| CD8_extended | im1 | raw | 0.274 | 0.000 | 0.286 | 0.058 | 1.000 |
| CD8_extended | ker0 | raw | 0.100 | 0.061 | 0.180 | 0.234 | 0.997 |
| CD8_extended | ker1 | raw | 0.230 | 0.000 | 0.308 | 0.042 | 0.499 |
| CD8_primary | cod0 | raw | 0.327 | 0.000 | 0.265 | 0.080 | 1.000 |
| CD8_primary | cod1 | raw | 0.312 | 0.000 | 0.356 | 0.019 | 0.200 |
| CD8_primary | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.276 |
| CD8_primary | cok1 | raw | 0.330 | 0.000 | 0.395 | 0.009 | 0.305 |
| CD8_primary | dom0 | raw | 0.305 | 0.000 | 0.278 | 0.066 | 0.256 |
| CD8_primary | dom1 | raw | 0.292 | 0.000 | 0.351 | 0.020 | 0.305 |
| CD8_primary | im0 | raw | 0.327 | 0.000 | 0.265 | 0.080 | 0.276 |
| CD8_primary | im1 | raw | 0.272 | 0.000 | 0.310 | 0.040 | 0.200 |
| CD8_primary | ker0 | raw | 0.124 | 0.021 | 0.222 | 0.143 | 0.453 |
| CD8_primary | ker1 | raw | 0.258 | 0.000 | 0.366 | 0.015 | 0.200 |
| CD8_strict | cod0 | raw | 0.336 | 0.000 | 0.284 | 0.060 | 1.000 |
| CD8_strict | cod1 | raw | 0.326 | 0.000 | 0.366 | 0.015 | 0.997 |
| CD8_strict | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.997 |
| CD8_strict | cok1 | raw | 0.387 | 0.000 | 0.444 | 0.003 | 0.997 |
| CD8_strict | dom0 | raw | 0.307 | 0.000 | 0.282 | 0.062 | 0.997 |
| CD8_strict | dom1 | raw | 0.291 | 0.000 | 0.345 | 0.022 | 0.997 |
| CD8_strict | im0 | raw | 0.336 | 0.000 | 0.284 | 0.060 | 0.997 |
| CD8_strict | im1 | raw | 0.261 | 0.000 | 0.297 | 0.049 | 1.000 |
| CD8_strict | ker0 | raw | 0.168 | 0.002 | 0.220 | 0.147 | 0.997 |
| CD8_strict | ker1 | raw | 0.304 | 0.000 | 0.399 | 0.008 | 0.997 |
| Fibroblast | cod0 | raw | 0.067 | 0.212 | 0.029 | 0.852 | 1.000 |
| Fibroblast | cod1 | raw | 0.051 | 0.345 | 0.036 | 0.818 | 0.239 |
| Fibroblast | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.285 |
| Fibroblast | cok1 | raw | 0.029 | 0.585 | -0.068 | 0.656 | 0.292 |
| Fibroblast | dom0 | raw | 0.031 | 0.558 | 0.012 | 0.943 | 0.403 |
| Fibroblast | dom1 | raw | 0.001 | 0.979 | 0.018 | 0.909 | 0.292 |
| Fibroblast | im0 | raw | 0.067 | 0.212 | 0.029 | 0.852 | 0.456 |
| Fibroblast | im1 | raw | 0.058 | 0.278 | 0.005 | 0.977 | 0.272 |
| Fibroblast | ker0 | raw | 0.003 | 0.963 | -0.105 | 0.491 | 0.239 |
| Fibroblast | ker1 | raw | 0.017 | 0.749 | -0.042 | 0.785 | 0.730 |
| Macrophage | cod0 | raw | 0.082 | 0.126 | 0.133 | 0.381 | 1.000 |
| Macrophage | cod1 | raw | 0.086 | 0.109 | 0.122 | 0.422 | 0.834 |
| Macrophage | cok0 | raw | 0.000 | 1.000 | 0.000 | 1.000 | 0.676 |
| Macrophage | cok1 | raw | 0.237 | 0.000 | 0.198 | 0.192 | 0.239 |
| Macrophage | dom0 | raw | 0.027 | 0.616 | 0.066 | 0.667 | 0.385 |
| Macrophage | dom1 | raw | 0.020 | 0.712 | 0.025 | 0.875 | 0.239 |
| Macrophage | im0 | raw | 0.082 | 0.126 | 0.133 | 0.381 | 0.239 |
| Macrophage | im1 | raw | 0.044 | 0.416 | 0.109 | 0.473 | 0.524 |
| Macrophage | ker0 | raw | 0.156 | 0.004 | 0.066 | 0.667 | 1.000 |
| Macrophage | ker1 | raw | 0.210 | 0.000 | 0.209 | 0.168 | 0.550 |

## Read

- CD8 aggregate definitions (primary/strict/extended), `/med_nn`: median per-patient delta = **+0.319**, 5/9 tests at p < 0.05.
- Control species (CD4, Bcell, Macrophage, Fibroblast), same variant: median delta = **+0.024**, 0/12 at p < 0.05.
- After BH across each variant's family, 0 of 240 tests survive q < 0.05.

**Caveat that matters.** The effect attenuates under `geom_adj`, and the same CD8-specific pattern appears in `dom0`/`dom1` -- the *monochromatic* complex, which contains no cross-species mixing information at all. That is consistent with 'responder tissue is more dispersed' rather than 'T cells infiltrate'. Regression adjustment cannot separate these, because density, dispersion and arrangement are collinear here. Only a null that holds cell positions fixed and shuffles labels can: see `perm_null.py` / `output/PERMUTATION_NULL.md`.

Script: `day2_length_features.py` (reads output/stats/*.parquet, no recomputation).
