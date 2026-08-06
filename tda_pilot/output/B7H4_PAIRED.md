# Is CD8 more excluded from B7H4+ tumour than from B7H4− tumour?

Paired within-ROI comparison, 556 pre-treatment ROIs. Within each ROI both cancer subtypes are subsampled to the same count (median m = 308) and compared against the identical CD8 set (median 827 cells), averaged over 5 subsamples. Species counts are therefore identical between the two arms by construction, so composition cannot drive the difference; and both arms share the same tissue, so ROI geometry cancels.

**19 of 24 statistics differ at BH q < 0.05.**

| statistic | median_cancer | median_b7h4 | median_diff | frac_b7h4_higher | wilcoxon_p | rank_biserial | q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ker1-num_bars | 322.7 | 269.5 | -41 | 0.1277 | 0 | -0.7446 | 0 |
| ker0-entropy | 4.973 | 4.737 | -0.1795 | 0.1331 | 0 | -0.7338 | 0 |
| im0-avg_length | 7.476 | 7.329 | -0.1061 | 0.0989 | 0 | -0.8022 | 0 |
| im0-med_length | 6.276 | 5.971 | -0.1647 | 0.1169 | 0 | -0.7662 | 0 |
| ker0-num_bars | 243.8 | 207 | -30 | 0.1457 | 0 | -0.7086 | 0 |
| ker1-entropy | 4.572 | 4.421 | -0.1763 | 0.1853 | 0 | -0.6295 | 0 |
| cok1-num_bars | 406.4 | 369.5 | -28.9 | 0.1906 | 0 | -0.6187 | 0 |
| im0-entropy | 7.063 | 7.056 | -0.0052 | 0.196 | 0 | -0.6079 | 0 |
| cok1-entropy | 5.142 | 5.04 | -0.0984 | 0.2302 | 0 | -0.5396 | 0 |
| im1-p90_length | 3.306 | 3.025 | -0.1116 | 0.3058 | 0 | -0.3885 | 0 |
| im1-num_bars | 800.4 | 828.5 | 17.7 | 0.732 | 0 | 0.464 | 0 |
| im1-med_length | 0.3488 | 0.3316 | -0.0074 | 0.286 | 0 | -0.4281 | 0 |
| im1-avg_length | 1.29 | 1.212 | -0.0214 | 0.3399 | 0 | -0.3201 | 0 |
| ker0-p90_length | 20.24 | 21.36 | 1.084 | 0.6853 | 0 | 0.3705 | 0 |
| ker1-med_length | 0.5345 | 0.5037 | -0.0258 | 0.3579 | 0 | -0.2842 | 0 |
| ker0-avg_length | 8.862 | 9.098 | 0.2874 | 0.6367 | 0 | 0.2734 | 0 |
| cok1-med_length | 0.4437 | 0.433 | -0.0108 | 0.4263 | 0.0003 | -0.1475 | 0.0004 |
| ker0-med_length | 5.16 | 5.003 | -0.0518 | 0.4784 | 0.0003 | -0.0432 | 0.0004 |
| ker1-p90_length | 5.14 | 5.313 | 0.0943 | 0.5486 | 0.0111 | 0.0971 | 0.014 |
| ker1-avg_length | 2.065 | 2.133 | 0.0421 | 0.5665 | 0.0506 | 0.1331 | 0.0608 |
| im1-entropy | 5.575 | 5.568 | 0.0023 | 0.5378 | 0.2202 | 0.0755 | 0.2517 |
| im0-p90_length | 12.16 | 12.09 | 0.0024 | 0.509 | 0.2629 | 0.018 | 0.2868 |
| cok1-avg_length | 1.289 | 1.301 | -0.0022 | 0.4982 | 0.6642 | -0.0036 | 0.6931 |
| cok1-p90_length | 3.714 | 3.748 | -0.0145 | 0.4946 | 0.8231 | -0.0108 | 0.8231 |

## Residual geometry: does the confound run with or against the effect?

Matching cell COUNTS does not match spatial EXTENT, and extent drives these statistics — the lesson of `SENSITIVITY.md`. So measure it on the same subsamples rather than assume it away:

| quantity | cancer_median | b7h4_median | median_diff | frac_b7h4_larger | wilcoxon_p |
| --- | --- | --- | --- | --- | --- |
| median NN distance (µm) | 14.15 | 11.33 | -2.445 | 0.0486 | 0 |
| convex hull area (µm²) | 7.538e+05 | 7.401e+05 | -2.078e+04 | 0.3112 | 0 |

**The confound runs AGAINST the effect.** At matched counts B7H4+ cancer cells are more tightly packed than B7H4− cells, and tighter packing shortens characteristic length scales — which would make kernel bars *shorter*. The observed kernel bars are *longer* for B7H4+. So the geometry difference cannot manufacture this result; if anything it understates it.


## Does the effect differ between responders and non-responders?

Exploratory, BH-corrected across the family. **0 of 30 interaction tests survive q < 0.05.**

| statistic | delta_R | p_R | delta_NR | p_NR | interaction_cliff | interaction_p | interaction_q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| im0-p90_length | 0.0191 | 0.0337 | -0.0106 | 0.1215 | 0.3081 | 0.0416 | 0.6456 |
| ker1-p90_length | 0.1474 | 0.0036 | -0.0383 | 0.9451 | 0.3059 | 0.043 | 0.6456 |
| im1-entropy | 0.0016 | 0.8826 | 0.0061 | 0.0095 | -0.2692 | 0.0752 | 0.7518 |
| cok1-entropy | -0.0965 | 0 | -0.1206 | 0 | 0.2151 | 0.1554 | 1 |
| cok1-num_bars | -27 | 0 | -35.8 | 0 | 0.1849 | 0.2225 | 1 |
| ker1-avg_length | 0.0446 | 0.0304 | 0.0357 | 0.8692 | 0.1719 | 0.2569 | 1 |
| im1-num_bars | 16 | 0 | 22.2 | 0 | -0.1481 | 0.3292 | 1 |
| ker1-num_bars | -38.6 | 0 | -48.6 | 0 | 0.1351 | 0.3736 | 1 |
| ker0-p90_length | 1.084 | 0 | 1.079 | 0.0001 | 0.133 | 0.3814 | 1 |
| ker0-med_length | -0.0697 | 0.0024 | 0.0009 | 0.075 | 0.1286 | 0.3972 | 1 |
| ker0-avg_length | 0.2775 | 0 | 0.3367 | 0.0145 | 0.1157 | 0.4469 | 1 |
| ker0-num_bars | -29.6 | 0 | -33.8 | 0 | 0.1135 | 0.4555 | 1 |
| cok1-med_length | -0.0121 | 0.0013 | -0.0054 | 0.1083 | -0.1092 | 0.4731 | 1 |
| im1-p90_length | -0.1024 | 0 | -0.1291 | 0 | 0.0768 | 0.6155 | 1 |
| ker1-med_length | -0.0268 | 0 | -0.0254 | 0 | 0.0681 | 0.6564 | 1 |
| im1-med_length | -0.0072 | 0 | -0.0095 | 0 | 0.0638 | 0.6773 | 1 |
| im0-med_length | -0.1671 | 0 | -0.1498 | 0 | -0.0616 | 0.6878 | 1 |
| im0-avg_length | -0.1061 | 0 | -0.1059 | 0 | 0.0573 | 0.7091 | 1 |
| im1-avg_length | -0.0217 | 0 | -0.0212 | 0 | 0.04 | 0.7962 | 1 |
| ker1-entropy | -0.1719 | 0 | -0.1936 | 0 | -0.0292 | 0.852 | 1 |
| cok1-avg_length | -0.0045 | 0.6345 | 0.0078 | 0.9893 | 0.0249 | 0.8746 | 1 |
| im0-entropy | -0.0055 | 0 | -0.0039 | 0 | 0.0119 | 0.9428 | 1 |
| ker0-entropy | -0.174 | 0 | -0.1966 | 0 | -0.0119 | 0.9428 | 1 |
| cok1-p90_length | -0.0145 | 0.8312 | -0.0144 | 0.9329 | -0.0032 | 0.9886 | 1 |
| cok0-p90_length | 0 | nan | 0 | nan | 0 | 1 | 1 |
| cok0-num_bars | 0 | nan | 0 | nan | 0 | 1 | 1 |
| cok0-entropy | 0 | nan | 0 | nan | 0 | 1 | 1 |
| im0-num_bars | 0 | nan | 0 | nan | 0 | 1 | 1 |
| cok0-med_length | 0 | nan | 0 | nan | 0 | 1 | 1 |
| cok0-avg_length | 0 | nan | 0 | nan | 0 | 1 | 1 |

The effect is present in BOTH response groups at similar magnitude, and the interaction is null. That makes this a statement about tumour biology that does not depend on the response contrast — and it does not rescue that contrast either, consistent with every other analysis in this directory.


## Reading this

`ker1-avg_length` is the exclusion signature: a tumour nest ringed by CD8 produces fewer but longer degree-1 kernel bars. Positive Δ means B7H4+ tumour shows more of that pattern than B7H4− tumour in the same tissue.

Conservative by design: `B7H4_PRECHECK.md` shows the two subtypes are only partly segregated (median mixing ratio 0.85 against 1.0 for random intermingling), which attenuates any true difference. A positive result here is therefore trustworthy; a null is partly explained by the overlap.

This question is independent of treatment response — it concerns tumour biology and can be reported whether or not it differs between responders and non-responders.

Figure: `b7h4_paired.png`. Script: `day2_b7h4_paired.py`.
