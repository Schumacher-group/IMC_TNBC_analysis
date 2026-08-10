# CLAHE batch-correction evaluation -- stage 3 metrics (pre)


## Cohort

- cohort: `pre`
- cells 3,227,482, ROIs 593, patients 62

Batch x response, at ROI level -- read the metrics against this:

```
Response     Non-Responder  Responder  All
Stain_Batch                               
1                       19         54   73
2                       21         16   37
3                       20          2   22
4                      100        157  257
6                        0        174  174
7                        0         30   30
All                    160        433  593
```

**Batches with a single response group: [6, 7].** Between-batch differences there are partly real biology, so improved mixing is not unambiguously good. See the batch-4-only run for a balanced comparison.

## Subsampling

- graph metrics (iLISI, connectivity, Leiden, UMAP, PCR): 100,000 cells
- silhouette metrics (ASW-batch, ASW-cell-type): 10,000 cells, since silhouette is O(n^2)
- the identical cell indices are used for both tables, so every comparison is exactly paired

## Clustering for ARI / NMI

- Leiden (resolution 1.0, flavor=leidenalg), identical for both tables

## Bootstrap

- 1000 replicates, resampling **patients** with replacement (62 patients), not cells
- per-cell metrics are re-aggregated over the resampled patients; ARI, NMI and PCR are fully recomputed on the resampled cells
- graph connectivity is a whole-graph quantity and is reported as a point estimate with no interval

## Results

Every metric is oriented so that **higher is better**; ASW-batch and PCR are inverted at computation to achieve this, following scib's normalised convention. `difference` is corrected minus uncorrected, so a positive value favours CLAHE.

```
                                 metric      family  corrected  uncorrected  difference  ci_low  ci_high  bootstrap_bias  interval_excludes_zero
0                  iLISI (batch mixing)       batch     1.8091       1.0577      0.7515  0.6737   0.8198          0.0002                    True
1                    graph connectivity       batch     0.7544       0.5939      0.1605     NaN      NaN             NaN                   False
2  ASW-batch (1-|sil|, scib convention)       batch     0.9434       0.8687      0.0747  0.0544   0.0943          0.0001                    True
3          PCR (1 - R2 of batch on PCs)       batch     0.9854       0.8244      0.1610  0.0920   0.2093          0.0075                    True
4                         ASW-cell-type  biological     0.4645       0.4402      0.0243  0.0172   0.0297          0.0003                    True
5                 ARI (leiden vs Pixie)  biological     0.1238       0.0259      0.0979  0.0776   0.1218         -0.0059                    True
6                 NMI (leiden vs Pixie)  biological     0.2316       0.1042      0.1274  0.1222   0.1466         -0.0067                    True
```

`ci_low`/`ci_high` are basic (reverse-percentile) bootstrap intervals. `bootstrap_bias` is the mean bootstrap difference minus the observed one: resampling patients with replacement duplicates whole patients, which biases ARI and NMI downward, and a raw percentile interval can then exclude the point estimate. The raw percentiles are kept in the CSV as `pct_low`/`pct_high` so the correction is auditable.

## Reading these numbers

- batch metrics significantly improved by CLAHE: 3 of 4
- biological conservation metrics significantly REDUCED by CLAHE: 0 of 3

**Caveat, to be stated in the rebuttal**: the Pixie labels were derived from CLAHE-corrected data, so ASW-cell-type, ARI and NMI mildly favour the corrected table by construction. Conservation metrics that are merely unchanged should be read in that light. Re-clustering and re-annotating the uncorrected data from scratch would remove this asymmetry and is a much larger job, deliberately not attempted here.

**Patient identity is not used as a secondary conservation label**: 59 of 62 patients sit in a single stain batch, so patient and batch are nearly the same variable and the check would be uninformative.

## scib cross-check

scib not importable (`No module named 'scib'`). All metrics above are the native implementations, which is why they were written that way.

Written: `metrics_pre.csv`, `metric_panel_pre.png`.
