# CLAHE batch-correction evaluation -- stage 3 metrics (balanced_batches)


## Cohort

- cohort: `pre`, batches [1, 2, 3, 4]
- cells 2,177,493, ROIs 389, patients 51

Batch x response, at ROI level -- read the metrics against this:

```
Response     Non-Responder  Responder  All
Stain_Batch                               
1                       19         54   73
2                       21         16   37
3                       20          2   22
4                      100        157  257
All                    160        229  389
```

## Subsampling

- graph metrics (iLISI, connectivity, Leiden, UMAP, PCR): 100,000 cells
- silhouette metrics (ASW-batch, ASW-cell-type): 10,000 cells, since silhouette is O(n^2)
- the identical cell indices are used for both tables, so every comparison is exactly paired

## Clustering for ARI / NMI

- Leiden (resolution 1.0, flavor=leidenalg), identical for both tables

## Bootstrap

- 1000 replicates, resampling **patients** with replacement (51 patients), not cells
- per-cell metrics are re-aggregated over the resampled patients; ARI, NMI and PCR are fully recomputed on the resampled cells
- graph connectivity is a whole-graph quantity and is reported as a point estimate with no interval

## Results

Every metric is oriented so that **higher is better**; ASW-batch and PCR are inverted at computation to achieve this, following scib's normalised convention. `difference` is corrected minus uncorrected, so a positive value favours CLAHE.

```
                                 metric      family  corrected  uncorrected  difference  ci_low  ci_high  bootstrap_bias  interval_excludes_zero
0                  iLISI (batch mixing)       batch     1.3579       1.0345      0.3234  0.2430   0.3862          0.0022                    True
1                    graph connectivity       batch     0.7051       0.5579      0.1471     NaN      NaN             NaN                   False
2  ASW-batch (1-|sil|, scib convention)       batch     0.9519       0.8291      0.1228  0.1069   0.1392         -0.0005                    True
3          PCR (1 - R2 of batch on PCs)       batch     0.9901       0.8600      0.1301  0.0471   0.1844          0.0056                    True
4                         ASW-cell-type  biological     0.4579       0.4424      0.0154  0.0082   0.0218          0.0000                    True
5                 ARI (leiden vs Pixie)  biological     0.0840       0.0361      0.0479  0.0351   0.0724         -0.0070                    True
6                 NMI (leiden vs Pixie)  biological     0.2034       0.1321      0.0713  0.0670   0.0881         -0.0059                    True
```

`ci_low`/`ci_high` are basic (reverse-percentile) bootstrap intervals. `bootstrap_bias` is the mean bootstrap difference minus the observed one: resampling patients with replacement duplicates whole patients, which biases ARI and NMI downward, and a raw percentile interval can then exclude the point estimate. The raw percentiles are kept in the CSV as `pct_low`/`pct_high` so the correction is auditable.

## Reading these numbers

- batch metrics significantly improved by CLAHE: 3 of 4
- biological conservation metrics significantly REDUCED by CLAHE: 0 of 3

**Caveat, to be stated in the rebuttal**: the Pixie labels were derived from CLAHE-corrected data, so ASW-cell-type, ARI and NMI mildly favour the corrected table by construction. Conservation metrics that are merely unchanged should be read in that light. Re-clustering and re-annotating the uncorrected data from scratch would remove this asymmetry and is a much larger job, deliberately not attempted here.

**Patient identity is not used as a secondary conservation label**: 59 of 62 patients sit in a single stain batch, so patient and batch are nearly the same variable and the check would be uninformative.

Written: `metrics_balanced_batches.csv`, `metric_panel_balanced_batches.png`.
