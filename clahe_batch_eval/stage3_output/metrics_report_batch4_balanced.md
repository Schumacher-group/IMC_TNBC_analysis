# CLAHE batch-correction evaluation -- stage 3 metrics (batch4_balanced)


## Cohort

- cohort: `pre`, batches [4]
- cells 1,637,521, ROIs 257, patients 32

Batch x response, at ROI level -- read the metrics against this:

```
Response     Non-Responder  Responder  All
Stain_Batch                               
4                      100        157  257
All                    100        157  257
```

## Subsampling

- graph metrics (iLISI, connectivity, Leiden, UMAP, PCR): 100,000 cells
- silhouette metrics (ASW-batch, ASW-cell-type): 10,000 cells, since silhouette is O(n^2)
- the identical cell indices are used for both tables, so every comparison is exactly paired

## Clustering for ARI / NMI

- Leiden (resolution 1.0), identical for both tables

## Bootstrap

- 1000 replicates, resampling **patients** with replacement (32 patients), not cells
- per-cell metrics are re-aggregated over the resampled patients; ARI, NMI and PCR are fully recomputed on the resampled cells
- graph connectivity is a whole-graph quantity and is reported as a point estimate with no interval

## Results

Every metric is oriented so that **higher is better**; ASW-batch and PCR are inverted at computation to achieve this, following scib's normalised convention. `difference` is corrected minus uncorrected, so a positive value favours CLAHE.

```
                                 metric      family  corrected  uncorrected  difference  ci_low  ci_high  interval_excludes_zero
0                  iLISI (batch mixing)       batch     1.0000       1.0000      0.0000  0.0000   0.0000                   False
1                    graph connectivity       batch     0.7101       0.5938      0.1163     NaN      NaN                   False
2  ASW-batch (1-|sil|, scib convention)       batch        NaN          NaN         NaN     NaN      NaN                   False
3          PCR (1 - R2 of batch on PCs)       batch        NaN          NaN         NaN     NaN      NaN                   False
4                         ASW-cell-type  biological     0.4587       0.4380      0.0207  0.0133   0.0290                    True
5                 ARI (leiden vs Pixie)  biological     0.1300       0.0532      0.0767  0.0441   0.0894                    True
6                 NMI (leiden vs Pixie)  biological     0.2387       0.1480      0.0906  0.0660   0.1009                    True
```

## Reading these numbers

- batch metrics significantly improved by CLAHE: 0 of 4
- biological conservation metrics significantly REDUCED by CLAHE: 0 of 3

**Caveat, to be stated in the rebuttal**: the Pixie labels were derived from CLAHE-corrected data, so ASW-cell-type, ARI and NMI mildly favour the corrected table by construction. Conservation metrics that are merely unchanged should be read in that light. Re-clustering and re-annotating the uncorrected data from scratch would remove this asymmetry and is a much larger job, deliberately not attempted here.

**Patient identity is not used as a secondary conservation label**: 59 of 62 patients sit in a single stain batch, so patient and batch are nearly the same variable and the check would be uninformative.

Written: `metrics_batch4_balanced.csv`, `metric_panel_batch4_balanced.png`.
