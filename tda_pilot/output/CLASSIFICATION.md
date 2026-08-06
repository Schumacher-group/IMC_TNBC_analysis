# M2S2-style classification: can TDA statistics predict response?

Pre-treatment, patient-level. 62 patients (25 NR / 37 R), 1520 statistics across 8 pair definitions, reduced by per-cell-group PCA inside the cross-validation.

| model | n_features | cv_auc |
| --- | --- | --- |
| TDA (all 8 pair definitions, per-group PCA) | 1520 | 0.544 |
| TDA, CD8_primary only | 190 | 0.573 |
| Composition + geometry only | 7 | 0.559 |
| DummyClassifier (stratified) | 0 | 0.495 |

- Label-permutation null (200 draws): null mean AUC = 0.496, 95th pct = 0.630; observed = **0.544**, empirical **p = 0.275**.
- A multivariate classifier over the full topological feature set therefore **does not** predict response better than chance.
- Composition + geometry alone reaches AUC = 0.559, so TDA is not being outperformed by a trivial abundance model either — neither carries signal.

## Two checks worth stating

- **The null is centred where it should be.** Mean null AUC = 0.496 over 200 label
  permutations. If the per-cell-group PCA were leaking across folds — as it does in the
  original M2S2 script, which fits PCA before `cross_validate` — the null would sit above
  0.5. It does not, so the cross-validation is honest.
- **Power.** The null's 95th percentile is AUC = 0.630, so with 62 patients this design can
  only detect a classifier of roughly AUC ≳ 0.63. A weak-but-real signal (AUC 0.55–0.60)
  would not be distinguishable from chance here. "No signal" means "no strong signal".

## Why this matters for the rebuttal

This is the analysis the Day-1 report deferred and never ran, and it is the one a reviewer is most likely to ask for: not a single hand-picked statistic, but the whole feature set given to a classifier with the authors' own hyperparameters. It is also the form of evidence least sensitive to the confounding arguments elsewhere in this directory — a classifier is free to use composition, geometry, topology or any combination, and still cannot separate the groups.

Script: `day2_classify_response.py`. Figure: `classification_figure.png`.
