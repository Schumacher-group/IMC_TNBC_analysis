# Critical confound: the image-persistence responder signal is class imbalance

**Date:** 2026-05-30. Definition: `CD8_primary` (Tumour vs CD8 + Memory CD8), 808-ROI cohort.
Feature: `im_dim1_total_persistence`. Effect size = Cliff's δ (R vs NR; −ve = responders lower).

## The finding
The apparent "topological infiltration" signal (responders = lower image persistence,
size-adjusted δ = −0.138) is **driven almost entirely by class imbalance** (tumour:CD8 ratio),
not spatial arrangement.

Minority fraction `f = min(n_tumour,n_other)/(n_tumour+n_other)`:
- responders are more balanced: median f 0.266 (NR) vs 0.309 (R), Cliff δ = +0.211.
- the size-adjusted residual is strongly tied to f: Spearman ρ = −0.54
  (imbalanced ROIs → high residual persistence, trivially: a sparse minority lets the
  majority form large monochromatic loops).

## Subset test (drop imbalanced ROIs, recompute) — decisive
| keep min_frac ≥ | n (NR/R) | raw δ | size-adj δ |
|---|---|--:|--:|
| 0.00 (all) | 330/467 | −0.168 | **−0.138** |
| 0.10 | 284/451 | −0.105 | −0.058 |
| 0.15 | 244/423 | −0.033 | +0.035 |
| 0.20 | 203/379 | −0.004 | +0.123 |
| 0.25 | 179/313 | +0.006 | +0.107 |
| 0.30 | 136/243 | +0.001 | +0.103 |

On balanced ROIs (minority ≥ 15–20%, still 580–670 ROIs), the responder-lower signal
**vanishes and mildly reverses**. The full-cohort −0.138 is carried by ROIs where one species
is nearly absent.

Note: a linear regression adjustment (adding `minority_frac` as a covariate) only moved
δ −0.138 → −0.067, *understating* the problem, because the imbalance effect is strongly
nonlinear. The subset test is the honest one.

## Interpretation
- The real responder difference is **compositional** (responders have relatively more CD8 per
  tumour cell) — capturable by abundance/composition analysis or normalised cross-PCF, **no TDA needed**.
- Where chromatic TDA's segregation-vs-interspersion question is well-posed (both species present),
  there is **no** responder difference in image topology.
- The TDA *method* is valid (size-matched galleries do separate segregated vs interspersed at the
  extremes), but it does **not** add spatial-arrangement signal beyond composition for the responder contrast.

## Implication
No-go indicator for TDA strengthening the responder rebuttal on this axis. Decision pending
comparison with the cross-PCF results.

Figures: `imbalance_check.png`, `threshold_check.png`. Scripts: `day1_imbalance_check.py`,
`day1_threshold_check.py`.
