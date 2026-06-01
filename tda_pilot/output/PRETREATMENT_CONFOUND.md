# Decisive confound: the responder signal is treatment-stage sampling, not topology

**Date:** 2026-05-30. Definition `CD8_primary`, feature `im_dim1_total_persistence`.
Cliff's δ = R vs NR (−ve = responders lower).

## Pre-treatment composition is flat
Per-ROI CD8:tumour ratio, **pre-treatment only** (`composition_ratio.png`):
- per-ROI: median NR 0.77 vs R 0.64, Cliff δ = −0.077, MWU p = 0.15
- **per-patient (honest unit): median NR 0.70 vs R 0.65, δ = +0.008, p = 0.97 — no difference.**

## The image-persistence signal is a pre/post artifact
| subset | nROI (NR/R) | raw δ | size-adj δ (ROI) | size-adj δ (patient) |
|---|---|--:|--:|--:|
| ALL (pre+post) | 797 (330/467) | −0.168 | −0.138 | −0.225 (p=0.13) |
| **PRE only** | 582 (160/422) | −0.061 | **+0.056** | +0.146 (p=0.34) |
| POST only | 177 (**170/7**) | −0.605 | +0.055 | +0.143 (p=0.91) |

**Mechanism:** post-treatment resections are ~96% non-responder (170 vs 7) — responders have
little/no residual tumour to image. Those residual-disease NR ROIs are tumour-imbalanced →
high image persistence → they inflate the NR average in the pooled cohort, manufacturing the
"responders lower" signal. On PRE-treatment biopsies the size-adjusted signal is null/slightly
reversed.

## Conclusion
- The −0.138 signal was **confounded by treatment stage**, not biology.
- Pre-treatment, neither **composition** (δ=+0.008, p=0.97) nor **size-adjusted topology**
  (δ=+0.146, p=0.34) distinguishes R from NR.
- **No-go** for TDA on the responder contrast. Honest statement: at baseline, tumour–CD8 spatial
  organisation does not distinguish responders from non-responders by composition or chromatic topology
  in this cohort.

Caveats: pre-treatment per-patient comparison is modestly powered (25 NR / 37 R patients) — rules out a
strong effect, not a tiny one. **Any pooled pre+post analysis (incl. cross-PCF) inherits the same
treatment-stage confound** and should be redone pre-treatment-only.

Script: `day1_composition_check.py` (+ inline pre/post decomposition).
