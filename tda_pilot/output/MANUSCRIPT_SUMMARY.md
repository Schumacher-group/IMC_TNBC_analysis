# Bottom line for the revision

`REVIEW.md` is the full audit trail, including two retractions and the reasoning behind each.
This file is the distilled version: what can be claimed, with which numbers, and what the
figures should be. Everything below is pre-treatment, patient-level, RCB endpoint
(Responder = RCB 0–1), 25 NR / 37 R patients, 593 ROIs.

---

## 1. What can be claimed

### A. B7H4+ tumour differs from B7H4− tumour in its relationship with CD8  — *positive*

Paired within-ROI, both cancer subtypes subsampled to identical counts and compared against the
same CD8 cells (556 ROIs, median 391 cells per arm).

| | |
|---|---|
| degree-1 kernel bar **count** | 322.7 → 269.5, rank-biserial −0.74, p < 1e-4 |
| degree-0 kernel bar count | 243.8 → 207.0, −0.71, p < 1e-4 |
| independent confirmation (Dowker) | same direction, ρ = 0.68 between methods |
| survives matching on cell clustering | yes, both methods (p < 0.002) |

**Claim it as a difference in the CD8–tumour topological relationship, not as "more exclusion".**
The *count* signal is robust across two structurally independent constructions. The *length*
signal — which is what would license the word "exclusion" — appears in the chromatic
construction but not the Dowker one, and B7H4+ cells are separately measurable as more tightly
clustered (median nearest-neighbour 11.3 vs 14.2 µm). So the architecture differs; whether CD8
*accessibility* differs beyond that is not established.

**No response interaction.** Present in both groups at similar magnitude (count: −38.6 in R,
−48.6 in NR, both p < 1e-4), 0 of 24 interaction tests survive BH correction. This is a
statement about tumour biology that does not depend on the response contrast — which is the
cleanest form it could take.

### B. Baseline tumour–CD8 architecture does not distinguish response  — *null, and calibrated*

Eight independent framings, all null:

| framing | result |
|---|---|
| supervised classifier, 1520 statistics, 8 pair definitions | AUC 0.544, permutation p = 0.275 |
| pairwise permutation null, colour-symmetric | δ = +0.17, p = 0.27 |
| pairwise permutation null, CD8-directional | δ = +0.38, p = 0.012 raw — see caveat |
| triple with fibroblast | AUC 0.338 |
| triple with collagen-high cells | AUC 0.514 |
| triple with real collagen fibres (pre-CLAHE) | AUC 0.453 |
| triple with real collagen fibres (post-CLAHE) | AUC 0.539 |
| direct geometric barrier test | δ = −0.07, p = 0.67 |

**The negative is calibrated, which is what makes it worth reporting.** The instrument
demonstrably separates the two Fig 4B panels — 18th vs 90th percentile of the cohort, ~9 SD of
their respective nulls, in the mechanistically predicted direction — and 42–45% of
ROI × statistic comparisons exceed |z| = 2 against a composition-matched null. So the absence
of a group difference reflects the data, not an insensitive measurement.

**Caveat that must travel with it.** With 25 NR / 37 R the classifier detects only AUC ≳ 0.63,
and rank tests only |δ| ≳ 0.42. Say *"no large difference"*, never *"no difference"*.

**The one unresolved thread**, and it should be stated rather than buried: the CD8-directional
statistic gives δ = +0.38 (p = 0.012) unadjusted and holds at +0.27–0.38 under geometry
*matching*, but collapses under one particular regression adjustment. Across reasonable
specifications it spans +0.01 to +0.40. Two things keep it short of a finding regardless of
normalisation: its direction is **opposite** to the exclusion hypothesis (responders higher),
and the macrophage control sits at p = 0.060, which undercuts a CD8-specific reading.

### C. The collagen barrier hypothesis  — *null, and now properly tested*

Earlier collagen tests used cell-based proxies (fibroblasts *make* collagen; collagen-high cells
*sit in* it), and both need a cell to place a point on, so both were blind to dense acellular
matrix — which is what a physical barrier is. The fibre segmentation is not. Tested with real
segmented collagen, under both image-processing variants, by two independent methods
(three-species topology and a direct geometric barrier test): null throughout.

A robust finding did come out of it, independent of response: **CD8 cells reach tumour along
routes with less intervening collagen than distance-matched comparable leukocytes** — 82% of
ROIs, median z = −3.74, equally in both groups (R −4.74, NR −3.51). This is about where CD8 are
found, not about whether collagen impedes them.

### D. Composition is flat

Pre-treatment CD8:tumour ratio, per patient: δ = +0.008, p = 0.97. Worth stating because it
forecloses the simplest alternative explanation.

---

## 2. Proposed figures

### Main figure 1 — B7H4  (`b7h4_paired.png`)
- **a** Design: within each ROI, both cancer subtypes subsampled to equal counts, same CD8 set.
- **b** Paired difference in degree-1 kernel bar count — the robust statistic.
- **c** Effect sizes across statistics, BH-corrected.

### Main figure 2 — the calibrated null  (`negative_result_figure.png` + `fig_counterexamples.png`)
- **a** Mechanistic plane (bar count vs bar length), Fig 4B panels highlighted — shows the
  readout separates them (18th vs 90th percentile).
- **b** Per-patient inference: distributions overlap, δ and p annotated.
- **c** **Counterexample gallery** — responder/non-responder pairs matched on the exclusion
  score across the whole spectrum. This is the panel that makes the null legible: the same
  architecture occurs in both groups at every level.

### Supplementary
| panel | asset | shows |
|---|---|---|
| SI-1 | `b7h4_paired.png` (interaction section) | B7H4 effect equal in R and NR, 0/24 at BH q<0.05 |
| SI-2 | `b7h4_precheck.png` | the two subtypes are spatially separable (mixing ratio 0.85) — the design check |
| SI-3 | `classification_figure.png` | classifier AUC 0.544 against a null centred at 0.496 |
| SI-4 | `barrier_test_nonclahe.png` | direct collagen barrier test, null |
| SI-5 | `geometry_supplement.png` | why the per-ROI null needs a cross-ROI geometry check |
| SI-6 | `raw_statistics_figure.png` | raw statistics are dominated by composition (ρ = −0.69) and cell count (ρ = +0.83) |

### Which null to lead with
**The classifier.** One number, no adjustment choices, a permutation null centred where it
should be, and it is the framing least vulnerable to the confounding arguments — a classifier is
free to exploit composition, geometry or topology and still cannot separate the groups. Put the
z-score/relative-deviation sensitivity analysis in SI; it is the most technically interesting
part and the least suitable for a main figure.

---

## 3. Code to move to `main`

**Tier 1 — needed to reproduce the figures**
```
cohort.py  pair_defs.py  build_pilot_data.py  roi_geometry.py
perm_null.py  day2_permutation_test.py
day2_b7h4_precheck.py  day2_b7h4_paired.py
day2_classify_response.py
fig_counterexamples.py  day2_negative_result_figure.py
```

**Tier 2 — supporting negatives cited in the text**
```
collagen_points.py  triple_defs.py  run_triples.py  day2_barrier_test.py
day2_sensitivity.py  day2_geometry_supplement.py  day2_raw_statistics_figure.py
```

**Tier 3 — robustness work, cite as available on request or keep on the branch**
```
dowker.py  day2_dowker.py  day2_method_comparison.py
day2_length_features.py  day2_interpret_directional.py  day2_directional_controls.py
```

Not for `main`: `M2S2_demo/` (upstream clone), `data/rois/`, `output/stats/`,
`output/diagrams/`, `fibre_skeletons*/` — all regenerable and large.

---

## 4. Two things to settle before submission

1. **Cohort consistency.** Figures currently mix the 813-ROI manuscript cohort and the 808-ROI
   revision cohort. Two of the five excluded ROIs are *duplicate acquisitions of the same
   tissue* — that is the part a reviewer objects to on principle rather than magnitude, and the
   current methods note folds it into "five ROIs". Recommend moving everything except the GNN to
   808 and naming the duplicates explicitly.
2. **Endpoint consistency.** Fig 1 as drawn uses a composite endpoint (RCB response *and* not
   `Extreme_NR`, i.e. early death); every analysis here uses RCB alone. On the RCB endpoint the
   808 cohort is 37 R / 26 NR patients and 478 / 330 ROIs. State which endpoint each figure uses.
