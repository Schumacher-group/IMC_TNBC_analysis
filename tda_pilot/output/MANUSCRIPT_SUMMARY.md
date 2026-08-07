# Bottom line for the revision

`REVIEW.md` is the full audit trail, including two retractions and the reasoning behind each.
This file is the distilled version: what can be claimed, with which numbers, and what the
figures should be. Everything below is pre-treatment, patient-level, RCB endpoint
(Responder = RCB 0–1), 25 NR / 37 R patients, 593 ROIs.

---

## 0. Definitions for the methods section

**Kernel bars.** The chromatic Delaunay–Čech filtration is built on the tumour+CD8 point cloud
with each cell coloured by type. The *kernel* diagram of the inclusion (monochromatic subcomplex
→ full complex) records features present in the single-colour structure that die once the other
colour is added. A tumour nest ringed by CD8 therefore appears as a long-lived kernel bar:
the ring persists in the tumour-only structure and is destroyed when CD8 cells are introduced.
Exclusion is thus a *joint* signature — **fewer but longer** kernel bars.

**Per-ROI null.** Each ROI is compared against itself. Cell positions are held fixed and the
tumour/CD8 labels are reshuffled 99 times, preserving the number of each type, and the
statistics are recomputed. This holds cell count, density, sampled area, tissue shape and
composition exactly constant, so any deviation reflects spatial arrangement alone.

**Relative deviation.** For statistic *s* in one ROI,

>  rd(*s*) = ( observed(*s*) − mean of the 99 shuffled values ) / mean of the 99 shuffled values

A z-score was avoided deliberately: dividing by the null *standard deviation* introduces a
denominator that shrinks as an ROI gains cells or density, making z incomparable between ROIs
(ρ with density = +0.31, against +0.13 for relative deviation).

**Exclusion score.** The two halves of the signature combined into one number:

>  exclusion = ( rd(kernel degree-1 mean bar length) − rd(kernel degree-1 bar count) ) / √2

Positive = longer and fewer bars than that ROI's own reshuffled null, i.e. more excluded.
Dimensionless, and by construction independent of the ROI's cell numbers and composition.

---

## 1. What can be claimed

### A. B7H4+ tumour differs from B7H4− tumour in its relationship with CD8  — *positive*

Paired within-ROI, both cancer subtypes subsampled to identical counts and compared against the
same CD8 cells (556 ROIs, median 391 cells per arm).

| | |
|---|---|
| degree-1 kernel bar **count** | 322.7 → 269.5, rank-biserial −0.74, p < 1e-4 |
| degree-0 kernel bar count | 243.8 → 207.0, −0.71, p < 1e-4 |
| independent confirmation (Dowker) — *internal check, NOT for the manuscript* | same direction, ρ = 0.68 between methods |
| survives matching on cell clustering | yes, both methods (p < 0.002) |

**Claim it as a difference in the CD8–tumour topological relationship, not as "more exclusion".**
The *count* signal is robust across two structurally independent constructions. The *length*
signal — which is what would license the word "exclusion" — appears in the chromatic
construction but not the Dowker one, and B7H4+ cells are separately measurable as more tightly
clustered (median nearest-neighbour 11.3 vs 14.2 µm). So the architecture differs; whether CD8
*accessibility* differs beyond that is not established.

Note this caution applies to **A only**. The figure shows count and length together because they
are the two halves of one signature, but the word "exclusion" is better supported by A2, where
both halves hold independently.

Dowker is not reported in the revision (one new method is enough), so its role here is internal:
it is why we are confident in the count half and cautious about the length half. If reviewers
ask for a second construction, `dowker.py` and `day2_method_comparison.py` are on the branch.

**No response interaction.** Present in both groups at similar magnitude (count: −38.6 in R,
−48.6 in NR, both p < 1e-4), 0 of 24 interaction tests survive BH correction. This is a
statement about tumour biology that does not depend on the response contrast — which is the
cleanest form it could take.

### A2. B7H4+ burden predicts exclusion ROI-wide — a *field* effect  — *positive, and the strongest result here*

A different hypothesis from A, and the paired design cannot see it: if B7H4+ cells suppress
infiltration across the whole microenvironment — helping neighbouring B7H4− cells evade CD8 —
then both arms of a within-ROI paired comparison sit in the same suppressed field and it
cancels. Tested at ROI level instead (`B7H4_ROI_LEVEL.md`).

| | ρ / δ | p | n |
|---|--:|--:|--:|
| **per-patient** Spearman(exclusion, B7H4+ fraction) | **+0.393** | **0.0016** | 62 patients |
| ROI-level, controlling CD8 + B7H4− counts + geometry | +0.208 | 3e-7 | 593 |
| …**also** controlling the CD8:tumour ratio | **+0.240** | 3e-9 | 593 |
| matched on (CD8, B7H4−) counts, high vs low B7H4+ | +0.288 | <1e-4 | 204 vs 202 |

**It is B7H4-specific, not "more tumour per T cell".** Exclusion is uncorrelated with the
CD8:tumour ratio (ρ = +0.017, p = 0.69), and controlling for that ratio *raises* the B7H4
association rather than lowering it.

**It does not depend on the method-dependent statistic.** Present on the count component alone
(ROI ρ = +0.302; per-patient +0.334, p = 0.008), which is the half Dowker independently
confirms, as well as on length.

**Present in both response groups** (per-patient ρ: R +0.408 p = 0.012; NR +0.369 p = 0.069 at
n = 25) and B7H4+ fraction itself does not differ by response (δ = +0.083, p = 0.59) — so this
neither creates nor is created by the response contrast.

**Caveat:** correlational. It cannot distinguish B7H4 driving exclusion from excluded tumours
upregulating B7H4, or a common cause.

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

## 2. Proposed figures  (with the script that makes each)

Dowker is left out of the revision entirely — one new method (chromatic TDA) is enough. The
code stays on the branch if reviewers press for a second construction.

### Main figure — B7H4  (two panels)
| panel | script | output | shows |
|---|---|---|---|
| **a** | `fig_b7h4_main.py --what signature` | `fig_b7h4_signature_dim1.png` (or `_dim0`) | within-ROI paired result as ONE plane: Δbar count on x, Δbar length on y. Exclusion quadrant is upper-left. |
| **b** | `fig_b7h4_main.py --what field` | `fig_b7h4_field.png` | field effect, one point per patient, ρ = +0.393, p = 0.0016 |

Panel **a** is generated at both degrees — pick per narrative:
- **degree 1** (loops; mechanistically "a ring around a nest"): Δcount −41 (87% of ROIs, p = 2e-72),
  Δlength +0.09 µm (55%, p = 0.011)
- **degree 0** (components): Δcount −30 (85%, p = 6e-69), Δlength +1.08 µm (69%, p = 2e-14)

Degree 1 matches the verbal story; degree 0 is markedly stronger on the length axis. Degree 1 is
the honest default, with degree 0 quotable as "the same pattern is stronger in degree 0".

### Main figure — the calibrated null
| panel | script | output |
|---|---|---|
| a–c | `day2_negative_result_figure.py` | `negative_result_figure.png` |

Three panels: mechanistic plane with the Fig 4B ROIs highlighted (18th vs 90th percentile);
per-ROI distributions overlapping; per-patient inference. **This is the one you flagged for
tweaking — the script above is what to edit.**

### Supplementary
| panel | script | output |
|---|---|---|
| SI-1 counterexample gallery | `fig_counterexamples.py` | `fig_counterexamples.png` |
| SI-2 B7H4 × response interaction (null) | `day2_b7h4_paired.py --reuse` | `B7H4_PAIRED.md`, `b7h4_paired.png` |
| SI-3 field effect, controls and matched design | `day2_b7h4_roi_level.py` | `b7h4_roi_level.png` |
| SI-4 B7H4 subtypes are spatially separable | `day2_b7h4_precheck.py` | `b7h4_precheck.png` |
| SI-5 classifier AUC vs permutation null | `day2_classify_response.py` | `classification_figure.png` |
| SI-6 collagen barrier test | `day2_barrier_test.py` | `barrier_test_nonclahe.png` |
| SI-7 why a cross-ROI geometry check is needed | `day2_geometry_supplement.py` | `geometry_supplement.png` |
| SI-8 raw statistics track composition and count | `day2_raw_statistics_figure.py` | `raw_statistics_figure.png` |
| SI-9 B7H4+ cells are more tightly clustered | *to make* — numbers in `METHOD_COMPARISON.md` | NN 11.3 vs 14.2 µm |

### Which null to lead with
**The classifier.** One number, no adjustment choices, permutation null centred at 0.496, and
it is the framing least vulnerable to the confounding arguments — a classifier may exploit
composition, geometry or topology and still cannot separate the groups. The
z-versus-relative-deviation sensitivity work is the most interesting part technically and the
least suitable for a main figure; keep it in SI or the response to reviewers.

## 3. Code to move to `main`

**Tier 1 — needed to reproduce the figures**
```
cohort.py  pair_defs.py  build_pilot_data.py  roi_geometry.py
perm_null.py  day2_permutation_test.py
day2_b7h4_precheck.py  day2_b7h4_paired.py  day2_b7h4_roi_level.py
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
