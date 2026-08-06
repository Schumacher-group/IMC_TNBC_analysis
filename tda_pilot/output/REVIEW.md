# Critical review of the TDA pilot, and the test that settles it

**Date:** 2026-08-05. Reviews the Day-0/Day-1 pilot (`REPORT.md`, `DAY1_REPORT.md`,
`DAY1_DIAGNOSTICS.md`, `IMBALANCE_CONFOUND.md`, `PRETREATMENT_CONFOUND.md`) and reports the
follow-up analyses run to check it (`LENGTH_FEATURES.md`, `PERMUTATION_NULL.md`).

---

> **Verdict confirmed, with one caveat on wording.** The directional inclusion was
> subsequently run (UPDATE 1) and initially looked positive, but does not survive cross-ROI
> geometry adjustment (UPDATE 2) — so the conclusion below stands, and has now survived a
> deliberate attempt to break it. The one correction: statements here that the per-ROI
> permutation null "holds geometry fixed" are true of the null mean but not of the z-score's
> scale across ROIs. See UPDATE 2 for why that matters.

## Verdict

**The pilot's conclusion stands, and is now on much firmer ground than it was.**

The original claim — *at baseline, tumour–CD8 spatial organisation does not distinguish
responders from non-responders in this cohort* — survives every check run here. But it was
originally reached from 12 of 456 computed statistics, via a single feature that turns out to
be a poor choice, with tissue geometry never controlled. Those gaps are now closed, and the
answer did not change.

The important new fact is not the null. It is this: **tumour and CD8 cells are strongly
non-randomly arranged in both responders and non-responders** — 42% of ROI × statistic
z-values exceed |z| = 2 against a composition- and geometry-matched null, with the most
deviant statistics at median |z| ≈ 17–18. The assay has plenty of resolving power. What it
does not resolve is *response*. That is a much stronger statement than "we found nothing".

---

## What holds up

- **The imbalance confound (`IMBALANCE_CONFOUND.md`) is real and correctly diagnosed.**
  `im_dim1_total_persistence` tracks minority fraction (ρ = −0.54) and the effect vanishes on
  balanced ROIs. Preferring a subset test to a linear covariate adjustment was the right call,
  and the note says why in terms that are still correct.
- **The treatment-stage confound (`PRETREATMENT_CONFOUND.md`) is real and decisive.**
  Post-treatment resections are 170 NR vs 7 R; pooling manufactures the "responders lower"
  effect. Catching this was the single most valuable thing the pilot did.
- **The engineering is sound.** Cohort definition, per-ROI failure isolation, resume,
  validation, and the size-matched residual galleries all check out. `Patient_ID` is the
  correct patient key and was used correctly (it merges core/resection pairs — LEAP001+LEAP002
  are one patient).
- **The wider re-scan agrees.** Pre-treatment, across 8 definitions × 10 features × 4 analysis
  levels (320 comparisons), nothing reaches nominal significance (min p = 0.056).

## What was wrong or missing

1. **The headline feature was a poor choice.** Median degree-1 image bar lifetime is 0.33 µm —
   far below one cell diameter — and only ~3% of bars clear the 5 µm threshold, yet
   `im_dim1_total_persistence` sums all of them. Since it factorises as
   `n_bars × mean lifetime` and only the count term was ever adjusted for, the feature largely
   re-encodes the very cell-number and imbalance confounds the pilot was fighting.
   The feature dict was also internally inconsistent (`total_persistence` over all bars,
   `n_features` over bars > 5 µm). Fixed: `extract_features.py` now emits `n_bars_all` and
   `n_bars_gt5um`.
2. **ROI area and density were never controlled** — only cell counts. Hull area spans 54×
   (3.7e4–2.0e6 µm²) and pre-treatment responder ROIs are significantly less dense (per-ROI
   Cliff δ = −0.28, p < 1e-3). This cut both ways: geometry was uncontrolled when the signal
   was claimed *and* when it was dismissed. Fixed: `roi_geometry.py`.
3. **444 of 456 computed statistics were never examined**, and the planned Day-2 hierarchical
   model was never run — so the go/no-go came from univariate descriptives that `DAY1_REPORT.md`
   itself flags as insufficient.
4. **No null model.** Confounding was handled entirely by regression and subsetting across
   ROIs, which cannot work when composition, density, area and dispersion are collinear with
   each other and with response. Fixed: `perm_null.py`.
5. **A metadata bug dropped 11 ROIs.** `Leap084a_1..11` parse to `LEAP084A`, absent from the
   metadata, so they lost their response label — this is why the pooled cohort is 797 and not
   808. They are the second tissue block of the LEAP084 core (`Leap084_b_12..21` is the first;
   only the separator differs), and `results_analysis/spatial_stats/conditions_pretreat.json`
   independently lists 433 pre-treatment responder ROIs against this table's 422 — exactly
   these 11. `roi_counts.py` now carries the alias, **not yet applied** (see Open decisions).
6. **Arithmetic in `PRETREATMENT_CONFOUND.md` did not reconcile** (797 ≠ 582 + 177): 38
   `mid-early`/`mid-late` ROIs, all responders, are in the pooled row but dropped from the split.
7. **Most p-values were per-ROI** (n ≈ 580–800) despite 63 patients, and pre-treatment ROIs per
   patient differ systematically (NR 6.4, R 11.4), so ROI pooling over-weights responders.

**One thing the pilot warned about that turned out not to be a problem:** the note that "any
pooled pre+post analysis (incl. cross-PCF) inherits the same confound". Checked —
`export_for_spoox.py` already stratifies correctly (`conditions_pretreat.json` is R vs NR
*within pre-treatment*; `conditions_nonresponder.json` is pre vs post *within non-responders*),
and neither config contains an out-of-cohort ROI. The warning applies to future pooled work only.

---

## The candidate signal, and why it did not survive

Decomposing `total_persistence` into its count and length terms surfaced something the pilot
never looked at. Pre-treatment, per-patient, the bar-*length* statistics are CD8-specific and
responder-higher across all three CD8 aggregate definitions (Cliff δ = +0.27 … +0.41,
p = 0.006 … 0.08), with the CD4, B-cell, fibroblast and NK/CD8 controls flat (macrophage was
intermediate, δ ≈ +0.18…+0.20, all n.s.) — and, unlike the signal the pilot killed, it
*strengthened* on balanced ROIs and was uncorrelated with minority fraction (ρ = 0.013). Note
also that no test survived BH correction across the 240-comparison scan even before any null
model was applied. Full table in `LENGTH_FEATURES.md`.

Two things argued against reading that as infiltration biology: it attenuated under geometry
adjustment, and the same CD8-specific pattern appeared in `dom0`/`dom1` — the *monochromatic*
complex, which contains no cross-species mixing information at all. That pointed at tissue
dispersion rather than arrangement, but regression could not separate the two.

**The permutation null separates them.** Holding each ROI's cell positions, counts, density,
area and dispersion exactly fixed and shuffling only the labels, the primary endpoint
(`ker1-avg_length` z, per-patient median) gives:

| | δ | p |
|---|---|---|
| descriptive, `/med_nn` (`LENGTH_FEATURES.md`) | +0.336 | 0.026 |
| **conditioned on the per-ROI null (`PERMUTATION_NULL.md`)** | **+0.174** | **0.25** |

Bootstrap CI [−0.140, +0.472]; 0 of 50 retained statistics survive BH. The signal was
substantially tissue geometry, exactly as the `dom0` diagnostic suggested.

**Honest caveats**, all recorded in `PERMUTATION_NULL.md`: the size-free companion estimator
gives δ = +0.269, p = 0.075 — borderline, so this is a soft null rather than a clean one; and
in a post-hoc scan the effect rises monotonically to δ = +0.34 (p = 0.04) on well-balanced
ROIs. Neither is evidence for an effect. Both are the places to look first if this is revisited.

---

## Statements you can defend

Usable again — the directional result that briefly contradicted this did not survive geometry
adjustment (UPDATE 2). One wording change from the original draft: the claim is that the
comparison is not explained by composition *or* by sampled-region geometry, and the honest
basis for the geometry half is the explicit residualisation, not the permutation null alone.

For the rebuttal, pre-treatment, per-patient, 25 NR vs 37 R:

> Baseline tumour–CD8 spatial organisation does not distinguish responders from
> non-responders in this cohort by composition (Cliff δ = +0.008, p = 0.97), and we found no
> robust evidence that it does so by chromatic topology. Each region was compared against its
> own label-permutation null, which holds cell positions, counts and composition fixed. The
> colour-symmetric formulation is null throughout (δ = +0.17, p = 0.27), as are CD4 and
> fibroblast comparators. A CD8-directional formulation gives δ = +0.38 (p = 0.011)
> unadjusted, but this is not separable from between-region differences in sampled area and
> cell density: across reasonable geometry adjustments the estimate spans δ = +0.01 to +0.40,
> and we therefore do not claim it. Notably, its direction is opposite to the
> immune-exclusion hypothesis — responders, not non-responders, score higher on the
> exclusion-like axis — so no analysis supports reduced cytotoxic T-cell infiltration in
> non-responders at baseline. The measurement is not insensitive: applied to
> individual regions selected for contrasting architecture, it separates a finely intermixed
> region from one with CD8 confined to the stroma around compact tumour nests by more than
> eight standard deviations of their respective nulls, with the excluded region's mean
> degree-1 kernel lifetime 48% above its own composition-matched null versus 8% for the
> intermixed one. Tumour and CD8 cells are strongly non-randomly arranged in both response
> groups (42–45% of region × statistic comparisons exceed |z| = 2). The absence of a robust
> group-level difference therefore reflects the data rather than a limitation of the method.

**Do not** state this as "no effect". With 25/37 patients, ~80% power arrives only at
|δ| ≈ 0.42; a moderate effect is not excluded. Say "no large difference".

---

## UPDATE 2 (2026-08-05, latest): the directional result does NOT survive — retracted

**UPDATE 1 below is wrong and is retained only for the record.** The directional effect is
carried by ROI geometry and collapses when that is adjusted for:

| estimator | raw | after adjusting for CD8 fraction, log cells, log density, log hull area |
|---|--:|--:|
| `ker1-avg_length` z | +0.379 (p = 0.012) | **+0.008 (p = 0.97)** |
| size-free companion | +0.356 (p = 0.019) | **−0.044 (p = 0.77)** |

**The error was mine, in the design of the null, and it is worth stating precisely because it
is easy to repeat.** A per-ROI label permutation holds composition, counts, density and area
fixed *within* each ROI, so the null **mean** is unbiased — that part was right. But
z = (obs − mean)/sd, and the null **sd** shrinks as an ROI gains cells or density. The
arrangement-to-z mapping is therefore itself geometry-dependent: **z is comparable within an
ROI, not across ROIs of different size and density.** Responder ROIs are systematically denser
and larger (Cliff δ = −0.28 on density), and that alone manufactured the effect. Every
statement of the form "conditioning on the per-ROI null removes geometry" in this document and
in `PERMUTATION_NULL*.md` was too strong.

Diagnostics confirming the mechanism: ρ(z, density) = +0.31, ρ(z, cells) = +0.21,
ρ(z, hull area) = +0.20 for the directional statistic. Composition itself is clean —
ρ(z, CD8 fraction) = +0.007, and per-patient CD8 fraction does not differ by response
(δ = +0.023, p = 0.89) — so this is a density/size artefact, not a composition one. (Note the
*symmetric* statistic is worse on composition: ρ(z, CD8 fraction) = −0.45.)

`day2_permutation_test.py` now runs this cross-ROI geometry adjustment as a required check and
will not report an endpoint as significant unless it survives. Both inclusions now report
**"Survives geometry adjustment: False"**.

**The control species confirm it** (`DIRECTIONAL_CONTROLS.md`, 3 × 593 ROIs × 100 six-packs):
z correlates with ROI density for *every* species (ρ = 0.18–0.31) and with cell count
(ρ = 0.21–0.39), so the geometry dependence is a property of the method, not of CD8.
**0 of 4 species survive adjustment on both estimators.** CD4 and fibroblast are null raw and
adjusted. Macrophage is the one wrinkle — raw δ = +0.304 (p = 0.045) barely moves under
z-adjustment (+0.295, p = 0.051) but does not hold on the size-free estimator (+0.165,
p = 0.275); across 4 species × 2 estimators that is 1 borderline result in 8 tests, i.e. what
chance produces. No claim, and it was not the hypothesis under test.

**Net position: the pilot's original no-go stands.** It stands for better reasons than the
pilot had, and it survived a genuine attempt to break it with the directional statistic.

**What the directional statistic is still good for.** It is a real and sensitive arrangement
readout — `directional_gallery.png` shows its extremes separate large-scale tumour/CD8
territorial segregation from diffuse tumour with sparsely intermixed CD8. It is also cleaner
on composition than the symmetric one. It just cannot be compared across ROIs of different
geometry without adjustment, and after adjustment there is no response difference.

---

## UPDATE 5 (2026-08-06): the multivariate analysis, finally run — AUC 0.544, p = 0.275

The Day-1 report deferred a multivariate model and it was never run. It has now been run, using
the M2S2 authors' own machinery on the 456-statistic table (`CLASSIFICATION.md`).

Pre-treatment, patient-level, 62 patients (25 NR / 37 R), 1,520 statistics across 8 pair
definitions after their `discard_feature` filter, reduced by per-cell-group PCA **inside** the
cross-validation:

| model | features | CV AUC |
|---|--:|--:|
| TDA, all 8 pair definitions | 1520 | 0.544 |
| TDA, CD8_primary only | 190 | 0.573 |
| Composition + geometry only | 7 | 0.559 |
| DummyClassifier | 0 | 0.495 |

Label-permutation null (200 draws): mean AUC 0.496, 95th pct 0.630, observed 0.544,
**empirical p = 0.275**.

**Why this is the most useful negative in this directory.** Every other result here depends on
an adjustment choice — which covariates, matched or regressed, z or relative deviation. A
classifier is free to exploit composition, geometry, topology or any combination, and it still
cannot separate the groups. There is nothing to adjust and nothing to argue about.

Two supporting checks: the null is centred at 0.496, which confirms the cross-validation is not
leaking (the original M2S2 script fits PCA before `cross_validate`, which would push the null
above 0.5 — we moved PCA inside the pipeline); and the null's 95th percentile is 0.630, so this
design detects only AUC ≳ 0.63 — "no signal" here means "no strong signal".

**Note the M2S2 pipeline answers a different question from ours.** They fit one classifier per
patient, using that patient's ROIs as samples, to predict each ROI's `sample_type` — a
within-patient question where patient identity is controlled by construction. We need a
between-patient question (predict the patient's response), so the unit is the patient. Their
machinery transfers; their study design does not.

---

## UPDATE 4 (2026-08-06): the retraction in UPDATE 2 was itself over-confident — indeterminate

UPDATE 2 concluded the directional effect was a geometry artefact, on the strength of a
four-covariate regression adjustment that drove δ from +0.379 to +0.008. That adjustment does
not hold up, and the honest answer is that **this design cannot resolve the question either
way**. Full detail: `SENSITIVITY.md`.

**Regression adjustment is unstable**, because `density = cells / area` makes the covariates
near-collinear (condition number **448**):

| adjustment | δ (z) | p |
|---|--:|--:|
| none (pre-specified endpoint) | +0.379 | 0.012 |
| log_density + log_cells | +0.343 | 0.023 |
| log_area only | +0.250 | 0.099 |
| log_cells + log_area (minimal, non-redundant) | +0.226 | 0.136 |
| all four (near-collinear) | **+0.008** | **0.966** |

The four covariates jointly explain only ~9% of the variance of z, yet extinguish a moderate
group difference. A 9%-of-variance block doing that is a signature of over-adjustment under
collinearity, not of confounding correctly removed.

**Matching does not collapse.** Restricting to a common band of area *and* density needs no
model and cannot suffer collinearity:

| area+density band | n patients (NR/R) | δ (z) | p |
|---|---|--:|--:|
| all ROIs | 25/37 | +0.379 | 0.012 |
| 10–90 pct | 24/34 | +0.270 | 0.084 |
| 20–80 pct | 20/29 | +0.376 | 0.027 |
| IQR | 14/25 | +0.326 | 0.098 |

δ stays at +0.27…+0.38 in every band. A genuine geometry artefact shrinks toward zero under
matching — that is exactly what the Day-1 signal did in `IMBALANCE_CONFOUND.md` (−0.138 → +0.035).
This one does not.

**So: indeterminate.** Across reasonable specifications δ spans +0.01 … +0.40 (p from 0.009 to
0.97). Matching, the more trustworthy arbiter, leaves the effect intact but underpowered.
Neither "real" nor "artefact" is established. Resolving it needs more patients, or ROI sampling
matched on area and density *by design* — not further reanalysis of these data.

**What is stable, and matters more than the significance: the direction.** Responders score
*higher* on the exclusion-like axis under every single specification. The immune-exclusion
hypothesis predicts the opposite. So **no specification supports "reduced infiltration of
cytotoxic T cells in non-responders"**; the open question is only whether there is a weak
effect running the other way.

`day2_permutation_test.py`'s built-in geometry check uses the collinear four-covariate set and
therefore reports "Survives geometry adjustment: False" — read `SENSITIVITY.md` instead of that
flag.

---

## UPDATE 3 (2026-08-06): the Fig 4B panels — the statistic IS a valid exclusion readout

The obvious worry about a negative result is that the measurement is blind. It is not, and the
cleanest demonstration is the manuscript's own illustrative figure.

**The two panels of Fig 4B are `Leap008_5` (responder, finely intermixed) and `Leap046_1`
(non-responder, compact tumour nests with CD8 confined to the surrounding stroma).** Both are
pre-treatment members of the 593-ROI cohort, so both already carry label-permutation nulls from
the main run and need no special-case computation. (Their extents, 977 × 971 µm and
1060 × 1075 µm, match the axes printed on the figure panels.) Scored against those nulls:

| statistic (directional) | Leap008_5 (R) | Leap046_1 (NR) | gap |
|---|--:|--:|--:|
| `ker1-avg_length`, geometry-adjusted z | −1.05 (47th pct) | **+7.68 (94th pct)** | **+8.7** |
| `ker0-avg_length`, geometry-adjusted z | −4.91 (25th pct) | **+6.73 (86th pct)** | **+11.6** |
| `ker1-num_bars`, geometry-adjusted z | +4.42 (79th pct) | −5.42 (13th pct) | −9.9 |

Raw observed values make it concrete: the excluded ROI's mean kernel bar is **48% longer than
its own composition-matched null** (3.24 vs 2.19 µm), the intermixed one only 8% longer
(1.73 vs 1.60 µm).

**Four things follow.**

1. **Face validity: the statistic measures what the eye calls exclusion**, and in the direction
   the geometry predicts. A tumour nest ringed by CD8 puts a long-lived loop in the CD8-only
   complex which the tumour then fills — so exclusion produces *fewer but much longer* kernel
   bars. That is precisely the signature: `ker1-num_bars` down, `ker1-avg_length` up. The
   intermixed ROI shows the reverse — many short-lived loops from fine-grained mixing.
2. **It is not geometry.** The two ROIs are near-identical in density (0.0053 vs 0.0051
   cells/µm²), hull area (0.94 vs 1.02 × 10⁶ µm²) and median NN spacing (10.3 vs 9.7 µm), and
   the gap *widens* under geometry adjustment (+7.5 → +8.7). Contrast this with the cohort-level
   effect in UPDATE 2, which vanished under the same adjustment.
3. **It corroborates the retraction independently.** The cohort trend (now retracted) had
   *responders* higher. On a pair chosen for genuine visual exclusion contrast, the sign
   **reverses** and the effect is an order of magnitude larger. A real exclusion contrast and
   the spurious cohort trend point in opposite directions — which is what you expect if the
   latter was a geometry artefact rather than a diluted version of the former.
4. **The cohort null is calibrated, not blind.** The method detects exclusion decisively when
   it is present. So the absence of a response difference across 593 ROIs / 62 patients means
   there is no *consistent* response-associated difference — not that the instrument cannot
   see one. This is the strongest form the negative result can take, and it is what makes the
   rebuttal paragraph defensible.

The honest framing for the manuscript: **Fig 4B is a fair illustration of a real and
quantifiable difference** — the panels are separated by ~9 sd of their own nulls, in the
direction the figure's caption asserts, and this is now measured rather than asserted from
visual inspection. What does not hold is the generalisation from those two panels to a
systematic baseline difference between response groups. Both statements can be made together,
and making them together is stronger than either alone: the illustration is genuine, and the
cohort-level claim it implies is not supported.

If the figure caption ("Reduced infiltration of cytotoxic T-cells in non-responders") is meant
as a claim about the groups rather than about the panels, it needs rewording. A defensible
version: *"Example regions illustrating the range of tumour–CD8 architectures observed, from
finely intermixed (responder, Leap008_5) to CD8 confined to peritumoural stroma
(non-responder, Leap046_1); this architectural range is not systematically associated with
response at baseline (see [rebuttal statement])."*

Script: `day2_interpret_directional.py --rois Leap008_5 Leap046_1`. Figure:
`directional_named_pair.png`. Full table: `DIRECTIONAL_INTERPRETATION.md`.

---

## UPDATE 1 (2026-08-05): RETRACTED — see UPDATE 2 above

The recommendation below was written before the directional run. This section reported that
run as a positive finding; **it does not survive cross-ROI geometry adjustment and should not
be relied on.** Retained unedited for the record.

**On the identical 593-ROI pre-treatment cohort, matched patient set:**

| inclusion | primary endpoint `ker1-avg_length` z, per-patient | δ | permutation p |
|---|---|--:|--:|
| `symmetric` — `KChromaticInclusion(filt, 1)`, what the entire pilot used | | +0.165 | 0.27 |
| **`other_only` — `SubChromaticInclusion(filt, [[CD8]])`, directional** | | **+0.379** | **0.011** |

Bootstrap 95% CI [+0.094, +0.643]; leave-one-patient-out +0.356 … +0.435; size-free companion
δ = +0.356, p = 0.019; nominally significant at every rung of the balance ladder from 0.00 to
0.25. This was a **single pre-specified endpoint**, fixed in writing and committed to git
(896eda7) before the run, so it does not need multiplicity correction — unlike the exploratory
scans elsewhere in this document. Full detail: `PERMUTATION_NULL_OTHER_ONLY.md`.

The dilution argument was therefore correct and quantitatively material: the symmetric
inclusion, whose domain is the monochromatic subcomplex of *both* colours, roughly halves this
effect and hides it. **Every negative result in this document above was computed with that
diluted statistic.**

### But this is not immune exclusion

The pre-specified direction matters. Exclusion in non-responders predicts **NR > R** on this
endpoint (δ < 0), because a tumour nest ringed by CD8 puts a long-lived loop in the CD8-only
complex which the tumour then fills. **The observed δ is positive — responders are higher.**
Whatever this is, it is the opposite of the exclusion hypothesis, and must not be written up
as support for it. `day2_permutation_test.py` refuses to phrase it that way.

### What is required before this is a finding

Per the decision rule fixed before the run, and not yet done:

1. **Control species under the same directional null** (Macrophage, Fibroblast, CD4). This is
   the load-bearing check: if the controls move too, the effect is general tissue architecture,
   not CD8 biology. ~3 h per species. **Nothing should be claimed until this is done.**
2. **Mechanistic interpretation.** "Responders have longer-lived CD8-derived kernel features"
   is a statement about barcodes, not about biology. The kernel bar's birth and death both
   depend on the filtration in ways that need working through on real ROIs — ideally with a
   size-matched point-cloud gallery sorted on this statistic, as `day1_residual_gallery.py`
   did for the old one — before any mechanism is asserted.
3. **The power limit is unchanged.** 25 NR / 37 R patients; δ = +0.379 sits just below the
   ~0.42 that this design detects with 80% power, so it is a borderline-powered result on a
   single endpoint. It should be treated as a hypothesis worth testing, not an established one.

---

## Superseded recommendation (kept for the record): one more run, then stop

**Stop chasing the length statistics.** Nothing there survives a null model, and further
digging in that direction is how false positives are manufactured.

**But run the directional inclusion once.** Every result above — Day-1 and Day-2 alike — uses
`KChromaticInclusion(filt, 1)`, whose domain is the monochromatic subcomplex of *both* colours.
That six-pack is colour-symmetric, and with tumour outnumbering CD8 ~2:1 its kernel is
tumour-dominated. It is not blind to exclusion, but it is diluted.

The geometry matters here. Immune exclusion means a compact tumour nest ringed by CD8. Those
CD8 cells form an annulus — a degree-1 loop in the CD8-only complex — which dies as soon as the
tumour cells filling the nest are added, i.e. it appears as a **long degree-1 kernel bar in the
CD8-directional six-pack**. `SubChromaticInclusion(filt, [[0]])` isolates exactly that, and it
has never been run. It is the most obvious reviewer objection to the negative result above
("you measured the wrong statistic"), and one overnight run removes it permanently.

Implemented as `perm_null.py --inclusion other_only` (smoke-tested; the `symmetric` path is
unchanged and still reproduces `perm_null_CD8_primary.parquet` exactly). **Not yet run.**

Pre-specify before running, so this stays one test and not a forking path:
- **Endpoint:** `ker1-avg_length` z, per-patient median, pre-treatment, `CD8_primary`.
- **Direction:** exclusion in non-responders predicts **NR > R**, i.e. Cliff's **δ < 0**.
  Worth noting up front: every trend measured under the symmetric inclusion runs the *other*
  way (δ = +0.17 conditioned on the null, +0.34 descriptively), so the signal chased above is
  not diluted exclusion. That lowers the prior on a positive result but does not make the test
  redundant.
- **Cost:** 582 ROIs × 100 six-packs ≈ 3.3 h on 10 workers, as for the symmetric run.
- **Read a null as "no large effect"**, per the power limit below — not as "no effect".

After that, stop. The remaining constraint is **patients**: 25 NR / 37 R, ~80% power only at
|δ| ≈ 0.42. That is not something a better statistic can fix.

## Open decisions for you

1. **Apply the `LEAP084A` alias?** It adds 11 pre-treatment responder ROIs to patient 44 and
   makes the labelled cohort 808 rather than 797, bringing this pipeline into line with the
   SpOOx configs. It re-cuts the cohort under every published Day-1 number, so it should be
   adopted deliberately and everything re-run together, or not at all.
2. **Regenerate the Day-1 reports?** `DAY1_DIAGNOSTICS.md` was produced before the
   `n_features` → `n_bars_gt5um` rename. Values are unchanged; only the column name differs.

Scripts added: `roi_geometry.py`, `day2_length_features.py`, `perm_null.py`,
`day2_permutation_test.py`. Compute: the permutation null was 582 ROIs × 100 six-packs,
198 min on 10 workers.
