# Field effect: are ROIs with more B7H4+ tumour less infiltrated overall?

593 pre-treatment ROIs, 62 patients. Outcome is the per-ROI exclusion score against that ROI's own label-permutation null (higher = more excluded). This asks a different question from the paired analysis in `B7H4_PAIRED.md`, which subsamples both subtypes within an ROI and so cannot see a field effect acting on the whole microenvironment.

## Continuous

| analysis | rho | p | n |
| --- | --- | --- | --- |
| raw Spearman(exclusion, B7H4+ fraction) | 0.291 | 0 | 593 |
| partial, controlling CD8 + B7H4- counts + geometry | 0.2081 | 0 | 593 |
| raw Spearman(exclusion, B7H4+ COUNT) | 0.0913 | 0.0262 | 593 |
| per-patient Spearman(exclusion, B7H4+ fraction) | 0.3928 | 0.0016 | 62 |

## Matched on CD8 and B7H4− counts

| comparison | median_high | median_low | cliffs_delta | p | n_high | n_low |
| --- | --- | --- | --- | --- | --- | --- |
| exclusion score, high vs low B7H4+ within (CD8, B7H4-) bins | 0.4073 | 0.2875 | 0.2883 | 0 | 204 | 202 |
| CD8:tumour ratio, same bins  [sanity: should differ by construction] | 0.5955 | 0.839 | -0.1537 | 0.0074 | 204 | 202 |

The CD8:tumour row is a sanity check, not a result: adding B7H4+ cells necessarily lowers the ratio when CD8 and B7H4− are held fixed, so it must differ. The question is whether the *arrangement* score differs beyond that.

Figure: `b7h4_roi_level.png`. Script: `day2_b7h4_roi_level.py`.
