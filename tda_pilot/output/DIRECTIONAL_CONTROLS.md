# Directional null, control species — is the geometry artefact generic?

Pre-treatment, per-patient Cliff's δ (R vs NR) on `ker1-avg_length`, 25 NR / 37 R patients. `_adj` = after residualising on non-tumour fraction, log(cells), log(density), log(hull area).

| definition | n_ROI | z_raw_d | z_raw_p | z_adj_d | z_adj_p | reldev_raw_d | reldev_raw_p | reldev_adj_d | reldev_adj_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 593 | 0.380 | 0.012 | 0.008 | 0.966 | 0.356 | 0.019 | -0.044 | 0.774 |
| CD4 | 591 | 0.150 | 0.322 | -0.090 | 0.556 | 0.131 | 0.389 | 0.005 | 0.977 |
| Fibroblast | 593 | -0.144 | 0.344 | -0.068 | 0.656 | -0.163 | 0.282 | -0.157 | 0.301 |
| Macrophage | 593 | 0.304 | 0.044 | 0.295 | 0.051 | 0.284 | 0.060 | 0.165 | 0.276 |

## 1. The geometry dependence of z is a property of the method, not of CD8

| definition | rho_z_density | rho_z_cells |
| --- | --- | --- |
| CD8_primary | 0.313 | 0.214 |
| CD4 | 0.303 | 0.385 |
| Fibroblast | 0.249 | 0.336 |
| Macrophage | 0.177 | 0.227 |

z correlates with ROI density (ρ = 0.18–0.31) and cell count for **every** species. This is the retraction's mechanism confirmed: z = (obs − mean)/sd, the null sd shrinks as an ROI gains cells or density, so z is comparable within an ROI but not across ROIs of differing geometry. Nothing about CD8 is special here.

## 2. Does anything survive adjustment?

**0 of 4 species survive on both estimators.**

- **CD8_primary** — the only species whose raw effect was significant, and it collapses completely (+0.379 → +0.008 on z; +0.356 → −0.044 on reldev). Retracted.
- **CD4, Fibroblast** — null raw and null adjusted. No signal to explain.
- **Macrophage** — raw δ = +0.304 (p = 0.045) barely moves under z-adjustment (+0.295, p = 0.051) but does not hold on the size-free estimator (+0.165, p = 0.275). Across 4 species × 2 estimators this is 1 borderline result in 8 tests, which is what chance produces. **No claim.** If anyone wants to pursue it, it needs its own pre-registered endpoint and a fresh cohort — it is not evidence as it stands, and it was not the hypothesis under test.

## Conclusion

The controls corroborate the retraction. The directional statistic's z-score carries a generic geometry dependence, the CD8 effect is fully explained by it, and no species shows a response difference that survives adjustment on both estimators. The pilot's original no-go stands.

Script: `day2_directional_controls.py`. Inputs: `perm_null_<definition>_other_only.parquet` (593 ROIs × 100 six-packs each).
