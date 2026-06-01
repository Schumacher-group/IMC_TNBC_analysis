# Phase 3 validation — full TDA pilot

Cohort: 808 ROIs. Stat columns expected: 456 (6 diagrams x 2 dims x 38 stats).

## Per-definition checks

| definition | ok | degenerate | failed | completeness_808 | finite_all_ok | nondegen_ok | ker1-num_bars | im1-num_bars |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD8_primary | 808 | 0 | 0 | True | True | True | 1.0 | 1.0 |
| CD8_strict | 808 | 0 | 0 | True | True | True | 1.0 | 1.0 |
| CD8_extended | 808 | 0 | 0 | True | True | True | 1.0 | 1.0 |
| CD4 | 803 | 5 | 0 | True | True | True | 1.0 | 1.0 |
| Bcell | 808 | 0 | 0 | True | True | True | 1.0 | 1.0 |
| Macrophage | 806 | 2 | 0 | True | True | True | 1.0 | 1.0 |
| Fibroblast | 807 | 1 | 0 | True | True | True | 1.0 | 1.0 |

### Non-degeneracy detail (fraction of ok ROIs with >0 bars)

| definition | ker0-num_bars | ker1-num_bars | im0-num_bars | im1-num_bars |
| --- | --- | --- | --- | --- |
| CD8_primary | 1.0 | 1.0 | 1.0 | 1.0 |
| CD8_strict | 1.0 | 1.0 | 1.0 | 1.0 |
| CD8_extended | 1.0 | 1.0 | 1.0 | 1.0 |
| CD4 | 1.0 | 1.0 | 1.0 | 1.0 |
| Bcell | 1.0 | 1.0 | 1.0 | 1.0 |
| Macrophage | 1.0 | 1.0 | 1.0 | 1.0 |
| Fibroblast | 1.0 | 1.0 | 1.0 | 1.0 |

## Summary

- Completeness (sum==808, unique roi_id) for all 7: **True**
- Finiteness (456 finite stats per ok row) for all 7: **True**
- Non-degeneracy (ker/im dim0&1 >0 in >=90% ok ROIs) for all 7: **True**
- Degenerate ROIs (a species had <3 cells): CD4=5, Macrophage=2, Fibroblast=1.
- Spot-check diagram overlay saved: `output/spotcheck_diagrams.png` (responder Leap066_11, non-responder Leap001_10).
