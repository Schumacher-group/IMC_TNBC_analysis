# CLAHE batch-correction evaluation -- stage 2 extraction

Host `roux`, masks `/mnt/data/Delta_Tissue/IMC/segmentation/sept2024_release/deepcell_output`.


## Channels

- markers in cell table: 38
- channels on disk: 38 processed, 37 non_processed
- derived, dropped from both: ['Carboplatin']
- **shared channel list: 37**

## Verification: does this extraction reproduce the published table?

```
          fov                 definition  min_channel_r  worst_channel_rel_err worst_channel
0  Leap001_10  mean over mask (sum/area)            1.0           1.108464e-07     Alpha-SMA
1  Leap001_10     sum / cell_size column            1.0           1.108464e-07     Alpha-SMA
2   Leap001_8  mean over mask (sum/area)            1.0           1.474047e-07     Alpha-SMA
3   Leap001_8     sum / cell_size column            1.0           1.474047e-07     Alpha-SMA
4   Leap001_9  mean over mask (sum/area)            1.0           1.306679e-07     Alpha-SMA
5   Leap001_9     sum / cell_size column            1.0           1.306679e-07     Alpha-SMA
```

Per candidate definition (worst channel, worst FOV):

```
                           min_r  worst_rel_err
definition                                     
mean over mask (sum/area)    1.0   1.474047e-07
sum / cell_size column       1.0   1.474047e-07
```

**Verified**: `mean over mask (sum/area)` reproduces the published intensities (min per-channel r 1.000000, worst-channel relative error 1.47e-07). The uncorrected table will be built with the identical procedure, so CLAHE is the only difference between them.

## Extraction

- FOVs to process: 829
- cells: 4,241,908
- reused from cache: 829
- newly extracted: 0
- cells with no matching mask object (NaN): 0

## Assembly

- **asserted**: identical cell IDs, order, labels, Pixie labels, coordinates and channel order between the two tables
- corrected: transformed (nan->0, log1p, q0.95 per-channel cap, scaled)
- uncorrected: transformed (nan->0, log1p, q0.95 per-channel cap, scaled)
- cells whose raw intensities differ between the tables: 4,241,908 of 4,241,908 (100.0%)

## Cohort flags recorded in .obs

```
           all  revision cohort  revision + pre-treatment
cells  4241908          4183318                   3227482
ROIs       829              808                       593
```

Nothing is filtered here. Stage 3 selects the cohort from these flags.

Written to `stage2_output`.
