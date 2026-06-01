#!/usr/bin/env python
"""One pass over the 8 GB cell table -> per-ROI point clouds for the full pilot.

For every revision-cohort ROI, writes data/rois/<fov>.csv with columns
(x, y, celltype) where celltype is the original `cell_meta_cluster` string,
restricted to the labels we ever need (pair_defs.ALL_LABELS). The per-definition
species aggregation/relabelling happens later, per worker, in run_full_pilot.py.

Coordinates: centroid-1 -> x, centroid-0 -> y (micrometres; Hyperion = 1 px/um).
No qc_pass filter (per user: use all cells). Non-cohort ROIs are skipped.
"""

from pathlib import Path

import pandas as pd

import cohort
import pair_defs

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CELL_TABLE = ROOT / "CellTable_CleanCohort" / "updated_cell_table.csv"
USECOLS = ["fov", "cell_meta_cluster", "centroid-0", "centroid-1"]
CHUNK = 500_000


def main() -> None:
    out_dir = HERE / "data" / "rois"
    out_dir.mkdir(parents=True, exist_ok=True)

    keep_labels = pair_defs.ALL_LABELS
    parts = []
    n_rows = 0
    for chunk in pd.read_csv(CELL_TABLE, usecols=USECOLS, chunksize=CHUNK):
        n_rows += len(chunk)
        sub = chunk[chunk["cell_meta_cluster"].isin(keep_labels)]
        # cohort filter (vectorised via the fov string)
        sub = sub[sub["fov"].map(cohort.is_cohort_member)]
        if len(sub):
            parts.append(sub)
    df = pd.concat(parts, ignore_index=True)
    print(f"scanned {n_rows:,} rows; kept {len(df):,} cohort cells of {len(keep_labels)} relevant labels")

    written = 0
    for fov, g in df.groupby("fov"):
        out = pd.DataFrame(
            {
                "x": g["centroid-1"].astype(float),
                "y": g["centroid-0"].astype(float),
                "celltype": g["cell_meta_cluster"].astype(str),
            }
        )
        out.to_csv(out_dir / f"{fov}.csv", index=False)
        written += 1
    print(f"wrote {written} per-ROI CSVs to {out_dir}")
    if written != 808:
        print(f"WARNING: expected 808 cohort ROIs, wrote {written}")


if __name__ == "__main__":
    main()
