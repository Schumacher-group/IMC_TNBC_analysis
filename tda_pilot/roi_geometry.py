#!/usr/bin/env python
"""Per-ROI geometry: the size/density covariates the Day-1 analysis never controlled.

The Day-1 "size adjustment" (`day1_diagnostics.py`) regressed features on cell COUNTS
only. But persistence lifetimes are lengths: they scale with the spacing between cells,
i.e. with density, and the imaged area varies ~54x across the cohort (convex hull
3.7e4 - 2.0e6 um^2). Counts alone therefore do not control for scale.

Writes one row per revision-cohort ROI with the geometric covariates needed to build
scale-free features and to residualise on tissue geometry:
  hull_area   convex hull of all cells in the ROI (um^2; Hyperion 1 px = 1 um)
  bbox_area   axis-aligned bounding box (um^2) -- cruder, kept as a robustness check
  med_nn      median nearest-neighbour distance over all cells (um); the natural
              length unit of the filtration, so lifetimes / med_nn are dimensionless
  density     n_all / hull_area (cells per um^2)

Counts are over ALL cells in data/rois/<fov>.csv, i.e. the union of every label in
pair_defs.ALL_LABELS -- the geometry of the imaged tissue, not of one species pair.

Output: output/roi_geometry.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree

import cohort

HERE = Path(__file__).resolve().parent
ROIS_DIR = HERE / "data" / "rois"
OUT = HERE / "output" / "roi_geometry.csv"


def geometry_for(fov: str) -> dict:
    df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
    pts = df[["x", "y"]].to_numpy(float)
    n = len(pts)

    # ConvexHull needs 3 non-collinear points; .volume is the area in 2D.
    try:
        hull_area = float(ConvexHull(pts).volume)
    except Exception:  # noqa: BLE001 - degenerate ROI, fall back to the bbox
        hull_area = np.nan
    bbox_area = float(
        (pts[:, 0].max() - pts[:, 0].min()) * (pts[:, 1].max() - pts[:, 1].min())
    )

    # k=2 because the first neighbour of a point is itself.
    dists, _ = cKDTree(pts).query(pts, k=2)
    med_nn = float(np.median(dists[:, 1]))

    return {
        "roi_id": fov,
        "n_all": n,
        "hull_area": hull_area,
        "bbox_area": bbox_area,
        "med_nn": med_nn,
        "density": n / hull_area if hull_area else np.nan,
    }


def main() -> None:
    rois = sorted(p.stem for p in ROIS_DIR.glob("*.csv") if cohort.is_cohort_member(p.stem))
    rows = [geometry_for(fov) for fov in rois]
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)

    print(f"wrote {len(df)} rows to {OUT}")
    if len(df) != 808:
        print(f"WARNING: expected 808 cohort ROIs, got {len(df)}")
    print(df[["n_all", "hull_area", "bbox_area", "med_nn", "density"]].describe().round(4).to_string())


if __name__ == "__main__":
    main()
