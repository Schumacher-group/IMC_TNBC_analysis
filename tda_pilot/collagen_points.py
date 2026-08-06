"""Turn coarse collagen fibre masks into a point cloud fit for chromatic TDA.

THE PROBLEM WITH THE OBJECT TABLE
---------------------------------
`fiber_object_table.csv` gives one centroid per fibre object, and that will not do. The
segmentation is bundle/blob level, not strand level: across 677 ROIs the single largest object
holds a median 20.7% of that ROI's collagen area (p90 44.5%), and the largest 10% of objects
hold 55.3%. For the Fig 4B panels the largest object is 31-32% of the collagen. Representing a
third of an ROI's collagen by one point is meaningless, and 61 points per ROI (median) is far
too few to carry degree-1 topology alongside ~1000-3000 cells. 15% of ROIs have under 30
objects at all.

WHY OCCUPANCY SAMPLING RATHER THAN SKELETONISATION
--------------------------------------------------
Skeletonising would place points on each blob's medial axis and discard its WIDTH. But for a
barrier hypothesis the width is the point -- a thick band of collagen is a thicker barrier than
a thin one, and a medial axis represents both by the same one-dimensional curve. Since the
segmentation is blob-level anyway, a skeleton would be describing the medial axis of an
arbitrary merge of nearby strands rather than any real fibre. Sampling the occupied area keeps
the representation honest about what the segmentation actually resolves: where collagen is,
not how it is threaded.

DENSITY IS A DESIGN CHOICE, SO IT IS MADE EXPLICITLY
----------------------------------------------------
The whole of this directory's history says these statistics are density-sensitive. So collagen
points are placed on a regular grid of fixed spacing intersected with the mask, which makes the
collagen point density a parameter we set rather than an artefact of segmentation. The default
spacing is chosen to match cell density: cohort median cell density is ~0.0051 cells/um^2, i.e.
about one cell per 196 um^2, so a 14 um grid puts collagen points at comparable density to the
cells they are being analysed alongside. `--spacings` runs a sensitivity check.

A regular grid also avoids the sampling noise of random placement, so results are deterministic.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

__all__ = ["load_mask", "sample_collagen_points", "mask_summary", "load_skeleton_points"]


def load_mask(path: str | Path) -> np.ndarray:
    """Labelled fibre mask as a 2-D integer array, indexed [row, col] = [y, x]."""
    return np.array(Image.open(path))


def mask_summary(mask: np.ndarray) -> dict:
    labels = np.bincount(mask.ravel().astype(int))
    sizes = labels[1:][labels[1:] > 0]
    return {
        "shape": mask.shape,
        "n_objects": int(len(sizes)),
        "coverage": float((mask > 0).mean()),
        "area_px": int(sizes.sum()) if len(sizes) else 0,
        "largest_share": float(sizes.max() / sizes.sum()) if len(sizes) else np.nan,
    }


SKELETON_DIR = Path(__file__).resolve().parent.parent / "fibre_skeletons"


def load_skeleton_points(fov: str, spacing: float = 20.0,
                         skeleton_dir: Path | None = None) -> np.ndarray | None:
    """Fibre-skeleton points for one ROI, thinned to `spacing` micrometres.

    The finer segmentation gives ~12k-98k skeleton pixels per ROI at 1 um resolution. Used
    raw they would outnumber the cells 20:1 and dominate the point cloud, so they are thinned.

    Thinning is by grid bucketing -- one point per `spacing`-sized cell -- rather than by
    arc length along a strand. That is deliberate: `strand_id` labels CONNECTED COMPONENTS,
    and a single component holds 72-92% of the skeleton, so the network is a branched graph
    and arc-length traversal is not defined on it. Grid bucketing is O(n), deterministic, and
    makes collagen point density a design parameter rather than a segmentation artefact --
    at 20 um it yields roughly as many points as there are cells.
    """
    d = (skeleton_dir or SKELETON_DIR) / f"{fov}_fibre_skeleton.csv"
    if not d.exists():
        return None
    pts = pd.read_csv(d, usecols=["x", "y"]).to_numpy(float)
    if not len(pts):
        return None
    key = np.floor(pts / spacing).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    return pts[np.sort(idx)]


def sample_collagen_points(mask: np.ndarray, spacing: float = 14.0,
                           offset: float = 0.0) -> np.ndarray:
    """Points on a `spacing`-pixel grid that fall inside the collagen mask.

    Returns an (n, 2) array of (x, y) to match the cell tables, where x is the column
    index and y the row index -- the same convention as `centroid-1` / `centroid-0`.
    1 px = 1 um for this cohort, so `spacing` is in micrometres.
    """
    rows, cols = mask.shape
    ys = np.arange(offset, rows, spacing)
    xs = np.arange(offset, cols, spacing)
    gx, gy = np.meshgrid(xs, ys)
    gx, gy = gx.ravel(), gy.ravel()
    ri = np.clip(np.round(gy).astype(int), 0, rows - 1)
    ci = np.clip(np.round(gx).astype(int), 0, cols - 1)
    keep = mask[ri, ci] > 0
    return np.column_stack([gx[keep], gy[keep]])
