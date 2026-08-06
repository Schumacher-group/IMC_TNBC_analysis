#!/usr/bin/env python
"""Fine-grained collagen fibre point-cloud extraction for topological analysis.

Context: `results_analysis/collagen/fiber_segmentation.ipynb` (ark-analysis
fiber_segmentation, Frangi filter + watershed) produces `fiber_object_table.csv`
with only ~2-180 (median 61) fused "fibre" objects per ROI -- each one an
ellipse fit to a merged blob of collagen signal, not an individual strand.
That's too coarse for network-topology TDA (persistent homology on the actual
fibre mesh: crossings, loops, branch structure).

This script instead re-segments the raw collagen channel per FOV with a
simpler, higher-resolution pipeline matching what Fig 6 describes (Otsu
threshold + edge/ridge extraction), then skeletonizes the mask to a 1px-wide
centerline and labels connected skeleton components as individual strands.
Output is one point per skeleton pixel (with a strand id), which is orders of
magnitude denser than the object table and preserves the mesh's branch/loop
topology.

Usage:
    python extract_fibre_skeleton.py \
        --img-root /mnt/data/Delta_Tissue/IMC/Img_Denoised/contrast_adj \
        --channel Collage-Type_I \
        --out-dir ./fibre_skeletons

Output: one CSV per FOV at <out-dir>/<fov>_fibre_skeleton.csv with columns
    fov, x, y, strand_id
(x = column/centroid-1, y = row/centroid-0, matching the convention used in
tda_pilot/adapter.py for cell point clouds -- 1px = 1um for Hyperion IMC).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.morphology import remove_small_objects, skeletonize


def segment_fov(img: np.ndarray, min_object_size: int = 20, gaussian_sigma: float = 0.0) -> np.ndarray:
    """Otsu-threshold + skeletonize a single-channel collagen image.

    Args:
        img: raw single-channel intensity image (e.g. Collage-Type_I.tiff)
        min_object_size: prune connected components smaller than this many
            pixels before skeletonizing, to avoid speckle noise producing
            spurious isolated skeleton points
        gaussian_sigma: optional light blur before thresholding (0 = off).
            IMC data is often already denoised upstream (Img_Denoised); only
            add blur here if the raw mask looks speckled.

    Returns:
        Boolean skeleton mask, same shape as img.
    """
    work = img.astype(float)
    if gaussian_sigma > 0:
        work = ndi.gaussian_filter(work, sigma=gaussian_sigma)

    mask = work > threshold_otsu(work)
    mask = remove_small_objects(mask, min_size=min_object_size)
    return skeletonize(mask)


def skeleton_to_points(skeleton: np.ndarray, fov: str) -> pd.DataFrame:
    """Label connected skeleton components (strands) and list their pixel coords."""
    labeled = label(skeleton, connectivity=2)
    rows, cols = np.nonzero(labeled)
    strand_id = labeled[rows, cols]
    return pd.DataFrame({
        "fov": fov,
        "x": cols.astype(float),
        "y": rows.astype(float),
        "strand_id": strand_id,
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--img-root", required=True, help="Directory containing one subfolder per FOV")
    ap.add_argument("--channel", default="Collage-Type_I", help="Channel filename (without .tiff)")
    ap.add_argument("--out-dir", required=True, help="Where to write per-FOV skeleton CSVs")
    ap.add_argument("--min-object-size", type=int, default=20, help="Small-object pruning threshold (px)")
    ap.add_argument("--gaussian-sigma", type=float, default=0.0, help="Optional pre-threshold blur radius")
    ap.add_argument("--fovs", nargs="*", default=None, help="Restrict to these FOV names (default: all)")
    ap.add_argument("--save-overlay-png", action="store_true", help="Save a QC overlay PNG per FOV")
    args = ap.parse_args()

    img_root = Path(args.img_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fov_dirs = sorted(d for d in img_root.iterdir() if d.is_dir())
    if args.fovs:
        fov_dirs = [d for d in fov_dirs if d.name in set(args.fovs)]

    summary = []
    for fov_dir in fov_dirs:
        channel_path = fov_dir / f"{args.channel}.tiff"
        if not channel_path.exists():
            print(f"skip {fov_dir.name}: no {args.channel}.tiff")
            continue

        img = tifffile.imread(channel_path)
        skeleton = segment_fov(img, min_object_size=args.min_object_size, gaussian_sigma=args.gaussian_sigma)
        points = skeleton_to_points(skeleton, fov_dir.name)
        points.to_csv(out_dir / f"{fov_dir.name}_fibre_skeleton.csv", index=False)

        n_strands = points["strand_id"].nunique()
        summary.append({"fov": fov_dir.name, "n_points": len(points), "n_strands": n_strands})
        print(f"{fov_dir.name}: {len(points)} skeleton points, {n_strands} connected strands")

        if args.save_overlay_png:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.imshow(img, cmap="gray")
            ax.imshow(np.ma.masked_where(~skeleton, skeleton), cmap="autumn", alpha=0.8)
            ax.set_title(fov_dir.name)
            ax.axis("off")
            fig.tight_layout()
            fig.savefig(out_dir / f"{fov_dir.name}_overlay.png", dpi=120)
            plt.close(fig)

    pd.DataFrame(summary).to_csv(out_dir / "_summary.csv", index=False)


if __name__ == "__main__":
    main()
