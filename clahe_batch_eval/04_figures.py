#!/usr/bin/env python
"""Stage 4: marker intensity distributions and example regions.

Two deliverables, both requested by Reviewer 2:

1. Per-batch marker intensity distributions before and after CLAHE, on a
   common axis, plus a between-batch dispersion summary for EVERY marker (not
   only the plotted ones) so the choice of markers to display cannot be
   accused of cherry-picking.

2. Example regions showing that cell type spatial structure survives
   correction: uncorrected composite, corrected composite, and the cell type
   map side by side, for ROIs drawn from different batches.

    python 04_figures.py                       # pre-treatment cohort
    python 04_figures.py --rois Leap001_8 ...  # pick the example ROIs by hand

On the example regions: because segmentation and cell type labels are SHARED
between the two tables by construction, the point these panels make is that
the corrected image preserves the tissue architecture visible in the
uncorrected one. They are not evidence that two independent annotations agree,
and the caption says so.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

IN_DIR = Path("./stage2_output")
OUT_DIR = Path("./stage4_output")

IMG_ROOT = Path("/mnt/data/Delta_Tissue/IMC/Img_Denoised/sept2024_release")
PROCESSED_DIR = IMG_ROOT / "processed"
NON_PROCESSED_DIR = IMG_ROOT / "non_processed"
MASK_DIR = Path("/mnt/data/Delta_Tissue/IMC/segmentation/sept2024_release/deepcell_output")
MASK_SUFFIX = "_whole_cell.tiff"

BATCH_KEY = "Stain_Batch"
LABEL_KEY = "Pixie"
ROI_KEY = "acquisition_ID"

# Chosen to span marker classes, as the spec asks: a lineage marker with clean
# bimodal expression, a structural marker, a functional marker, and one known
# to be dim. The dispersion table covers all 37 regardless, so these are for
# display only and the full ranking is auditable.
DISPLAY_MARKERS = ["Pan-keratin", "Collage-Type_I", "Ki-67", "PD-L2"]

# RGB composite for the example regions.
COMPOSITE_CHANNELS = {"r": "Pan-keratin", "g": "Alpha-SMA", "b": "DNA1"}

N_EXAMPLE_ROIS = 4
COMPOSITE_PERCENTILE = 99.0   # display stretch only, applied identically to both
SEED = 0


def log(msg, t0):
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


# --------------------------------------------------------------------------
# marker distributions
# --------------------------------------------------------------------------

def dispersion_table(cor, unc, batches, markers):
    """Between-batch dispersion per marker, before and after correction.

    Two summaries, because they answer slightly different questions:

    sd_of_batch_medians  -- how far apart the batches' typical intensities sit.
        Simple and directly interpretable, but blind to shape differences.
    mean_pairwise_emd    -- mean Earth Mover's distance between every pair of
        batch distributions. Catches differences in spread and shape that
        matching medians would hide.

    Lower is better for both. Reported for every marker so the display
    selection cannot bias the conclusion.
    """
    from scipy.stats import wasserstein_distance

    uniq = np.unique(batches)
    rows = []
    for j, marker in enumerate(markers):
        row = {"marker": marker}
        for name, X in [("uncorrected", unc), ("corrected", cor)]:
            per_batch = [X[batches == b, j] for b in uniq]
            medians = [np.median(v) for v in per_batch if len(v)]
            row[f"sd_of_batch_medians_{name}"] = float(np.std(medians))
            emds = [wasserstein_distance(per_batch[a], per_batch[b])
                    for a in range(len(uniq)) for b in range(a + 1, len(uniq))
                    if len(per_batch[a]) and len(per_batch[b])]
            row[f"mean_pairwise_emd_{name}"] = float(np.mean(emds)) if emds else np.nan
        row["sd_reduction"] = (row["sd_of_batch_medians_uncorrected"]
                               - row["sd_of_batch_medians_corrected"])
        row["emd_reduction"] = (row["mean_pairwise_emd_uncorrected"]
                                - row["mean_pairwise_emd_corrected"])
        rows.append(row)
    return pd.DataFrame(rows).sort_values("emd_reduction", ascending=False)


def plot_distributions(cor, unc, batches, markers, display, out_path):
    """Per-batch densities before and after, one row per marker, common axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    uniq = np.unique(batches)
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(uniq)))
    present = [m for m in display if m in markers]
    fig, axes = plt.subplots(len(present), 2, figsize=(11, 2.6 * len(present)),
                             squeeze=False)

    for i, marker in enumerate(present):
        j = markers.index(marker)
        # One x-range per marker, shared by both panels, so the eye compares
        # spread rather than axis scaling.
        both = np.concatenate([unc[:, j], cor[:, j]])
        lo, hi = np.percentile(both, [0.5, 99.5])
        grid = np.linspace(lo, hi, 200)

        for col, (name, X) in enumerate([("uncorrected (denoised only)", unc),
                                         ("corrected (CLAHE)", cor)]):
            ax = axes[i][col]
            for c, b in zip(colours, uniq):
                v = X[batches == b, j]
                if len(v) < 10 or np.allclose(v, v[0]):
                    continue
                ax.plot(grid, gaussian_kde(v)(grid), color=c, lw=1.4,
                        label=f"batch {b}" if i == 0 and col == 0 else None)
            ax.set_xlim(lo, hi)
            ax.set_yticks([])
            if i == 0:
                ax.set_title(name, fontweight="bold")
            if col == 0:
                ax.set_ylabel(marker, fontweight="bold")
        axes[i][1].sharey(axes[i][0])

    axes[0][0].legend(fontsize=8, frameon=False)
    fig.suptitle("Marker intensity distributions by staining batch",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------
# example regions
# --------------------------------------------------------------------------

def composite(fov, img_dir, tifffile, shape=None):
    """RGB composite from three channels, percentile-stretched per channel."""
    planes = []
    for key in ("r", "g", "b"):
        path = img_dir / fov / f"{COMPOSITE_CHANNELS[key]}.tiff"
        img = tifffile.imread(path).astype(np.float32)
        top = np.percentile(img, COMPOSITE_PERCENTILE)
        planes.append(np.clip(img / top, 0, 1) if top > 0 else np.zeros_like(img))
    return np.dstack(planes)


def plot_examples(rois, obs, tifffile, out_path):
    """Uncorrected | corrected | cell type map, one row per ROI."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    types = sorted(obs[LABEL_KEY].astype(str).unique())
    palette = {t: plt.cm.tab20(i % 20) for i, t in enumerate(types)}

    fig, axes = plt.subplots(len(rois), 3, figsize=(13, 4.3 * len(rois)),
                             squeeze=False)
    for i, fov in enumerate(rois):
        sub = obs[obs[ROI_KEY] == fov]
        batch = sub[BATCH_KEY].iloc[0]

        for col, (img_dir, title) in enumerate(
                [(NON_PROCESSED_DIR, "uncorrected (denoised only)"),
                 (PROCESSED_DIR, "corrected (CLAHE)")]):
            ax = axes[i][col]
            ax.imshow(composite(fov, img_dir, tifffile))
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(title, fontweight="bold")
            if col == 0:
                ax.set_ylabel(f"{fov}\nbatch {batch}", fontweight="bold", fontsize=9)

        ax = axes[i][2]
        ax.scatter(sub["centroid-1"], sub["centroid-0"], s=1.5,
                   c=[palette[t] for t in sub[LABEL_KEY].astype(str)], linewidths=0)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        if i == 0:
            ax.set_title("cell types (shared by both)", fontweight="bold")

    handles = [plt.Line2D([], [], marker="o", ls="", color=palette[t], label=t)
               for t in types]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Tissue architecture before and after CLAHE\n"
                 "segmentation and cell types are shared by construction; the point "
                 "is that correction preserves the structure visible in the raw image",
                 fontweight="bold", fontsize=10)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def pick_rois(obs, n, rng):
    """One ROI from each of n different batches, preferring cell-rich ROIs.

    Different batches matter here: the figure is meant to show that correction
    preserves architecture across the batch range, not in one favourable field.
    """
    # observed=True matters: with a categorical grouper pandas otherwise emits a
    # zero-count row for every batch x ROI combination, and nlargest then pads
    # the candidate pool with ROIs that belong to a different batch entirely --
    # which produced an example figure labelled with the wrong batch.
    counts = (obs.groupby([BATCH_KEY, ROI_KEY], observed=True)
              .size().rename("n").reset_index())
    counts = counts[counts["n"] > 0]

    chosen = []
    for batch in sorted(counts[BATCH_KEY].unique())[:n]:
        pool = counts[counts[BATCH_KEY] == batch].nlargest(10, "n")
        pick = pool.iloc[rng.integers(len(pool))]
        assert pick[BATCH_KEY] == batch, "ROI selected from the wrong batch"
        chosen.append(pick[ROI_KEY])
    return chosen


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", type=Path, default=IN_DIR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--cohort", choices=["pre", "revision", "all"], default="pre")
    ap.add_argument("--markers", nargs="*", default=DISPLAY_MARKERS)
    ap.add_argument("--rois", nargs="*", default=None,
                    help="Example ROIs; default picks one per batch")
    ap.add_argument("--n-rois", type=int, default=N_EXAMPLE_ROIS)
    ap.add_argument("--skip-examples", action="store_true",
                    help="Distributions only; no image access needed")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()

    import anndata as ad
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "reconcile", Path(__file__).resolve().parent / "01_reconcile.py")
    reconcile = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reconcile)

    t0 = time.time()
    tag = args.tag or args.cohort
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rep = reconcile.Report()
    rep(f"# CLAHE batch-correction evaluation -- stage 4 figures ({tag})\n")

    cor = ad.read_h5ad(args.in_dir / "adata_corrected.h5ad")
    unc = ad.read_h5ad(args.in_dir / "adata_uncorrected.h5ad")

    keep = np.ones(cor.n_obs, dtype=bool)
    if args.cohort in ("pre", "revision"):
        keep &= cor.obs["in_revision_cohort"].to_numpy()
    if args.cohort == "pre":
        keep &= cor.obs["is_pre_treatment"].to_numpy()

    obs = cor.obs[keep].copy()
    obs["centroid-0"] = cor.obsm["spatial"][keep, 0]
    obs["centroid-1"] = cor.obsm["spatial"][keep, 1]
    markers = list(cor.var_names)
    batches = obs[BATCH_KEY].to_numpy()

    # Distributions use the transformed values -- the same ones the stage 3
    # metrics saw -- so the figure and the numbers describe the same data.
    Xc, Xu = cor.X[keep], unc.X[keep]

    rep.h("Cohort")
    rep(f"- `{args.cohort}`: {int(keep.sum()):,} cells, "
        f"{obs[ROI_KEY].nunique()} ROIs, {len(np.unique(batches))} batches")
    rep("- distributions are of the transformed intensities (log1p, per-channel "
        "0.95-quantile cap), identical to what stage 3 measured")

    rep.h("Between-batch dispersion, all markers")
    disp = dispersion_table(Xc, Xu, batches, markers)
    disp.to_csv(out_dir / f"marker_dispersion_{tag}.csv", index=False)
    log("dispersion table done", t0)

    improved_emd = int((disp["emd_reduction"] > 0).sum())
    improved_sd = int((disp["sd_reduction"] > 0).sum())
    rep(f"- markers whose between-batch EMD FELL after CLAHE: "
        f"{improved_emd} of {len(disp)}")
    rep(f"- markers whose SD of per-batch medians fell: {improved_sd} of {len(disp)}")
    rep("\nLower dispersion is better. Computed for every marker, so the choice of "
        "which to plot cannot bias the conclusion.\n")
    rep("Five most improved:\n")
    rep.table(disp.head(5).round(4))
    rep("\nFive least improved (negative = CLAHE made between-batch agreement worse):\n")
    rep.table(disp.tail(5).round(4))
    if improved_emd < len(disp):
        rep(f"\n**{len(disp) - improved_emd} marker(s) got worse, not better.** "
            "Listed above; this belongs in the rebuttal rather than being dropped.")

    plot_distributions(Xc, Xu, batches, markers, args.markers,
                       out_dir / f"marker_distributions_{tag}.png")
    rep(f"\nPlotted: {[m for m in args.markers if m in markers]}")
    log("distribution figure done", t0)

    if not args.skip_examples:
        import tifffile

        rois = args.rois or pick_rois(obs, args.n_rois, rng)
        rep.h("Example regions")
        for fov in rois:
            sub = obs[obs[ROI_KEY] == fov]
            rep(f"- `{fov}` -- batch {sub[BATCH_KEY].iloc[0]}, {len(sub):,} cells")
        missing = [f for f in rois if not (NON_PROCESSED_DIR / f).exists()]
        if missing:
            rep.fail(f"ROIs with no image directory: {missing}")
        else:
            plot_examples(rois, obs, tifffile,
                          out_dir / f"example_regions_{tag}.png")
            rep(f"\nComposite channels: R={COMPOSITE_CHANNELS['r']}, "
                f"G={COMPOSITE_CHANNELS['g']}, B={COMPOSITE_CHANNELS['b']}, each "
                f"stretched to its {COMPOSITE_PERCENTILE}th percentile. The stretch "
                f"is a display choice applied identically to both columns; it does "
                f"not equalise them, since CLAHE is local and a percentile stretch "
                f"is global.")
            log("example figure done", t0)

    rep.write(out_dir / f"figures_report_{tag}.md")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
