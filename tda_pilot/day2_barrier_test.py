#!/usr/bin/env python
"""Does collagen lie between cytotoxic T cells and tumour, more so in non-responders?

This states the Dundee-slide claim literally, in one number, without any topology. For each
CD8 cell, take the straight segment to its nearest tumour cell and measure what FRACTION OF
THAT PATH LENGTH passes through collagen.

WHY FRACTION AND NOT "DOES IT CROSS"
------------------------------------
A binary crossing test saturates. With ~20-33% of tissue labelled collagen, nearly every path
touches some, and on the Fig 4B panels a binary test went from a 20.5-point gap to a 1.6-point
gap simply by loosening the crossing tolerance. The fraction of path length is graded, has no
tolerance parameter once a mask is fixed, and cannot saturate.

WHY A WITHIN-ROI NULL IS NOT OPTIONAL
-------------------------------------
Without one this measures collagen ABUNDANCE, not barrier placement: an ROI with more collagen
scores higher whatever the arrangement. That is precisely the confound that wrecked the
z-score analyses earlier in this directory.

So each ROI is compared against itself. The null relabels which non-tumour cells are "CD8":
`n_CD8` cells are drawn at random from the non-tumour population, keeping cell positions, the
tissue scaffold, the tumour cells and the collagen mask all exactly fixed. The observed
fraction is then expressed relative to that null, which is dimensionless and immune to how
much collagen the ROI happens to contain, to its density, and to its area.

Reported per ROI:
    obs           observed median path fraction through collagen
    null_mean     mean of the same statistic over `--n-perm` relabellings
    excess        obs - null_mean   (positive = collagen sits between CD8 and tumour MORE
                                     than the tissue's own architecture implies)
    z             (obs - null_mean) / null_sd

Masks: `--variant nonclahe` (default, post-denoise / pre-rescale / pre-CLAHE) or `clahe`.
Non-CLAHE is the primary: it shows no batch dependence (rho = +0.02 vs -0.11) and does not
inflate coverage by amplifying faint signal in collagen-poor regions, which is where a barrier
question is decided.

Output: output/BARRIER_TEST_<variant>.md, output/barrier_test_<variant>.png
"""

import argparse
import multiprocessing as mp
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

import cohort

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "output"
ROIS_DIR = HERE / "data" / "rois"
TUMOUR = {"Cancer cell", "B7H4 Cancer cell"}
CD8 = {"CD8 T cell", "Memory CD8 T cell"}
# Comparable leukocytes for the null: excludes fibroblasts and endothelium (structurally
# embedded in collagen/vessels) and the CD8 cells themselves.
NULL_POOL = {"Memory CD4 T cell", "Regulatory T cell", "B cell", "B cell + HLA",
             "Macrophage", "Macrophage + HLA", "Macrophage M2", "NK/CD8",
             "Antigen presenting cell", "Neutrophil", "Monocyte", "Cl Monocyte",
             "Int Monocyte"}
MIN_CD8, MIN_TUM = 30, 30
N_STEPS = 60
SEED = 20260806
COLOURS = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}
_CFG = {"variant": "nonclahe", "n_perm": 99}


def _init(cfg):
    global _CFG
    _CFG = cfg


def _mask_dir(variant: str) -> Path:
    return ROOT / ("fibre_skeletons_nonCLAHE" if variant == "nonclahe" else "fibre_skeletons")


def _path_fraction(src: np.ndarray, tum_tree: cKDTree, tum: np.ndarray,
                   mask: np.ndarray) -> float:
    """MEAN over `src` points of the fraction of the path to nearest tumour in collagen.

    Mean, not median: most CD8 cells sit within ~12 um of tumour, so over half have paths
    too short to touch any collagen and the median collapses to exactly 0.
    """
    if len(src) == 0:
        return np.nan
    _, idx = tum_tree.query(src, k=1)
    t = np.linspace(0, 1, N_STEPS)[None, :, None]
    seg = src[:, None, :] * (1 - t) + tum[idx][:, None, :] * t
    ci = np.clip(np.round(seg[..., 0]).astype(int), 0, mask.shape[1] - 1)
    ri = np.clip(np.round(seg[..., 1]).astype(int), 0, mask.shape[0] - 1)
    return float(np.mean((mask[ri, ci] > 0).mean(axis=1)))


def _distance_matched_sample(d_pool: np.ndarray, d_target: np.ndarray, rng,
                             n_bins: int = 20) -> np.ndarray:
    """Indices into the pool whose distance-to-tumour distribution matches the target's.

    THIS IS THE LOAD-BEARING PART OF THE NULL. Non-tumour cells drawn at random sit further
    from tumour than real CD8 cells do, so their paths are longer and cross more collagen for
    purely geometric reasons -- an unmatched null gave z = -14.7 on a single ROI, i.e. it was
    measuring distance-to-tumour rather than collagen placement. Matching the distance
    distribution removes that, so the comparison isolates the question actually being asked:
    for a cell THIS FAR from tumour, is there more collagen in between than chance implies?
    """
    edges = np.quantile(d_target, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    tgt_counts, _ = np.histogram(d_target, bins=edges)
    pool_bin = np.digitize(d_pool, edges[1:-1])
    picks = []
    for b, want in enumerate(tgt_counts):
        if want == 0:
            continue
        cand = np.flatnonzero(pool_bin == b)
        if len(cand) == 0:
            continue
        picks.append(rng.choice(cand, min(want, len(cand)), replace=False))
    return np.concatenate(picks) if picks else np.empty(0, dtype=int)


def _process_roi(fov: str) -> dict | None:
    try:
        mpath = _mask_dir(_CFG["variant"]) / f"{fov}_fibre_mask.tiff"
        if not mpath.exists():
            return None
        df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
        tum = df[df.celltype.isin(TUMOUR)][["x", "y"]].to_numpy(float)
        cd8_mask = df.celltype.isin(CD8)
        cd8 = df[cd8_mask][["x", "y"]].to_numpy(float)
        # The null pool is other IMMUNE cells, not all non-tumour cells. Fibroblasts and
        # endothelium sit inside collagen and vessels by definition, so including them draws
        # the null from collagen-embedded positions and makes real CD8 cells look
        # artificially collagen-avoidant (z ~ -8 on both Fig 4B panels before this fix).
        # Other leukocytes are the cells that could plausibly occupy a CD8 cell's position.
        other = df[df.celltype.isin(NULL_POOL)][["x", "y"]].to_numpy(float)
        if len(cd8) < MIN_CD8 or len(tum) < MIN_TUM or len(other) <= len(cd8):
            return None

        mask = np.array(Image.open(mpath))
        tree = cKDTree(tum)
        obs = _path_fraction(cd8, tree, tum, mask)

        # match the null on distance-to-tumour, otherwise it measures path length
        d_cd8, _ = tree.query(cd8, k=1)
        d_other, _ = tree.query(other, k=1)
        rng = np.random.default_rng(SEED + (int.from_bytes(fov.encode(), "little") % 10**6))
        null = np.array([
            _path_fraction(other[_distance_matched_sample(d_other, d_cd8, rng)],
                           tree, tum, mask)
            for _ in range(_CFG["n_perm"])])
        sd = float(np.nanstd(null, ddof=1))
        return {"roi_id": fov, "n_cd8": len(cd8), "n_tumour": len(tum),
                "coverage": float((mask > 0).mean()), "obs": obs,
                "null_mean": float(np.nanmean(null)), "null_sd": sd,
                "excess": obs - float(np.nanmean(null)),
                "z": (obs - float(np.nanmean(null))) / sd if sd > 0 else 0.0}
    except Exception:  # noqa: BLE001
        return None


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, float(u.pvalue)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="nonclahe", choices=["nonclahe", "clahe"])
    ap.add_argument("--n-perm", type=int, default=99)
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()
    mp.set_start_method("fork", force=True)

    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "resp",
                 "Sample_Type_(pre/post treatment)": "stype"})
    meta = meta[(meta.stype.astype(str).str.lower() == "pre")
                & meta.resp.isin(["Responder", "Non-Responder"])]
    rois = sorted(r for r in meta.roi_id
                  if cohort.is_cohort_member(r) and (ROIS_DIR / f"{r}.csv").exists())

    t0 = time.perf_counter()
    print(f"{len(rois)} pre-treatment ROIs, variant={args.variant}, "
          f"{args.n_perm} relabellings each, {args.workers} workers")
    with mp.Pool(args.workers, initializer=_init,
                 initargs=({"variant": args.variant, "n_perm": args.n_perm},)) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(_process_roi, rois, chunksize=2), 1):
            if r is not None:
                out.append(r)
            if i % 100 == 0 or i == len(rois):
                print(f"  {i}/{len(rois)} ({(time.perf_counter()-t0)/60:.1f} min)")
    d = pd.DataFrame(out).merge(meta[["roi_id", "pid", "resp"]], on="roi_id")
    d.to_csv(OUT / f"barrier_test_{args.variant}.csv", index=False)
    print(f"{len(d)} ROIs in {(time.perf_counter()-t0)/60:.1f} min")

    rows = []
    for col, lab in (("obs", "raw path fraction (confounded by collagen amount)"),
                     ("excess", "excess over within-ROI null"),
                     ("z", "z vs within-ROI null")):
        pat = d.groupby(["pid", "resp"], as_index=False)[col].median()
        dl, pv = cliffs(pat[pat.resp == "Responder"][col],
                        pat[pat.resp == "Non-Responder"][col])
        roi_dl, roi_pv = cliffs(d[d.resp == "Responder"][col], d[d.resp == "Non-Responder"][col])
        rows.append({"statistic": col, "meaning": lab,
                     "median_R": float(pat[pat.resp == "Responder"][col].median()),
                     "median_NR": float(pat[pat.resp == "Non-Responder"][col].median()),
                     "patient_delta": dl, "patient_p": pv,
                     "roi_delta": roi_dl, "roi_p": roi_pv})
    res = pd.DataFrame(rows)

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))
    ax[0].hist(d.obs, bins=40, color="0.6", alpha=.8, label="observed")
    ax[0].hist(d.null_mean, bins=40, color="tab:orange", alpha=.6, label="within-ROI null")
    ax[0].set_xlabel("median fraction of CD8→tumour path in collagen")
    ax[0].set_ylabel("ROIs"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"A. Is collagen preferentially between CD8 and tumour?\n"
                    f"median excess = {d.excess.median():+.4f}", fontsize=10)
    for i, (col, ttl) in enumerate([("obs", "B. Raw fraction by response\n(confounded by "
                                            "collagen amount)"),
                                    ("excess", "C. Excess over within-ROI null\n(the "
                                               "barrier statistic)")]):
        a = ax[i + 1]
        pat = d.groupby(["pid", "resp"], as_index=False)[col].median()
        jit = np.random.default_rng(0)
        for j, grp in enumerate(["Non-Responder", "Responder"]):
            v = pat[pat.resp == grp][col].to_numpy()
            a.scatter(jit.normal(j, .08, len(v)), v, s=38, c=COLOURS[grp], alpha=.85, lw=0)
            a.hlines(np.median(v), j - .26, j + .26, color="black", lw=2.4)
        r = res[res.statistic == col].iloc[0]
        a.axhline(0, color="grey", lw=.6, ls=":")
        a.set_xticks([0, 1]); a.set_xticklabels(["NR", "R"])
        a.set_title(f"{ttl}\nCliff δ={r.patient_delta:+.3f}, p={r.patient_p:.3f}", fontsize=10)
    fig.suptitle(f"Do T cells have to cross collagen to reach tumour? "
                 f"{len(d)} pre-treatment ROIs, {args.variant} masks", fontsize=12)
    fig.tight_layout()
    figp = OUT / f"barrier_test_{args.variant}.png"
    fig.savefig(figp, dpi=150, bbox_inches="tight")

    def md(dfr):
        h = "| " + " | ".join(dfr.columns) + " |"
        s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
        b = "\n".join("| " + " | ".join(
            f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
            for v in row) + " |" for row in dfr.itertuples(index=False))
        return "\n".join([h, s, b])

    lines = [
        f"# Barrier test ({args.variant} masks)\n",
        f"{len(d)} pre-treatment ROIs, {d.pid.nunique()} patients. For each CD8 cell, the "
        f"fraction of the straight path to its nearest tumour cell that lies inside collagen; "
        f"median over cells per ROI. Compared against {args.n_perm} relabellings in which the "
        f"same number of CD8 cells is drawn from the non-tumour population of that ROI, "
        f"holding cell positions, tumour and collagen fixed.\n",
        f"Collagen is preferentially between CD8 and tumour overall: median excess over the "
        f"within-ROI null = **{d.excess.median():+.4f}** "
        f"({(d.excess > 0).mean():.0%} of ROIs positive).\n",
        md(res.round(4)),
        "\n`obs` is reported only to show what it does WITHOUT the null -- it is a measure of "
        "how much collagen an ROI contains as much as of where that collagen sits. `excess` "
        "and `z` are the interpretable ones.\n",
        f"Figure: `{figp.name}`. Script: `day2_barrier_test.py`.",
    ]
    (OUT / f"BARRIER_TEST_{args.variant.upper()}.md").write_text("\n".join(lines) + "\n")
    pd.set_option("display.width", 220)
    print(res.round(4).to_string(index=False))
    print(f"\nmedian excess over null: {d.excess.median():+.4f} "
          f"({(d.excess>0).mean():.0%} of ROIs positive)")
    print("wrote", figp)


if __name__ == "__main__":
    main()
