#!/usr/bin/env python
"""Is CD8 more excluded from B7H4+ tumour than from B7H4- tumour? A paired, matched design.

The naive approach -- run (Cancer, CD8) and (B7H4-Cancer, CD8) as two separate pair
definitions and compare the statistics -- does not work, for a reason established earlier in
this directory: bar LENGTH tracks the minority fraction at rho = -0.69. Cancer cells outnumber
B7H4+ cancer cells (median split 69% / 31%), so the CD8 fraction differs systematically
between the two comparisons and the "more exclusion" signal would be manufactured by
composition alone.

Two design choices fix that, and together they make this a much better-powered question than
the responder contrast:

  MATCHED   Within each ROI, subsample both cancer subtypes to the SAME number of cells,
            m = min(n_Cancer, n_B7H4), and use the identical CD8 set for both. The two
            six-packs then have identical species counts by construction, so composition
            cannot differ. Averaged over `--repeats` independent subsamples to remove
            sampling noise. This needs no permutation null -- the matching does the job
            directly, and far more cheaply.

  PAIRED    Both statistics come from the same ROI: same tissue, same area, same density,
            same CD8 cells. Every between-ROI confound that dogged the response analysis
            (geometry, sampling, patient) cancels in the within-ROI difference.

The comparison is also well powered: n = 557 ROIs, not 62 patients. Note this question is
independent of response -- it asks about tumour biology, and can be reported whether or not it
differs between responders and non-responders.

CAVEAT, established by `day2_b7h4_precheck.py`: the two subtypes are only partly segregated
(median mixing ratio 0.85 vs 1.0 for random intermingling). Partial intermingling attenuates
any true difference, so this design is conservative -- a positive result is trustworthy, a
null is partly explained by the overlap.

Output: output/B7H4_PAIRED.md, output/b7h4_paired.png
"""

import argparse
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import cohort

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))

from chalc import chromatic  # noqa: E402
from chalc.sixpack import KChromaticInclusion  # noqa: E402
from utils.generate_stats import canonicalize_labels, get_stats_from_barcodes_dict  # noqa: E402

OUT = HERE / "output"
ROIS_DIR = HERE / "data" / "rois"
CANCER, B7H4 = "Cancer cell", "B7H4 Cancer cell"
CD8 = {"CD8 T cell", "Memory CD8 T cell"}
FILTRATION = "delaunay_cech"
MIN_EACH, MIN_CD8 = 30, 30
SEED = 20260806

# Exclusion-relevant degree-1 and degree-0 summaries; kernel bars are the exclusion signature
# (a tumour nest ringed by CD8 gives fewer but longer kernel bars).
KEEP = tuple(f"{dg}{dim}-{st}"
             for dg in ("ker", "im", "cok")
             for dim in (0, 1)
             for st in ("avg_length", "med_length", "p90_length", "num_bars", "entropy"))
_REPEATS = 5


def _init(repeats: int) -> None:
    global _REPEATS
    _REPEATS = repeats


def _sixpack_stats(pts: np.ndarray, is_tumour: np.ndarray) -> dict:
    """Two-colour six-pack, same code path as the pair runs (KChromaticInclusion k=1)."""
    colours = canonicalize_labels(is_tumour.astype(int))
    filt = getattr(chromatic, FILTRATION)(pts.transpose(), colours)
    dgms = KChromaticInclusion(filt, 1).sixpack()
    barcodes = {n: dgms.get_matrix(n, [0, 1]) for n in dgms}
    allstats = get_stats_from_barcodes_dict(barcodes)
    return {k: float(allstats[k]) for k in KEEP if k in allstats}


def _process_roi(fov: str) -> dict | None:
    try:
        df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
        canc = df[df.celltype == CANCER][["x", "y"]].to_numpy(float)
        b7h4 = df[df.celltype == B7H4][["x", "y"]].to_numpy(float)
        cd8 = df[df.celltype.isin(CD8)][["x", "y"]].to_numpy(float)
        m = min(len(canc), len(b7h4))
        if m < MIN_EACH or len(cd8) < MIN_CD8:
            return None

        rng = np.random.default_rng(SEED + (int.from_bytes(fov.encode(), "little") % 10**6))
        acc: dict[str, list[float]] = {}
        for _ in range(_REPEATS):
            for tag, arr in (("cancer", canc), ("b7h4", b7h4)):
                sel = arr[rng.choice(len(arr), m, replace=False)]
                pts = np.vstack([sel, cd8])
                is_t = np.concatenate([np.ones(m, bool), np.zeros(len(cd8), bool)])
                for k, v in _sixpack_stats(pts, is_t).items():
                    acc.setdefault(f"{tag}__{k}", []).append(v)
        row = {"roi_id": fov, "n_matched": m, "n_cd8": len(cd8),
               "n_cancer_total": len(canc), "n_b7h4_total": len(b7h4)}
        row |= {k: float(np.mean(v)) for k, v in acc.items()}
        return row
    except Exception:  # noqa: BLE001
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()
    mp.set_start_method("fork", force=True)

    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "response",
                 "Sample_Type_(pre/post treatment)": "sample_type"})
    meta = meta[(meta.sample_type.astype(str).str.lower() == "pre")
                & meta.response.isin(["Responder", "Non-Responder"])]
    rois = sorted(r for r in meta.roi_id
                  if cohort.is_cohort_member(r) and (ROIS_DIR / f"{r}.csv").exists())

    t0 = time.perf_counter()
    print(f"{len(rois)} candidate ROIs, {args.repeats} matched subsamples each, "
          f"{args.workers} workers")
    with mp.Pool(args.workers, initializer=_init, initargs=(args.repeats,)) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(_process_roi, rois, chunksize=2), 1):
            if r is not None:
                out.append(r)
            if i % 100 == 0 or i == len(rois):
                print(f"  {i}/{len(rois)}  ({(time.perf_counter()-t0)/60:.1f} min)")
    d = pd.DataFrame(out).merge(meta[["roi_id", "pid", "response"]], on="roi_id")
    d.to_csv(OUT / "b7h4_paired.csv", index=False)
    print(f"\n{len(d)} ROIs analysed in {(time.perf_counter()-t0)/60:.1f} min")

    rows = []
    for stat in KEEP:
        a, b = f"b7h4__{stat}", f"cancer__{stat}"
        if a not in d.columns or b not in d.columns:
            continue
        diff = (d[a] - d[b]).to_numpy(float)
        diff = diff[np.isfinite(diff)]
        if len(diff) < 20 or np.allclose(diff, 0):
            continue
        stat_w, p = wilcoxon(diff)
        # rank-biserial correlation: signed effect size for a paired rank test
        pos = float((diff > 0).mean())
        rows.append({"statistic": stat, "median_cancer": float(d[b].median()),
                     "median_b7h4": float(d[a].median()), "median_diff": float(np.median(diff)),
                     "frac_b7h4_higher": pos, "wilcoxon_p": float(p),
                     "rank_biserial": 2 * pos - 1})
    res = pd.DataFrame(rows)
    # Benjamini-Hochberg across the family
    p = res.wilcoxon_p.to_numpy()
    order = np.argsort(p)
    adj = np.minimum.accumulate((p[order] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    q = np.empty(len(p)); q[order] = np.minimum(adj, 1.0)
    res["q"] = q
    res = res.sort_values("wilcoxon_p").reset_index(drop=True)

    key = "ker1-avg_length"
    kr = res[res.statistic == key].iloc[0] if (res.statistic == key).any() else None

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    if kr is not None:
        diff = (d[f"b7h4__{key}"] - d[f"cancer__{key}"]).to_numpy(float)
        ax[0].hist(diff, bins=45, color="tab:purple", alpha=0.8)
        ax[0].axvline(0, color="black", lw=1.5, ls="--")
        ax[0].axvline(np.median(diff), color="crimson", lw=2,
                      label=f"median Δ = {np.median(diff):+.3f}")
        ax[0].set_xlabel(f"Δ {key}  (B7H4+ − B7H4−), µm")
        ax[0].set_ylabel("ROIs")
        ax[0].set_title(f"A. Paired difference, matched cell counts\n"
                        f"{kr.frac_b7h4_higher:.0%} of ROIs higher for B7H4+, "
                        f"Wilcoxon p={kr.wilcoxon_p:.2g}", fontsize=10)
        ax[0].legend(fontsize=8)

        lim = np.percentile(np.r_[d[f"cancer__{key}"], d[f"b7h4__{key}"]], [1, 99])
        ax[1].scatter(d[f"cancer__{key}"], d[f"b7h4__{key}"], s=11, c="tab:purple",
                      alpha=0.45, linewidths=0)
        ax[1].plot(lim, lim, "k--", lw=1.2, label="equality")
        ax[1].set_xlim(*lim); ax[1].set_ylim(*lim)
        ax[1].set_xlabel(f"{key}, B7H4− cancer vs CD8")
        ax[1].set_ylabel(f"{key}, B7H4+ cancer vs CD8")
        ax[1].set_title("B. Per ROI, same tissue, matched counts\n(above the line = more "
                        "exclusion from B7H4+)", fontsize=10)
        ax[1].legend(fontsize=8)

    top = res.head(10).iloc[::-1]
    colours = ["tab:purple" if q < 0.05 else "0.7" for q in top.q]
    ax[2].barh(top.statistic, top.rank_biserial, color=colours)
    ax[2].axvline(0, color="black", lw=1)
    ax[2].set_xlabel("rank-biserial effect size  (+ve = higher for B7H4+)")
    ax[2].set_title("C. Strongest paired differences\n(purple = BH q < 0.05)", fontsize=10)

    fig.suptitle("Is CD8 more excluded from B7H4+ than B7H4− tumour? "
                 f"Paired within-ROI, cell counts matched, {len(d)} pre-treatment ROIs",
                 fontsize=12)
    fig.tight_layout()
    figp = OUT / "b7h4_paired.png"
    fig.savefig(figp, dpi=150, bbox_inches="tight")

    def md(dfr):
        h = "| " + " | ".join(dfr.columns) + " |"
        s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
        b = "\n".join("| " + " | ".join(
            f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
            for v in row) + " |" for row in dfr.itertuples(index=False))
        return "\n".join([h, s, b])

    n_sig = int((res.q < 0.05).sum())
    lines = [
        "# Is CD8 more excluded from B7H4+ tumour than from B7H4− tumour?\n",
        f"Paired within-ROI comparison, {len(d)} pre-treatment ROIs. Within each ROI both "
        f"cancer subtypes are subsampled to the same count "
        f"(median m = {int(d.n_matched.median())}) and compared against the identical CD8 set "
        f"(median {int(d.n_cd8.median())} cells), averaged over {args.repeats} subsamples. "
        f"Species counts are therefore identical between the two arms by construction, so "
        f"composition cannot drive the difference; and both arms share the same tissue, so "
        f"ROI geometry cancels.\n",
        f"**{n_sig} of {len(res)} statistics differ at BH q < 0.05.**\n",
        md(res.round(4)),
        "\n## Reading this\n",
        "`ker1-avg_length` is the exclusion signature: a tumour nest ringed by CD8 produces "
        "fewer but longer degree-1 kernel bars. Positive Δ means B7H4+ tumour shows more of "
        "that pattern than B7H4− tumour in the same tissue.\n",
        "Conservative by design: `B7H4_PRECHECK.md` shows the two subtypes are only partly "
        "segregated (median mixing ratio 0.85 against 1.0 for random intermingling), which "
        "attenuates any true difference. A positive result here is therefore trustworthy; a "
        "null is partly explained by the overlap.\n",
        "This question is independent of treatment response — it concerns tumour biology and "
        "can be reported whether or not it differs between responders and non-responders.\n",
        f"Figure: `{figp.name}`. Script: `day2_b7h4_paired.py`.",
    ]
    (OUT / "B7H4_PAIRED.md").write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 220)
    print(res.head(12).round(4).to_string(index=False))
    print(f"\nBH q<0.05: {n_sig}/{len(res)}")
    print("wrote", figp, "and", OUT / "B7H4_PAIRED.md")


if __name__ == "__main__":
    main()
