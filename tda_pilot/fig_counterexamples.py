#!/usr/bin/env python
"""Counterexample gallery: the architecture spectrum, with both response groups at every level.

The honest visual accompaniment to a null result. Fig 4B shows one intermixed responder and
one excluded non-responder, which invites the reading that architecture tracks response. This
figure shows the same axis across the whole cohort and finds responders AND non-responders at
every point along it -- which is what "no consistent difference" looks like, as opposed to
merely asserting it.

Axis: an exclusion score built from the CD8-directional Dowker-style signature measured against
each ROI's own label-permutation null. Immune exclusion means CD8 ringing compact tumour nests,
which produces FEWER but LONGER degree-1 kernel bars, so

    exclusion = (relative deviation of ker1-avg_length − relative deviation of ker1-num_bars) / sqrt(2)

Relative deviation, (observed − null mean) / null mean, rather than a z-score: it keeps the
within-ROI pairing but avoids dividing by the null SD, which shrinks with ROI size and density
and so makes z incomparable across ROIs (rho(z, density) = +0.31 against +0.13 for relative
deviation). See REVIEW.md UPDATE 2/4 for how that was learned the hard way.

ROIs are drawn from a mid-size band so the panels are comparable, and matched in pairs on the
exclusion score so each column contrasts a responder and a non-responder with near-identical
architecture.

Output: output/fig_counterexamples.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pair_defs

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
ROIS_DIR = HERE / "data" / "rois"
CD8 = pair_defs.DEFINITIONS["CD8_primary"]
COLOURS = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}
PANELS = {"Leap008_5": "Fig 4B R", "Leap046_1": "Fig 4B NR"}


def load(defn: str = "CD8_primary") -> pd.DataFrame:
    d = pd.read_parquet(OUT / f"perm_null_{defn}_other_only.parquet")
    d = d[d["status"] == "ok"].copy()
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "response"},
    )[["roi_id", "pid", "response"]]
    d = d.merge(meta, on="roi_id").query("response in @COLOURS")
    for s in ("ker1-avg_length", "ker1-num_bars"):
        m = d[f"{s}__null_mean"]
        d[f"rd_{s}"] = np.where(m != 0, (d[f"{s}__obs"] - m) / m, np.nan)
    d["exclusion"] = (d["rd_ker1-avg_length"] - d["rd_ker1-num_bars"]) / np.sqrt(2)
    d["n_pair"] = d.n_tumour + d.n_other
    return d.dropna(subset=["exclusion"])


def draw(ax, fov: str, title: str, colour: str) -> None:
    df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
    t = df[df.celltype.isin(pair_defs.TUMOUR)]
    s = df[df.celltype.isin(CD8)]
    ax.scatter(t.x, t.y, s=1.4, c="0.75", alpha=0.85, linewidths=0)
    ax.scatter(s.x, s.y, s=1.6, c=colour, alpha=0.9, linewidths=0)
    ax.set_title(title, fontsize=8, color=colour)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(colour)
        sp.set_linewidth(1.6)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cols", type=int, default=5)
    args = ap.parse_args()
    d = load()

    lo, hi = d.n_pair.quantile([0.25, 0.75])
    band = d[(d.n_pair >= lo) & (d.n_pair <= hi)].copy()
    R = band[band.response == "Responder"].sort_values("exclusion")
    N = band[band.response == "Non-Responder"].sort_values("exclusion")

    # pick responder ROIs spanning the score range, then match each with the closest
    # non-responder -- so every column is a like-for-like architectural comparison
    targets = np.quantile(band.exclusion, np.linspace(0.05, 0.95, args.n_cols))
    pairs, used = [], set()
    for t in targets:
        r = R.iloc[(R.exclusion - t).abs().argsort()]
        r = r[~r.roi_id.isin(used)].iloc[0]
        n = N.iloc[(N.exclusion - r.exclusion).abs().argsort()]
        n = n[~n.roi_id.isin(used)].iloc[0]
        used |= {r.roi_id, n.roi_id}
        pairs.append((r, n))

    fig, axes = plt.subplots(2, args.n_cols, figsize=(2.7 * args.n_cols, 6.0))
    for j, (r, n) in enumerate(pairs):
        for i, row in enumerate((r, n)):
            pct = 100 * (d.exclusion < row.exclusion).mean()
            draw(axes[i, j], row.roi_id,
                 f"{row.roi_id}\nexclusion {row.exclusion:+.2f}  ({pct:.0f}th pct)",
                 COLOURS[row.response])
    axes[0, 0].set_ylabel("Responder", fontsize=11, color=COLOURS["Responder"])
    axes[1, 0].set_ylabel("Non-Responder", fontsize=11, color=COLOURS["Non-Responder"])
    fig.suptitle(
        "The same architectural spectrum occurs in both response groups\n"
        "grey = tumour, coloured = CD8.  Columns are matched on the exclusion score; "
        "left = intermixed, right = excluded.\n"
        f"Each column pairs a responder and a non-responder with near-identical "
        f"architecture ({len(band)} size-matched ROIs).",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    path = OUT / "fig_counterexamples.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")

    print(f"{len(d)} ROIs scored; {len(band)} in the size band "
          f"({lo:.0f}-{hi:.0f} tumour+CD8 cells)")
    print(f"exclusion score: median {d.exclusion.median():+.3f}, "
          f"range {d.exclusion.min():+.2f} to {d.exclusion.max():+.2f}")
    for grp in COLOURS:
        v = d[d.response == grp].exclusion
        print(f"  {grp:15s} median {v.median():+.3f}  "
              f"IQR [{v.quantile(.25):+.3f}, {v.quantile(.75):+.3f}]  n={len(v)}")
    for fov, lab in PANELS.items():
        if fov in set(d.roi_id):
            r = d.set_index("roi_id").loc[fov]
            print(f"  {lab} ({fov}): exclusion {r.exclusion:+.3f} "
                  f"({100*(d.exclusion < r.exclusion).mean():.0f}th pct)")
    print("wrote", path)


if __name__ == "__main__":
    main()
