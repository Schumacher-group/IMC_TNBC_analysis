#!/usr/bin/env python
"""Step 1 of the presentation: the RAW statistics, and how they depend on ROI geometry.

No z-scores, no permutation nulls, no adjustment. These are the numbers that come straight
out of the persistence diagrams -- the same quantities the M2S2 pipeline feeds to its
classifier -- so this is the natural starting point before any correction is introduced.

Figure 1 (raw_statistics_figure.png)
  A  the mechanistic plane: degree-1 kernel bar COUNT vs mean bar LENGTH, per ROI,
     coloured by response, with the two Fig 4B panels highlighted
  B  the same axes, but coloured by ROI cell density, to show the geometry gradient
  C  raw mean bar length vs density
  D  raw bar count vs number of cells

Figure 2 (raw_statistics_geometry.png)
  Spearman correlations of every retained raw statistic against each geometry variable,
  as a heatmap -- a compact statement of how much of the raw signal is geometry.

Inclusion is the CD8-directional six-pack (`--inclusion other_only`, the default) since that
is the one relevant to an exclusion hypothesis; `--inclusion symmetric` gives the other.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
COLOURS = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}
PANELS = {"Leap008_5": "Fig 4B R\n(intermixed)", "Leap046_1": "Fig 4B NR\n(excluded)"}
X_STAT, Y_STAT = "ker1-num_bars", "ker1-avg_length"
GEOM = {"cell count": "n_pair", "cell density": "density",
        "hull area": "hull_area", "median NN dist": "med_nn", "CD8 fraction": "ofrac"}
STATS = ["ker1-avg_length", "ker1-med_length", "ker1-p90_length", "ker1-num_bars",
         "ker0-avg_length", "ker0-num_bars", "im1-avg_length", "cok1-avg_length"]


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def load(defn: str, inclusion: str) -> pd.DataFrame:
    tag = "" if inclusion == "symmetric" else f"_{inclusion}"
    geom = pd.read_csv(OUT / "roi_geometry.csv")
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response"},
    )[["roi_id", "patient_id", "response"]]
    d = pd.read_parquet(OUT / f"perm_null_{defn}{tag}.parquet")
    d = d[d["status"] == "ok"].merge(meta, on="roi_id").merge(geom, on="roi_id")
    d = d[d.response.isin(COLOURS)].copy()
    for s in STATS:                       # RAW observed values, straight from the diagrams
        d[s] = d[f"{s}__obs"]
    d["n_pair"] = d.n_tumour + d.n_other
    d["ofrac"] = d.n_other / d.n_pair
    return d


def mark_panels(ax, d, xcol, ycol, dx=10, dy=6):
    for fov, lab in PANELS.items():
        if fov not in set(d.roi_id):
            continue
        r = d.set_index("roi_id").loc[fov]
        ax.scatter([r[xcol]], [r[ycol]], s=175, facecolors="none", edgecolors="black",
                   linewidths=2.0, zorder=6)
        ax.annotate(lab, (r[xcol], r[ycol]), textcoords="offset points", xytext=(dx, dy),
                    fontsize=8, fontweight="bold", color=COLOURS[r.response],
                    arrowprops=dict(arrowstyle="-", lw=0.8, color="black"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--definition", default="CD8_primary")
    ap.add_argument("--inclusion", default="other_only", choices=["other_only", "symmetric"])
    args = ap.parse_args()
    d = load(args.definition, args.inclusion)

    # ---------------- Figure 1 ----------------
    fig = plt.figure(figsize=(17.5, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.24,
                          top=0.88, bottom=0.07, left=0.055, right=0.985)

    # A handful of tiny/degenerate ROIs give extreme mean lengths; clip so the bulk is legible.
    ylim = (0, float(np.percentile(d[Y_STAT], 99.3)))

    axA = fig.add_subplot(gs[0, 0])
    for grp, c in COLOURS.items():
        m = d.response == grp
        axA.scatter(d[X_STAT][m], d[Y_STAT][m], s=14, c=c, alpha=0.5, linewidths=0,
                    label=f"{grp} (n={int(m.sum())})")
    mark_panels(axA, d, X_STAT, Y_STAT)
    axA.set_ylim(*ylim)
    dl, pv = cliffs(d[d.response == "Responder"][Y_STAT],
                    d[d.response == "Non-Responder"][Y_STAT])
    axA.set_xlabel(f"{X_STAT}  (raw count of degree-1 kernel bars)")
    axA.set_ylabel(f"{Y_STAT}  (raw mean bar lifetime, µm)")
    axA.set_title(f"A. Raw statistics, coloured by response\n"
                  f"mean-length Cliff δ={dl:+.3f} (p={pv:.1g}), per ROI", fontsize=10)
    axA.legend(fontsize=8, loc="upper right")

    axB = fig.add_subplot(gs[0, 1])
    sc = axB.scatter(d[X_STAT], d[Y_STAT], s=14, c=d.ofrac, cmap="viridis",
                     alpha=0.8, linewidths=0)
    plt.colorbar(sc, ax=axB, label="CD8 fraction of the tumour+CD8 pair")
    mark_panels(axB, d, X_STAT, Y_STAT)
    axB.set_ylim(*ylim)
    axB.set_xlabel(f"{X_STAT}  (raw)")
    axB.set_ylabel(f"{Y_STAT}  (raw)")
    axB.set_title("B. Same axes, coloured by CD8 fraction\n"
                  "the cloud is organised by COMPOSITION, not by response", fontsize=10)

    axC = fig.add_subplot(gs[1, 0])
    for grp, c in COLOURS.items():
        m = d.response == grp
        axC.scatter(d.ofrac[m], d[Y_STAT][m], s=14, c=c, alpha=0.5, linewidths=0, label=grp)
    rho, p = spearmanr(d.ofrac, d[Y_STAT])
    mark_panels(axC, d.assign(_x=d.ofrac), "_x", Y_STAT, dx=12, dy=14)
    axC.set_ylim(*ylim)
    axC.set_xlabel("CD8 fraction  =  CD8 / (tumour + CD8)")
    axC.set_ylabel(f"{Y_STAT}  (raw, µm)")
    axC.set_title(f"C. Raw mean bar length is mostly CD8 FRACTION\n"
                  f"Spearman ρ = {rho:+.3f} (p={p:.1g}) — sparse CD8 ⇒ long bars",
                  fontsize=10)
    axC.legend(fontsize=8)

    axD = fig.add_subplot(gs[1, 1])
    for grp, c in COLOURS.items():
        m = d.response == grp
        axD.scatter(d.n_pair[m], d[X_STAT][m], s=14, c=c, alpha=0.5, linewidths=0, label=grp)
    rho2, p2 = spearmanr(d.n_pair, d[X_STAT])
    mark_panels(axD, d.assign(_x=d.n_pair), "_x", X_STAT)
    axD.set_xlabel("tumour + CD8 cells in ROI")
    axD.set_ylabel(f"{X_STAT}  (raw)")
    axD.set_title(f"D. Raw bar count vs number of cells\nSpearman ρ = {rho2:+.3f} "
                  f"(p={p2:.1g}) — near-deterministic", fontsize=10)
    axD.legend(fontsize=8)

    fig.suptitle(
        "Raw chromatic-TDA statistics, before any null model or adjustment  —  "
        f"{args.definition}, {args.inclusion} inclusion, pre-treatment (n={len(d)} ROIs)\n"
        "Panels C and D are the problem: bar COUNT is essentially cell number (ρ=0.83) and "
        "bar LENGTH is essentially CD8 fraction (ρ=−0.69).\n"
        "So the raw 'exclusion-like' corner (few, long bars) is mostly just 'few CD8 cells' "
        "— a composition statement, not a spatial-arrangement one.",
        fontsize=12, y=0.975)
    p1 = OUT / "raw_statistics_figure.png"
    fig.savefig(p1, dpi=150, bbox_inches="tight")

    # ---------------- Figure 2: correlation heatmap ----------------
    rows = []
    for s in STATS:
        rows.append({"statistic": s,
                     **{g: spearmanr(d[col], d[s])[0] for g, col in GEOM.items()}})
    corr = pd.DataFrame(rows).set_index("statistic")

    fig2, ax2 = plt.subplots(figsize=(7.6, 5.4))
    im = ax2.imshow(corr.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax2.set_xticks(range(len(corr.columns)))
    ax2.set_xticklabels(corr.columns, rotation=25, ha="right", fontsize=9)
    ax2.set_yticks(range(len(corr.index)))
    ax2.set_yticklabels(corr.index, fontsize=9)
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            v = corr.iat[i, j]
            ax2.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8,
                     color="white" if abs(v) > 0.55 else "black")
    plt.colorbar(im, ax=ax2, label="Spearman ρ")
    ax2.set_title("Raw TDA statistics vs ROI geometry\n"
                  "(pre-treatment; every raw statistic is geometry-coupled)", fontsize=11)
    fig2.tight_layout()
    p2f = OUT / "raw_statistics_geometry.png"
    fig2.savefig(p2f, dpi=150, bbox_inches="tight")

    pd.set_option("display.width", 200)
    print("Spearman ρ of raw statistics against geometry:")
    print(corr.round(3).to_string())
    print(f"\nFig 4B panels, raw values:")
    print(d[d.roi_id.isin(PANELS)][["roi_id", "response", X_STAT, Y_STAT, "n_pair",
                                    "density", "hull_area"]].round(3).to_string(index=False))
    print("\nwrote", p1, "\nwrote", p2f)


if __name__ == "__main__":
    main()
