#!/usr/bin/env python
"""Supplementary figure: why the z-scores need a second, cross-ROI geometry adjustment.

There are TWO separate corrections in this analysis and they do different jobs.

STEP 1 — the per-ROI label-permutation null (this is the shuffling step).
    For one ROI, keep every cell exactly where it is and randomly reassign which cells are
    "tumour" and which are "CD8", preserving the exact number of each. Recompute the full
    chromatic filtration and six-pack. Repeat 99x. That gives a null distribution for that
    ROI, and z = (observed − null mean) / null sd.
    Because the shuffle reuses the same points and the same counts, it controls WITHIN that
    ROI for: cell number, cell density, tissue shape, gaps in the tissue, and the CD8:tumour
    ratio. z therefore answers: "is this ROI's arrangement unusual for its own cells?"

STEP 2 — the cross-ROI regression (this figure; no shuffling involved).
    z is a ratio, and its denominator is not a constant. The null sd shrinks as an ROI gains
    cells or density, so a large dense ROI converts a modest arrangement deviation into a
    large z, while a small sparse ROI does not. z is a signal-to-noise ratio whose SCALE
    depends on ROI geometry — comparable within an ROI, not across ROIs.
    That matters here because responder and non-responder ROIs differ systematically in
    density, so a group difference in raw z can be manufactured by geometry alone. The fix is
    to regress z across ROIs on composition and geometry and keep the residual.

Panels:
    A  ROI density by response — is there a group imbalance to worry about? (yes)
    B  raw z vs density, coloured by response — the coupling, with Spearman rho
    C  adjusted z vs density — coupling removed
    D  rho(z, density) before/after for all four cell types — the coupling is a property of
       the method, not of CD8

Output: output/geometry_supplement.png
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
STAT = "ker1-avg_length"
DEFS = ["CD8_primary", "CD4", "Fibroblast", "Macrophage"]
PANELS = {"Leap008_5": "Fig 4B R", "Leap046_1": "Fig 4B NR"}


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def load(defn: str) -> pd.DataFrame:
    geom = pd.read_csv(OUT / "roi_geometry.csv")
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response"},
    )[["roi_id", "patient_id", "response"]]
    d = pd.read_parquet(OUT / f"perm_null_{defn}_other_only.parquet")
    d = d[d["status"] == "ok"].merge(meta, on="roi_id").merge(geom, on="roi_id")
    d = d[d.response.isin(COLOURS)].copy()
    d["z"] = d[f"{STAT}__z"]
    n_pair = d.n_tumour + d.n_other
    X = np.column_stack([np.ones(len(d)), d.n_other / n_pair, np.log(n_pair),
                         np.log(d.density), np.log(d.hull_area)])
    beta, *_ = np.linalg.lstsq(X, d.z.to_numpy(float), rcond=None)
    d["z_adj"] = d.z - X @ beta
    return d


def binned(ax, x, y, c):
    q = pd.qcut(x, 8, duplicates="drop")
    g = pd.DataFrame({"x": x, "y": y, "b": q}).groupby("b", observed=True)
    ax.plot(g["x"].median(), g["y"].median(), c=c, lw=2.4, marker="o", ms=4, zorder=4)


def scatter_panel(ax, d, ycol, title, ylab):
    for grp, c in COLOURS.items():
        m = d.response == grp
        ax.scatter(d.density[m] * 1e3, d[ycol][m], s=12, c=c, alpha=0.42, linewidths=0,
                   label=grp)
        binned(ax, d.density[m].to_numpy() * 1e3, d[ycol][m].to_numpy(), c)
    rho, p = spearmanr(d.density, d[ycol])
    ax.axhline(0, color="grey", lw=0.6, ls=":")
    ax.set_ylim(*np.percentile(d[ycol], [1, 99]) + np.array([-2, 2]))
    ax.set_xlabel("ROI cell density (cells per 1000 µm²)")
    ax.set_ylabel(ylab)
    ax.set_title(f"{title}\nSpearman ρ = {rho:+.3f} (p={p:.1g})", fontsize=10)
    for fov, lab in PANELS.items():
        if fov in set(d.roi_id):
            r = d.set_index("roi_id").loc[fov]
            ax.scatter([r.density * 1e3], [r[ycol]], s=150, facecolors="none",
                       edgecolors="black", linewidths=1.8, zorder=6)
            ax.annotate(lab, (r.density * 1e3, r[ycol]), textcoords="offset points",
                        xytext=(9, 4), fontsize=8, fontweight="bold")
    return rho


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--definition", default="CD8_primary")
    args = ap.parse_args()
    d = load(args.definition)

    fig = plt.figure(figsize=(18.5, 5.4))
    gs = fig.add_gridspec(1, 4, wspace=0.30, top=0.74, bottom=0.16, left=0.045, right=0.99)

    # A: density by response -- is there an imbalance?
    axA = fig.add_subplot(gs[0, 0])
    bins = np.linspace(*np.percentile(d.density * 1e3, [0.5, 99.5]), 40)
    for grp, c in COLOURS.items():
        v = d[d.response == grp].density * 1e3
        axA.hist(v, bins=bins, color=c, alpha=0.45, density=True, label=f"{grp} (n={len(v)})")
        axA.axvline(v.median(), color=c, lw=2, ls="--")
    dd, pp = cliffs(d[d.response == "Responder"].density,
                    d[d.response == "Non-Responder"].density)
    pat = d.groupby(["patient_id", "response"], as_index=False)["density"].median()
    dp, ppat = cliffs(pat[pat.response == "Responder"].density,
                      pat[pat.response == "Non-Responder"].density)
    axA.set_xlabel("ROI cell density (cells per 1000 µm²)")
    axA.set_ylabel("density")
    axA.set_title(f"A. Is there a geometry imbalance between groups?\n"
                  f"per ROI Cliff δ={dd:+.3f} (p={pp:.1g}); "
                  f"per patient δ={dp:+.3f} (p={ppat:.2f})", fontsize=10)
    axA.legend(fontsize=8)

    axB = fig.add_subplot(gs[0, 1])
    rho_raw = scatter_panel(axB, d, "z", "B. Raw z is coupled to density",
                            f"{STAT}  (raw z)")
    axB.legend(fontsize=8, loc="upper left")

    axC = fig.add_subplot(gs[0, 2])
    scatter_panel(axC, d, "z_adj", "C. After cross-ROI adjustment", f"{STAT}  (adjusted z)")

    # D: coupling across all species
    axD = fig.add_subplot(gs[0, 3])
    rows = []
    for defn in DEFS:
        f = OUT / f"perm_null_{defn}_other_only.parquet"
        if not f.exists():
            continue
        dd_ = load(defn)
        rows.append({"defn": defn,
                     "raw": spearmanr(dd_.density, dd_.z)[0],
                     "adj": spearmanr(dd_.density, dd_.z_adj)[0]})
    r = pd.DataFrame(rows)
    x = np.arange(len(r))
    axD.bar(x - 0.19, r["raw"], 0.38, label="raw z", color="0.45")
    axD.bar(x + 0.19, r["adj"], 0.38, label="adjusted z", color="tab:blue")
    axD.axhline(0, color="black", lw=0.8)
    axD.set_xticks(x)
    axD.set_xticklabels(r.defn, rotation=20, ha="right", fontsize=8)
    axD.set_ylabel("Spearman ρ(statistic, ROI density)")
    axD.set_title("D. The coupling is a property of the method,\nnot of CD8 — every cell "
                  "type shows it", fontsize=10)
    axD.legend(fontsize=8)

    fig.suptitle(
        "Supplementary: the per-ROI permutation null controls arrangement WITHIN an ROI; a "
        "second cross-ROI adjustment is needed to compare BETWEEN ROIs\n"
        "z = (observed − mean of that ROI's label-shuffled null) / sd of that null. The null "
        "sd shrinks in larger, denser ROIs, so raw z is a signal-to-noise ratio whose scale "
        "depends on ROI geometry.",
        fontsize=11, y=0.985)
    path = OUT / "geometry_supplement.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(r.round(3).to_string(index=False))
    print(f"\ndensity by response: per-ROI δ={dd:+.3f} p={pp:.2g}; "
          f"per-patient δ={dp:+.3f} p={ppat:.3f}")
    print(f"rho(raw z, density) = {rho_raw:+.3f}")
    print("wrote", path)


if __name__ == "__main__":
    main()
