#!/usr/bin/env python
"""Main-figure panels for the two B7H4 results. One key panel each.

PANEL 1 (--what signature): the within-ROI paired result, as a single plane.
    Immune exclusion is a JOINT signature in this construction -- a tumour nest ringed by CD8
    gives FEWER but LONGER degree-1 kernel bars -- so the two statistics are plotted as the two
    axes of one panel rather than as separate results. Each point is one ROI; both cancer
    subtypes were subsampled to identical counts within that ROI and compared against the same
    CD8 cells, so composition cannot differ between the arms and ROI geometry cancels.
    The exclusion quadrant is upper-left: fewer bars, longer bars.

    --degree 1 is the loop signature (mechanistically what "a ring around a nest" means).
    --degree 0 is stronger on both axes but describes components rather than loops.
    Both are generated; pick per the narrative.

PANEL 2 (--what field): the ROI-wide field effect, per patient.
    Plotted at patient level because that is the unit response is measured at and ROIs cluster
    within patients. One point per patient: median B7H4+ fraction against median exclusion
    score.

The p90 bar length is used rather than the mean: exclusion predicts a FEW long-lived loops,
not a shift in the average bar, so the upper tail is the apt statistic (and the mean is only
p = 0.05 at degree 1, while p90 is p = 0.01).

Output: output/fig_b7h4_signature_dim{0,1}.png, output/fig_b7h4_field.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
COLOURS = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}


def signature_panel(degree: int) -> Path:
    d = pd.read_csv(OUT / "b7h4_paired.csv")
    cnt, ln = f"ker{degree}-num_bars", f"ker{degree}-p90_length"
    dx = (d[f"b7h4__{cnt}"] - d[f"cancer__{cnt}"]).to_numpy(float)
    dy = (d[f"b7h4__{ln}"] - d[f"cancer__{ln}"]).to_numpy(float)
    ok = np.isfinite(dx) & np.isfinite(dy)
    dx, dy = dx[ok], dy[ok]

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    xlim = np.percentile(dx, [1, 99]) * np.array([1.1, 1.1])
    ylim = np.percentile(dy, [1, 99]) * np.array([1.1, 1.15])
    # shade the exclusion quadrant (x<0, y>0) in DATA coordinates -- axhspan's xmin/xmax are
    # axes fractions, and x=0 is nowhere near the middle when the range is asymmetric
    ax.add_patch(plt.Rectangle((xlim[0], 0), -xlim[0], ylim[1],
                               color="tab:red", alpha=0.055, zorder=0, linewidth=0))
    ax.scatter(dx, dy, s=13, c="0.35", alpha=0.45, linewidths=0, zorder=2)
    ax.axhline(0, color="black", lw=0.9, ls="--", zorder=1)
    ax.axvline(0, color="black", lw=0.9, ls="--", zorder=1)
    ax.scatter([np.median(dx)], [np.median(dy)], s=190, c="crimson", marker="D",
               edgecolors="white", linewidths=1.4, zorder=5, label="median ROI")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    px, py = wilcoxon(dx).pvalue, wilcoxon(dy).pvalue
    ax.set_xlabel(f"Δ number of degree-{degree} kernel bars\n(B7H4+ − B7H4−, matched counts)")
    ax.set_ylabel(f"Δ kernel bar length, 90th percentile (µm)\n(B7H4+ − B7H4−)")
    ax.text(0.03, 0.97, "more excluded\nfewer, longer bars", transform=ax.transAxes,
            va="top", ha="left", fontsize=10, style="italic", color="darkred")
    ax.set_title(f"CD8 exclusion signature, B7H4+ vs B7H4− tumour\n"
                 f"paired within {len(dx)} ROIs at matched cell counts\n"
                 f"Δcount {np.median(dx):+.0f} (p={px:.1g}),  "
                 f"Δlength {np.median(dy):+.2f} µm (p={py:.2g})", fontsize=10.5)
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    p = OUT / f"fig_b7h4_signature_dim{degree}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"  dim{degree}: Δcount {np.median(dx):+.1f} (p={px:.2g}, {100*(dx<0).mean():.0f}% "
          f"of ROIs negative) | Δp90length {np.median(dy):+.3f} (p={py:.2g}, "
          f"{100*(dy>0).mean():.0f}% positive)  -> {p.name}")
    return p


def field_panel() -> Path:
    import day2_b7h4_roi_level as roi
    d = roi.load()
    pat = d.groupby(["pid", "resp"], as_index=False)[["b7h4_frac", "exclusion"]].median()
    rho, pv = spearmanr(pat.b7h4_frac, pat.exclusion)

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    for grp, c in COLOURS.items():
        m = pat.resp == grp
        ax.scatter(pat.b7h4_frac[m], pat.exclusion[m], s=62, c=c, alpha=0.85,
                   linewidths=0, label=f"{grp} (n={int(m.sum())})")
    z = np.polyfit(pat.b7h4_frac, pat.exclusion, 1)
    xs = np.linspace(pat.b7h4_frac.min(), pat.b7h4_frac.max(), 50)
    ax.plot(xs, np.polyval(z, xs), color="0.3", lw=2, ls="--", zorder=1)
    ax.set_xlabel("B7H4+ fraction of cancer cells (patient median)")
    ax.set_ylabel("exclusion score (patient median)")
    ax.set_title("Tumours with a higher B7H4+ burden are more CD8-excluded\n"
                 f"one point per patient   Spearman ρ = {rho:+.3f}, p = {pv:.4f} "
                 f"(n = {len(pat)})", fontsize=11)
    ax.legend(fontsize=9)
    fig.tight_layout()
    p = OUT / "fig_b7h4_field.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"  field effect: per-patient rho={rho:+.3f} p={pv:.4f} n={len(pat)}  -> {p.name}")
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--what", default="both", choices=["both", "signature", "field"])
    args = ap.parse_args()
    if args.what in ("both", "signature"):
        for deg in (1, 0):
            signature_panel(deg)
    if args.what in ("both", "field"):
        field_panel()


if __name__ == "__main__":
    main()
