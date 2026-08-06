#!/usr/bin/env python
"""Figure for the negative result: a calibrated null, not a blind measurement.

The point of the figure is to make two claims visible at once:

  (a) the readout WORKS -- the two Fig 4B panels, chosen by eye as intermixed vs excluded,
      land at opposite ends of the axis, in the mechanistically predicted direction;
  (b) the RESPONSE GROUPS OVERLAP -- so the architectural difference Fig 4B illustrates does
      not generalise to a systematic baseline difference between responders and non-responders.

Axes. Immune exclusion means CD8 ringing compact tumour nests. In the CD8-directional six-pack
that produces FEWER but LONGER degree-1 kernel bars (one big long-lived loop per nest, rather
than many short-lived loops from fine-grained intermixing). So the mechanistic plane is
`ker1-num_bars` vs `ker1-avg_length`, and the exclusion direction is up-and-left.

Values plotted are GEOMETRY-ADJUSTED z-scores:
  z    = (observed - mean of the ROI's own label-permutation null) / sd of that null
  adj  = residual of z on non-tumour fraction, log(cells), log(density), log(hull area)
The adjustment is not cosmetic. Raw bar counts and lengths scale with ROI size and density,
and z itself retains a geometry coupling (rho ~ 0.2-0.3, see REVIEW.md UPDATE 2) because the
null sd shrinks as an ROI gains cells. Plotting unadjusted values would show separation that
is partly the artefact this analysis retracted. `--raw` emits the unadjusted version for
comparison.

Panel C is the actual inference. Response is a patient-level outcome and ROIs cluster within
patients (6-11 each), so the per-ROI test in panel B is anticonservative -- it is shown
precisely to make that point, not as a result.

Output: output/negative_result_figure.png (+ _raw.png with --raw)
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import perm_null

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"

# Matches the manuscript slides and the Cell Counts palette.
COLOURS = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}
PANELS = {"Leap008_5": "Fig 4B responder\n(intermixed)",
          "Leap046_1": "Fig 4B non-responder\n(excluded)"}
X_STAT, Y_STAT = "ker1-num_bars", "ker1-avg_length"
N_LABEL_PERM = 20_000
SEED = 20260806


def cliffs(a, b) -> tuple[float, float]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def label_perm_p(pat: pd.DataFrame, col: str, rng) -> float:
    lab, v = pat.response.to_numpy(), pat[col].to_numpy(float)
    obs = cliffs(v[lab == "Responder"], v[lab == "Non-Responder"])[0]
    null = np.array([cliffs(v[(s := rng.permutation(lab)) == "Responder"],
                            v[s == "Non-Responder"])[0] for _ in range(N_LABEL_PERM)])
    return float((np.abs(null) >= abs(obs)).mean())


def load(defn: str, adjust: bool) -> pd.DataFrame:
    geom = pd.read_csv(OUT / "roi_geometry.csv")
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response"},
    )[["roi_id", "patient_id", "response"]]
    d = pd.read_parquet(OUT / f"perm_null_{defn}_other_only.parquet")
    d = d[d["status"] == "ok"].merge(meta, on="roi_id").merge(geom, on="roi_id")
    d = d[d.response.isin(COLOURS)].copy()

    # Minimal, NON-REDUNDANT geometry basis. density = cells/area, so adding log(density)
    # to log(cells) and log(area) makes the design near-collinear (condition number ~450)
    # and over-adjusts: see day2_sensitivity.py / SENSITIVITY.md.
    n_pair = d.n_tumour + d.n_other
    X = np.column_stack([np.ones(len(d)), np.log(n_pair), np.log(d.hull_area)])
    for stat in (X_STAT, Y_STAT):
        y = d[f"{stat}__z"].to_numpy(float)
        if adjust:
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            y = y - X @ beta
        d[stat] = y
    # exclusion score: fewer bars AND longer bars, projected onto that diagonal
    d["exclusion"] = (d[Y_STAT] - d[X_STAT]) / np.sqrt(2)
    return d


def robust_lim(v, pad=0.06):
    lo, hi = np.percentile(v, [1, 99])
    m = (hi - lo) * pad
    return lo - m, hi + m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--definition", default="CD8_primary")
    ap.add_argument("--raw", action="store_true", help="plot unadjusted z (for comparison)")
    args = ap.parse_args()
    adjust = not args.raw
    d = load(args.definition, adjust)
    rng = np.random.default_rng(SEED)
    unit = "geometry-adjusted z" if adjust else "raw z"

    fig = plt.figure(figsize=(16.5, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 0.85], wspace=0.28,
                          top=0.72, bottom=0.13, left=0.05, right=0.985)

    # ---- Panel A: mechanistic plane ------------------------------------------
    axA = fig.add_subplot(gs[0, 0])
    for grp, c in COLOURS.items():
        m = d.response == grp
        axA.scatter(d[X_STAT][m], d[Y_STAT][m], s=13, c=c, alpha=0.45, linewidths=0,
                    label=f"{grp} (n={int(m.sum())} ROIs)")
    axA.axhline(0, color="grey", lw=0.6, ls=":")
    axA.axvline(0, color="grey", lw=0.6, ls=":")
    xl, yl = robust_lim(d[X_STAT]), robust_lim(d[Y_STAT])
    axA.set_xlim(*xl)
    axA.set_ylim(*yl)
    # exclusion direction: fewer bars (left) + longer bars (up)
    axA.annotate("", xy=(xl[0] + 0.16 * (xl[1] - xl[0]), yl[0] + 0.86 * (yl[1] - yl[0])),
                 xytext=(xl[0] + 0.52 * (xl[1] - xl[0]), yl[0] + 0.50 * (yl[1] - yl[0])),
                 arrowprops=dict(arrowstyle="-|>", lw=2.0, color="0.25"))
    axA.text(xl[0] + 0.13 * (xl[1] - xl[0]), yl[0] + 0.90 * (yl[1] - yl[0]),
             "more excluded\n(fewer, longer loops)", fontsize=8.5, color="0.2",
             ha="left", va="bottom", style="italic")
    axA.text(xl[0] + 0.66 * (xl[1] - xl[0]), yl[0] + 0.24 * (yl[1] - yl[0]),
             "more intermixed\n(many, short loops)", fontsize=8.5, color="0.2",
             ha="left", va="top", style="italic")
    for fov, lab in PANELS.items():
        if fov not in set(d.roi_id):
            continue
        r = d.set_index("roi_id").loc[fov]
        axA.scatter([r[X_STAT]], [r[Y_STAT]], s=190, facecolors="none",
                    edgecolors="black", linewidths=2.0, zorder=5)
        axA.annotate(lab, (r[X_STAT], r[Y_STAT]), textcoords="offset points",
                     xytext=(12, -6), fontsize=8.5, fontweight="bold",
                     color=COLOURS[r.response],
                     arrowprops=dict(arrowstyle="-", lw=0.8, color="black"))
    axA.set_xlabel(f"kernel degree-1 bar COUNT  ({unit})")
    axA.set_ylabel(f"kernel degree-1 mean bar LENGTH  ({unit})")
    axA.set_title("A. Tumour–CD8 architecture, per ROI\n"
                  "(axes clipped to 1–99th pct; a few extreme ROIs lie outside)",
                  fontsize=10)
    axA.legend(loc="lower left", fontsize=8, framealpha=0.9)

    # ---- Panel B: per-ROI exclusion score -------------------------------------
    axB = fig.add_subplot(gs[0, 1])
    lim = robust_lim(d.exclusion)
    bins = np.linspace(*lim, 46)
    for grp, c in COLOURS.items():
        v = d[d.response == grp].exclusion
        axB.hist(v, bins=bins, color=c, alpha=0.45, density=True, label=grp)
        axB.axvline(v.median(), color=c, lw=2, ls="--")
    for fov in PANELS:
        if fov not in set(d.roi_id):
            continue
        r = d.set_index("roi_id").loc[fov]
        axB.axvline(r.exclusion, color="black", lw=1.6)
        axB.annotate(fov, (r.exclusion, axB.get_ylim()[1] * 0.96),
                     rotation=90, fontsize=8, ha="right", va="top", fontweight="bold")
    dr, pr = cliffs(d[d.response == "Responder"].exclusion,
                    d[d.response == "Non-Responder"].exclusion)
    axB.set_xlim(*lim)
    axB.set_xlabel(f"exclusion score  ({unit})\n← intermixed        excluded →")
    axB.set_ylabel("density")
    axB.set_title(f"B. Per ROI (n={len(d)}): distributions overlap\n"
                  f"Cliff δ={dr:+.3f}, p={pr:.3f} — but ROIs cluster\n"
                  f"within patients, so this p is anticonservative", fontsize=10)
    axB.legend(fontsize=8)

    # ---- Panel C: per-patient (the actual inference) --------------------------
    axC = fig.add_subplot(gs[0, 2])
    pat = d.groupby(["patient_id", "response"], as_index=False)["exclusion"].median()
    jit = np.random.default_rng(0)
    for i, grp in enumerate(["Non-Responder", "Responder"]):
        v = pat[pat.response == grp].exclusion.to_numpy()
        axC.scatter(jit.normal(i, 0.075, len(v)), v, s=42, c=COLOURS[grp], alpha=0.85,
                    linewidths=0)
        axC.hlines(np.median(v), i - 0.26, i + 0.26, color="black", lw=2.4, zorder=4)
    dp, _ = cliffs(pat[pat.response == "Responder"].exclusion,
                   pat[pat.response == "Non-Responder"].exclusion)
    pp = label_perm_p(pat, "exclusion", rng)
    axC.axhline(0, color="grey", lw=0.6, ls=":")
    axC.set_xticks([0, 1])
    axC.set_xticklabels([f"NR\n(n={int((pat.response=='Non-Responder').sum())})",
                         f"R\n(n={int((pat.response=='Responder').sum())})"])
    axC.set_ylabel(f"patient-median exclusion score ({unit})")
    axC.set_title(f"C. Per patient — the inference\nCliff δ={dp:+.3f}, "
                  f"permutation p={pp:.3f}\nδ spans +0.01…+0.40 across geometry\n"
                  f"adjustments — not robust either way", fontsize=9.5)

    fig.suptitle(
        "Tumour–CD8 architecture is measurable and varies widely; no robust evidence that it "
        "distinguishes response at baseline\n"
        "Chromatic TDA, CD8-directional six-pack, pre-treatment; each ROI scored against its "
        "own label-permutation null (cell positions and composition held fixed). Any residual "
        "trend runs OPPOSITE to the immune-exclusion hypothesis.",
        fontsize=11.5, y=0.985)
    path = OUT / ("negative_result_figure.png" if adjust else "negative_result_figure_raw.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")

    print(f"panel B per-ROI     : Cliff δ={dr:+.3f} p={pr:.4f}")
    print(f"panel C per-patient : Cliff δ={dp:+.3f} permutation p={pp:.4f}")
    for fov in PANELS:
        if fov in set(d.roi_id):
            r = d.set_index("roi_id").loc[fov]
            pct = 100 * (d.exclusion < r.exclusion).mean()
            print(f"  {fov:12s} exclusion={r.exclusion:+7.2f}  ({pct:.0f}th pct)")
    print("wrote", path)


if __name__ == "__main__":
    main()
