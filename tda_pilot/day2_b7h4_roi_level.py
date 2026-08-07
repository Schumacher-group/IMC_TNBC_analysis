#!/usr/bin/env python
"""Field effect: are ROIs with more B7H4+ tumour less infiltrated overall?

This is a DIFFERENT hypothesis from `day2_b7h4_paired.py`, and the paired design cannot test
it. That one subsamples both cancer subtypes to equal counts within the same ROI and compares
them against the same CD8 cells -- a cell-intrinsic question: is CD8 arranged differently
around B7H4+ than around B7H4- tumour in the same tissue? If B7H4+ cells instead suppress
infiltration across the whole microenvironment -- helping neighbouring B7H4- cells evade CD8 --
then both arms of the paired design sit in the same suppressed field and the effect cancels.

So ask it at ROI level: holding the number of CD8 cells and the number of B7H4- cancer cells
roughly fixed, are ROIs carrying more B7H4+ cancer less infiltrated?

OUTCOME. The per-ROI exclusion score used elsewhere in this directory: tumour (both subtypes
aggregated) against CD8, scored relative to that ROI's own label-permutation null, as

    exclusion = (relative deviation of ker1-avg_length − relative deviation of ker1-num_bars)/sqrt(2)

Relative deviation rather than z, because z divides by the null SD, which shrinks with ROI size
and density (rho = +0.31 against +0.13). Conditioning on each ROI's own null means the score
already accounts for that ROI's cell counts and composition, so it measures arrangement.

CONFOUND TO RESPECT. More B7H4+ cancer means more total tumour, hence a higher tumour:CD8
ratio. That is not a nuisance to be removed -- it is partly the mechanism under test -- but it
must not be the *only* thing driving a result. Two complementary analyses:

  CONTINUOUS  partial Spearman of exclusion against B7H4+ fraction of cancer, controlling for
              log CD8 count, log cancer count, density and hull area.
  MATCHED     ROIs binned on (n_CD8, n_B7H4-negative) by quantile, then within each bin the
              top vs bottom third of B7H4+ burden compared. Assumption-free, costs power.

Both at ROI level and aggregated per patient, since response is patient-level and ROIs cluster.

Output: output/B7H4_ROI_LEVEL.md, output/b7h4_roi_level.png
"""

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


def partial_spearman(x, y, covars):
    """Spearman correlation of x and y after linear regression on ranked covariates."""
    def resid(v):
        R = np.column_stack([np.ones(len(v))] + [pd.Series(c).rank().to_numpy() for c in covars])
        b, *_ = np.linalg.lstsq(R, pd.Series(v).rank().to_numpy(), rcond=None)
        return pd.Series(v).rank().to_numpy() - R @ b
    return spearmanr(resid(x), resid(y))


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, float(u.pvalue)


def load() -> pd.DataFrame:
    d = pd.read_parquet(OUT / "perm_null_CD8_primary_other_only.parquet")
    d = d[d["status"] == "ok"].copy()
    for s in ("ker1-avg_length", "ker1-num_bars"):
        m = d[f"{s}__null_mean"]
        d[f"rd_{s}"] = np.where(m != 0, (d[f"{s}__obs"] - m) / m, np.nan)
    d["exclusion"] = (d["rd_ker1-avg_length"] - d["rd_ker1-num_bars"]) / np.sqrt(2)

    counts = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "resp"})
    counts["n_b7h4"] = counts["tumour"] - 0  # placeholder, recomputed below
    geom = pd.read_csv(OUT / "roi_geometry.csv")

    # exact per-subtype counts from the point clouds
    rows = []
    for fov in d.roi_id:
        c = pd.read_csv(HERE / "data" / "rois" / f"{fov}.csv", usecols=["celltype"])
        rows.append({"roi_id": fov,
                     "n_b7h4": int((c.celltype == "B7H4 Cancer cell").sum()),
                     "n_cancer_neg": int((c.celltype == "Cancer cell").sum())})
    d = (d.merge(pd.DataFrame(rows), on="roi_id")
           .merge(counts[["roi_id", "pid", "resp", "cd8_primary"]], on="roi_id")
           .merge(geom, on="roi_id"))
    d = d.rename(columns={"cd8_primary": "n_cd8"})
    d["n_tumour_tot"] = d.n_b7h4 + d.n_cancer_neg
    d["b7h4_frac"] = d.n_b7h4 / d.n_tumour_tot
    d["cd8_ratio"] = d.n_cd8 / d.n_tumour_tot
    return d.dropna(subset=["exclusion"])


def md(dfr):
    h = "| " + " | ".join(dfr.columns) + " |"
    s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    b = "\n".join("| " + " | ".join(
        f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
        for v in row) + " |" for row in dfr.itertuples(index=False))
    return "\n".join([h, s, b])


def main() -> None:
    d = load()
    covars = [np.log(d.n_cd8 + 1), np.log(d.n_cancer_neg + 1),
              np.log(d.density), np.log(d.hull_area)]

    rows = [
        {"analysis": "raw Spearman(exclusion, B7H4+ fraction)",
         "rho": spearmanr(d.b7h4_frac, d.exclusion)[0],
         "p": spearmanr(d.b7h4_frac, d.exclusion)[1], "n": len(d)},
        {"analysis": "partial, controlling CD8 + B7H4- counts + geometry",
         "rho": partial_spearman(d.b7h4_frac, d.exclusion, covars)[0],
         "p": partial_spearman(d.b7h4_frac, d.exclusion, covars)[1], "n": len(d)},
        {"analysis": "raw Spearman(exclusion, B7H4+ COUNT)",
         "rho": spearmanr(d.n_b7h4, d.exclusion)[0],
         "p": spearmanr(d.n_b7h4, d.exclusion)[1], "n": len(d)},
    ]
    pat = d.groupby("pid", as_index=False)[["b7h4_frac", "exclusion"]].median()
    rows.append({"analysis": "per-patient Spearman(exclusion, B7H4+ fraction)",
                 "rho": spearmanr(pat.b7h4_frac, pat.exclusion)[0],
                 "p": spearmanr(pat.b7h4_frac, pat.exclusion)[1], "n": len(pat)})
    cont = pd.DataFrame(rows)

    # matched: bin on CD8 and B7H4-negative counts, contrast B7H4+ burden within bin
    d["bin"] = (pd.qcut(d.n_cd8, 4, labels=False, duplicates="drop").astype(str) + "_" +
                pd.qcut(d.n_cancer_neg, 4, labels=False, duplicates="drop").astype(str))
    hi, lo = [], []
    for _, g in d.groupby("bin"):
        if len(g) < 6:
            continue
        t = g.b7h4_frac.quantile([1 / 3, 2 / 3]).to_numpy()
        hi.append(g[g.b7h4_frac >= t[1]])
        lo.append(g[g.b7h4_frac <= t[0]])
    hi, lo = pd.concat(hi), pd.concat(lo)
    dl, pv = cliffs(hi.exclusion, lo.exclusion)
    dl_r, pv_r = cliffs(hi.cd8_ratio, lo.cd8_ratio)
    matched = pd.DataFrame([
        {"comparison": "exclusion score, high vs low B7H4+ within (CD8, B7H4-) bins",
         "median_high": hi.exclusion.median(), "median_low": lo.exclusion.median(),
         "cliffs_delta": dl, "p": pv, "n_high": len(hi), "n_low": len(lo)},
        {"comparison": "CD8:tumour ratio, same bins  [sanity: should differ by construction]",
         "median_high": hi.cd8_ratio.median(), "median_low": lo.cd8_ratio.median(),
         "cliffs_delta": dl_r, "p": pv_r, "n_high": len(hi), "n_low": len(lo)},
    ])

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    for grp, c in COLOURS.items():
        m = d.resp == grp
        ax[0].scatter(d.b7h4_frac[m], d.exclusion[m], s=11, c=c, alpha=0.5, lw=0, label=grp)
    ax[0].set_ylim(*np.percentile(d.exclusion, [1, 99]))
    ax[0].set_xlabel("B7H4+ fraction of cancer cells")
    ax[0].set_ylabel("exclusion score (vs own null)")
    ax[0].set_title(f"A. Raw association\nSpearman ρ = {cont.rho.iloc[0]:+.3f} "
                    f"(p = {cont.p.iloc[0]:.2g})", fontsize=10)
    ax[0].legend(fontsize=8)

    resid_y = partial_spearman(d.b7h4_frac, d.exclusion, covars)
    R = np.column_stack([np.ones(len(d))] + [pd.Series(c).rank().to_numpy() for c in covars])
    b, *_ = np.linalg.lstsq(R, d.exclusion.rank().to_numpy(), rcond=None)
    ry = d.exclusion.rank().to_numpy() - R @ b
    b2, *_ = np.linalg.lstsq(R, d.b7h4_frac.rank().to_numpy(), rcond=None)
    rx = d.b7h4_frac.rank().to_numpy() - R @ b2
    ax[1].scatter(rx, ry, s=11, c="tab:purple", alpha=0.45, lw=0)
    ax[1].set_xlabel("B7H4+ fraction (residual rank)")
    ax[1].set_ylabel("exclusion score (residual rank)")
    ax[1].set_title(f"B. After controlling CD8 + B7H4− counts, geometry\n"
                    f"partial ρ = {resid_y[0]:+.3f} (p = {resid_y[1]:.2g})", fontsize=10)

    jit = np.random.default_rng(0)
    for i, (lab, g) in enumerate([("low B7H4+", lo), ("high B7H4+", hi)]):
        v = g.exclusion.to_numpy()
        ax[2].scatter(jit.normal(i, 0.08, len(v)), v, s=10, c="0.4", alpha=0.45, lw=0)
        ax[2].hlines(np.median(v), i - 0.26, i + 0.26, color="crimson", lw=2.4)
    ax[2].set_ylim(*np.percentile(d.exclusion, [1, 99]))
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["low B7H4+", "high B7H4+"])
    ax[2].set_ylabel("exclusion score")
    ax[2].set_title(f"C. Matched on CD8 and B7H4− counts\n"
                    f"Cliff δ = {dl:+.3f} (p = {pv:.2g})", fontsize=10)
    fig.suptitle("Field effect: does a higher B7H4+ burden reduce infiltration ROI-wide?",
                 fontsize=12)
    fig.tight_layout()
    figp = OUT / "b7h4_roi_level.png"
    fig.savefig(figp, dpi=150, bbox_inches="tight")

    lines = [
        "# Field effect: are ROIs with more B7H4+ tumour less infiltrated overall?\n",
        f"{len(d)} pre-treatment ROIs, {d.pid.nunique()} patients. Outcome is the per-ROI "
        "exclusion score against that ROI's own label-permutation null (higher = more "
        "excluded). This asks a different question from the paired analysis in "
        "`B7H4_PAIRED.md`, which subsamples both subtypes within an ROI and so cannot see a "
        "field effect acting on the whole microenvironment.\n",
        "## Continuous\n", md(cont.round(4)),
        "\n## Matched on CD8 and B7H4− counts\n", md(matched.round(4)),
        "\nThe CD8:tumour row is a sanity check, not a result: adding B7H4+ cells necessarily "
        "lowers the ratio when CD8 and B7H4− are held fixed, so it must differ. The question "
        "is whether the *arrangement* score differs beyond that.\n",
        f"Figure: `{figp.name}`. Script: `day2_b7h4_roi_level.py`.",
    ]
    (OUT / "B7H4_ROI_LEVEL.md").write_text("\n".join(lines) + "\n")
    pd.set_option("display.width", 220)
    print(cont.round(4).to_string(index=False))
    print()
    print(matched.round(4).to_string(index=False))
    print("\nwrote", figp)


if __name__ == "__main__":
    main()
