#!/usr/bin/env python
"""Pre-check: are B7H4+ and B7H4- cancer cells spatially separated enough to compare?

The proposed question is whether CD8 cells are more excluded from B7H4+ tumour than from
B7H4- tumour. That comparison only has a chance of working if the two cancer subtypes occupy
distinguishable territory. If they are finely intermingled inside the same nests, then any
CD8 cell near one is near the other by construction, both statistics measure the same tissue,
and a null result would be uninformative -- a design failure rather than a biological finding.

So measure the mixing first, cheaply, before spending compute on six-packs.

Mixing index, per ROI, over cancer cells only:
    observed  = mean fraction of a cancer cell's k nearest cancer neighbours that carry the
                OTHER B7H4 status
    expected  = the same quantity after shuffling B7H4 status among those same cancer cells,
                which preserves the counts of each subtype exactly
    ratio     = observed / expected
      ratio ~ 1  -> intermingled at random; the two subtypes are not spatially distinguishable
                    and the paired CD8 comparison is not worth running
      ratio < 1  -> segregated into distinct territory; the comparison is meaningful
      ratio > 1  -> actively interleaved, more alternating than chance

This uses only nearest-neighbour counting -- no persistence -- so it runs in seconds and is
independent of every modelling choice elsewhere in this directory.

Output: output/B7H4_PRECHECK.md, output/b7h4_precheck.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import cohort

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
ROIS_DIR = HERE / "data" / "rois"
CANCER, B7H4 = "Cancer cell", "B7H4 Cancer cell"
K = 10
N_SHUFFLE = 20
MIN_EACH = 30          # need enough of both subtypes for the index to mean anything
SEED = 20260806


def mixing_for_roi(fov: str, rng) -> dict | None:
    df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
    sub = df[df["celltype"].isin([CANCER, B7H4])]
    n_c = int((sub["celltype"] == CANCER).sum())
    n_b = int((sub["celltype"] == B7H4).sum())
    if min(n_c, n_b) < MIN_EACH:
        return None

    pts = sub[["x", "y"]].to_numpy(float)
    is_b = (sub["celltype"] == B7H4).to_numpy()
    k = min(K, len(pts) - 1)
    # k+1 because the first neighbour of a point is itself
    _, idx = cKDTree(pts).query(pts, k=k + 1)
    nbr = idx[:, 1:]

    def hetero_frac(labels):
        return float(np.mean(labels[nbr] != labels[:, None]))

    obs = hetero_frac(is_b)
    exp = float(np.mean([hetero_frac(rng.permutation(is_b)) for _ in range(N_SHUFFLE)]))
    return {"roi_id": fov, "n_cancer": n_c, "n_b7h4": n_b,
            "b7h4_frac": n_b / (n_b + n_c), "observed": obs, "expected": exp,
            "ratio": obs / exp if exp > 0 else np.nan}


def main() -> None:
    rng = np.random.default_rng(SEED)
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Response": "response",
                 "Sample_Type_(pre/post treatment)": "sample_type"})
    meta = meta[(meta.sample_type.astype(str).str.lower() == "pre")
                & meta.response.isin(["Responder", "Non-Responder"])]
    rois = [r for r in meta.roi_id if cohort.is_cohort_member(r)
            and (ROIS_DIR / f"{r}.csv").exists()]

    rows = [m for r in rois if (m := mixing_for_roi(r, rng)) is not None]
    d = pd.DataFrame(rows).merge(meta[["roi_id", "response"]], on="roi_id")
    d.to_csv(OUT / "b7h4_precheck.csv", index=False)

    med = d.ratio.median()
    frac_segregated = float((d.ratio < 0.9).mean())
    verdict = (
        "**intermingled** — the two cancer subtypes are not spatially distinguishable, so a "
        "paired CD8 comparison would be measuring the same tissue twice"
        if med > 0.95 else
        "**spatially segregated** — the subtypes occupy distinguishable territory, so a "
        "paired CD8 comparison is meaningful"
    )

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    ax[0].hist(d.ratio, bins=40, color="tab:purple", alpha=0.75)
    ax[0].axvline(1.0, color="black", lw=1.5, ls="--", label="random intermingling")
    ax[0].axvline(med, color="crimson", lw=2, label=f"median = {med:.3f}")
    ax[0].set_xlabel("mixing ratio  (observed / label-shuffled)")
    ax[0].set_ylabel("ROIs")
    ax[0].set_title(f"A. Are B7H4+ and B7H4− cancer cells mixed?\n{len(d)} pre-treatment ROIs",
                    fontsize=10)
    ax[0].legend(fontsize=8)

    ax[1].scatter(d.b7h4_frac, d.ratio, s=12, c="tab:purple", alpha=0.5, linewidths=0)
    ax[1].axhline(1.0, color="black", lw=1.2, ls="--")
    ax[1].set_xlabel("B7H4+ fraction of cancer cells")
    ax[1].set_ylabel("mixing ratio")
    ax[1].set_title("B. Mixing vs subtype balance\n(is the index just tracking abundance?)",
                    fontsize=10)

    cols = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}
    for i, (grp, c) in enumerate(cols.items()):
        v = d[d.response == grp].ratio.to_numpy()
        ax[2].scatter(np.random.default_rng(0).normal(i, 0.07, len(v)), v, s=14, c=c,
                      alpha=0.55, linewidths=0)
        ax[2].hlines(np.median(v), i - 0.25, i + 0.25, color="black", lw=2.2)
    ax[2].axhline(1.0, color="black", lw=1.2, ls="--")
    ax[2].set_xticks([0, 1])
    ax[2].set_xticklabels(["Responder", "Non-Responder"])
    ax[2].set_ylabel("mixing ratio")
    ax[2].set_title("C. By response\n(descriptive only)", fontsize=10)

    fig.suptitle("Pre-check: do B7H4+ and B7H4− cancer cells occupy distinguishable territory?",
                 fontsize=12)
    fig.tight_layout()
    figp = OUT / "b7h4_precheck.png"
    fig.savefig(figp, dpi=150, bbox_inches="tight")

    lines = [
        "# Pre-check: are B7H4+ and B7H4− cancer cells spatially separated?\n",
        f"{len(d)} pre-treatment ROIs with at least {MIN_EACH} cells of each subtype "
        f"(of {len(rois)} eligible). k = {K} nearest cancer neighbours, expectation from "
        f"{N_SHUFFLE} shuffles of B7H4 status among the same cells.\n",
        f"- Median mixing ratio (observed / shuffled) = **{med:.3f}**",
        f"- Interquartile range: {d.ratio.quantile(.25):.3f} – {d.ratio.quantile(.75):.3f}",
        f"- ROIs with ratio < 0.9 (clearly segregated): **{frac_segregated:.0%}**",
        f"- Median B7H4+ share of cancer cells: {d.b7h4_frac.median():.1%}\n",
        f"**Verdict: {verdict}.**\n",
        "Why this matters: it is the design check for the proposed "
        "'is CD8 more excluded from B7H4+ tumour?' comparison. A ratio near 1 would mean any "
        "CD8 cell near one subtype is near the other by construction, so both statistics "
        "would describe the same tissue and a null result would be uninformative.\n",
        "Method note: this is nearest-neighbour counting only, no persistent homology, so it "
        "is independent of every modelling choice elsewhere in this directory.\n",
        f"Figure: `{figp.name}`. Script: `day2_b7h4_precheck.py`.",
    ]
    (OUT / "B7H4_PRECHECK.md").write_text("\n".join(lines) + "\n")

    print(f"ROIs analysed: {len(d)} of {len(rois)} eligible")
    print(f"median mixing ratio = {med:.3f}  (1.0 = random intermingling)")
    print(f"IQR {d.ratio.quantile(.25):.3f}-{d.ratio.quantile(.75):.3f}; "
          f"{frac_segregated:.0%} of ROIs below 0.9")
    print(f"median B7H4+ share of cancer cells: {d.b7h4_frac.median():.1%}")
    print("wrote", figp, "and", OUT / "B7H4_PRECHECK.md")


if __name__ == "__main__":
    main()
