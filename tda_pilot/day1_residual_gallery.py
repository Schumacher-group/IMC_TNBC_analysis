#!/usr/bin/env python
"""Size-controlled version of the segregation gallery.

The raw image-persistence sort is dominated by cell count (highest = big dense
ROIs, lowest = tiny sparse ones). To isolate topology we:
  1. restrict to a mid-size band (cohort IQR of tumour+CD8 count), and
  2. sort by the size-adjusted residual of im_dim1_total_persistence
     (residual after OLS on n_tumour + n_other).
High residual = MORE image persistence than expected for its size -> expected
SEGREGATED; low residual -> expected INTERSPERSED. Cell counts are comparable
across the two galleries, so any visual difference reflects spatial arrangement.
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
ROIS = HERE / "data" / "rois"
OUT = HERE / "output"
N = 8
FEAT = "im_dim1_total_persistence"

ap = argparse.ArgumentParser()
ap.add_argument("--definition", default="CD8_primary", choices=list(pair_defs.DEFINITIONS))
DEFN = ap.parse_args().definition
SPECIES = pair_defs.DEFINITIONS[DEFN]
SPECIES_LABEL = {"CD8_primary": "CD8", "CD8_strict": "CD8", "CD8_extended": "CD8",
                 "CD8_NK_only": "NK/CD8", "CD4": "CD4", "Bcell": "B cell",
                 "Macrophage": "Macrophage", "Fibroblast": "Fibroblast"}.get(DEFN, DEFN)
SUFFIX = "" if DEFN == "CD8_primary" else f"_{DEFN}"

df = pd.read_parquet(OUT / "per_roi_summary.parquet")
d = df[(df.pair_definition == DEFN) & (df.status == "ok")].copy()
d["n_pair"] = d["n_tumour"] + d["n_other"]

# size-adjusted residual (OLS on n_tumour + n_other, fit on all ok ROIs)
X = np.column_stack([np.ones(len(d)), d["n_tumour"], d["n_other"]])
beta, *_ = np.linalg.lstsq(X, d[FEAT].to_numpy(), rcond=None)
d["resid"] = d[FEAT].to_numpy() - X @ beta

# restrict to mid-size band so the two galleries have comparable cell counts
lo, hi = d["n_pair"].quantile([0.25, 0.75])
band = d[(d["n_pair"] >= lo) & (d["n_pair"] <= hi)].sort_values("resid")
print(f"mid-size band: {lo:.0f}-{hi:.0f} cells, {len(band)} ROIs")

low = band.head(N)            # most interspersed-for-their-size
high = band.tail(N).iloc[::-1]  # most segregated-for-their-size


def plot_gallery(sel, title, fname):
    ncol = 4
    nrow = int(np.ceil(len(sel) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, row in zip(axes, sel.itertuples(index=False)):
        pc = pd.read_csv(ROIS / f"{row.roi_id}.csv")
        t = pc[pc["celltype"].isin(pair_defs.TUMOUR)]
        c = pc[pc["celltype"].isin(SPECIES)]
        ax.scatter(t["x"], t["y"], s=5, c="tab:red", alpha=0.6, label=f"Tumour {len(t)}")
        ax.scatter(c["x"], c["y"], s=5, c="tab:blue", alpha=0.6, label=f"{SPECIES_LABEL} {len(c)}")
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        resp = {"Responder": "R", "Non-Responder": "NR"}.get(row.response, "?")
        ax.set_title(f"{row.roi_id} [{resp}]  n={row.n_pair}\n"
                     f"im_pers={getattr(row, FEAT):.0f}  resid={row.resid:+.0f}", fontsize=8)
        ax.legend(fontsize=6, loc="upper right", markerscale=2)
    for ax in axes[len(sel):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")


plot_gallery(high, f"{DEFN}: HIGH size-adj residual (top {N}, mid-size) — expected SEGREGATED "
             f"(red=Tumour, blue={SPECIES_LABEL})", f"gallery_resid_high{SUFFIX}.png")
plot_gallery(low, f"{DEFN}: LOW size-adj residual (bottom {N}, mid-size) — expected INTERSPERSED "
             f"(red=Tumour, blue={SPECIES_LABEL})", f"gallery_resid_low{SUFFIX}.png")

print("\nHIGH residual (segregated?) response mix:", high.response.value_counts().to_dict())
print("LOW residual (interspersed?) response mix:", low.response.value_counts().to_dict())
print("\nHIGH residual ROIs:\n", high[["roi_id", "response", "n_pair", FEAT, "resid"]].to_string(index=False))
print("\nLOW residual ROIs:\n", low[["roi_id", "response", "n_pair", FEAT, "resid"]].to_string(index=False))
