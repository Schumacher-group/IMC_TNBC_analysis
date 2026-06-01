#!/usr/bin/env python
"""Histogram of CD8_primary image degree-1 total persistence by response, plus
point-cloud galleries for the highest- and lowest-persistence ROIs.

Hypothesis: high image total persistence -> tumour/CD8 spatially SEGREGATED;
low image total persistence -> INTERSPERSED (more infiltration).
Point clouds use the same style as the Leap066_11 sanity overlay.
"""

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
N = 8  # top/bottom count
CD8 = pair_defs.DEFINITIONS["CD8_primary"]
FEAT = "im_dim1_total_persistence"

df = pd.read_parquet(OUT / "per_roi_summary.parquet")
d = df[(df.pair_definition == "CD8_primary") & (df.status == "ok")].copy()

# ---- histogram coloured by response -----------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
bins = np.linspace(0, d[FEAT].quantile(0.99), 40)
colours = {"Non-Responder": "tab:blue", "Responder": "tab:orange"}
for grp, c in colours.items():
    v = d[d.response == grp][FEAT]
    axes[0].hist(v, bins=bins, alpha=0.5, label=f"{grp} (n={len(v)})", color=c)
    axes[0].axvline(v.median(), color=c, lw=2, ls="--")
    axes[1].hist(v, bins=bins, density=True, histtype="step", lw=2, label=grp, color=c)
axes[0].set_title("counts"); axes[1].set_title("density (group-normalised)")
for ax in axes:
    ax.set_xlabel("im_dim1_total_persistence"); ax.legend()
fig.suptitle("CD8_primary image degree-1 total persistence by response (808 ROIs)")
fig.tight_layout(); fig.savefig(OUT / "hist_im_persistence.png", dpi=120, bbox_inches="tight")

# ---- select highest / lowest ROIs -------------------------------------------
d = d.sort_values(FEAT)
low = d.head(N)   # lowest persistence -> expected interspersed
high = d.tail(N).iloc[::-1]  # highest persistence -> expected segregated


def plot_gallery(sel: pd.DataFrame, title: str, fname: str) -> None:
    ncol = 4
    nrow = int(np.ceil(len(sel) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, row in zip(axes, sel.itertuples(index=False)):
        pc = pd.read_csv(ROIS / f"{row.roi_id}.csv")
        is_t = pc["celltype"].isin(pair_defs.TUMOUR)
        is_c = pc["celltype"].isin(CD8)
        t = pc[is_t]; c = pc[is_c]
        ax.scatter(t["x"], t["y"], s=4, c="tab:red", alpha=0.6, label=f"Tumour {len(t)}")
        ax.scatter(c["x"], c["y"], s=4, c="tab:blue", alpha=0.6, label=f"CD8 {len(c)}")
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        resp = {"Responder": "R", "Non-Responder": "NR"}.get(row.response, "?")
        ax.set_title(f"{row.roi_id} [{resp}]\nim_pers={getattr(row, FEAT):.0f}", fontsize=9)
        ax.legend(fontsize=6, loc="upper right", markerscale=2)
    for ax in axes[len(sel):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")


plot_gallery(high, f"HIGHEST image persistence (top {N}) — expected SEGREGATED",
             "gallery_high_persistence.png")
plot_gallery(low, f"LOWEST image persistence (bottom {N}) — expected INTERSPERSED",
             "gallery_low_persistence.png")

print("highest im_persistence ROIs:")
print(high[["roi_id", "response", "n_tumour", "n_other", FEAT]].to_string(index=False))
print("\nlowest im_persistence ROIs:")
print(low[["roi_id", "response", "n_tumour", "n_other", FEAT]].to_string(index=False))
print("\nwrote hist_im_persistence.png, gallery_high_persistence.png, gallery_low_persistence.png")
