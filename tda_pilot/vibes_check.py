#!/usr/bin/env python
"""End-of-Day-1 qualitative check (NOT a model).

Stripplots of ker_dim1_total_persistence for CD8_primary, responder vs
non-responder: (left) one point per ROI, (right) one point per patient (mean).
Saves output/vibes_ker_dim1.png. Prints group medians for a quick read.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
df = pd.read_parquet(HERE / "output" / "per_roi_summary.parquet")
d = df[(df["pair_definition"] == "CD8_primary") & (df["status"] == "ok")].copy()
d = d[d["response"].isin(["Responder", "Non-Responder"])]
metric = "ker_dim1_total_persistence"

order = ["Non-Responder", "Responder"]
rng = np.random.default_rng(0)
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=False)

# left: per-ROI
ax = axes[0]
for i, grp in enumerate(order):
    y = d[d["response"] == grp][metric].to_numpy()
    ax.scatter(rng.normal(i, 0.06, len(y)), y, s=10, alpha=0.4)
    ax.hlines(np.median(y), i - 0.25, i + 0.25, color="k", lw=2)
ax.set_xticks([0, 1]); ax.set_xticklabels(order)
ax.set_ylabel(metric); ax.set_title(f"per-ROI (n={len(d)})")

# right: per-patient mean
pat = d.groupby(["patient_id", "response"])[metric].mean().reset_index()
ax = axes[1]
for i, grp in enumerate(order):
    y = pat[pat["response"] == grp][metric].to_numpy()
    ax.scatter(rng.normal(i, 0.06, len(y)), y, s=30, alpha=0.7)
    ax.hlines(np.median(y), i - 0.25, i + 0.25, color="k", lw=2)
ax.set_xticks([0, 1]); ax.set_xticklabels(order)
ax.set_ylabel(f"patient-mean {metric}"); ax.set_title(f"per-patient (n={pat['patient_id'].nunique()})")

fig.suptitle("CD8_primary kernel degree-1 total persistence — vibes check (no model)")
fig.tight_layout()
fig.savefig(HERE / "output" / "vibes_ker_dim1.png", dpi=120, bbox_inches="tight")

for grp in order:
    roi_med = d[d["response"] == grp][metric].median()
    pat_med = pat[pat["response"] == grp][metric].median()
    print(f"{grp:14s}  per-ROI median={roi_med:10.1f}   per-patient median={pat_med:10.1f}")
print("wrote output/vibes_ker_dim1.png")
