#!/usr/bin/env python
"""Does the CD8_primary responder signal survive when we DROP imbalanced ROIs?

Subset (not regression) test: restrict to ROIs whose minority fraction
  f = min(n_tumour, n_other)/(n_tumour+n_other)
is >= threshold, then recompute the responder-vs-NR Cliff's delta on
im_dim1_total_persistence, both raw and size-adjusted (OLS on n_tumour + n_other,
refit WITHIN each subset). Reports how the effect moves with threshold.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
FEAT = "im_dim1_total_persistence"

df = pd.read_parquet(OUT / "per_roi_summary.parquet")
d0 = df[(df.pair_definition == "CD8_primary") & (df.status == "ok")
        & df.response.isin(["Responder", "Non-Responder"])].copy()
d0["minority_frac"] = np.minimum(d0.n_tumour, d0.n_other) / (d0.n_tumour + d0.n_other)


def cliffs(r, nr):
    if len(r) == 0 or len(nr) == 0:
        return np.nan
    u = mannwhitneyu(r, nr).statistic
    return 2 * u / (len(r) * len(nr)) - 1


def size_adj_resid(d):
    X = np.column_stack([np.ones(len(d)), d.n_tumour, d.n_other])
    beta, *_ = np.linalg.lstsq(X, d[FEAT].to_numpy(), rcond=None)
    return d[FEAT].to_numpy() - X @ beta


rows = []
for thr in [0.00, 0.10, 0.15, 0.20, 0.25, 0.30]:
    d = d0[d0.minority_frac >= thr].copy()
    rmask = (d.response == "Responder").to_numpy()
    raw = cliffs(d[FEAT].to_numpy()[rmask], d[FEAT].to_numpy()[~rmask])
    res = size_adj_resid(d)
    adj = cliffs(res[rmask], res[~rmask])
    mf = cliffs(d.minority_frac.to_numpy()[rmask], d.minority_frac.to_numpy()[~rmask])
    rows.append({"min_frac>=": thr, "n_NR": int((~rmask).sum()), "n_R": int(rmask.sum()),
                 "raw_delta": raw, "size_adj_delta": adj, "remaining_minorityfrac_delta": mf})

tab = pd.DataFrame(rows)
print(tab.round(3).to_string(index=False))

fig, ax = plt.subplots(figsize=(7, 4.8))
ax.plot(tab["min_frac>="], tab["raw_delta"], "o-", label="raw δ")
ax.plot(tab["min_frac>="], tab["size_adj_delta"], "s-", label="size-adjusted δ")
ax.axhline(0, color="grey", lw=0.6, ls=":")
for _, r in tab.iterrows():
    ax.annotate(f"n={r.n_NR+r.n_R:.0f}", (r["min_frac>="], r["size_adj_delta"]),
                textcoords="offset points", xytext=(0, -12), fontsize=8, ha="center")
ax.set_xlabel("minority-fraction threshold (keep ROIs with min(n_t,n_o)/(n_t+n_o) >= x)")
ax.set_ylabel("Responder vs NR Cliff's δ  (im_dim1_total_persistence)")
ax.set_title("CD8_primary: responder signal vs balance threshold")
ax.legend()
fig.tight_layout(); fig.savefig(OUT / "threshold_check.png", dpi=120, bbox_inches="tight")
print("\nwrote output/threshold_check.png")
