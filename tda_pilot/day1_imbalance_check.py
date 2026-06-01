#!/usr/bin/env python
"""Is high image persistence explained by class IMBALANCE rather than segregation?

For CD8_primary, plot im_dim1_total_persistence vs the minority fraction
  f = min(n_tumour, n_other) / (n_tumour + n_other)      (0 = fully imbalanced, 0.5 = balanced)
coloured by response, both raw and as the size-adjusted residual (the quantity behind the
-0.138 delta). Then test whether the responder-vs-NR persistence signal survives adding the
minority fraction as a covariate (if it collapses, the "signal" is imbalance).
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
FEAT = "im_dim1_total_persistence"
DEFN = "CD8_primary"

df = pd.read_parquet(OUT / "per_roi_summary.parquet")
d = df[(df.pair_definition == DEFN) & (df.status == "ok")
       & df.response.isin(["Responder", "Non-Responder"])].copy()
d["minority_frac"] = np.minimum(d.n_tumour, d.n_other) / (d.n_tumour + d.n_other)


def cliffs(r, nr):
    r, nr = np.asarray(r), np.asarray(nr)
    u = mannwhitneyu(r, nr).statistic
    return 2 * u / (len(r) * len(nr)) - 1


def resid(cols):
    X = np.column_stack([np.ones(len(d))] + [d[c].to_numpy(float) for c in cols])
    beta, *_ = np.linalg.lstsq(X, d[FEAT].to_numpy(), rcond=None)
    return d[FEAT].to_numpy() - X @ beta


d["resid_size"] = resid(["n_tumour", "n_other"])
r_mask = (d.response == "Responder").to_numpy()


def binned(x, y, q=10):
    bins = pd.qcut(x, q, duplicates="drop")
    g = pd.DataFrame({"x": x, "y": y, "b": bins}).groupby("b", observed=True)
    return g["x"].median(), g["y"].median()


fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
for col, (yv, lab) in zip(ax, [(d[FEAT], "raw im_dim1_total_persistence"),
                               (d["resid_size"], "size-adjusted residual")]):
    for grp, c in [("Non-Responder", "tab:blue"), ("Responder", "tab:orange")]:
        m = d.response == grp
        col.scatter(d.minority_frac[m], yv[m], s=10, c=c, alpha=0.4, label=grp)
        bx, by = binned(d.minority_frac[m].to_numpy(), yv[m].to_numpy())
        col.plot(bx, by, c=c, lw=2.5, marker="o", ms=4)
    rho, p = spearmanr(d.minority_frac, yv)
    col.set_xlabel("minority fraction  min(n_t,n_o)/(n_t+n_o)")
    col.set_ylabel(lab)
    col.set_title(f"{lab}\nSpearman ρ(minority_frac) = {rho:+.3f}")
    col.legend()
fig.suptitle(f"{DEFN}: image persistence vs class imbalance (lines = per-response binned medians)")
fig.tight_layout(); fig.savefig(OUT / "imbalance_check.png", dpi=120, bbox_inches="tight")

# ---- quantitative summary ----------------------------------------------------
print(f"minority_frac: median NR={d[~r_mask].minority_frac.median():.3f} "
      f"R={d[r_mask].minority_frac.median():.3f}  Cliff δ(R vs NR)={cliffs(d[r_mask].minority_frac, d[~r_mask].minority_frac):+.3f}")
print(f"Spearman ρ(persistence, minority_frac)        = {spearmanr(d.minority_frac, d[FEAT])[0]:+.3f}")
print(f"Spearman ρ(size-adj residual, minority_frac)  = {spearmanr(d.minority_frac, d.resid_size)[0]:+.3f}")

print("\nResponder-vs-NR Cliff δ on image persistence under different adjustments:")
for name, cols in [("raw (none)", []),
                   ("size (n_t + n_o)", ["n_tumour", "n_other"]),
                   ("imbalance only (minority_frac)", ["minority_frac"]),
                   ("size + imbalance", ["n_tumour", "n_other", "minority_frac"])]:
    rr = resid(cols) if cols else d[FEAT].to_numpy()
    print(f"  {name:34s} δ = {cliffs(rr[r_mask], rr[~r_mask]):+.3f}")
print("\nwrote output/imbalance_check.png")
