#!/usr/bin/env python
"""Sanity-check the size adjustment behind CD8_NK_only's -0.132 size-adjusted delta.

The raw->size-adjusted swing for CD8_NK_only (-0.026 -> -0.132) is much larger than
for the other CD8 definitions, so we verify the OLS adjustment (im_dim1_total_persistence
~ n_tumour + n_other) is well-behaved:
  - marginal scatters vs n_tumour and vs n_NK/CD8 with OLS lines,
  - observed-vs-fitted and residual-vs-fitted (curvature / heteroscedasticity),
  - leverage (hat values) + Cook's distance to flag influential ROIs,
  - robustness of the size-adjusted Cliff's delta under (i) OLS, (ii) OLS after
    dropping high-leverage ROIs, (iii) Huber robust regression, (iv) log-log fit.
Reference: CD8_primary (known well-behaved) run through the same checks.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.linear_model import HuberRegressor

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
FEAT = "im_dim1_total_persistence"
df = pd.read_parquet(OUT / "per_roi_summary.parquet")


def cliffs(r, nr):
    r, nr = np.asarray(r), np.asarray(nr)
    u = mannwhitneyu(r, nr).statistic
    return 2 * u / (len(r) * len(nr)) - 1


def get(defn):
    d = df[(df.pair_definition == defn) & (df.status == "ok")
           & df.response.isin(["Responder", "Non-Responder"])].copy()
    return d


def ols_resid(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta, beta


def delta_from_resid(d, resid):
    r = resid[(d.response == "Responder").to_numpy()]
    nr = resid[(d.response == "Non-Responder").to_numpy()]
    return cliffs(r, nr)


def robustness(defn):
    d = get(defn)
    y = d[FEAT].to_numpy()
    nt, no = d["n_tumour"].to_numpy(float), d["n_other"].to_numpy(float)
    X = np.column_stack([np.ones(len(d)), nt, no])
    # (i) OLS
    resid, beta = ols_resid(X, y)
    d_ols = delta_from_resid(d, resid)
    # leverage + Cook's distance
    H = X @ np.linalg.pinv(X.T @ X) @ X.T
    hat = np.diag(H)
    p = X.shape[1]
    mse = (resid ** 2).sum() / (len(d) - p)
    cook = resid ** 2 / (p * mse) * hat / (1 - hat) ** 2
    keep = cook < (4 / len(d))            # standard Cook's cutoff
    keepl = hat < (3 * p / len(d))        # high-leverage cutoff
    # (ii) OLS after dropping high-influence ROIs
    m = keep & keepl
    resid2, _ = ols_resid(X[m], y[m])
    d_trim = delta_from_resid(d[m], resid2)
    # (iii) Huber robust regression
    hub = HuberRegressor().fit(np.column_stack([nt, no]), y)
    resid_h = y - hub.predict(np.column_stack([nt, no]))
    d_hub = delta_from_resid(d, resid_h)
    # (iv) log-log (guard against zeros)
    pos = (y > 0)
    dl = d[pos]
    Xl = np.column_stack([np.ones(pos.sum()), np.log(nt[pos] + 1), np.log(no[pos] + 1)])
    resid_l, _ = ols_resid(Xl, np.log(y[pos]))
    d_log = delta_from_resid(dl, resid_l)
    return {"definition": defn, "n": len(d), "raw_delta": cliffs(
                d[d.response == "Responder"][FEAT], d[d.response == "Non-Responder"][FEAT]),
            "ols": d_ols, "ols_drop_influential": d_trim, "n_dropped": int((~m).sum()),
            "huber": d_hub, "loglog": d_log,
            "max_cook": float(cook.max()), "max_hat": float(hat.max()),
            "beta_n_tumour": float(beta[1]), "beta_n_other": float(beta[2])}, d, resid, hat, cook, beta


# ---- diagnostic figure for CD8_NK_only --------------------------------------
res_nk, d, resid, hat, cook, beta = robustness("CD8_NK_only")
fitted = d[FEAT].to_numpy() - resid
colors = np.where(d.response == "Responder", "tab:orange", "tab:blue")
infl = cook >= (4 / len(d))

fig, ax = plt.subplots(2, 2, figsize=(13, 10))
# (a) vs n_NK/CD8
xo = d["n_other"].to_numpy()
ax[0, 0].scatter(xo, d[FEAT], s=10, c=colors, alpha=0.5)
b1 = np.polyfit(xo, d[FEAT], 1)
xs = np.linspace(xo.min(), xo.max(), 50)
ax[0, 0].plot(xs, np.polyval(b1, xs), "k-", lw=2, label=f"OLS slope={b1[0]:.2f}")
ax[0, 0].set_xlabel("n_NK/CD8 (n_other)"); ax[0, 0].set_ylabel(FEAT); ax[0, 0].legend()
ax[0, 0].set_title("(a) image persistence vs NK/CD8 count")
# (b) vs n_tumour
xt = d["n_tumour"].to_numpy()
ax[0, 1].scatter(xt, d[FEAT], s=10, c=colors, alpha=0.5)
b2 = np.polyfit(xt, d[FEAT], 1)
xs2 = np.linspace(xt.min(), xt.max(), 50)
ax[0, 1].plot(xs2, np.polyval(b2, xs2), "k-", lw=2, label=f"OLS slope={b2[0]:.2f}")
ax[0, 1].set_xlabel("n_tumour"); ax[0, 1].set_ylabel(FEAT); ax[0, 1].legend()
ax[0, 1].set_title("(b) image persistence vs tumour count")
# (c) observed vs fitted (the actual 2-predictor OLS), highlight influential
ax[1, 0].scatter(fitted, d[FEAT], s=10, c=colors, alpha=0.5)
ax[1, 0].scatter(fitted[infl], d[FEAT].to_numpy()[infl], s=60, facecolors="none",
                 edgecolors="red", label=f"Cook's D > 4/n ({infl.sum()})")
lim = [min(fitted.min(), d[FEAT].min()), max(fitted.max(), d[FEAT].max())]
ax[1, 0].plot(lim, lim, "k--", lw=1)
ax[1, 0].set_xlabel("OLS fitted (n_tumour + n_other)"); ax[1, 0].set_ylabel("observed " + FEAT)
ax[1, 0].legend(); ax[1, 0].set_title("(c) observed vs fitted")
# (d) residual vs fitted
ax[1, 1].scatter(fitted, resid, s=10, c=colors, alpha=0.5)
ax[1, 1].axhline(0, color="grey", ls=":")
ax[1, 1].set_xlabel("OLS fitted"); ax[1, 1].set_ylabel("residual")
ax[1, 1].set_title("(d) residual vs fitted (curvature / heteroscedasticity)")
fig.suptitle("CD8_NK_only — size-adjustment regression diagnostics "
             "(orange=Responder, blue=Non-Responder)")
fig.tight_layout(); fig.savefig(OUT / "nk_ols_diagnostics.png", dpi=120, bbox_inches="tight")

# ---- robustness table: NK vs primary ----------------------------------------
rows = [res_nk, robustness("CD8_primary")[0]]
tab = pd.DataFrame(rows)[["definition", "n", "raw_delta", "ols", "ols_drop_influential",
                          "n_dropped", "huber", "loglog", "max_cook", "max_hat",
                          "beta_n_tumour", "beta_n_other"]]
pd.set_option("display.width", 200, "display.max_columns", 30)
print(tab.round(3).to_string(index=False))
print(f"\nCook's cutoff 4/n = {4/res_nk['n']:.4f}; high-leverage cutoff 3p/n = {9/res_nk['n']:.4f}")
print("wrote output/nk_ols_diagnostics.png")
