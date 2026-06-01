#!/usr/bin/env python
"""Day-1 descriptive diagnostics (NO models).

(2) Scan all 7 pair definitions: R-vs-NR separation per focal feature, to check
    (a) whether CD8_extended (adds 134k-cell NK/CD8) differs from CD8_primary, and
    (b) whether any control pair separates MORE than CD8_primary (red flag: signal
    would not be CD8-specific).
(3) Pressure-test the image-persistence hint: does responders-lower im_dim1
    total persistence survive a per-ROI size adjustment (linear regression of the
    feature on cell counts), or was it a cell-count confound?

Effect size = Cliff's delta (descriptive rank statistic; +ve = responders higher).
Outputs: output/diag_separation.png, output/diag_image_residual.png,
         output/DAY1_DIAGNOSTICS.md
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
df = pd.read_parquet(HERE / "output" / "per_roi_summary.parquet")
df = df[(df["status"] == "ok") & df["response"].isin(["Responder", "Non-Responder"])].copy()
DEFS = ["CD8_primary", "CD8_strict", "CD8_extended", "CD8_NK_only",
        "CD4", "Bcell", "Macrophage", "Fibroblast"]


def cliffs_delta(r, nr):
    r, nr = np.asarray(r), np.asarray(nr)
    if len(r) == 0 or len(nr) == 0:
        return np.nan
    u = mannwhitneyu(r, nr, alternative="two-sided").statistic
    return 2 * u / (len(r) * len(nr)) - 1  # +ve => responders tend higher


def sep_row(d, feat):
    r = d[d.response == "Responder"][feat].dropna()
    nr = d[d.response == "Non-Responder"][feat].dropna()
    return {"NR_med": nr.median(), "R_med": r.median(),
            "ratio": (r.median() / nr.median()) if nr.median() else np.nan,
            "cliffs_delta": cliffs_delta(r, nr)}


# ---- (2) cross-definition separation table -----------------------------------
FEATS = ["ker_dim1_total_persistence", "im_dim1_total_persistence",
         "cok_dim1_total_persistence", "im_dim1_n_features", "ker_dim1_n_features"]
rows = []
for dfn in DEFS:
    d = df[df.pair_definition == dfn]
    for f in FEATS:
        rows.append({"definition": dfn, "feature": f, **sep_row(d, f)})
sep = pd.DataFrame(rows)
sep.to_csv(HERE / "output" / "diag_separation.csv", index=False)

# stripplot grid: headline feature im_dim1_total_persistence across 7 definitions
rng = np.random.default_rng(0)
fig, axes = plt.subplots(1, len(DEFS), figsize=(2.9 * len(DEFS), 4.2), sharey=True)
for ax, dfn in zip(axes, DEFS):
    d = df[df.pair_definition == dfn]
    for i, grp in enumerate(["Non-Responder", "Responder"]):
        y = d[d.response == grp]["im_dim1_total_persistence"].dropna().to_numpy()
        ax.scatter(rng.normal(i, 0.07, len(y)), y, s=7, alpha=0.35)
        ax.hlines(np.median(y), i - 0.25, i + 0.25, color="k", lw=2)
    dlt = sep[(sep.definition == dfn) & (sep.feature == "im_dim1_total_persistence")]["cliffs_delta"].iloc[0]
    ax.set_xticks([0, 1]); ax.set_xticklabels(["NR", "R"])
    ax.set_title(f"{dfn}\nCliff δ={dlt:+.3f}", fontsize=9)
axes[0].set_ylabel("im_dim1_total_persistence")
fig.suptitle("Image degree-1 total persistence, R vs NR, across 7 pair definitions (per-ROI)")
fig.tight_layout()
fig.savefig(HERE / "output" / "diag_separation.png", dpi=120, bbox_inches="tight")

# ---- (3) size-adjusted image-persistence residual ----------------------------
def size_adjust(d, feat):
    """Residual of feat after OLS on [1, n_tumour, n_other]."""
    sub = d.dropna(subset=[feat, "n_tumour", "n_other"])
    X = np.column_stack([np.ones(len(sub)), sub["n_tumour"], sub["n_other"]])
    y = sub[feat].to_numpy()
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    out = sub[["response"]].copy()
    out["resid"] = resid
    return out, beta

resid_rows = []
RESID_DEFS = ["CD8_primary", "CD8_strict", "CD8_extended", "CD8_NK_only"]
fig2, axes2 = plt.subplots(1, len(RESID_DEFS), figsize=(4.3 * len(RESID_DEFS), 4.5))
for ax, dfn in zip(axes2, RESID_DEFS):
    d = df[df.pair_definition == dfn]
    radj, beta = size_adjust(d, "im_dim1_total_persistence")
    r = radj[radj.response == "Responder"]["resid"]; nr = radj[radj.response == "Non-Responder"]["resid"]
    raw = sep[(sep.definition == dfn) & (sep.feature == "im_dim1_total_persistence")]["cliffs_delta"].iloc[0]
    adj = cliffs_delta(r, nr)
    resid_rows.append({"definition": dfn, "raw_cliffs_delta": raw, "size_adj_cliffs_delta": adj,
                       "resid_med_NR": nr.median(), "resid_med_R": r.median()})
    for i, (grp, yv) in enumerate([("Non-Responder", nr), ("Responder", r)]):
        ax.scatter(rng.normal(i, 0.07, len(yv)), yv, s=7, alpha=0.35)
        ax.hlines(np.median(yv), i - 0.25, i + 0.25, color="k", lw=2)
    ax.axhline(0, color="grey", lw=0.6, ls=":")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["NR", "R"])
    ax.set_title(f"{dfn}\nraw δ={raw:+.3f} -> size-adj δ={adj:+.3f}", fontsize=9)
axes2[0].set_ylabel("im_dim1_total_persistence residual\n(after OLS on n_tumour + n_other)")
fig2.suptitle("Image total persistence — does the responders-lower hint survive size adjustment?")
fig2.tight_layout()
fig2.savefig(HERE / "output" / "diag_image_residual.png", dpi=120, bbox_inches="tight")
resid = pd.DataFrame(resid_rows)

# ---- write markdown ----------------------------------------------------------
def md(dfr):
    h = "| " + " | ".join(dfr.columns) + " |"
    s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    b = "\n".join("| " + " | ".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in row) + " |"
                  for row in dfr.itertuples(index=False))
    return "\n".join([h, s, b])

# rank definitions by |delta| for the two total-persistence features
piv = sep.pivot(index="definition", columns="feature", values="cliffs_delta").reindex(DEFS)
lines = ["# Day-1 descriptive diagnostics (no models)\n",
         "Effect size = Cliff's delta on per-ROI values (+ve = responders higher; "
         "|δ|≈0.1 small, 0.3 medium, 0.5 large).\n",
         "## (2) Separation across all 7 definitions (Cliff's δ, R vs NR)\n",
         md(piv.reset_index().round(3)),
         "\n## Per-definition x feature detail\n", md(sep.round(3)),
         "\n## (3) Image total persistence: raw vs size-adjusted (CD8 definitions)\n", md(resid.round(3)),
         "\nFigures: `output/diag_separation.png`, `output/diag_image_residual.png`.\n"]
(HERE / "output" / "DAY1_DIAGNOSTICS.md").write_text("\n".join(lines) + "\n")

# compact stdout
print("== Cliff's delta (R vs NR), per definition ==")
print(piv.round(3).to_string())
print("\n== image persistence raw vs size-adjusted ==")
print(resid.round(3).to_string(index=False))
print("\nmax |delta| im_dim1_total_persistence by def:")
print(sep[sep.feature == "im_dim1_total_persistence"].set_index("definition")["cliffs_delta"].round(3).to_string())
