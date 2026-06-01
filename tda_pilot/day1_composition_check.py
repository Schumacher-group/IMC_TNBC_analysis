#!/usr/bin/env python
"""Is the responder difference simply compositional? Per-ROI CD8:tumour ratio.

Restricted to PRE-treatment ROIs (resection/post composition reflects chemo, not a
predictor). Shows: (a) per-ROI histogram of log2(CD8/tumour) by response, and
(b) per-patient median log2 ratio by response (the honest unit; multiple ROIs/patient).
Tests: Mann-Whitney + Cliff's delta at ROI and patient level.
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

df = pd.read_parquet(OUT / "per_roi_summary.parquet")
d = df[(df.pair_definition == "CD8_primary") & (df.status == "ok")
       & df.response.isin(["Responder", "Non-Responder"])].copy()

# sample type (pre/post) from per_roi_counts
meta = pd.read_csv(HERE / "per_roi_counts.csv")[["fov", "Sample_Type_(pre/post treatment)"]]
meta = meta.rename(columns={"fov": "roi_id", "Sample_Type_(pre/post treatment)": "sample_type"})
d = d.merge(meta, on="roi_id", how="left")
print("sample_type values:", d.sample_type.value_counts(dropna=False).to_dict())

pre = d[d.sample_type.astype(str).str.lower() == "pre"].copy()
pre["ratio"] = pre.n_other / pre.n_tumour            # CD8 : tumour
pre["log2ratio"] = np.log2(pre.ratio.replace(0, np.nan))
pre = pre.dropna(subset=["log2ratio"])
print(f"pre-treatment ROIs: {len(pre)}  "
      f"(NR={int((pre.response=='Non-Responder').sum())}, R={int((pre.response=='Responder').sum())})")


def cliffs(r, nr):
    r, nr = np.asarray(r), np.asarray(nr)
    u = mannwhitneyu(r, nr).statistic
    return 2 * u / (len(r) * len(nr)) - 1, mannwhitneyu(r, nr, alternative="two-sided").pvalue


cols = {"Non-Responder": "tab:blue", "Responder": "tab:orange"}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))

# (a) per-ROI histogram of log2 ratio
bins = np.linspace(pre.log2ratio.min(), pre.log2ratio.max(), 35)
for grp, c in cols.items():
    v = pre[pre.response == grp]["log2ratio"]
    ax[0].hist(v, bins=bins, alpha=0.5, color=c, label=f"{grp} (n={len(v)})", density=True)
    ax[0].axvline(v.median(), color=c, lw=2, ls="--")
ax[0].set_xlabel("log2(CD8 / tumour)  per ROI"); ax[0].set_ylabel("density")
d_roi, p_roi = cliffs(pre[pre.response == "Responder"].log2ratio,
                      pre[pre.response == "Non-Responder"].log2ratio)
ax[0].set_title(f"(a) per-ROI  (Cliff δ={d_roi:+.3f}, MWU p={p_roi:.1e})")
ax[0].legend()

# (b) per-patient median log2 ratio (one point per patient)
pat = pre.groupby(["patient_id", "response"], as_index=False)["log2ratio"].median()
rng = np.random.default_rng(0)
for i, grp in enumerate(["Non-Responder", "Responder"]):
    v = pat[pat.response == grp]["log2ratio"]
    ax[1].scatter(rng.normal(i, 0.06, len(v)), v, s=28, c=cols[grp], alpha=0.7,
                  label=f"{grp} (n={len(v)} patients)")
    ax[1].hlines(v.median(), i - 0.25, i + 0.25, color="k", lw=2)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["NR", "R"])
ax[1].set_ylabel("patient-median log2(CD8 / tumour)")
d_pat, p_pat = cliffs(pat[pat.response == "Responder"].log2ratio,
                      pat[pat.response == "Non-Responder"].log2ratio)
ax[1].set_title(f"(b) per-patient  (Cliff δ={d_pat:+.3f}, MWU p={p_pat:.3f})")
ax[1].legend()
fig.suptitle("Pre-treatment CD8:tumour composition by response")
fig.tight_layout(); fig.savefig(OUT / "composition_ratio.png", dpi=120, bbox_inches="tight")

print("\n-- per ROI --")
print(f"  median CD8:tumour  NR={2**pre[pre.response=='Non-Responder'].log2ratio.median():.3f}  "
      f"R={2**pre[pre.response=='Responder'].log2ratio.median():.3f}")
print(f"  Cliff δ={d_roi:+.3f}  MWU p={p_roi:.2e}")
print("-- per patient --")
print(f"  median  NR={2**pat[pat.response=='Non-Responder'].log2ratio.median():.3f}  "
      f"R={2**pat[pat.response=='Responder'].log2ratio.median():.3f}")
print(f"  Cliff δ={d_pat:+.3f}  MWU p={p_pat:.3f}  "
      f"(n_pat NR={int((pat.response=='Non-Responder').sum())}, R={int((pat.response=='Responder').sum())})")
print("\nwrote output/composition_ratio.png")
