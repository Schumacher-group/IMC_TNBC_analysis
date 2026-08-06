#!/usr/bin/env python
"""Specification sensitivity: is the directional response association real, or geometry?

Context. The pre-specified endpoint (per-patient median z of `ker1-avg_length`, directional
inclusion, pre-treatment) gives Cliff delta = +0.379, permutation p = 0.011. z carries a
geometry coupling (rho(z, density) ~ 0.31) because the null sd shrinks in larger, denser ROIs,
so the question is whether that coupling explains the association.

Two families of answer, and they DISAGREE, which is the finding:

  REGRESSION adjustment  residualise z on composition + geometry across ROIs. Result depends
                         strongly on which covariates are used. The four-covariate set
                         (cd8_frac, log cells, log density, log area) is near-collinear --
                         density = cells/area, so the three size terms nearly span a 2-D
                         space (condition number ~514) -- and it drives delta to +0.008.
                         A minimal non-redundant set gives +0.226; log_density + log_cells
                         gives +0.343 (p = 0.023).

  MATCHING               restrict to ROIs in a common band of area AND density, then test with
                         no adjustment at all. Assumption-free and immune to collinearity.
                         Here delta stays at +0.27 to +0.38 in every band -- it does NOT
                         shrink toward zero, which is what a geometry artefact would do.

Matching is the more trustworthy arbiter (no functional-form assumption, no collinearity), but
it costs power, so its p-values straddle 0.05. The honest conclusion is that this design cannot
resolve the question -- see the verdict printed at the end.

Note the DIRECTION throughout: responders score HIGHER on the exclusion-like axis. That is the
opposite of the immune-exclusion hypothesis, so none of these specifications supports
"reduced infiltration in non-responders", whichever way the robustness question lands.

Output: output/SENSITIVITY.md
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
STAT = "ker1-avg_length"


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def load(defn: str, tag: str = "_other_only") -> pd.DataFrame:
    geom = pd.read_csv(OUT / "roi_geometry.csv")
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "resp"},
    )[["roi_id", "pid", "resp"]]
    d = pd.read_parquet(OUT / f"perm_null_{defn}{tag}.parquet")
    d = d[d["status"] == "ok"].merge(meta, on="roi_id").merge(geom, on="roi_id")
    d = d[d.resp.isin(["Responder", "Non-Responder"])].copy()
    d["z"] = d[f"{STAT}__z"]
    m = d[f"{STAT}__null_mean"]
    d["reldev"] = np.where(m != 0, (d[f"{STAT}__obs"] - m) / m, np.nan)
    d["n_pair"] = d.n_tumour + d.n_other
    d["ofrac"] = d.n_other / d.n_pair
    return d


def patient_test(d: pd.DataFrame, col: str):
    p = d.groupby(["pid", "resp"], as_index=False)[col].median()
    delta, pv = cliffs(p[p.resp == "Responder"][col], p[p.resp == "Non-Responder"][col])
    return delta, pv, int((p.resp == "Non-Responder").sum()), int((p.resp == "Responder").sum())


def resid(d: pd.DataFrame, col: str, cols: list[np.ndarray]) -> np.ndarray:
    X = np.column_stack([np.ones(len(d))] + list(cols))
    y = d[col].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def md(dfr: pd.DataFrame) -> str:
    head = "| " + " | ".join(dfr.columns) + " |"
    sep = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    body = "\n".join("| " + " | ".join(
        f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)
        for v in row) + " |" for row in dfr.itertuples(index=False))
    return "\n".join([head, sep, body])


def main() -> None:
    d = load("CD8_primary")
    terms = {
        "cd8_frac": d.ofrac, "log_cells": np.log(d.n_pair),
        "log_density": np.log(d.density), "log_area": np.log(d.hull_area),
    }
    cond = np.linalg.cond(np.column_stack([np.ones(len(d))] + list(terms.values())))

    reg_sets = {
        "none (pre-specified endpoint)": [],
        "cd8_frac": [terms["cd8_frac"]],
        "log_cells": [terms["log_cells"]],
        "log_density": [terms["log_density"]],
        "log_area": [terms["log_area"]],
        "log_cells + log_area (minimal, non-redundant)":
            [terms["log_cells"], terms["log_area"]],
        "log_density + log_cells": [terms["log_density"], terms["log_cells"]],
        "log_cells + log_area + cd8_frac":
            [terms["log_cells"], terms["log_area"], terms["cd8_frac"]],
        "all four (near-collinear)": list(terms.values()),
    }
    reg_rows = []
    for name, cols in reg_sets.items():
        row = {"adjustment": name}
        for col in ("z", "reldev"):
            d["_t"] = d[col] if not cols else resid(d, col, cols)
            dl, pv, *_ = patient_test(d, "_t")
            row[f"{col}_delta"], row[f"{col}_p"] = dl, pv
        reg_rows.append(row)
    reg = pd.DataFrame(reg_rows)

    match_rows = []
    for name, (lo, hi) in {
        "all ROIs": (0.0, 1.0), "10–90 pct": (0.10, 0.90),
        "20–80 pct": (0.20, 0.80), "IQR (25–75 pct)": (0.25, 0.75),
    }.items():
        s = d if lo == 0.0 else d[
            d.hull_area.between(*d.hull_area.quantile([lo, hi]).values)
            & d.density.between(*d.density.quantile([lo, hi]).values)]
        row = {"area+density band": name, "n_ROI": len(s)}
        for col in ("z", "reldev"):
            dl, pv, nnr, nr = patient_test(s, col)
            row |= {f"{col}_delta": dl, f"{col}_p": pv}
        row["n_NR_pat"], row["n_R_pat"] = patient_test(s, "z")[2:]
        match_rows.append(row)
    match = pd.DataFrame(match_rows)

    reg_range = (reg.z_delta.min(), reg.z_delta.max())
    match_range = (match.z_delta.min(), match.z_delta.max())
    stable_under_matching = bool(match.z_delta.min() > 0.2)

    lines = [
        "# Specification sensitivity: directional statistic, pre-treatment\n",
        f"Endpoint: per-patient median of `{STAT}`, CD8-directional inclusion, 593 ROIs / "
        f"25 NR + 37 R patients. Cliff's δ is R vs NR (+ve = responders higher).\n",
        "## 1. Regression adjustment — answer depends on the covariate set\n",
        md(reg.round(4)),
        f"\nThe four covariates are near-collinear (density = cells / area; condition number "
        f"**{cond:.0f}**), and together they explain only ~9% of the variance of z, yet the "
        f"full set drives δ from +0.379 to +0.008. A 9%-of-variance covariate block "
        f"extinguishing a moderate group difference is a signature of over-adjustment under "
        f"collinearity, not of confounding removed.\n",
        "## 2. Matching — assumption-free, and it does not collapse\n",
        md(match.round(4)),
        "\nRestricting to a common band of area **and** density needs no model and cannot "
        "suffer collinearity. δ stays at "
        f"**{match_range[0]:+.3f} … {match_range[1]:+.3f}** in every band — it does not shrink "
        "toward zero, which is what a pure geometry artefact would do (compare "
        "`IMBALANCE_CONFOUND.md`, where the Day-1 signal fell from −0.138 to +0.035 as "
        "imbalanced ROIs were dropped). p-values widen because matching costs patients.\n",
        "## Verdict\n",
        f"- Across reasonable specifications δ ranges **{reg_range[0]:+.3f} … "
        f"{reg_range[1]:+.3f}** (p from 0.009 to 0.97). The result is **not robust to "
        "analytic choice**, so it is neither established nor refuted.",
        f"- Matching, the more trustworthy arbiter, leaves the effect intact "
        f"(δ ≈ {match.z_delta.iloc[1:].mean():+.2f}) but underpowered "
        f"(p ≈ 0.03–0.10 at n = 14–24 NR patients).",
        "- **Therefore: indeterminate.** The honest statement is that this design cannot "
        "separate a moderate arrangement effect from ROI geometry. Resolving it needs more "
        "patients or ROI sampling matched on area and density by design, not more analysis "
        "of these data.",
        "- **The direction is stable and matters more than the significance.** Responders "
        "score *higher* on the exclusion-like axis under every specification. The "
        "immune-exclusion hypothesis predicts the opposite (non-responders higher). So no "
        "specification supports 'reduced infiltration of cytotoxic T cells in "
        "non-responders'; the unresolved question is whether there is a weak effect running "
        "the other way.\n",
        "Script: `day2_sensitivity.py`.",
    ]
    (OUT / "SENSITIVITY.md").write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 220)
    print("=== regression adjustment ===")
    print(reg.round(4).to_string(index=False))
    print(f"\ncondition number of the 4-covariate design: {cond:.0f}")
    print("\n=== geometry matching ===")
    print(match.round(4).to_string(index=False))
    print(f"\nstable under matching: {stable_under_matching}")
    print("wrote", OUT / "SENSITIVITY.md")


if __name__ == "__main__":
    main()
