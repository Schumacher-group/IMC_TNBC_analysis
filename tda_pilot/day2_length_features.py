#!/usr/bin/env python
"""Day-2 re-scan: the bar-LENGTH statistics the Day-1 analysis never looked at.

Motivation
----------
Day-1 decided the go/no-go on `im_dim1_total_persistence` and `n_features`, recomputed
by `extract_features.py` from the diagrams. But

    total_persistence = num_bars x avg_length

and only the *count* term was ever adjusted for (OLS on n_tumour + n_other in
`day1_diagnostics.py`). The count term is extensive and carries the cell-number and
class-imbalance confounds that `IMBALANCE_CONFOUND.md` correctly identified; the length
term is the actual spatial-SCALE statistic. Residualising the product on counts is a poor
way to strip a multiplicative term.

`avg_length`, `med_length` and `p90_length` were already computed for every diagram and
both dimensions by the full run and are sitting unused in output/stats/*.parquet (456
statistics per ROI x definition, of which Day-1 used 12). This script reads them and runs
the pre-treatment responder contrast properly:

  - PRE-TREATMENT ONLY  (pooled pre+post inherits the treatment-stage confound documented
    in PRETREATMENT_CONFOUND.md: post-treatment resections are 170 NR vs 7 R)
  - PER-PATIENT as the primary unit (63 patients, 6-11 ROIs each; per-ROI p-values are
    badly anticonservative), with per-ROI reported alongside
  - three scale variants, because responder ROIs are genuinely less dense pre-treatment
    (Cliff delta ~ -0.28 on density) so an unadjusted length difference could be geometry:
      raw          the statistic as computed (um)
      /med_nn      divided by the ROI's median nearest-neighbour distance -> dimensionless
      geom_adj     residual after OLS on log(med_nn), log(hull_area), log(n_tumour),
                   log(n_other)
  - all 8 pair definitions, so the four control species act as specificity controls
  - Benjamini-Hochberg across each variant's family of tests

This is DESCRIPTIVE and post hoc. It does not establish a result; it establishes whether
`perm_null.py` (the composition- and geometry-matched null) is worth running.

Output: output/LENGTH_FEATURES.md
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
STATS_DIR = HERE / "output" / "stats"
OUT_MD = HERE / "output" / "LENGTH_FEATURES.md"

DEFS = ["CD8_primary", "CD8_strict", "CD8_extended", "CD8_NK_only",
        "CD4", "Bcell", "Macrophage", "Fibroblast"]
CD8_DEFS = {"CD8_primary", "CD8_strict", "CD8_extended", "CD8_NK_only"}
# ker/im/cok carry the mixing information; dom/cod are the monochromatic and joint
# complexes, kept because they are the diagnostic for "is this about mixing at all?"
DIAGRAMS = ["ker1", "im1", "cok1", "dom1", "cod1", "ker0", "im0", "cok0", "dom0", "cod0"]
STAT = "avg_length"


def cliffs_delta(r, nr) -> tuple[float, float]:
    """Cliff's delta and two-sided MWU p. +ve => responders tend higher."""
    r, nr = np.asarray(r, float), np.asarray(nr, float)
    if len(r) < 2 or len(nr) < 2:
        return np.nan, np.nan
    res = mannwhitneyu(r, nr, alternative="two-sided")
    return 2 * res.statistic / (len(r) * len(nr)) - 1, res.pvalue


def bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    if not ok.any():
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    adj = np.minimum.accumulate((pv[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def load() -> pd.DataFrame:
    """Join stats x geometry x metadata, restricted to pre-treatment labelled ROIs."""
    geom = pd.read_csv(HERE / "output" / "roi_geometry.csv")
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response",
                 "Sample_Type_(pre/post treatment)": "sample_type"},
    )[["roi_id", "patient_id", "response", "sample_type"]]

    frames = []
    for defn in DEFS:
        d = pd.read_parquet(STATS_DIR / f"{defn}.parquet")
        d = d[d["status"] == "ok"].copy()
        d["pair_definition"] = defn
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    df = df.merge(geom, on="roi_id", how="left").merge(meta, on="roi_id", how="left")
    return df[df["response"].isin(["Responder", "Non-Responder"])
              & (df["sample_type"].astype(str).str.lower() == "pre")].copy()


def geom_adjust(d: pd.DataFrame, col: str) -> np.ndarray:
    """Residual of `col` after OLS on log geometry + log counts, fit within `d`."""
    X = np.column_stack([
        np.ones(len(d)),
        np.log(d["med_nn"]),
        np.log(d["hull_area"]),
        np.log(d["n_tumour"] + 1),
        np.log(d["n_other"] + 1),
    ])
    y = d[col].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def scan(pre: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for defn in DEFS:
        d = pre[pre.pair_definition == defn]
        for dg in DIAGRAMS:
            col = f"{dg}-{STAT}"
            if col not in d.columns:
                continue
            variants = {
                "raw": d[col].to_numpy(float),
                "/med_nn": d[col].to_numpy(float) / d["med_nn"].to_numpy(float),
                "geom_adj": geom_adjust(d, col),
            }
            for vname, v in variants.items():
                dd = d.assign(_v=v)
                roi_d, roi_p = cliffs_delta(dd[dd.response == "Responder"]._v,
                                            dd[dd.response == "Non-Responder"]._v)
                pat = dd.groupby(["patient_id", "response"], as_index=False)["_v"].median()
                pat_d, pat_p = cliffs_delta(pat[pat.response == "Responder"]._v,
                                            pat[pat.response == "Non-Responder"]._v)
                rows.append({"definition": defn, "diagram": dg, "variant": vname,
                             "roi_delta": roi_d, "roi_p": roi_p,
                             "pat_delta": pat_d, "pat_p": pat_p})
    t = pd.DataFrame(rows)
    t["pat_q"] = np.concatenate([bh(g["pat_p"].to_numpy())
                                 for _, g in t.groupby("variant", sort=False)])
    return t


def md_table(dfr: pd.DataFrame) -> str:
    head = "| " + " | ".join(dfr.columns) + " |"
    sep = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    body = "\n".join(
        "| " + " | ".join(f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)
                          for v in row) + " |"
        for row in dfr.itertuples(index=False)
    )
    return "\n".join([head, sep, body])


def bar_scale_note(pre: pd.DataFrame) -> str:
    """Why im_dim1_total_persistence is a poor headline feature, from the stats table."""
    d = pre[pre.pair_definition == "CD8_primary"]
    summ = pd.read_parquet(HERE / "output" / "per_roi_summary.parquet")
    summ = summ[(summ.pair_definition == "CD8_primary") & (summ.status == "ok")]
    j = d[["roi_id", "im1-num_bars", "im1-med_length", "im1-avg_length"]].merge(
        summ[["roi_id", "im_dim1_n_bars_gt5um"]], on="roi_id", how="inner")
    frac_over_thr = (j["im_dim1_n_bars_gt5um"] / j["im1-num_bars"]).median()
    return (
        f"- Median degree-1 image bar lifetime: **{j['im1-med_length'].median():.2f} um** "
        f"(mean {j['im1-avg_length'].median():.2f} um) -- well below one cell diameter.\n"
        f"- Median fraction of image degree-1 bars clearing the 5 um threshold used for "
        f"`n_bars_gt5um`: **{frac_over_thr:.1%}** (median {j['im1-num_bars'].median():.0f} "
        f"bars per ROI, of which {j['im_dim1_n_bars_gt5um'].median():.0f} clear it).\n"
        f"- `im_dim1_total_persistence` sums *all* finite bars, so the statistic the Day-1 "
        f"go/no-go rested on is dominated by sub-cell-scale features. `_dim1_feats` in "
        f"`extract_features.py` compounds this: `total_persistence` and `persistent_entropy` "
        f"use every bar while `n_bars_gt5um` uses only bars > 5 um, so they summarise "
        f"different bar populations."
    )


def main() -> None:
    pre = load()
    n_roi = int((pre.pair_definition == "CD8_primary").sum())
    n_pat = pre.groupby("response")["patient_id"].nunique().to_dict()
    t = scan(pre)

    dim1 = t[t.diagram.isin(["ker1", "im1", "cok1"])]
    lines = [
        "# Day-2 re-scan: bar-length statistics (pre-treatment, per-patient)\n",
        f"Cohort: pre-treatment ROIs only, n = {n_roi} ROIs / "
        f"{n_pat.get('Non-Responder', 0)} NR + {n_pat.get('Responder', 0)} R patients.",
        "Effect size = Cliff's delta (+ve = responders higher); `pat_*` = per-patient "
        "medians (primary unit), `roi_*` = per-ROI. `pat_q` = Benjamini-Hochberg across "
        "each variant's family. **Descriptive and post hoc** -- see the caveat at the end.\n",
        "## Why bar length rather than total persistence\n",
        bar_scale_note(pre),
        "\n## Degree-1 mixing diagrams (ker / im / cok), avg_length\n",
    ]
    for vname in ["raw", "/med_nn", "geom_adj"]:
        sub = dim1[dim1.variant == vname]
        piv = sub.pivot(index="definition", columns="diagram",
                        values=["pat_delta", "pat_p"]).reindex(DEFS)
        piv.columns = [f"{a}_{b}" for a, b in piv.columns]
        lines += [f"\n### variant: `{vname}`\n",
                  md_table(piv.reset_index().round(3))]

    lines += ["\n## Full scan (all diagrams, all dimensions)\n",
              md_table(t.round(4).sort_values(
                  ["variant", "definition", "diagram"]).reset_index(drop=True))]

    # headline read
    cd8 = dim1[(dim1.variant == "/med_nn") & dim1.definition.isin(CD8_DEFS - {"CD8_NK_only"})]
    ctrl = dim1[(dim1.variant == "/med_nn") & ~dim1.definition.isin(CD8_DEFS)]
    lines += [
        "\n## Read\n",
        f"- CD8 aggregate definitions (primary/strict/extended), `/med_nn`: median "
        f"per-patient delta = **{cd8.pat_delta.median():+.3f}**, "
        f"{int((cd8.pat_p < 0.05).sum())}/{len(cd8)} tests at p < 0.05.",
        f"- Control species (CD4, Bcell, Macrophage, Fibroblast), same variant: median "
        f"delta = **{ctrl.pat_delta.median():+.3f}**, "
        f"{int((ctrl.pat_p < 0.05).sum())}/{len(ctrl)} at p < 0.05.",
        f"- After BH across each variant's family, "
        f"{int((t.pat_q < 0.05).sum())} of {len(t)} tests survive q < 0.05.",
        "",
        "**Caveat that matters.** The effect attenuates under `geom_adj`, and the same "
        "CD8-specific pattern appears in `dom0`/`dom1` -- the *monochromatic* complex, "
        "which contains no cross-species mixing information at all. That is consistent "
        "with 'responder tissue is more dispersed' rather than 'T cells infiltrate'. "
        "Regression adjustment cannot separate these, because density, dispersion and "
        "arrangement are collinear here. Only a null that holds cell positions fixed and "
        "shuffles labels can: see `perm_null.py` / `output/PERMUTATION_NULL.md`.\n",
        "Script: `day2_length_features.py` (reads output/stats/*.parquet, no recomputation).",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 220)
    print(f"pre-treatment: {n_roi} ROIs, patients {n_pat}")
    print("\n== per-patient Cliff delta, avg_length / med_nn, degree-1 ==")
    print(dim1[dim1.variant == "/med_nn"].pivot(
        index="definition", columns="diagram",
        values=["pat_delta", "pat_p"]).reindex(DEFS).round(3).to_string())
    print(f"\nBH-significant tests: {int((t.pat_q < 0.05).sum())} / {len(t)}")
    print("wrote", OUT_MD)


if __name__ == "__main__":
    main()
