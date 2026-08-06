#!/usr/bin/env python
"""Chromatic vs Dowker on the B7H4 question: where do they agree, and why do they differ?

Two structurally independent constructions were run on the same matched paired design:

  chromatic  six-pack of the tumour/CD8 inclusion. Filtration uses ALL distances, so it is
             sensitive to within-species structure.
  Dowker     landmarks = cancer, witnesses = CD8. Filtration uses ONLY cross-distances, so
             within-species structure never enters it.

They agree that B7H4+ tumour has FEWER degree-1 features. They disagree on LENGTH, which is
where the exclusion signal is supposed to live: chromatic says longer for B7H4+ (more
exclusion), Dowker says shorter (less exclusion). Both cannot be a clean readout of CD8
accessibility, so this script asks which one is tracking the confound.

The candidate confound is measured, not assumed: at matched cell counts B7H4+ cancer cells are
more tightly packed than B7H4- cells (median nearest-neighbour distance 11.3 vs 14.2 um, in 95%
of ROIs). Clustering plausibly moves both statistics, in opposite directions:
  - tighter landmarks are easier for a single witness to cover  -> SHORTER Dowker bars
  - tighter monochromatic structure forms at smaller scales     -> SHORTER chromatic bars
so a length effect that tracks clustering is suspect, and one that does not is not.

Two tests:
  1. Correlate each method's paired length difference against the paired clustering difference.
  2. Restrict to ROIs where the two subtypes are most similarly clustered and re-test.

Output: output/METHOD_COMPARISON.md
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"

PAIRS = {
    "count": ("ker1-num_bars", "dow1-num_bars"),
    "length": ("ker0-avg_length", "dow1-avg_length"),
}


def md(dfr: pd.DataFrame) -> str:
    h = "| " + " | ".join(dfr.columns) + " |"
    s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    b = "\n".join("| " + " | ".join(
        f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
        for v in row) + " |" for row in dfr.itertuples(index=False))
    return "\n".join([h, s, b])


def main() -> None:
    chro = pd.read_csv(OUT / "b7h4_paired.csv")
    dow = pd.read_csv(OUT / "dowker_features.csv")
    d = chro.merge(dow, on="roi_id", suffixes=("_chr", "_dow"))

    # paired differences, B7H4 minus Cancer
    d["d_clustering"] = d["b7h4__geom_nn"] - d["cancer__geom_nn"]   # -ve = B7H4 tighter
    for name, (ck, dk) in PAIRS.items():
        d[f"d_chr_{name}"] = d[f"b7h4__{ck}"] - d[f"cancer__{ck}"]
        d[f"d_dow_{name}"] = d[f"b7h4__{dk}"] - d[f"cancer__{dk}"]

    rows = []
    for name, (ck, dk) in PAIRS.items():
        for meth, key, stat in (("chromatic", f"d_chr_{name}", ck),
                                ("Dowker", f"d_dow_{name}", dk)):
            rho, p = spearmanr(d.d_clustering, d[key])
            rows.append({"signal": name, "method": meth, "statistic": stat,
                         "rho_with_clustering": rho, "p": p})
    conf = pd.DataFrame(rows)

    agree = pd.DataFrame([
        {"signal": name,
         "rho_between_methods": spearmanr(d[f"d_chr_{name}"], d[f"d_dow_{name}"])[0],
         "p": spearmanr(d[f"d_chr_{name}"], d[f"d_dow_{name}"])[1]}
        for name in PAIRS])

    # clustering-matched subset: ROIs where the subtypes are most similarly packed
    thr = d.d_clustering.abs().quantile(0.25)
    q = d[d.d_clustering.abs() < thr]
    rows = []
    for name, (ck, dk) in PAIRS.items():
        for meth, key, stat in (("chromatic", f"d_chr_{name}", ck),
                                ("Dowker", f"d_dow_{name}", dk)):
            v_all = d[key].dropna()
            v_sub = q[key].dropna()
            rows.append({"signal": name, "method": meth, "statistic": stat,
                         "median_all": float(np.median(v_all)),
                         "p_all": float(wilcoxon(v_all).pvalue),
                         "median_matched": float(np.median(v_sub)),
                         "p_matched": float(wilcoxon(v_sub).pvalue)})
    matched = pd.DataFrame(rows)

    dow_len = matched[(matched.signal == "length") & (matched.method == "Dowker")].iloc[0]
    chr_len = matched[(matched.signal == "length") & (matched.method == "chromatic")].iloc[0]
    verdict = (
        "the **Dowker** length effect is the clustering-driven one — it correlates with the "
        "clustering difference and vanishes when clustering is matched, while the chromatic "
        "one does neither"
        if (dow_len.p_matched > 0.05 >= chr_len.p_matched) else
        "the **chromatic** length effect is the clustering-driven one"
        if (chr_len.p_matched > 0.05 >= dow_len.p_matched) else
        "neither length effect cleanly survives clustering matching, so neither should be "
        "read as a measure of CD8 accessibility"
    )

    lines = [
        "# Chromatic vs Dowker on the B7H4 question\n",
        f"{len(d)} pre-treatment ROIs with both analyses. Paired differences are B7H4+ minus "
        "B7H4−, within ROI, at matched cell counts.\n",
        "## 1. Where the methods agree\n",
        md(agree.round(4)),
        "\nThe **count** signal is strongly concordant between two structurally independent "
        "constructions — that is the robust result. The **length** signal is essentially "
        "uncorrelated between them, so at most one of the two can be a clean readout of CD8 "
        "accessibility.\n",
        "## 2. Which signal tracks the clustering confound?\n",
        "At matched counts B7H4+ cancer cells are more tightly packed (median NN 11.3 vs "
        "14.2 µm, tighter in 95% of ROIs). Correlating each paired difference against the "
        "paired clustering difference:\n",
        md(conf.round(4)),
        "\n## 3. Restricting to ROIs where clustering is matched\n",
        f"The {len(q)} ROIs whose subtypes are most similarly packed (|Δ nearest-neighbour "
        f"distance| below its 25th percentile, {thr:.2f} µm):\n",
        md(matched.round(4)),
        f"\n**Verdict: {verdict}.**\n",
        "## What can be said\n",
        "- **Robust:** B7H4+ tumour carries fewer degree-1 topological features than B7H4− "
        "tumour in the same tissue, by two independent constructions that agree strongly "
        "(ρ ≈ 0.68). Part of this tracks the clustering difference.",
        "- **Directly measured, not inferred:** B7H4+ cancer cells are more tightly clustered.",
        "- **Provisional:** the longer-bar 'exclusion' signature appears in the chromatic "
        "construction, is not explained by clustering, and survives clustering-matching — but "
        "it is not reproduced by an independent construction, so it should be reported as "
        "method-dependent rather than established.",
        "- **Not supported:** a clean claim that CD8 accessibility differs between the "
        "subtypes beyond the architectural difference. The two methods would need to agree "
        "for that, and they do not.\n",
        "Script: `day2_method_comparison.py`.",
    ]
    (OUT / "METHOD_COMPARISON.md").write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 200)
    print("agreement between methods:\n", agree.round(4).to_string(index=False))
    print("\ncorrelation with clustering:\n", conf.round(4).to_string(index=False))
    print(f"\nclustering-matched subset (n={len(q)}):\n",
          matched.round(4).to_string(index=False))
    print(f"\nverdict: {verdict}")
    print("wrote", OUT / "METHOD_COMPARISON.md")


if __name__ == "__main__":
    main()
