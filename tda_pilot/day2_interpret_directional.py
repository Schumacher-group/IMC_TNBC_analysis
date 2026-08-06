#!/usr/bin/env python
"""What does the directional statistic actually see? Face validity for ker1-avg_length z.

The directional permutation null (`perm_null.py --inclusion other_only`) gives a significant
responder-vs-non-responder difference (delta = +0.379, p = 0.011) whose sign is OPPOSITE to
the immune-exclusion hypothesis. Before that can be interpreted at all, we need to know what
a high or low z-score looks like on real tissue.

Two modes:

  --gallery   Size-matched point-cloud gallery of the extreme ROIs by directional z, in the
              style of `day1_residual_gallery.py`. ROIs are restricted to the cohort
              interquartile band of tumour+CD8 count so the contrast is spatial arrangement
              rather than cell number. If high z looks interspersed and low z looks excluded
              (or the reverse), that fixes the meaning of the statistic.

  --rois A B  Score named ROIs against their own label-permutation nulls and place them in
              the cohort distribution. Use this to check whether the statistic reproduces a
              difference that is obvious by eye -- e.g. the two example panels of Fig 4B.
              Recomputes the six-pack for those ROIs directly, so it works for any ROI,
              including ones outside the pre-treatment set the null was run on.

Output: output/directional_gallery.png, output/DIRECTIONAL_INTERPRETATION.md
"""

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pair_defs
import perm_null

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
ROIS_DIR = HERE / "data" / "rois"
STAT = perm_null.PRIMARY  # ker1-avg_length
N_PANEL = 6


def load_scores(defname: str = "CD8_primary") -> pd.DataFrame:
    d = pd.read_parquet(OUT / f"perm_null_{defname}_other_only.parquet")
    d = d[d["status"] == "ok"].copy()
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response"},
    )[["roi_id", "patient_id", "response"]]
    d = d.merge(meta, on="roi_id", how="left")
    d["n_pair"] = d.n_tumour + d.n_other
    d["z"] = d[f"{STAT}__z"]
    return d


def draw(ax, fov: str, title: str, species: set[str]) -> None:
    df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
    is_t = df["celltype"].isin(pair_defs.TUMOUR)
    is_s = df["celltype"].isin(species)
    t, s = df[is_t], df[is_s]
    ax.scatter(t.x, t.y, s=1.5, c="tab:red", alpha=0.55, linewidths=0, label="Tumour")
    ax.scatter(s.x, s.y, s=1.5, c="tab:blue", alpha=0.85, linewidths=0, label="CD8")
    ax.set_title(title, fontsize=8)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])


def gallery(d: pd.DataFrame, species: set[str]) -> Path:
    lo, hi = d.n_pair.quantile([0.25, 0.75])
    band = d[(d.n_pair >= lo) & (d.n_pair <= hi)].sort_values("z")
    low, high = band.head(N_PANEL), band.tail(N_PANEL).iloc[::-1]

    fig, axes = plt.subplots(2, N_PANEL, figsize=(3.0 * N_PANEL, 6.4))
    for ax, (_, r) in zip(axes[0], high.iterrows()):
        draw(ax, r.roi_id, f"{r.roi_id} ({r.response[0]})\nz={r.z:+.1f}  n={r.n_pair:.0f}",
             species)
    for ax, (_, r) in zip(axes[1], low.iterrows()):
        draw(ax, r.roi_id, f"{r.roi_id} ({r.response[0]})\nz={r.z:+.1f}  n={r.n_pair:.0f}",
             species)
    axes[0, 0].set_ylabel(f"HIGHEST {STAT} z", fontsize=9)
    axes[1, 0].set_ylabel(f"LOWEST {STAT} z", fontsize=9)
    fig.suptitle(
        f"Directional (CD8-only domain) {STAT} z-score extremes — size-matched band "
        f"{lo:.0f}–{hi:.0f} tumour+CD8 cells (n={len(band)}).  red = tumour, blue = CD8",
        fontsize=11)
    fig.tight_layout()
    path = OUT / "directional_gallery.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    return path


STATS_REPORT = ("ker1-avg_length", "ker0-avg_length", "ker1-med_length",
                "ker1-p90_length", "cok1-avg_length", "im1-avg_length", "ker1-num_bars")


def geom_adjusted(d: pd.DataFrame, col: str) -> np.ndarray:
    """Residual of `col` on composition + geometry, fit across the cohort."""
    n_pair = d.n_tumour + d.n_other
    X = np.column_stack([
        np.ones(len(d)), d.n_other / n_pair, np.log(n_pair),
        np.log(d.density), np.log(d.hull_area),
    ])
    y = d[col].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def pair_report(rois: list[str], defname: str) -> tuple[str, pd.DataFrame]:
    """Compare named ROIs using the stored cohort scores, under both inclusions.

    Reports raw z, geometry-adjusted z and cohort percentile for each statistic. The
    adjustment matters: a raw z gap between two ROIs of different size or density is
    partly geometry (see the retraction in REVIEW.md). If the gap survives adjustment,
    it is arrangement.
    """
    geom = pd.read_csv(OUT / "roi_geometry.csv")
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Response": "response"})[["roi_id", "response"]]
    out, blocks = [], []
    for tag, name in (("_other_only", "directional"), ("", "symmetric")):
        f = OUT / f"perm_null_{defname}{tag}.parquet"
        if not f.exists():
            continue
        d = pd.read_parquet(f)
        d = d[d["status"] == "ok"].merge(geom, on="roi_id").merge(meta, on="roi_id")
        if not set(rois).issubset(set(d.roi_id)):
            missing = sorted(set(rois) - set(d.roi_id))
            blocks.append(f"\n**{name}**: not scored for {missing} "
                          f"(outside the pre-treatment set the null was run on).\n")
            continue
        s = d.set_index("roi_id")
        for stat in STATS_REPORT:
            col = f"{stat}__z"
            if col not in d.columns:
                continue
            d["_adj"] = geom_adjusted(d, col)
            sa = d.set_index("roi_id")
            row = {"inclusion": name, "statistic": stat}
            for r in rois:
                row[f"z[{r}]"] = s.loc[r, col]
                row[f"adj[{r}]"] = sa.loc[r, "_adj"]
                row[f"pct[{r}]"] = 100 * (d["_adj"] < sa.loc[r, "_adj"]).mean()
            row["adj_gap"] = row[f"adj[{rois[1]}]"] - row[f"adj[{rois[0]}]"]
            out.append(row)
    return "\n".join(blocks), pd.DataFrame(out)


def pair_figure(rois: list[str], species: set[str], labels: list[str]) -> Path:
    fig, axes = plt.subplots(1, len(rois), figsize=(6.2 * len(rois), 6.0))
    for ax, fov, lab in zip(np.atleast_1d(axes), rois, labels):
        draw(ax, fov, lab, species)
    fig.suptitle("Named ROI comparison — red = tumour, blue = CD8 (CD8_primary)", fontsize=12)
    fig.tight_layout()
    path = OUT / "directional_named_pair.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    return path


def score_rois(rois: list[str], defname: str, nperm: int, seed: int) -> pd.DataFrame:
    """Recompute obs + null for arbitrary ROIs under the directional inclusion."""
    perm_null._init(pair_defs.DEFINITIONS[defname], defname, nperm, seed, "other_only")
    rows = []
    for fov in rois:
        r = perm_null._process_roi(fov)
        if r.get("status") != "ok":
            rows.append({"roi_id": fov, "status": r.get("status"), "error": r.get("error", "")})
            continue
        rows.append({
            "roi_id": fov, "status": "ok",
            "n_tumour": r["n_tumour"], "n_other": r["n_other"],
            "obs": r[f"{STAT}__obs"], "null_mean": r[f"{STAT}__null_mean"],
            "null_sd": r[f"{STAT}__null_sd"], "z": r[f"{STAT}__z"],
            "null_rank": r[f"{STAT}__null_rank"],
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--definition", default="CD8_primary")
    ap.add_argument("--gallery", action="store_true")
    ap.add_argument("--rois", nargs="*", default=[],
                    help="ROI ids to score individually, e.g. the two Fig 4B panels")
    ap.add_argument("--n-perm", type=int, default=99)
    ap.add_argument("--seed", type=int, default=20260805)
    args = ap.parse_args()

    d = load_scores(args.definition)
    species = pair_defs.DEFINITIONS[args.definition]
    lines = [f"# Directional statistic: what does `{STAT}` z actually measure?\n",
             f"Cohort reference: {len(d)} pre-treatment ROIs, "
             f"z median {d.z.median():+.2f}, IQR "
             f"[{d.z.quantile(.25):+.2f}, {d.z.quantile(.75):+.2f}], "
             f"range [{d.z.min():+.1f}, {d.z.max():+.1f}].\n"]

    if args.gallery:
        p = gallery(d, species)
        lines += [f"## Size-matched extremes\n\nFigure: `{p.name}`. Top row = highest z, "
                  f"bottom row = lowest z, within the interquartile band of tumour+CD8 count "
                  f"so the contrast is arrangement rather than cell number.\n"]
        print("wrote", p)

    if len(args.rois) == 2:
        # Both ROIs are usually already scored in the cohort run; prefer those numbers so
        # the comparison uses exactly the same seeds and nulls as everything else.
        note, tab = pair_report(args.rois, args.definition)
        if len(tab):
            g = pd.read_csv(OUT / "roi_geometry.csv")
            meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
                columns={"fov": "roi_id", "Response": "response",
                         "Sample_Type_(pre/post treatment)": "sample_type"})
            gg = g[g.roi_id.isin(args.rois)].merge(
                meta[["roi_id", "response", "sample_type", "tumour", "cd8_primary"]],
                on="roi_id")
            gg["cd8_frac"] = gg.cd8_primary / (gg.tumour + gg.cd8_primary)
            labels = [f"{r} ({gg.set_index('roi_id').loc[r,'response']})" for r in args.rois]
            fig = pair_figure(args.rois, species, labels)
            pd.set_option("display.width", 240)
            print("\n=== geometry / composition ===")
            print(gg[["roi_id", "response", "tumour", "cd8_primary", "cd8_frac",
                      "hull_area", "density", "med_nn"]].round(4).to_string(index=False))
            print("\n=== z-scores, raw and geometry-adjusted ===")
            print(tab.round(2).to_string(index=False))
            lines += ["## Named pair\n", "### Geometry / composition\n",
                      "| " + " | ".join(gg.columns) + " |",
                      "| " + " | ".join("---" for _ in gg.columns) + " |"]
            lines += ["| " + " | ".join(
                f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
                for v in row) + " |" for row in gg.itertuples(index=False)]
            lines += ["\n### z-scores (raw, geometry-adjusted, cohort percentile)\n",
                      "| " + " | ".join(tab.columns) + " |",
                      "| " + " | ".join("---" for _ in tab.columns) + " |"]
            lines += ["| " + " | ".join(
                f"{v:.2f}" if isinstance(v, (float, np.floating)) else str(v)
                for v in row) + " |" for row in tab.round(2).itertuples(index=False)]
            lines += [note, f"\nFigure: `{fig.name}`.\n"]
            print("wrote", fig)
        args.rois = []  # handled

    if args.rois:
        sc = score_rois(args.rois, args.definition, args.n_perm, args.seed)
        # locate each in the cohort distribution
        ok = sc[sc.status == "ok"].copy()
        if len(ok):
            ok["cohort_pctile"] = [float((d.z < z).mean() * 100) for z in ok["z"]]
        meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(columns={"fov": "roi_id"})
        ok = ok.merge(meta[["roi_id", "Response", "Sample_Type_(pre/post treatment)"]],
                      on="roi_id", how="left")
        pd.set_option("display.width", 220)
        print("\n=== named-ROI scores (directional inclusion) ===")
        print(ok.round(3).to_string(index=False))
        if len(sc[sc.status != "ok"]):
            print("\nfailed/degenerate:")
            print(sc[sc.status != "ok"].to_string(index=False))
        lines += ["## Named ROIs\n",
                  "| " + " | ".join(ok.columns) + " |",
                  "| " + " | ".join("---" for _ in ok.columns) + " |"]
        lines += ["| " + " | ".join(
            f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v) for v in row) + " |"
            for row in ok.round(3).itertuples(index=False)]
        if len(ok) == 2:
            a, b = ok.iloc[0], ok.iloc[1]
            lines += [f"\n**Difference:** z {a.roi_id} = {a.z:+.2f} vs {b.roi_id} = {b.z:+.2f} "
                      f"(gap {abs(a.z - b.z):.2f} sd of the per-ROI null); cohort percentiles "
                      f"{a.cohort_pctile:.0f} vs {b.cohort_pctile:.0f}.\n"]

    (OUT / "DIRECTIONAL_INTERPRETATION.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUT / "DIRECTIONAL_INTERPRETATION.md")


if __name__ == "__main__":
    main()
