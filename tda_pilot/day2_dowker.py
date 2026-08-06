#!/usr/bin/env python
"""Dowker PH: independent confirmation of the B7H4 result, plus the response contrast.

Motivation. Every result in this directory so far rests on chromatic six-packs, whose
statistics are entangled with within-species structure -- bar count tracks cell number at
rho = 0.83, bar length tracks composition at rho = -0.69, and z-scores stay coupled to ROI
density. We had to argue around that repeatedly, and once got it wrong (see REVIEW.md
UPDATE 2/4). Dowker PH uses ONLY cross-distances between the two species, so it is
structurally immune to that whole family of confounds. That makes it a genuinely independent
check rather than a re-summary of the same information. See `dowker.py` for the construction
and its validation.

Landmarks are cancer cells, witnesses are CD8 cells, so a long degree-1 bar is a loop of
cancer cells that no single CD8 cell is close to: a region of tumour CD8 does not reach.

Two questions, one pass:

  PRIMARY   Is CD8 less able to reach B7H4+ tumour than B7H4- tumour? Same matched paired
            design as `day2_b7h4_paired.py` -- within each ROI both cancer subtypes are
            subsampled to the same count and witnessed by the identical CD8 set, averaged
            over repeats. Prediction if the chromatic finding is real: LONGER degree-1
            Dowker bars for B7H4+.

  SECONDARY Does whole-tumour Dowker structure differ by response? Near-free once the
            machinery exists. Expect null -- five independent framings already are -- but a
            better-matched instrument yields a larger standardised effect, so it is not
            strictly true that the descriptor cannot matter. Reported with that expectation
            set, and it would need external validation either way.

Landmarks are capped (`--cap`) for tractability; the cap is not confound control -- the paired
design already handles between-ROI variation in cell number, since both arms of a pair share
the same count.

Output: output/DOWKER.md, output/dowker_b7h4.png
"""

import argparse
import multiprocessing as mp
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

import cohort
import dowker

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
ROIS_DIR = HERE / "data" / "rois"
CANCER, B7H4 = "Cancer cell", "B7H4 Cancer cell"
CD8 = {"CD8 T cell", "Memory CD8 T cell"}
MIN_EACH, MIN_CD8 = 30, 30
SEED = 20260806
COLOURS = {"Responder": "darkgreen", "Non-Responder": "#8B1A1A"}

_CFG = {"cap": 300, "cutoff": 100.0, "repeats": 5}


def _init(cfg):
    global _CFG
    _CFG = cfg


def _process_roi(fov: str) -> dict | None:
    try:
        df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
        canc = df[df.celltype == CANCER][["x", "y"]].to_numpy(float)
        b7h4 = df[df.celltype == B7H4][["x", "y"]].to_numpy(float)
        cd8 = df[df.celltype.isin(CD8)][["x", "y"]].to_numpy(float)
        tum = np.vstack([canc, b7h4]) if len(canc) or len(b7h4) else np.zeros((0, 2))
        if min(len(canc), len(b7h4)) < MIN_EACH or len(cd8) < MIN_CD8:
            return None

        cap, cut, reps = _CFG["cap"], _CFG["cutoff"], _CFG["repeats"]
        m = min(len(canc), len(b7h4), cap)
        rng = np.random.default_rng(SEED + (int.from_bytes(fov.encode(), "little") % 10**6))
        acc: dict[str, list[float]] = {}
        for _ in range(reps):
            for tag, arr in (("cancer", canc), ("b7h4", b7h4)):
                sel = arr[rng.choice(len(arr), m, replace=False)]
                for k, v in dowker.dowker_features(sel, cd8, cut).items():
                    acc.setdefault(f"{tag}__{k}", []).append(v)
            # whole tumour vs CD8, for the response contrast
            mt = min(len(tum), cap)
            selt = tum[rng.choice(len(tum), mt, replace=False)]
            for k, v in dowker.dowker_features(selt, cd8, cut).items():
                acc.setdefault(f"tumour__{k}", []).append(v)
        row = {"roi_id": fov, "m": m, "n_cd8": len(cd8),
               "n_cancer": len(canc), "n_b7h4": len(b7h4)}
        row |= {k: float(np.mean(v)) for k, v in acc.items()}
        return row
    except Exception:  # noqa: BLE001
        return None


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, float(u.pvalue)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    q = np.empty(len(p))
    q[o] = np.minimum(adj, 1.0)
    return q


def md(dfr):
    h = "| " + " | ".join(dfr.columns) + " |"
    s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    b = "\n".join("| " + " | ".join(
        f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
        for v in row) + " |" for row in dfr.itertuples(index=False))
    return "\n".join([h, s, b])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=300)
    ap.add_argument("--cutoff", type=float, default=100.0)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--reuse", action="store_true")
    args = ap.parse_args()
    mp.set_start_method("fork", force=True)

    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "response",
                 "Sample_Type_(pre/post treatment)": "sample_type"})
    meta = meta[(meta.sample_type.astype(str).str.lower() == "pre")
                & meta.response.isin(["Responder", "Non-Responder"])]
    rois = sorted(r for r in meta.roi_id
                  if cohort.is_cohort_member(r) and (ROIS_DIR / f"{r}.csv").exists())

    cache = OUT / "dowker_features.csv"
    if args.reuse and cache.exists():
        d = pd.read_csv(cache)
        print(f"reusing {cache.name}: {len(d)} ROIs")
    else:
        cfg = {"cap": args.cap, "cutoff": args.cutoff, "repeats": args.repeats}
        t0 = time.perf_counter()
        print(f"{len(rois)} ROIs, cap={args.cap} landmarks, cutoff={args.cutoff} um, "
              f"{args.repeats} repeats, {args.workers} workers")
        with mp.Pool(args.workers, initializer=_init, initargs=(cfg,)) as pool:
            out = []
            for i, r in enumerate(pool.imap_unordered(_process_roi, rois, chunksize=2), 1):
                if r is not None:
                    out.append(r)
                if i % 100 == 0 or i == len(rois):
                    print(f"  {i}/{len(rois)} ({(time.perf_counter()-t0)/60:.1f} min)")
        d = pd.DataFrame(out).merge(meta[["roi_id", "pid", "response"]], on="roi_id")
        d.to_csv(cache, index=False)
        print(f"{len(d)} ROIs in {(time.perf_counter()-t0)/60:.1f} min")

    stats = sorted({c.split("__", 1)[1] for c in d.columns if "__" in c})

    # ---- PRIMARY: paired B7H4 vs Cancer ------------------------------------
    rows = []
    for s in stats:
        a, b = f"b7h4__{s}", f"cancer__{s}"
        if a not in d or b not in d:
            continue
        diff = (d[a] - d[b]).replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
        if len(diff) < 20 or np.allclose(diff, 0):
            continue
        rows.append({"statistic": s, "median_cancer": float(d[b].median()),
                     "median_b7h4": float(d[a].median()),
                     "median_diff": float(np.median(diff)),
                     "frac_b7h4_higher": float((diff > 0).mean()),
                     "wilcoxon_p": float(wilcoxon(diff).pvalue),
                     "rank_biserial": 2 * float((diff > 0).mean()) - 1})
    prim = pd.DataFrame(rows)
    prim["q"] = bh(prim.wilcoxon_p.to_numpy())
    prim = prim.sort_values("wilcoxon_p").reset_index(drop=True)

    # ---- SECONDARY: response contrast on whole-tumour Dowker ----------------
    rows = []
    for s in stats:
        c = f"tumour__{s}"
        if c not in d:
            continue
        pat = d.groupby(["pid", "response"], as_index=False)[c].median()
        dl, pv = cliffs(pat[pat.response == "Responder"][c],
                        pat[pat.response == "Non-Responder"][c])
        if np.isfinite(pv):
            rows.append({"statistic": s, "cliff_R_vs_NR": dl, "p": pv})
    sec = pd.DataFrame(rows)
    sec["q"] = bh(sec.p.to_numpy())
    sec = sec.sort_values("p").reset_index(drop=True)

    key = "dow1-avg_length"
    n_prim = int((prim.q < 0.05).sum())
    n_sec = int((sec.q < 0.05).sum())

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    if (prim.statistic == key).any():
        diff = (d[f"b7h4__{key}"] - d[f"cancer__{key}"]).dropna()
        ax[0].hist(diff, bins=45, color="teal", alpha=0.8)
        ax[0].axvline(0, color="black", lw=1.5, ls="--")
        ax[0].axvline(diff.median(), color="crimson", lw=2,
                      label=f"median Δ = {diff.median():+.3f}")
        kr = prim[prim.statistic == key].iloc[0]
        ax[0].set_xlabel(f"Δ {key} (B7H4+ − B7H4−), µm")
        ax[0].set_ylabel("ROIs")
        ax[0].set_title(f"A. Dowker: paired B7H4 difference\n{kr.frac_b7h4_higher:.0%} of "
                        f"ROIs higher for B7H4+, p={kr.wilcoxon_p:.2g}", fontsize=10)
        ax[0].legend(fontsize=8)
        lim = np.percentile(np.r_[d[f"cancer__{key}"], d[f"b7h4__{key}"]], [1, 99])
        ax[1].scatter(d[f"cancer__{key}"], d[f"b7h4__{key}"], s=11, c="teal",
                      alpha=0.45, linewidths=0)
        ax[1].plot(lim, lim, "k--", lw=1.2)
        ax[1].set_xlim(*lim); ax[1].set_ylim(*lim)
        ax[1].set_xlabel("B7H4− cancer, witnessed by CD8")
        ax[1].set_ylabel("B7H4+ cancer, witnessed by CD8")
        ax[1].set_title("B. Same tissue, matched counts\n(above line = CD8 reaches B7H4+ less)",
                        fontsize=10)
    pat = d.groupby(["pid", "response"], as_index=False)[f"tumour__{key}"].median()
    jit = np.random.default_rng(0)
    for i, grp in enumerate(["Non-Responder", "Responder"]):
        v = pat[pat.response == grp][f"tumour__{key}"].to_numpy()
        ax[2].scatter(jit.normal(i, 0.08, len(v)), v, s=38, c=COLOURS[grp], alpha=0.85,
                      linewidths=0)
        ax[2].hlines(np.median(v), i - 0.26, i + 0.26, color="black", lw=2.4)
    sr = sec[sec.statistic == key]
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["NR", "R"])
    ax[2].set_ylabel(f"patient-median {key}")
    ax[2].set_title("C. Response contrast (secondary)\n"
                    + (f"Cliff δ={sr.cliff_R_vs_NR.iloc[0]:+.3f}, p={sr.p.iloc[0]:.3f}"
                       if len(sr) else ""), fontsize=10)
    fig.suptitle(f"Dowker persistent homology — cancer witnessed by CD8, {len(d)} "
                 f"pre-treatment ROIs (cap {args.cap} landmarks, cutoff {args.cutoff:.0f} µm)",
                 fontsize=12)
    fig.tight_layout()
    figp = OUT / "dowker_b7h4.png"
    fig.savefig(figp, dpi=150, bbox_inches="tight")

    lines = [
        "# Dowker persistent homology: independent check on the B7H4 finding\n",
        f"{len(d)} pre-treatment ROIs. Landmarks = cancer cells (cap {args.cap}), witnesses = "
        f"CD8 cells, cutoff {args.cutoff:.0f} µm, {args.repeats} matched subsamples averaged. "
        "A long degree-1 bar is a loop of cancer cells that no single CD8 cell is close to — "
        "a region of tumour CD8 does not reach.\n",
        "**Why this is an independent check, not a re-summary:** the Dowker filtration uses "
        "only cross-distances between the two species. Within-species structure — the source "
        "of every confound in the chromatic analyses (bar count ~ cell number ρ = 0.83, bar "
        "length ~ composition ρ = −0.69, z coupled to density) — never enters it.\n",
        f"## Primary: paired B7H4+ vs B7H4−  ({n_prim}/{len(prim)} at BH q < 0.05)\n",
        md(prim.round(4)),
        f"\n## Secondary: response contrast on whole-tumour Dowker  "
        f"({n_sec}/{len(sec)} at BH q < 0.05)\n",
        md(sec.round(4)),
        f"\nFigure: `{figp.name}`. Scripts: `dowker.py` (construction + validation), "
        "`day2_dowker.py`.",
    ]
    (OUT / "DOWKER.md").write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 220)
    print("\n=== PRIMARY: paired B7H4 vs Cancer (Dowker) ===")
    print(prim.round(4).to_string(index=False))
    print(f"\n=== SECONDARY: response contrast ({n_sec}/{len(sec)} q<0.05) ===")
    print(sec.head(8).round(4).to_string(index=False))
    print("\nwrote", OUT / "DOWKER.md")


if __name__ == "__main__":
    main()
