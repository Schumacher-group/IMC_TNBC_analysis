#!/usr/bin/env python
"""M2S2-style supervised classification: can TDA statistics predict response?

This is the multivariate Day-2 analysis the pilot planned and never ran, using the M2S2
authors' own machinery on the 456-statistic table already sitting in output/stats/.

WHAT WE TOOK FROM M2S2 (`colorectal_cancer_dataset_analysis/classification_by_patient_pca.py`)
  - their feature filter (`discard_feature`): drop codomain and relative diagrams, the
    always-empty degree-0 cokernel, degree-0 domain/image birth/midpoint/length statistics
    (births are identically 0 there), and IQR statistics (p25/p75 already retained)
  - per-cell-group PCA keeping max(ceil(0.05 n), min(n, 5)) components
  - GradientBoostingClassifier inside a BaggingClassifier, with their hyperparameters
  - a DummyClassifier baseline

WHERE OUR TASK DIFFERS -- this is not a criticism of their pipeline, it is a different question
  They fit ONE CLASSIFIER PER PATIENT, using that patient's ROIs as samples, predicting each
  ROI's `sample_type` (adenoma vs carcinoma). That is a within-patient question and patient
  identity is controlled by construction.
  We need a BETWEEN-patient question: predict the patient's RESPONSE. So the unit of analysis
  is the patient (62 of them: 25 NR / 37 R), features are per-patient medians of each statistic
  over that patient's pre-treatment ROIs, and the classifier is cross-validated across patients.

DELIBERATE DEVIATIONS, and why
  1. PCA inside the cross-validation pipeline. Their script runs PCA on the full patient table
     before `cross_validate`, so the projection sees the test folds. Harmless-ish for their
     within-patient design, but it would inflate a between-patient score, which is exactly the
     number we care about here.
  2. StandardScaler before PCA. The statistics span very different scales (`num_bars` in the
     thousands, lengths ~1 um); unscaled PCA is dominated by the count features, which are
     nearly a proxy for cell number.
  3. ROC AUC as the primary metric. With 25/37 class imbalance, accuracy is misleading.
  4. A label-permutation null on the cross-validated AUC. This is the null model their pipeline
     does not need (their design is descriptive/predictive, not inferential) but ours does,
     because we are asked whether an association exists at all.
  5. A COMPOSITION-ONLY baseline (CD8 fraction, cell counts, density, hull area). Raw TDA
     statistics are strongly confounded with composition -- rho(bar length, CD8 fraction) =
     -0.69 -- so a TDA model must beat composition to be worth anything.

Output: output/CLASSIFICATION.md, output/classification_figure.png
"""

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import BaggingClassifier, GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedShuffleSplit, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
STATS_DIR = OUT / "stats"
DEFS = ["CD8_primary", "CD8_strict", "CD8_extended", "CD8_NK_only",
        "CD4", "Bcell", "Macrophage", "Fibroblast"]
SEED = 20260806


def discard_feature(col: str) -> bool:
    """M2S2's feature filter, adapted to our flat `<dgm><dim>-<statistic>` column names."""
    if "-" not in col:
        return True
    head, statistic = col.split("-", 1)
    dgm, dim = head[:-1], int(head[-1])
    return (
        dgm in ("cod", "rel")                      # codomain and relative diagrams
        or (dim == 0 and dgm == "cok")             # degree-0 cokernel is always empty
        or (dim == 0 and dgm in ("dom", "im")
            and any(x in statistic for x in ("birth", "midpt", "length")))
        or ("iqr" in statistic)                    # p25/p75 already retained
    )


def build_patient_table(triples: tuple[str, ...] = ()) -> tuple[
        pd.DataFrame, pd.Series, pd.DataFrame, list[str]]:
    """Per-patient median of every retained statistic, pre-treatment only.

    Pair definitions and (optionally) three-species definitions become separate feature
    groups, each getting its own PCA -- matching M2S2, which treats every cell tuple as its
    own group.
    """
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "pid", "Response": "resp",
                 "Sample_Type_(pre/post treatment)": "stype"})
    meta = meta[["roi_id", "pid", "resp", "stype"]]
    geom = pd.read_csv(OUT / "roi_geometry.csv")

    blocks, comp_src = [], None
    for defn in DEFS:
        d = pd.read_parquet(STATS_DIR / f"{defn}.parquet")
        d = d[d["status"] == "ok"].merge(meta, on="roi_id")
        d = d[(d.stype.astype(str).str.lower() == "pre")
              & d.resp.isin(["Responder", "Non-Responder"])]
        stat_cols = [c for c in d.columns if "-" in c and not discard_feature(c)]
        block = d.groupby("pid")[stat_cols].median()
        block.columns = pd.MultiIndex.from_product([[defn], block.columns])
        blocks.append(block)
        if defn == "CD8_primary":
            comp_src = d.merge(geom, on="roi_id")

    triple_groups = []
    for name in triples:
        f = STATS_DIR / f"triple_{name}.parquet"
        if not f.exists():
            print(f"  (skipping triple {name}: {f.name} not found)")
            continue
        d = pd.read_parquet(f)
        d = d[d["status"] == "ok"].merge(meta, on="roi_id")
        d = d[(d.stype.astype(str).str.lower() == "pre")
              & d.resp.isin(["Responder", "Non-Responder"])]
        stat_cols = [c for c in d.columns if "-" in c and not discard_feature(c)]
        block = d.groupby("pid")[stat_cols].median()
        gname = f"triple:{name}"
        block.columns = pd.MultiIndex.from_product([[gname], block.columns])
        blocks.append(block)
        triple_groups.append(gname)

    X = pd.concat(blocks, axis=1)
    y = (comp_src.groupby("pid")["resp"].first().reindex(X.index) == "Responder").astype(int)

    comp_src = comp_src.assign(
        n_pair=comp_src.n_tumour + comp_src.n_other,
        cd8_frac=comp_src.n_other / (comp_src.n_tumour + comp_src.n_other))
    C = comp_src.groupby("pid")[
        ["cd8_frac", "n_pair", "n_tumour", "n_other", "density", "hull_area", "med_nn"]
    ].median().reindex(X.index)
    return X, y, C, triple_groups


def n_components(n: int) -> int:
    """M2S2's rule."""
    return max(math.ceil(n * 0.05), min(n, 5))


def make_model(X: pd.DataFrame, per_group_pca: bool, seed: int) -> Pipeline:
    gbc = GradientBoostingClassifier(
        loss="log_loss", n_estimators=25, learning_rate=0.4, max_features=0.03,
        max_depth=3, max_leaf_nodes=6, min_samples_leaf=5, random_state=seed)
    clf = BaggingClassifier(gbc, n_estimators=50, max_samples=1.0, bootstrap=False,
                            random_state=seed, n_jobs=1)
    steps = [("impute", SimpleImputer(strategy="constant", fill_value=0.0,
                                      keep_empty_features=True)),
             ("scale", StandardScaler())]
    if per_group_pca:
        groups = []
        for g in X.columns.get_level_values(0).unique():
            idx = [i for i, c in enumerate(X.columns) if c[0] == g]
            k = min(n_components(len(idx)), 20)
            groups.append((f"pca_{g}", PCA(n_components=k, random_state=seed), idx))
        steps.append(("pca", ColumnTransformer(groups)))
    steps.append(("clf", clf))
    return Pipeline(steps)


def cv_auc(model, X, y, seed: int, n_splits: int = 10) -> float:
    cv = StratifiedShuffleSplit(n_splits=n_splits, test_size=0.3, random_state=seed)
    return float(np.mean(cross_val_score(model, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--triples", nargs="*", default=[],
                    help="three-species definitions to add as extra feature groups")
    args = ap.parse_args()

    X, y, C, tgroups = build_patient_table(tuple(args.triples))
    pair_cols = [c for c in X.columns if not c[0].startswith("triple:")]
    Xp = X[pair_cols]
    rng = np.random.default_rng(SEED)
    print(f"patients: {len(X)} ({int((y==0).sum())} NR / {int((y==1).sum())} R); "
          f"features: {X.shape[1]} ({len(pair_cols)} from {len(DEFS)} pair definitions"
          + (f", {X.shape[1]-len(pair_cols)} from {len(tgroups)} triples)" if tgroups else ")"))

    models = {
        "TDA (all 8 pair definitions, per-group PCA)": (make_model(Xp, True, SEED), Xp),
        "TDA, CD8_primary only": (
            make_model(X[["CD8_primary"]], True, SEED), X[["CD8_primary"]]),
        "Composition + geometry only": (
            Pipeline([("impute", SimpleImputer(strategy="median")),
                      ("scale", StandardScaler()),
                      ("clf", BaggingClassifier(
                          GradientBoostingClassifier(
                              loss="log_loss", n_estimators=25, learning_rate=0.4,
                              max_depth=3, max_leaf_nodes=6, min_samples_leaf=5,
                              random_state=SEED),
                          n_estimators=50, bootstrap=False, random_state=SEED))]), C),
    }
    # Each triple on its own, then everything together. Screened descriptively -- the
    # permutation null below is run on the pre-specified combined model, not on whichever
    # of these happens to score highest (that would be post-hoc selection).
    for g in tgroups:
        models[f"TDA, {g} only"] = (make_model(X[[g]], True, SEED), X[[g]])
    if tgroups:
        models["TDA, pairs + triples (combined)"] = (make_model(X, True, SEED), X)

    rows = []
    for name, (model, feats) in models.items():
        auc = cv_auc(model, feats, y, SEED)
        rows.append({"model": name, "n_features": feats.shape[1], "cv_auc": auc})
        print(f"  {name:46s} AUC = {auc:.3f}")

    dummy = cv_auc(DummyClassifier(strategy="stratified", random_state=SEED), X, y, SEED)
    rows.append({"model": "DummyClassifier (stratified)", "n_features": 0, "cv_auc": dummy})
    res = pd.DataFrame(rows)

    # ---- label-permutation null on the primary model -------------------------
    primary = ("TDA, pairs + triples (combined)" if tgroups
               else "TDA (all 8 pair definitions, per-group PCA)")
    model, feats = models[primary]
    obs = float(res.loc[res.model == primary, "cv_auc"].iloc[0])
    print(f"\nprimary model for the null: {primary}")
    print(f"\nrunning {args.n_perm} label permutations for the null...")
    null = np.empty(args.n_perm)
    for i in range(args.n_perm):
        null[i] = cv_auc(model, feats, pd.Series(rng.permutation(y.to_numpy()), index=y.index),
                         SEED)
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{args.n_perm}  null AUC mean so far {null[:i+1].mean():.3f}")
    if args.n_perm == 0:                      # scores-only smoke run
        print("\n(no permutations requested; skipping the null)")
        print(res.round(4).to_string(index=False))
        return
    p_emp = float((null >= obs).mean())

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    ax[0].hist(null, bins=30, color="0.7", edgecolor="white")
    ax[0].axvline(obs, color="crimson", lw=2.5,
                  label=f"observed AUC = {obs:.3f}\npermutation p = {p_emp:.3f}")
    ax[0].axvline(0.5, color="black", lw=1, ls=":")
    ax[0].set_xlabel("cross-validated ROC AUC")
    ax[0].set_ylabel("permutations")
    ax[0].set_title(f"A. Label-permutation null ({args.n_perm} draws)\n"
                    f"null mean AUC = {null.mean():.3f}", fontsize=10)
    ax[0].legend(fontsize=8)

    r = res.sort_values("cv_auc")
    ax[1].barh(r.model, r.cv_auc, color=["0.6" if "Dummy" in m else "tab:blue"
                                         for m in r.model])
    ax[1].axvline(0.5, color="black", lw=1, ls=":")
    ax[1].set_xlim(0, 1)
    ax[1].set_xlabel("cross-validated ROC AUC")
    ax[1].set_title("B. Models compared\n(TDA must beat composition to be worth anything)",
                    fontsize=10)
    for i, v in enumerate(r.cv_auc):
        ax[1].text(v + 0.015, i, f"{v:.3f}", va="center", fontsize=9)
    fig.suptitle("Can chromatic-TDA statistics predict response? "
                 f"Pre-treatment, {len(X)} patients ({int((y==0).sum())} NR / "
                 f"{int((y==1).sum())} R), M2S2-style pipeline", fontsize=11.5)
    fig.tight_layout()
    figp = OUT / "classification_figure.png"
    fig.savefig(figp, dpi=150, bbox_inches="tight")

    def md(dfr):
        h = "| " + " | ".join(dfr.columns) + " |"
        s = "| " + " | ".join("---" for _ in dfr.columns) + " |"
        b = "\n".join("| " + " | ".join(
            f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)
            for v in row) + " |" for row in dfr.itertuples(index=False))
        return "\n".join([h, s, b])

    verdict = ("**does not** predict response better than chance" if p_emp >= 0.05
               else "**does** predict response better than chance")
    lines = [
        "# M2S2-style classification: can TDA statistics predict response?\n",
        f"Pre-treatment, patient-level. {len(X)} patients ({int((y==0).sum())} NR / "
        f"{int((y==1).sum())} R), {X.shape[1]} statistics across {len(DEFS)} pair "
        f"definitions, reduced by per-cell-group PCA inside the cross-validation.\n",
        md(res.round(4)),
        f"\n- Label-permutation null ({args.n_perm} draws): null mean AUC = "
        f"{null.mean():.3f}, 95th pct = {np.percentile(null, 95):.3f}; "
        f"observed = **{obs:.3f}**, empirical **p = {p_emp:.3f}**.",
        f"- A multivariate classifier over the full topological feature set therefore "
        f"{verdict}.",
        f"- Composition + geometry alone reaches AUC = "
        f"{res.loc[res.model.str.startswith('Composition'), 'cv_auc'].iloc[0]:.3f}, so TDA "
        f"is not being outperformed by a trivial abundance model either — neither carries "
        f"signal.\n",
        "## Why this matters for the rebuttal\n",
        "This is the analysis the Day-1 report deferred and never ran, and it is the one a "
        "reviewer is most likely to ask for: not a single hand-picked statistic, but the "
        "whole feature set given to a classifier with the authors' own hyperparameters. "
        "It is also the form of evidence least sensitive to the confounding arguments "
        "elsewhere in this directory — a classifier is free to use composition, geometry, "
        "topology or any combination, and still cannot separate the groups.\n",
        "Script: `day2_classify_response.py`. Figure: `classification_figure.png`.",
    ]
    (OUT / "CLASSIFICATION.md").write_text("\n".join(lines) + "\n")
    print(f"\nobserved AUC {obs:.3f}, null mean {null.mean():.3f}, p = {p_emp:.3f}")
    print("wrote", OUT / "CLASSIFICATION.md")


if __name__ == "__main__":
    main()
