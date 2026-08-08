#!/usr/bin/env python
"""Stage 2 of the CLAHE batch-correction evaluation: paired intensity extraction.

Builds two cell tables holding the SAME cells, in the same order, with the same
cell type labels and coordinates, differing only in whether CLAHE was applied
to the images the intensities came from:

    adata_corrected.h5ad     intensities from `processed/`     (CLAHE applied)
    adata_uncorrected.h5ad   intensities from `non_processed/` (denoised only)

Segmentation is shared by construction -- the same Mesmer whole-cell masks are
used for both -- so no segmentation difference can leak into the comparison.

Typical use:

    # 1. Pilot: verify the extraction reproduces the published cell table,
    #    then do a handful of ROIs end to end. Do this first.
    python 02_extract.py --verify-only
    python 02_extract.py --limit 5

    # 2. Full run. Resumable: per-FOV results are cached, so re-running skips
    #    what is already done.
    python 02_extract.py

Why the verification step exists
--------------------------------
The paired design assumes the two tables differ ONLY by CLAHE. That holds only
if this script extracts intensities the same way the original pipeline did. If
it used, say, a different per-cell normalisation, then the "corrected" table
here would differ from the published one in a second way, and any metric
difference could no longer be attributed to CLAHE. So --verify-only extracts
from `processed/` and checks it reproduces `updated_cell_table_Ki67.csv`,
testing the candidate definitions of "size normalized" against each other.
Nothing else runs until that passes.

Outputs (all under --out-dir):
    cache/<fov>.npz          per-FOV extracted means, so a rerun is cheap
    adata_corrected.h5ad     .X = published intensities, 37 shared markers
    adata_uncorrected.h5ad   .X = intensities extracted from non_processed/
    extract_report.md        what happened, including the verification result
"""

import argparse
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# CONFIG -- edit here, or override on the command line
# --------------------------------------------------------------------------

IMG_ROOT = Path("/mnt/data/Delta_Tissue/IMC/Img_Denoised/sept2024_release")
PROCESSED_DIR = IMG_ROOT / "processed"          # CLAHE applied
NON_PROCESSED_DIR = IMG_ROOT / "non_processed"  # denoised only

# Identified by stage 1 check 5a as the segmentation that produced the cell
# table (labels and centroids reproduced to 1e-13 px). Do not change without
# re-running that check: three other directories on this machine look
# plausible and are wrong by ~1000 px.
MASK_DIR = Path("/mnt/data/Delta_Tissue/IMC/segmentation/sept2024_release/deepcell_output")
MASK_SUFFIX = "_whole_cell.tiff"

CELL_TABLE = Path("/mnt/data/Delta_Tissue/IMC/CleanCohort/updated_cell_table_Ki67.csv")
METADATA = Path("/home/linus/CleanCohort_Metadata.csv")

OUT_DIR = Path("./stage2_output")

# Carboplatin is derived during preprocessing, not measured, so it exists only
# in `processed/`. Dropped from BOTH tables to keep them like-for-like. The
# pre-treatment cohort is all CORE, where the pipeline zeroes it anyway.
DERIVED_CHANNELS = ["Carboplatin"]

# Cohort bookkeeping. Everything is extracted; these only set .obs flags so the
# cohort choice stays a stage-3 decision (see stage 1 report, point 8).
EXCLUDED_PATIENTS = ["LEAP149", "LEAP150"]
EXCLUDED_ROIS = ["Leap005_2_1", "Leap005_2_2", "Leap010_7", "Leap094_7", "Leap095_13"]
PRE_TREATMENT_VALUES = ["pre"]

BATCH_COL = "Stain_Batch"
SAMPLE_TYPE_COL = "Sample_Type_(pre/post treatment)"

N_VERIFY_FOVS = 3
VERIFY_MIN_R = 0.999      # per-channel Pearson r required against the published table
VERIFY_MAX_MEDIAN_REL = 1e-3

QUANTILE = 0.95           # matches phenotyping_utils.normalise()


def log_mem(msg, t0):
    """Peak RSS so far, plus elapsed. ru_maxrss is KB on Linux, bytes on macOS."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    gb = peak / (1024 ** 2) if sys.platform != "darwin" else peak / (1024 ** 3)
    print(f"[{time.time() - t0:7.1f}s | peak RSS {gb:6.2f} GB] {msg}", flush=True)


# --------------------------------------------------------------------------
# channels
# --------------------------------------------------------------------------

def marker_columns(header):
    """Marker columns are those between `cell_size` and `label` in the table."""
    cols = list(header)
    return cols[cols.index("cell_size") + 1: cols.index("label")]


def resolve_channels(rep, table_markers):
    """The channel list both directories and the cell table agree on."""
    proc_dirs = [d for d in PROCESSED_DIR.iterdir() if d.is_dir()]
    raw_dirs = [d for d in NON_PROCESSED_DIR.iterdir() if d.is_dir()]
    proc = {f.stem for f in proc_dirs[0].glob("*.tiff")}
    raw = {f.stem for f in raw_dirs[0].glob("*.tiff")}

    shared = sorted((proc & raw & set(table_markers)) - set(DERIVED_CHANNELS))
    rep(f"- markers in cell table: {len(table_markers)}")
    rep(f"- channels on disk: {len(proc)} processed, {len(raw)} non_processed")
    rep(f"- derived, dropped from both: {DERIVED_CHANNELS}")
    rep(f"- **shared channel list: {len(shared)}**")

    missing = set(table_markers) - set(shared) - set(DERIVED_CHANNELS)
    if missing:
        raise SystemExit(f"Cell table markers with no TIFF in both directories: "
                         f"{sorted(missing)}. Resolve before extracting.")
    return shared


# --------------------------------------------------------------------------
# extraction
# --------------------------------------------------------------------------

def extract_fov(fov, labels, channels, img_dir, tifffile, ndi):
    """Mean intensity per (cell, channel) for one FOV.

    Args:
        labels: mask object IDs, IN THE CELL TABLE'S ROW ORDER. Indexing the
            output by these keeps both tables cell-for-cell aligned without
            relying on any sort order.

    Returns (n_cells, n_channels) float32, and the per-cell pixel area.
    """
    mask = tifffile.imread(MASK_DIR / f"{fov}{MASK_SUFFIX}")
    out = np.zeros((len(labels), len(channels)), dtype=np.float32)

    for j, ch in enumerate(channels):
        img = tifffile.imread(img_dir / fov / f"{ch}.tiff")
        if img.shape != mask.shape:
            raise SystemExit(f"{fov}/{ch}: image {img.shape} vs mask {mask.shape}")
        out[:, j] = ndi.mean(img, labels=mask, index=labels)

    area = np.asarray(ndi.sum(np.ones_like(mask, dtype=np.float32),
                              labels=mask, index=labels))
    return out, area


def verify_against_published(rep, table, channels, tifffile, ndi):
    """Check this script's extraction reproduces the published cell table.

    Tests two candidate readings of "size normalized": mean intensity over the
    mask (sum/area), and sum divided by the table's own `cell_size`. Whichever
    matches is the pipeline's definition; if neither does, extraction differs
    from the original in some other way and stage 2 must not proceed.
    """
    rep.h("Verification: does this extraction reproduce the published table?")

    fovs = list(dict.fromkeys(table["fov"]))[:N_VERIFY_FOVS]
    rows = []
    for fov in fovs:
        sub = table[table["fov"] == fov]
        labels = sub["label"].to_numpy()
        got, area = extract_fov(fov, labels, channels, PROCESSED_DIR, tifffile, ndi)
        published = sub[channels].to_numpy(dtype=np.float64)

        # Candidate definitions of the published quantity.
        candidates = {
            "mean over mask (sum/area)": got,
            "sum / cell_size column": got * area[:, None] / sub["cell_size"].to_numpy()[:, None],
        }
        for name, cand in candidates.items():
            with np.errstate(invalid="ignore", divide="ignore"):
                rs, rels = [], []
                for j in range(len(channels)):
                    rs.append(np.corrcoef(cand[:, j], published[:, j])[0, 1])
                    col = published[:, j]
                    denom = np.where(np.abs(col) > 0, np.abs(col), np.nan)
                    # Per channel, then worst channel. Correlation alone is
                    # affine-invariant, so a single rescaled channel would slip
                    # past r; and a median pooled over all channels would be
                    # outvoted by the channels that are fine.
                    rels.append(np.nanmedian(np.abs(cand[:, j] - col) / denom))
            rows.append({"fov": fov, "definition": name,
                         "min_channel_r": float(np.nanmin(rs)),
                         "worst_channel_rel_err": float(np.nanmax(rels)),
                         "worst_channel": channels[int(np.nanargmax(rels))]})

    df = pd.DataFrame(rows)
    rep.table(df)

    summary = df.groupby("definition").agg(
        min_r=("min_channel_r", "min"),
        worst_rel_err=("worst_channel_rel_err", "max"))
    rep("\nPer candidate definition (worst channel, worst FOV):\n")
    rep.table(summary)

    ok = summary[(summary["min_r"] >= VERIFY_MIN_R)
                 & (summary["worst_rel_err"] <= VERIFY_MAX_MEDIAN_REL)]
    if len(ok) == 0:
        rep.fail(
            "This script's extraction does not reproduce the published cell table "
            "under either candidate definition. The two tables would then differ by "
            "more than CLAHE and the paired design would be invalid. Do not run the "
            "full extraction; work out what the original pipeline did differently.")
        return None
    chosen = ok.index[0]
    rep(f"\n**Verified**: `{chosen}` reproduces the published intensities "
        f"(min per-channel r {ok.iloc[0]['min_r']:.6f}, worst-channel relative error "
        f"{ok.iloc[0]['worst_rel_err']:.2e}). The uncorrected table will be built "
        f"with the identical procedure, so CLAHE is the only difference between them.")
    return chosen


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

def build_obs(table, meta):
    """Cell-level metadata, shared verbatim by both tables."""
    obs = pd.DataFrame(index=pd.Index(
        [f"{f}_{l}" for f, l in zip(table["fov"], table["label"])], name="cell_id"))
    obs["acquisition_ID"] = table["fov"].values
    obs["label"] = table["label"].values
    obs["Pixie"] = pd.Categorical(table["cell_meta_cluster"].astype(str).values)
    obs["cell_size"] = table["cell_size"].values
    obs["LEAP_ID"] = table["fov"].str.extract(r"^(Leap\d+)")[0].str.upper().values

    joined = obs.reset_index().merge(meta, on="LEAP_ID", how="left").set_index("cell_id")
    joined = joined.loc[obs.index]

    joined["is_pre_treatment"] = joined[SAMPLE_TYPE_COL].isin(PRE_TREATMENT_VALUES)
    joined["in_revision_cohort"] = (~joined["LEAP_ID"].isin(EXCLUDED_PATIENTS)
                                    & ~joined["acquisition_ID"].isin(EXCLUDED_ROIS))
    return joined


def transform_like_pipeline(adata, sc, preprocessing):
    """Reproduce phenotyping_utils.generate_anndata_from_cell_table's chain.

    NaN -> 0, log1p, per-channel quantile normalisation capped at 1, then
    StandardScaler into layers['scaled']. Applied identically to both tables;
    each computes its own quantiles, which is what "same transformation" means
    for a procedure whose parameters are estimated from the data.
    """
    adata.layers["raw"] = adata.X.copy()
    adata.X = np.nan_to_num(adata.X, nan=0.0)
    sc.pp.log1p(adata)

    q = np.nanquantile(adata.X, QUANTILE, axis=0)
    adata.X = np.divide(adata.X, q, out=np.zeros_like(adata.X), where=q != 0)
    adata.X[adata.X > 1] = 1
    adata.uns["normalisation_quantiles"] = q

    adata.layers["scaled"] = preprocessing.StandardScaler().fit_transform(adata.X)
    return adata


def assert_paired(a, b):
    """The two tables must be the same cells in the same order. Assert it."""
    assert a.shape == b.shape, f"shape {a.shape} vs {b.shape}"
    assert list(a.obs_names) == list(b.obs_names), "cell IDs differ or are reordered"
    assert list(a.var_names) == list(b.var_names), "channel order differs"
    assert (a.obs["acquisition_ID"].values == b.obs["acquisition_ID"].values).all()
    assert (a.obs["label"].values == b.obs["label"].values).all()
    assert (a.obs["Pixie"].astype(str).values == b.obs["Pixie"].astype(str).values).all()
    assert np.array_equal(a.obsm["spatial"], b.obsm["spatial"])


# --------------------------------------------------------------------------

def main():
    global CELL_TABLE, METADATA, MASK_DIR

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--cell-table", type=Path, default=CELL_TABLE)
    ap.add_argument("--metadata", type=Path, default=METADATA)
    ap.add_argument("--mask-dir", type=Path, default=None)
    ap.add_argument("--limit", type=int, default=None,
                    help="Extract only the first N FOVs (pilot run)")
    ap.add_argument("--fovs", nargs="*", default=None, help="Restrict to these FOVs")
    ap.add_argument("--verify-only", action="store_true",
                    help="Run the reproduction check against the published table and stop")
    ap.add_argument("--skip-verify", action="store_true",
                    help="Skip the check. Only for reruns where it already passed.")
    args = ap.parse_args()

    CELL_TABLE, METADATA = args.cell_table, args.metadata
    if args.mask_dir:
        MASK_DIR = args.mask_dir

    import anndata as ad
    import scanpy as sc
    import tifffile
    from scipy import ndimage as ndi
    from sklearn import preprocessing

    # Reuse stage 1's report helper so both stages produce the same shape of output.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "reconcile", Path(__file__).resolve().parent / "01_reconcile.py")
    reconcile = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reconcile)

    t0 = time.time()
    out_dir = args.out_dir
    (out_dir / "cache").mkdir(parents=True, exist_ok=True)

    rep = reconcile.Report()
    rep("# CLAHE batch-correction evaluation -- stage 2 extraction\n")
    rep(f"Host `{os.uname().nodename}`, masks `{MASK_DIR}`.\n")
    rep.h("Channels")

    header = pd.read_csv(CELL_TABLE, nrows=0).columns
    table_markers = marker_columns(header)
    channels = resolve_channels(rep, table_markers)

    keep_cols = ["fov", "label", "cell_size", "centroid-0", "centroid-1",
                 "cell_meta_cluster"] + channels
    log_mem("reading cell table (selected columns)", t0)
    table = pd.read_csv(CELL_TABLE, usecols=keep_cols)
    table = table[keep_cols]
    log_mem(f"cell table loaded: {table.shape}", t0)

    if not args.skip_verify:
        chosen = verify_against_published(rep, table, channels, tifffile, ndi)
        if chosen is None:
            rep.write(out_dir / "extract_report.md")
            return 1
        if args.verify_only:
            rep.write(out_dir / "extract_report.md")
            print("\nVerification passed. Re-run without --verify-only to extract.")
            return 0

    fovs = list(dict.fromkeys(table["fov"]))
    if args.fovs:
        fovs = [f for f in fovs if f in set(args.fovs)]
    if args.limit:
        fovs = fovs[: args.limit]
    table = table[table["fov"].isin(set(fovs))].reset_index(drop=True)
    # Group rows by FOV while preserving the table's own ordering within each.
    table = table.sort_values("fov", kind="stable").reset_index(drop=True)

    rep.h("Extraction")
    rep(f"- FOVs to process: {len(fovs)}")
    rep(f"- cells: {len(table):,}")

    uncorrected = np.zeros((len(table), len(channels)), dtype=np.float32)
    offsets = {}
    start = 0
    for fov, sub in table.groupby("fov", sort=False):
        offsets[fov] = (start, start + len(sub))
        start += len(sub)

    n_cached = 0
    for i, (fov, (lo, hi)) in enumerate(offsets.items(), 1):
        cache = out_dir / "cache" / f"{fov}.npz"
        labels = table["label"].to_numpy()[lo:hi]

        if cache.exists():
            stored = np.load(cache, allow_pickle=False)
            if (stored["labels"] == labels).all() and stored["means"].shape[1] == len(channels):
                uncorrected[lo:hi] = stored["means"]
                n_cached += 1
                continue

        means, _ = extract_fov(fov, labels, channels, NON_PROCESSED_DIR, tifffile, ndi)
        uncorrected[lo:hi] = means
        np.savez_compressed(cache, means=means, labels=labels,
                            channels=np.array(channels))
        if i % 25 == 0 or i == len(offsets):
            log_mem(f"extracted {i}/{len(offsets)} FOVs", t0)

    rep(f"- reused from cache: {n_cached}")
    rep(f"- newly extracted: {len(offsets) - n_cached}")

    # A table row whose label is absent from the mask yields NaN here. The
    # transform would silently turn that into 0, so count it instead.
    n_nan = int(np.isnan(uncorrected).any(axis=1).sum())
    rep(f"- cells with no matching mask object (NaN): {n_nan:,}")
    if n_nan:
        bad = table.loc[np.isnan(uncorrected).any(axis=1), "fov"].value_counts()
        rep(f"\nAffected FOVs:\n\n```\n{bad.head(20).to_string()}\n```")
        rep.fail(f"{n_nan:,} cells in the cell table have no object in the mask. "
                 "They would become zeros downstream. Investigate before trusting "
                 "the extracted table.")
    log_mem("extraction complete", t0)

    rep.h("Assembly")
    meta = pd.read_csv(METADATA)
    obs = build_obs(table, meta)
    spatial = table[["centroid-0", "centroid-1"]].to_numpy()

    adata_unc = ad.AnnData(X=uncorrected.astype(np.float32), obs=obs.copy(),
                           var=pd.DataFrame(index=pd.Index(channels, name="marker")),
                           obsm={"spatial": spatial.copy()})
    adata_cor = ad.AnnData(X=table[channels].to_numpy(dtype=np.float32), obs=obs.copy(),
                           var=pd.DataFrame(index=pd.Index(channels, name="marker")),
                           obsm={"spatial": spatial.copy()})

    # Both AnnData objects now hold their own copy of the intensities; the
    # source frame and array are the largest things still alive.
    del table, uncorrected
    log_mem("assembled, source frame released", t0)

    assert_paired(adata_cor, adata_unc)
    rep("- **asserted**: identical cell IDs, order, labels, Pixie labels, "
        "coordinates and channel order between the two tables")

    for name, adata in [("corrected", adata_cor), ("uncorrected", adata_unc)]:
        transform_like_pipeline(adata, sc, preprocessing)
        rep(f"- {name}: transformed (nan->0, log1p, q{QUANTILE} per-channel cap, scaled)")

    assert_paired(adata_cor, adata_unc)

    n_diff = int((adata_cor.layers["raw"] != adata_unc.layers["raw"]).any(axis=1).sum())
    rep(f"- cells whose raw intensities differ between the tables: {n_diff:,} of "
        f"{adata_cor.n_obs:,} ({100 * n_diff / adata_cor.n_obs:.1f}%)")

    rep.h("Cohort flags recorded in .obs")
    flags = pd.DataFrame({
        "all": [adata_cor.n_obs, adata_cor.obs["acquisition_ID"].nunique()],
        "revision cohort": [int(adata_cor.obs["in_revision_cohort"].sum()),
                            adata_cor.obs.loc[adata_cor.obs["in_revision_cohort"],
                                              "acquisition_ID"].nunique()],
        "revision + pre-treatment": [
            int((adata_cor.obs["in_revision_cohort"] & adata_cor.obs["is_pre_treatment"]).sum()),
            adata_cor.obs.loc[adata_cor.obs["in_revision_cohort"]
                              & adata_cor.obs["is_pre_treatment"], "acquisition_ID"].nunique()],
    }, index=["cells", "ROIs"])
    rep.table(flags)
    rep("\nNothing is filtered here. Stage 3 selects the cohort from these flags.")

    adata_cor.write(out_dir / "adata_corrected.h5ad")
    adata_unc.write(out_dir / "adata_uncorrected.h5ad")
    log_mem("written", t0)
    rep(f"\nWritten to `{out_dir}`.")

    rep.write(out_dir / "extract_report.md")
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
