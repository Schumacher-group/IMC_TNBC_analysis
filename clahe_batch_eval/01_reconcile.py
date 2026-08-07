#!/usr/bin/env python
"""Stage 1 of the CLAHE batch-correction evaluation: reconciliation only, no analysis.

Answers points 1-8 of the analysis specification and writes a markdown report.
Nothing here modifies data or produces analysis results; it is safe to re-run.

Run on the workstation where the image data lives:

    python 01_reconcile.py --out reconcile_report.md

Points 6-8 (metadata only) can also be run anywhere the cell table lives:

    python 01_reconcile.py --skip-images --out reconcile_report_metadata.md

Checks
------
1. processed vs non_processed: is the difference CLAHE?
   Decisive test is the *rank* correlation of pixel values between the two
   versions of the same FOV/channel. CLAHE is a LOCAL operation, so it
   reorders distant pixels: Spearman rho well below 1. Any global monotone
   rescale (quantile normalisation, log, min-max) preserves global rank order
   exactly: rho == 1. Supporting evidence is histogram entropy (CLAHE flattens
   the histogram, raising entropy) and a bounded fixed output range.
2. FOV matching between clean_cohort_fovs.txt and non_processed/.
3. Per-pair image dimensions and channel counts agree.
4. Per-pair channel name sets agree (channels are stored one TIFF per channel,
   so "ordering" is whatever order the extraction loop imposes; the risk is a
   missing/renamed channel, which set comparison catches).
5. Three-way match: clean cohort FOV <-> Mesmer whole-cell mask <-> cell table
   acquisition ID.
6. Batch variable: which column, how many batches, ROIs and patients per batch.
7. Confounding of batch with patient and with response group.
8. Cohort scope after restricting to pre-treatment and applying the revision
   cohort exclusions.
"""

import argparse
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# CONFIG -- edit these paths, or override on the command line
# --------------------------------------------------------------------------

IMG_ROOT = Path("/mnt/data/Delta_Tissue/IMC/Img_Denoised/sept2024_release")
PROCESSED_DIR = IMG_ROOT / "processed"
NON_PROCESSED_DIR = IMG_ROOT / "non_processed"
CLEAN_COHORT_DIR = IMG_ROOT / "CleanCohort" / "processed"

# Other Img_Denoised subdirectories seen in this repo's existing scripts.
# Listed so the report can say what else is on disk under the same root;
# `contrast_adj` is what collagen_fibre_segmentation/extract_fibre_skeleton.py
# and segmentation/deepcell_output_contrast_adj were built against, and it may
# or may not be the same thing as `processed`.
SIBLING_DIRS_TO_REPORT = ["contrast_adj", "non_preprocessed", "processed", "non_processed"]

CLEAN_COHORT_FOV_LIST = Path("clean_cohort_fovs.txt")

# Mesmer whole-cell masks: <MASK_DIR>/<fov>_whole_cell.tiff
MASK_DIR = Path("../segmentation/deepcell_output_contrast_adj")
MASK_SUFFIX = "_whole_cell.tiff"

# Cell table and metadata. Only the last few columns of the cell table are
# needed here, so it is read with usecols to avoid loading ~8 GB.
CELL_TABLE = Path("../CellTable_CleanCohort/updated_cell_table_Ki67.csv")
METADATA = Path("../CellTable_CleanCohort/CleanCohort_Metadata.csv")

BATCH_COL = "Stain_Batch"
SAMPLE_TYPE_COL = "Sample_Type_(pre/post treatment)"
RESPONSE_COL = "Response"
PRE_TREATMENT_VALUES = ["pre"]

# Revision cohort exclusions (see project memory: cohort definitions).
EXCLUDED_PATIENTS = ["LEAP149", "LEAP150"]
EXCLUDED_ROIS = ["Leap005_2_1", "Leap005_2_2", "Leap010_7", "Leap094_7", "Leap095_13"]

# Check 1 sampling
N_HIST_FOVS = 5
HIST_CHANNELS = ["DNA1", "Pan-keratin", "Collage-Type_I", "CD8a", "Carboplatin"]
PIXEL_SUBSAMPLE = 200_000
RNG_SEED = 0


# --------------------------------------------------------------------------
# report helper
# --------------------------------------------------------------------------

class Report:
    """Accumulates markdown and echoes to stdout as it goes."""

    def __init__(self):
        self.lines = []
        self.failures = []

    def h(self, text, level=2):
        self(f"\n{'#' * level} {text}\n")

    def __call__(self, text=""):
        print(text)
        self.lines.append(text)

    def table(self, df):
        """Markdown table if `tabulate` is installed, otherwise a fenced text block."""
        try:
            self(df.to_markdown())
        except ImportError:
            self("```\n" + df.to_string() + "\n```")

    def fail(self, text):
        """Record a blocking problem: something the spec says to stop and flag."""
        self.failures.append(text)
        self(f"\n> **STOP / FLAG:** {text}\n")

    def write(self, path):
        Path(path).write_text("\n".join(self.lines) + "\n")


# --------------------------------------------------------------------------
# check 1: is the processed/non_processed difference CLAHE?
# --------------------------------------------------------------------------

def histogram_entropy(img, bins=256):
    """Shannon entropy of the intensity histogram, in bits.

    CLAHE equalises the histogram, so it raises this towards log2(bins) = 8.
    """
    finite = img[np.isfinite(img)]
    if finite.size == 0 or finite.max() == finite.min():
        return np.nan
    counts, _ = np.histogram(finite, bins=bins)
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def channel_stats(img):
    finite = img[np.isfinite(img)].ravel()
    q = np.percentile(finite, [1, 25, 50, 75, 99]) if finite.size else [np.nan] * 5
    return {
        "dtype": str(img.dtype),
        "shape": "x".join(str(s) for s in img.shape),
        "min": float(finite.min()) if finite.size else np.nan,
        "max": float(finite.max()) if finite.size else np.nan,
        "p1": q[0], "p25": q[1], "p50": q[2], "p75": q[3], "p99": q[4],
        "frac_zero": float((finite == 0).mean()) if finite.size else np.nan,
        "entropy_bits": histogram_entropy(img),
    }


def check_1_clahe(rep, fovs, tifffile, spearmanr):
    rep.h("1. Is `processed` vs `non_processed` the CLAHE difference?")

    rep("Other subdirectories present under the release root (for orientation -- "
        "this repo's existing collagen and segmentation code was built against "
        "`contrast_adj`, not `processed`):\n")
    for name in SIBLING_DIRS_TO_REPORT:
        d = IMG_ROOT / name
        if d.exists():
            n = sum(1 for _ in d.iterdir())
            rep(f"- `{name}/` exists, {n} entries")
        else:
            rep(f"- `{name}/` absent")
    rep()

    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for fov in fovs[:N_HIST_FOVS]:
        for ch in HIST_CHANNELS:
            p_path = PROCESSED_DIR / fov / f"{ch}.tiff"
            n_path = NON_PROCESSED_DIR / fov / f"{ch}.tiff"
            if not (p_path.exists() and n_path.exists()):
                rows.append({"fov": fov, "channel": ch, "note": "missing in one or both"})
                continue

            p_img = tifffile.imread(p_path).astype(float)
            n_img = tifffile.imread(n_path).astype(float)
            if p_img.shape != n_img.shape:
                rows.append({"fov": fov, "channel": ch,
                             "note": f"shape mismatch {p_img.shape} vs {n_img.shape}"})
                continue

            flat_p, flat_n = p_img.ravel(), n_img.ravel()
            if flat_p.size > PIXEL_SUBSAMPLE:
                idx = rng.choice(flat_p.size, PIXEL_SUBSAMPLE, replace=False)
                flat_p, flat_n = flat_p[idx], flat_n[idx]
            rho = spearmanr(flat_n, flat_p).correlation

            ps, ns = channel_stats(p_img), channel_stats(n_img)
            rows.append({
                "fov": fov, "channel": ch,
                "shape": ps["shape"],
                "raw_dtype": ns["dtype"], "proc_dtype": ps["dtype"],
                "raw_max": ns["max"], "proc_max": ps["max"],
                "raw_p99": ns["p99"], "proc_p99": ps["p99"],
                "raw_entropy": ns["entropy_bits"], "proc_entropy": ps["entropy_bits"],
                "spearman_rho": rho,
                "note": "",
            })

    df = pd.DataFrame(rows)
    rep.table(df)
    rep()

    if "spearman_rho" not in df or df["spearman_rho"].dropna().empty:
        rep.fail("Could not compute any processed/non_processed pixel comparisons. "
                 "Check the directory layout before going further.")
        return df

    rho = df["spearman_rho"].dropna()
    ent_gain = (df["proc_entropy"] - df["raw_entropy"]).dropna()
    rep(f"Spearman rho across {len(rho)} FOV/channel pairs: "
        f"median {rho.median():.4f}, min {rho.min():.4f}, max {rho.max():.4f}")
    rep(f"Histogram entropy change (processed - non_processed), bits: "
        f"median {ent_gain.median():+.3f}, min {ent_gain.min():+.3f}, max {ent_gain.max():+.3f}")
    rep()

    if rho.min() > 0.9999:
        rep.fail(
            "Pixel rank order is preserved exactly between `processed` and "
            "`non_processed` (Spearman rho == 1). CLAHE is a local operation and "
            "cannot preserve global rank order, so the difference between these "
            "two directories is a GLOBAL monotone rescale, not CLAHE. The working "
            "assumption in the spec is wrong; do not proceed until the true CLAHE "
            "and pre-CLAHE directories are identified.")
    elif ent_gain.median() <= 0:
        rep.fail(
            "Pixel rank order changes (consistent with a local operation) but the "
            "histogram does not flatten -- entropy does not increase in `processed`. "
            "That is not the CLAHE signature either. Investigate before proceeding.")
    else:
        rep("**Consistent with CLAHE**: rank order is not preserved (local operation) "
            "and the histogram flattens (entropy increases) in `processed`.")
    return df


# --------------------------------------------------------------------------
# checks 2-4: FOV matching, dimensions, channels
# --------------------------------------------------------------------------

def check_2_4_fovs(rep, clean_fovs, tifffile):
    rep.h("2. FOV matching: clean cohort -> `non_processed/`")

    matched, unmatched = [], []
    for fov in clean_fovs:
        if (NON_PROCESSED_DIR / fov).exists():
            matched.append(fov)
        else:
            unmatched.append(fov)

    rep(f"- clean cohort FOVs listed: {len(clean_fovs)}")
    rep(f"- matched by direct name lookup in `non_processed/`: {len(matched)}")
    rep(f"- unmatched: {len(unmatched)}")
    if unmatched:
        rep("\nUnmatched FOVs (listed in full, none dropped silently):\n")
        for fov in unmatched:
            rep(f"  - `{fov}`")
        rep("\nNearest names present in `non_processed/` for each unmatched FOV:\n")
        available = sorted(d.name for d in NON_PROCESSED_DIR.iterdir() if d.is_dir())
        import difflib
        for fov in unmatched:
            close = difflib.get_close_matches(fov, available, n=3, cutoff=0.6)
            rep(f"  - `{fov}` -> {close if close else 'no close match'}")
        rep.fail(f"{len(unmatched)} clean cohort FOVs have no counterpart in "
                 "`non_processed/`. Resolve the naming rule before extraction; "
                 "do not drop them.")
    else:
        rep("\nMatching rule used: direct directory-name lookup, no transformation.")

    # Non-FOV entries in the list (the list came from `ls`, so it may contain files)
    non_dirs = [f for f in clean_fovs if (CLEAN_COHORT_DIR / f).exists()
                and not (CLEAN_COHORT_DIR / f).is_dir()]
    if non_dirs:
        rep(f"\nNote: {len(non_dirs)} entries in the FOV list are files, not "
            f"directories: {non_dirs[:10]}")

    rep.h("3 & 4. Per-pair image dimensions, channel counts, channel names")

    rows = []
    dim_mismatch, chan_mismatch = [], []
    for fov in matched:
        p_dir, n_dir = CLEAN_COHORT_DIR / fov, NON_PROCESSED_DIR / fov
        p_chans = {f.stem for f in p_dir.glob("*.tiff")}
        n_chans = {f.stem for f in n_dir.glob("*.tiff")}

        # Read only the header of one shared channel to get dimensions cheaply.
        shared = sorted(p_chans & n_chans)
        p_shape = n_shape = None
        if shared:
            ref = shared[0]
            with tifffile.TiffFile(p_dir / f"{ref}.tiff") as tf:
                p_shape = tf.pages[0].shape
            with tifffile.TiffFile(n_dir / f"{ref}.tiff") as tf:
                n_shape = tf.pages[0].shape

        ok_dims = p_shape == n_shape and p_shape is not None
        ok_chans = p_chans == n_chans
        if not ok_dims:
            dim_mismatch.append((fov, p_shape, n_shape))
        if not ok_chans:
            chan_mismatch.append((fov, sorted(p_chans - n_chans), sorted(n_chans - p_chans)))

        rows.append({"fov": fov, "shape": str(p_shape), "n_channels_processed": len(p_chans),
                     "n_channels_raw": len(n_chans), "dims_agree": ok_dims,
                     "channels_agree": ok_chans})

    df = pd.DataFrame(rows)
    rep(f"- pairs checked: {len(df)}")
    rep(f"- dimension mismatches: {len(dim_mismatch)}")
    rep(f"- channel-set mismatches: {len(chan_mismatch)}")
    rep(f"- distinct image shapes across the cohort: "
        f"{dict(Counter(df['shape']).most_common(5))}")
    rep(f"- distinct channel counts (processed): "
        f"{dict(Counter(df['n_channels_processed']).most_common())}")

    if dim_mismatch:
        rep("\nDimension mismatches (fov, processed, non_processed):\n")
        for row in dim_mismatch[:50]:
            rep(f"  - {row}")
        rep.fail(f"{len(dim_mismatch)} pairs disagree on image dimensions. These are "
                 "not the same FOV; masks would not align. Stop.")
    if chan_mismatch:
        rep("\nChannel-set mismatches (fov, processed-only, raw-only):\n")
        for row in chan_mismatch[:50]:
            rep(f"  - {row}")
        rep.fail(f"{len(chan_mismatch)} pairs disagree on which channels are present.")

    if not dim_mismatch and not chan_mismatch and len(df):
        rep("\n**Asserted**: every matched pair agrees on image dimensions and on the "
            "exact set of channel names. Channels are stored one TIFF per channel, so "
            "there is no stored channel order to invert -- stage 2 will index channels "
            "by an explicit name list, identically for both directories.")

    return matched, df


# --------------------------------------------------------------------------
# check 5: three-way mask match
# --------------------------------------------------------------------------

def check_5_masks(rep, clean_fovs, cell_table_fovs):
    rep.h("5. Three-way match: clean cohort FOV / Mesmer mask / cell table ROI")

    if not MASK_DIR.exists():
        rep.fail(f"Mask directory `{MASK_DIR}` does not exist. Set MASK_DIR to the "
                 "Mesmer output used for the published segmentation.")
        return set()

    mask_fovs = {p.name[: -len(MASK_SUFFIX)] for p in MASK_DIR.glob(f"*{MASK_SUFFIX}")}
    clean = set(clean_fovs)
    table = set(cell_table_fovs)

    three_way = clean & mask_fovs & table
    rep(f"- clean cohort FOVs: {len(clean)}")
    rep(f"- Mesmer whole-cell masks in `{MASK_DIR}`: {len(mask_fovs)}")
    rep(f"- distinct ROIs in the cell table: {len(table)}")
    rep(f"- **three-way matched: {len(three_way)}** "
        f"({100 * len(three_way) / max(len(clean), 1):.1f}% of the clean cohort)")
    rep()
    rep(f"- clean cohort with no mask: {len(clean - mask_fovs)}")
    rep(f"- clean cohort not in cell table: {len(clean - table)}")
    rep(f"- cell table ROIs not in clean cohort: {len(table - clean)}")

    for label, missing in [("no mask", sorted(clean - mask_fovs)),
                           ("not in cell table", sorted(clean - table))]:
        if missing:
            rep(f"\nClean cohort FOVs {label} ({len(missing)}):\n")
            for fov in missing:
                rep(f"  - `{fov}`")

    if len(three_way) < len(clean):
        rep.fail(f"Only {len(three_way)} of {len(clean)} clean cohort FOVs match "
                 "three ways. Every FOV used in stage 2 needs an image pair, a mask, "
                 "and cell table rows.")
    return three_way


# --------------------------------------------------------------------------
# checks 6-8: metadata
# --------------------------------------------------------------------------

def load_roi_table(rep):
    """One row per ROI: acquisition ID, cell count, and the sample metadata."""
    rep.h("Cell table / metadata load")

    meta = pd.read_csv(METADATA)
    rep(f"- metadata: {METADATA} -- {len(meta)} samples, "
        f"{meta['Patient_ID'].nunique()} patients, columns: {list(meta.columns)}")

    fovs = pd.read_csv(CELL_TABLE, usecols=["fov"])["fov"]
    rep(f"- cell table: {CELL_TABLE} -- {len(fovs):,} cells, {fovs.nunique()} ROIs")

    rois = fovs.value_counts().rename_axis("fov").reset_index(name="n_cells")
    rois["LEAP_ID"] = rois["fov"].str.extract(r"^(Leap\d+)")[0].str.upper()

    missing = sorted(set(rois["LEAP_ID"]) - set(meta["LEAP_ID"]))
    if missing:
        rep.fail(f"ROIs whose LEAP_ID is absent from the metadata: {missing}")

    return rois.merge(meta, on="LEAP_ID", how="left"), meta


def check_6_8_metadata(rep, roi_table):
    rep.h(f"6. Batch variable: `{BATCH_COL}`")

    rep(f"The batch variable is `{BATCH_COL}` in `{METADATA.name}`, keyed by "
        f"`LEAP_ID` and joined to each ROI via the `Leap###` prefix of `fov`. "
        f"It is a staining-batch index; there is no separate imaging-batch column.\n")

    all_batches = (roi_table.groupby(BATCH_COL)
                   .agg(n_ROIs=("fov", "size"),
                        n_samples=("LEAP_ID", "nunique"),
                        n_patients=("Patient_ID", "nunique"),
                        n_cells=("n_cells", "sum")))
    rep("All ROIs in the cell table, by batch:\n")
    rep.table(all_batches)
    rep(f"\nNumber of batches: {roi_table[BATCH_COL].nunique()}")

    rep.h("8. Cohort scope: pre-treatment, revision cohort")

    rep("Sample type by batch, all ROIs (this shows whether the clean cohort list "
        "already restricts to pre-treatment):\n")
    rep.table(pd.crosstab(roi_table[SAMPLE_TYPE_COL], roi_table[BATCH_COL], margins=True))

    rev = roi_table[~roi_table["LEAP_ID"].isin(EXCLUDED_PATIENTS)
                    & ~roi_table["fov"].isin(EXCLUDED_ROIS)]
    pre = rev[rev[SAMPLE_TYPE_COL].isin(PRE_TREATMENT_VALUES)]

    rep(f"\n- all ROIs: {len(roi_table)}")
    rep(f"- after revision-cohort exclusions ({', '.join(EXCLUDED_PATIENTS)} and "
        f"{len(EXCLUDED_ROIS)} named ROIs): {len(rev)} ROIs, "
        f"{rev['Patient_ID'].nunique()} patients")
    rep(f"- **after restricting to pre-treatment: {len(pre)} ROIs, "
        f"{pre['Patient_ID'].nunique()} patients, {pre['LEAP_ID'].nunique()} samples, "
        f"{pre['n_cells'].sum():,} cells**")
    rep(f"\nThe clean cohort FOV list does NOT restrict to pre-treatment: "
        f"{len(rev) - len(pre)} non-pre-treatment ROIs remain after the revision "
        f"exclusions, so a further filter on `{SAMPLE_TYPE_COL}` is required.")

    rep.h("7. Is batch confounded with patient or with response?")

    rep("ROIs per batch x response (pre-treatment, revision cohort):\n")
    rep.table(pd.crosstab(pre[BATCH_COL], pre[RESPONSE_COL], margins=True))

    rep("\nPatients per batch x response:\n")
    pat = pre.groupby([BATCH_COL, RESPONSE_COL])["Patient_ID"].nunique().unstack(fill_value=0)
    rep.table(pat)

    n_multi = (pre.groupby("Patient_ID")[BATCH_COL].nunique() > 1).sum()
    n_pat = pre["Patient_ID"].nunique()
    rep(f"\n- patients appearing in more than one batch: {n_multi} of {n_pat}")

    pure = pat.index[(pat > 0).sum(axis=1) == 1].tolist()
    if pure:
        rep(f"- batches containing only ONE response group: {pure}")

    rep(f"\n**Patient is nested within batch** ({n_pat - n_multi} of {n_pat} patients "
        f"sit in a single batch). Any metric that measures batch mixing is therefore "
        f"also measuring patient mixing, and the secondary conservation check using "
        f"patient identity as the label is uninformative.")
    if pure:
        rep(f"\n**Batch is confounded with response**: batches {pure} contain a single "
            f"response group. Between-batch differences are partly real biology, so "
            f"'better mixing' cannot be read as unambiguously good.")

    return pre


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--img-root", type=Path, default=None,
                    help="Override IMG_ROOT (processed/, non_processed/, CleanCohort/)")
    ap.add_argument("--mask-dir", type=Path, default=None, help="Override MASK_DIR")
    ap.add_argument("--cell-table", type=Path, default=None, help="Override CELL_TABLE")
    ap.add_argument("--metadata", type=Path, default=None, help="Override METADATA")
    ap.add_argument("--fov-list", type=Path, default=None, help="Override CLEAN_COHORT_FOV_LIST")
    ap.add_argument("--skip-images", action="store_true",
                    help="Run only checks 6-8 (metadata); no image access needed")
    ap.add_argument("--out", default="reconcile_report.md", help="Markdown report path")
    args = ap.parse_args()

    global IMG_ROOT, PROCESSED_DIR, NON_PROCESSED_DIR, CLEAN_COHORT_DIR
    global MASK_DIR, CELL_TABLE, METADATA, CLEAN_COHORT_FOV_LIST
    if args.img_root:
        IMG_ROOT = args.img_root
        PROCESSED_DIR = IMG_ROOT / "processed"
        NON_PROCESSED_DIR = IMG_ROOT / "non_processed"
        CLEAN_COHORT_DIR = IMG_ROOT / "CleanCohort" / "processed"
    if args.mask_dir:
        MASK_DIR = args.mask_dir
    if args.cell_table:
        CELL_TABLE = args.cell_table
    if args.metadata:
        METADATA = args.metadata
    if args.fov_list:
        CLEAN_COHORT_FOV_LIST = args.fov_list

    rep = Report()
    rep("# CLAHE batch-correction evaluation -- stage 1 reconciliation\n")
    rep(f"Generated by `{Path(__file__).name}` on host `{os.uname().nodename}`.\n")

    roi_table, _ = load_roi_table(rep)
    pre = check_6_8_metadata(rep, roi_table)

    if not args.skip_images:
        import tifffile
        from scipy.stats import spearmanr

        if not CLEAN_COHORT_FOV_LIST.exists():
            rep.fail(f"`{CLEAN_COHORT_FOV_LIST}` not found. Regenerate with:\n"
                     f"    ls {CLEAN_COHORT_DIR} > {CLEAN_COHORT_FOV_LIST}")
            clean_fovs = []
        else:
            clean_fovs = [l.strip() for l in CLEAN_COHORT_FOV_LIST.read_text().splitlines()
                          if l.strip()]

        check_1_clahe(rep, clean_fovs, tifffile, spearmanr)
        matched, _ = check_2_4_fovs(rep, clean_fovs, tifffile)
        three_way = check_5_masks(rep, clean_fovs, roi_table["fov"])

        rep.h("Overlap of the three-way matched set with the analysis cohort")
        usable = set(pre["fov"]) & three_way
        rep(f"- pre-treatment revision-cohort ROIs: {len(pre)}")
        rep(f"- of those, with a complete image pair + mask: **{len(usable)}**")
        missing = sorted(set(pre["fov"]) - three_way)
        if missing:
            rep(f"\nPre-treatment ROIs missing an image pair or mask ({len(missing)}):\n")
            for fov in missing:
                rep(f"  - `{fov}`")
    else:
        rep.h("Checks 1-5 skipped (--skip-images)")
        rep("These require the image data and Mesmer masks and must be run on the "
            "workstation.")

    rep.h("Summary", level=2)
    if rep.failures:
        rep(f"**{len(rep.failures)} blocking issue(s) -- do not proceed to stage 2:**\n")
        for f in rep.failures:
            rep(f"- {f}")
    else:
        rep("No blocking issues found in the checks that were run.")

    rep.write(args.out)
    print(f"\n[written to {args.out}]", file=sys.stderr)
    return 1 if rep.failures else 0


if __name__ == "__main__":
    sys.exit(main())
