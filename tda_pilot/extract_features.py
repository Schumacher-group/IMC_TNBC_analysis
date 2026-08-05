#!/usr/bin/env python
"""Phase 4: focused per-(ROI, definition) feature table for Day-2 modelling.

Computes interpretable degree-1 summaries of the kernel / image / cokernel
diagrams, EXACTLY from the persistence diagrams (the >5 um lifetime threshold and
true max-lifetime are not recoverable from the 456 summary stats, so we use the
diagrams). CD8_primary diagrams are loaded from the saved .h5; the other six are
recomputed (same chalc filtration + KChromaticInclusion.sixpack as the full run).

Output: output/per_roi_summary.parquet  (one row per ROI x pair_definition).
"""

import multiprocessing as mp
import os
import sys
from pathlib import Path

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

import numpy as np
import pandas as pd

import cohort
import pair_defs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))
from chalc import chromatic  # noqa: E402
from chalc.sixpack import KChromaticInclusion, SixPack  # noqa: E402
from utils.generate_stats import canonicalize_labels  # noqa: E402

ROIS_DIR = HERE / "data" / "rois"
DIAG_DIR = HERE / "output" / "diagrams"
STATS_DIR = HERE / "output" / "stats"
FILTRATION = "delaunay_cech"
MIN_COUNT = 3
LIFETIME_THRESHOLD = 5.0  # micrometres

# Emitted per diagram as <diagram>_dim1_<name>. `n_bars_all` vs `n_bars_gt5um` makes the
# two different bar populations explicit; the old single `n_features` column (bars > 5 um)
# hid the distinction. See the note in _dim1_feats.
FEATURE_NAMES = ("total_persistence", "n_bars_all", "n_bars_gt5um",
                 "persistent_entropy", "max_lifetime")

_SPECIES: set[str] = set()
_DEFNAME = ""
_USE_H5 = False


def _init(species, defname, use_h5):
    global _SPECIES, _DEFNAME, _USE_H5
    _SPECIES, _DEFNAME, _USE_H5 = species, defname, use_h5


def _dim1_feats(bar, thr=LIFETIME_THRESHOLD) -> dict:
    """Degree-1 summaries of one diagram.

    NOTE on bar populations -- these features do NOT all summarise the same bars, and
    conflating them is what made `im_dim1_total_persistence` a misleading headline feature
    in the Day-1 analysis:

      * `total_persistence`, `persistent_entropy`, `max_lifetime`, `n_bars_all`
        are over EVERY finite positive-lifetime bar. Median degree-1 bar lifetime in this
        cohort is ~0.33 um -- far below one cell diameter -- so `total_persistence` is
        dominated by sub-cell-scale features and is close to a rescaled bar count.
      * `n_bars_gt5um` counts only bars longer than `thr` (5 um), i.e. ~3% of them.

    `total_persistence` also factorises as n_bars_all x mean lifetime; the count term
    carries the cell-number and class-imbalance confounds, the length term is the actual
    spatial-scale statistic. See `day2_length_features.py`.
    """
    bar = np.asarray(bar) if bar is not None else np.zeros((0, 2))
    out = {"total_persistence": 0.0, "n_bars_all": 0, "n_bars_gt5um": 0,
           "persistent_entropy": 0.0, "max_lifetime": 0.0}
    if len(bar) == 0:
        return out
    fin = bar[np.isfinite(bar).all(axis=1)]
    life = fin[:, 1] - fin[:, 0]
    life = life[life > 0]
    if len(life) == 0:
        return out
    p = life / life.sum()
    out["total_persistence"] = float(life.sum())
    out["n_bars_all"] = int(len(life))
    out["n_bars_gt5um"] = int((life > thr).sum())
    out["persistent_entropy"] = float(-(p * np.log(p)).sum())
    out["max_lifetime"] = float(life.max())
    return out


def _sixpack_for(fov: str):
    """Return (dgms_or_None, n_tumour, n_other, status)."""
    df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
    is_t = df["celltype"].isin(pair_defs.TUMOUR)
    is_s = df["celltype"].isin(_SPECIES)
    df = df[is_t | is_s].copy()
    df["label"] = np.where(df["celltype"].isin(pair_defs.TUMOUR), "Tumour", _DEFNAME)
    n_t = int((df["label"] == "Tumour").sum())
    n_o = int((df["label"] == _DEFNAME).sum())
    if n_t < MIN_COUNT or n_o < MIN_COUNT:
        return None, n_t, n_o, "degenerate"
    if _USE_H5:
        import h5py
        with h5py.File(DIAG_DIR / _DEFNAME / f"{fov}.h5", "r") as fh:
            return SixPack.from_file(fh), n_t, n_o, "ok"
    points = df[["x", "y"]].to_numpy().transpose()
    colours = canonicalize_labels(df["label"].astype("category").cat.codes.to_numpy())
    filt = getattr(chromatic, FILTRATION)(points, colours)
    return KChromaticInclusion(filt, 1).sixpack(), n_t, n_o, "ok"


def _worker(fov: str) -> dict:
    try:
        dgms, n_t, n_o, status = _sixpack_for(fov)
        row = {"roi_id": fov, "pair_definition": _DEFNAME, "n_tumour": n_t, "n_other": n_o,
               "status": status}
        if dgms is None:
            for dg in ("ker", "im", "cok"):
                for f in FEATURE_NAMES:
                    row[f"{dg}_dim1_{f}"] = np.nan
            return row
        for dg in ("ker", "im", "cok"):
            f = _dim1_feats(dgms.get_matrix(dg, [1])[0])
            for k, v in f.items():
                row[f"{dg}_dim1_{k}"] = v
        return row
    except Exception as e:  # noqa: BLE001
        return {"roi_id": fov, "pair_definition": _DEFNAME, "status": f"failed: {e!r}"}


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--definitions", nargs="*", default=list(pair_defs.DEFINITIONS),
                    help="subset of definitions to (re)compute; merged into existing summary")
    args = ap.parse_args()

    mp.set_start_method("fork", force=True)
    rois = sorted(p.stem for p in ROIS_DIR.glob("*.csv") if cohort.is_cohort_member(p.stem))
    all_rows = []
    for defname in args.definitions:
        species = pair_defs.DEFINITIONS[defname]
        use_h5 = defname in pair_defs.SAVE_DIAGRAMS and (DIAG_DIR / defname).exists()
        with mp.Pool(10, initializer=_init, initargs=(species, defname, use_h5)) as pool:
            rows = list(pool.imap_unordered(_worker, rois, chunksize=4))
        all_rows.extend(rows)
        print(f"[{defname}] {len(rows)} rows (use_h5={use_h5})")

    df = pd.DataFrame(all_rows)
    # attach patient_id + response
    meta = pd.read_csv(HERE / "per_roi_counts.csv")[["fov", "Patient_ID", "Response"]]
    meta = meta.rename(columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response"})
    df = df.merge(meta, on="roi_id", how="left")

    # tidy column order
    feat_cols = [f"{dg}_dim1_{f}" for dg in ("ker", "im", "cok") for f in FEATURE_NAMES]
    # cokernel: keep only the persistence/count features per the Phase-4 spec
    feat_cols = [c for c in feat_cols if not (c.startswith("cok_dim1_")
                 and c.endswith(("persistent_entropy", "max_lifetime")))]
    lead = ["roi_id", "patient_id", "response", "pair_definition", "n_tumour", "n_other", "status"]
    df = df[lead + [c for c in feat_cols if c in df.columns]]
    out = HERE / "output" / "per_roi_summary.parquet"

    # merge into any existing summary: replace the recomputed definitions, keep the rest
    if out.exists() and set(args.definitions) != set(pair_defs.DEFINITIONS):
        prev = pd.read_parquet(out)
        prev = prev[~prev["pair_definition"].isin(args.definitions)]
        df = pd.concat([prev, df], ignore_index=True)
    df = df.sort_values(["pair_definition", "roi_id"]).reset_index(drop=True)
    df.to_parquet(out, index=False)

    n_expected = len(rois) * len(pair_defs.DEFINITIONS)
    print(f"PHASE4_ROWS={len(df)} EXPECTED={n_expected} "
          f"OK={int((df['status']=='ok').sum())} DEGEN={int((df['status']=='degenerate').sum())} "
          f"FAILED={int(df['status'].astype(str).str.startswith('failed').sum())}")
    print("wrote", out)


if __name__ == "__main__":
    main()
