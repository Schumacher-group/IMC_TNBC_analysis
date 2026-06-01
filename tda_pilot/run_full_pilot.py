#!/usr/bin/env python
"""Day-1 full pilot: chromatic six-pack persistent statistics for the 808-ROI
revision cohort across 7 pair definitions, parallelised at the ROI level.

Why not loop M2S2's generate_stats: it processes ROIs serially and parallelises
only over the (<=3) label-combinations within one ROI, so --num-workers 10 would
use ~3 workers, and a single bad ROI crashes the whole run. This driver maps ROIs
across a 10-worker pool, calling M2S2's OWN computation primitives unchanged
(chalc chromatic filtration + KChromaticInclusion.sixpack, and utils'
get_stats_from_barcodes_dict / canonicalize_labels). It adds per-ROI failure
isolation, resume, and single-compute raw-diagram saving for the focal pair.

Output per definition: output/stats/<definition>.parquet  (one row per ROI, 456
stats columns + metadata) and output/stats/<definition>_failures.csv.
Raw .h5 six-pack diagrams saved to output/diagrams/CD8_primary/ (focal only).
"""

import argparse
import multiprocessing as mp
import os
import sys
import time
import traceback
from pathlib import Path

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

import numpy as np
import pandas as pd

import cohort
import pair_defs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))

from chalc import chromatic  # noqa: E402
from chalc.sixpack import KChromaticInclusion  # noqa: E402
from utils.generate_stats import canonicalize_labels, get_stats_from_barcodes_dict  # noqa: E402

ROIS_DIR = HERE / "data" / "rois"
STATS_DIR = HERE / "output" / "stats"
DIAG_DIR = HERE / "output" / "diagrams"
FILTRATION = "delaunay_cech"  # M2S2 default
MAXDIM = 1                    # degrees 0 and 1
MIN_COUNT = 3                 # M2S2 default: drop a species with < this many cells in an ROI

# module-level globals set per definition before the pool runs (inherited via fork)
_SPECIES: set[str] = set()
_DEFNAME: str = ""
_SAVE_DIAG: bool = False


def _init(species: set[str], defname: str, save_diag: bool) -> None:
    global _SPECIES, _DEFNAME, _SAVE_DIAG
    _SPECIES, _DEFNAME, _SAVE_DIAG = species, defname, save_diag


def _process_roi(fov: str) -> dict:
    """Compute the (Tumour, species) six-pack stats for one ROI. Never raises."""
    try:
        df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
        is_t = df["celltype"].isin(pair_defs.TUMOUR)
        is_s = df["celltype"].isin(_SPECIES)
        df = df[is_t | is_s].copy()
        df["label"] = np.where(df["celltype"].isin(pair_defs.TUMOUR), "Tumour", _DEFNAME)
        n_tumour = int((df["label"] == "Tumour").sum())
        n_other = int((df["label"] == _DEFNAME).sum())

        base = {"roi_id": fov, "n_tumour": n_tumour, "n_other": n_other,
                "status": "ok", "error": ""}
        if n_tumour < MIN_COUNT or n_other < MIN_COUNT:
            base["status"] = "degenerate"  # too few cells of one species; no diagram
            return base

        points = df[["x", "y"]].to_numpy().transpose()
        colours = canonicalize_labels(df["label"].astype("category").cat.codes.to_numpy())
        filt = getattr(chromatic, FILTRATION)(points, colours)
        dgms = KChromaticInclusion(filt, 1).sixpack()  # k=1 (2-label codomain)

        if _SAVE_DIAG:
            import h5py
            ddir = DIAG_DIR / _DEFNAME
            ddir.mkdir(parents=True, exist_ok=True)
            with h5py.File(ddir / f"{fov}.h5", "w") as fh:
                dgms.save(fh)

        barcodes = {name: dgms.get_matrix(name, list(range(MAXDIM + 1))) for name in dgms}
        stats = get_stats_from_barcodes_dict(barcodes)
        base.update({k: float(v) for k, v in stats.items()})
        return base
    except Exception:  # noqa: BLE001 - isolate per-ROI failures
        return {"roi_id": fov, "status": "failed", "error": traceback.format_exc(limit=3)}


def run_definition(defname: str, species: set[str], n_workers: int, resume: bool) -> dict:
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    out_parquet = STATS_DIR / f"{defname}.parquet"
    fail_csv = STATS_DIR / f"{defname}_failures.csv"
    save_diag = defname in pair_defs.SAVE_DIAGRAMS

    all_rois = sorted(p.stem for p in ROIS_DIR.glob("*.csv") if cohort.is_cohort_member(p.stem))
    done: set[str] = set()
    prev_rows: list[dict] = []
    if resume and out_parquet.exists():
        prev = pd.read_parquet(out_parquet)
        done = set(prev["roi_id"])
        prev_rows = prev.to_dict("records")
    todo = [r for r in all_rois if r not in done]

    t0 = time.perf_counter()
    print(f"[{defname}] {len(all_rois)} cohort ROIs, {len(todo)} to compute "
          f"({len(done)} resumed), workers={n_workers}, save_diagrams={save_diag}")
    rows: list[dict] = []
    if todo:
        with mp.Pool(n_workers, initializer=_init, initargs=(species, defname, save_diag)) as pool:
            for i, row in enumerate(pool.imap_unordered(_process_roi, todo, chunksize=4), 1):
                rows.append(row)
                if i % 100 == 0 or i == len(todo):
                    print(f"  [{defname}] {i}/{len(todo)}")
    elapsed = time.perf_counter() - t0

    all_rows = prev_rows + rows
    res = pd.DataFrame(all_rows)
    failed = res[res["status"] == "failed"]
    degen = res[res["status"] == "degenerate"]
    ok = res[res["status"] == "ok"]

    # write stats parquet (ok + degenerate rows; failures only to the failure log)
    res[res["status"] != "failed"].to_parquet(out_parquet, index=False)
    if len(failed):
        failed[["roi_id", "error"]].to_csv(fail_csv, index=False)
    elif fail_csv.exists():
        fail_csv.unlink()

    summary = {"definition": defname, "n_cohort": len(all_rois), "completed": len(ok),
               "degenerate": len(degen), "failed": len(failed), "elapsed_s": round(elapsed, 1)}
    print(f"[{defname}] completed={len(ok)} degenerate={len(degen)} failed={len(failed)} "
          f"in {elapsed:.1f}s")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--no-resume", action="store_true")
    ap.add_argument("--definitions", nargs="*", default=list(pair_defs.DEFINITIONS),
                    help="subset of definitions to run (default: all 7)")
    args = ap.parse_args()

    mp.set_start_method("fork", force=True)
    t0 = time.perf_counter()
    summaries = []
    for defname in args.definitions:
        species = pair_defs.DEFINITIONS[defname]
        s = run_definition(defname, species, args.workers, resume=not args.no_resume)
        summaries.append(s)
        # Phase 2b safety: stop if failures exceed 5% of the cohort for any definition
        if s["failed"] > 0.05 * s["n_cohort"]:
            print(f"\nSTOP: {defname} failures ({s['failed']}) exceed 5% of cohort. "
                  "Halting before the next definition (per Phase 2b).")
            break

    total = time.perf_counter() - t0
    print("\n===== FULL PILOT SUMMARY =====")
    sm = pd.DataFrame(summaries)
    print(sm.to_string(index=False))
    print(f"total wall clock: {total/60:.1f} min")
    sm.to_csv(STATS_DIR / "run_summary.csv", index=False)


if __name__ == "__main__":
    main()
