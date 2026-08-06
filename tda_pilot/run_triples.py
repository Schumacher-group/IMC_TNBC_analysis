#!/usr/bin/env python
"""Three-species chromatic six-packs: the analysis the tool was chosen for.

Every previous run in this directory used `KChromaticInclusion(filt, 1)` on exactly two
labels. M2S2's own `get_stats_k_chromatic_quotient` uses `KChromaticQuotient(filt, k)` with
k = n_labels - 1 whenever more than two labels are present, so for a triple that is
`KChromaticQuotient(filt, 2)`. This script uses that code path unchanged.

Scope note: this computes OBSERVED statistics only -- no per-ROI permutation null. That is
deliberate. The most robust instrument we have is the supervised classifier
(`day2_classify_response.py`, AUC 0.544, p = 0.275), which needs only observed statistics and
is insensitive to the adjustment arguments that dogged the z-score analyses. So triples are
screened by adding them as feature groups to that classifier. Only a triple that moves the
AUC would justify the ~3 h per definition that a permutation null costs.

Output: output/stats/triple_<name>.parquet, one row per ROI, 456 statistics + metadata,
in the same layout as the pair tables so the classifier can consume both.
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
import triple_defs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))

from chalc import chromatic  # noqa: E402
from chalc.sixpack import KChromaticQuotient  # noqa: E402
from utils.generate_stats import canonicalize_labels, get_stats_from_barcodes_dict  # noqa: E402

ROIS_DIR = HERE / "data" / "rois"
STATS_DIR = HERE / "output" / "stats"
FILTRATION = "delaunay_cech"
MAXDIM = 1
MIN_COUNT = 3          # M2S2 default: a species with fewer cells is dropped for that ROI

_SPEC: list[tuple[str, set[str]]] = []
_NAME = ""


def _init(spec, name) -> None:
    global _SPEC, _NAME
    _SPEC, _NAME = spec, name


def _process_roi(fov: str) -> dict:
    """Compute the three-species six-pack statistics for one ROI. Never raises."""
    try:
        df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
        # Assign each cell to its species; cells matching none are dropped.
        label = pd.Series(pd.NA, index=df.index, dtype="object")
        counts = {}
        for name, labels in _SPEC:
            m = df["celltype"].isin(labels)
            label[m] = name
            counts[f"n_{name}"] = int(m.sum())
        df = df.assign(label=label).dropna(subset=["label"])

        base = {"roi_id": fov, "triple": _NAME, "status": "ok", "error": "", **counts}
        if any(v < MIN_COUNT for v in counts.values()):
            base["status"] = "degenerate"
            return base

        points = df[["x", "y"]].to_numpy().transpose()
        # Categorical codes are alphabetical; canonicalize_labels maps them to 0..k-1.
        colours = canonicalize_labels(df["label"].astype("category").cat.codes.to_numpy())
        n_colours = len(set(colours.tolist()))
        if n_colours != len(_SPEC):
            base["status"] = "degenerate"
            return base

        filt = getattr(chromatic, FILTRATION)(points, colours)
        # k = n_labels - 1, exactly as M2S2's get_stats_k_chromatic_quotient does.
        dgms = KChromaticQuotient(filt, n_colours - 1).sixpack()
        barcodes = {nm: dgms.get_matrix(nm, list(range(MAXDIM + 1))) for nm in dgms}
        base.update({k: float(v) for k, v in get_stats_from_barcodes_dict(barcodes).items()})
        return base
    except Exception:  # noqa: BLE001 - isolate per-ROI failures
        return {"roi_id": fov, "triple": _NAME, "status": "failed",
                "error": traceback.format_exc(limit=3)}


def pretreatment_rois() -> list[str]:
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Sample_Type_(pre/post treatment)": "sample_type"})
    meta = meta[(meta["sample_type"].astype(str).str.lower() == "pre")
                & meta["Response"].isin(["Responder", "Non-Responder"])]
    have = {p.stem for p in ROIS_DIR.glob("*.csv")}
    return sorted(r for r in meta["roi_id"] if cohort.is_cohort_member(r) and r in have)


def run_triple(name: str, workers: int, resume: bool) -> None:
    spec = triple_defs.TRIPLES[name]
    out = STATS_DIR / f"triple_{name}.parquet"
    rois = pretreatment_rois()

    prev_rows, done = [], set()
    if resume and out.exists():
        prev = pd.read_parquet(out)
        prev_rows, done = prev.to_dict("records"), set(prev["roi_id"])
    todo = [r for r in rois if r not in done]
    sizes = pd.read_csv(HERE / "per_roi_counts.csv").set_index("fov")["n_cells_total"]
    todo.sort(key=lambda r: -float(sizes.get(r, 0)))   # big ROIs first, as in perm_null.py

    t0 = time.perf_counter()
    print(f"[{name}] {len(rois)} pre-treatment ROIs, {len(todo)} to compute "
          f"({len(done)} resumed), species={[s for s, _ in spec]}, workers={workers}")
    rows = []
    if todo:
        with mp.Pool(workers, initializer=_init, initargs=(spec, name)) as pool:
            for i, row in enumerate(pool.imap_unordered(_process_roi, todo, chunksize=2), 1):
                rows.append(row)
                if i % 50 == 0 or i == len(todo):
                    el = time.perf_counter() - t0
                    print(f"  [{name}] {i}/{len(todo)}  ({el/60:.1f} min, "
                          f"~{el/i*(len(todo)-i)/60:.0f} min left)")

    res = pd.DataFrame(prev_rows + rows)
    res[res["status"] != "failed"].to_parquet(out, index=False)
    n_ok = int((res.status == "ok").sum())
    print(f"[{name}] ok={n_ok} degenerate={int((res.status=='degenerate').sum())} "
          f"failed={int((res.status=='failed').sum())} in "
          f"{(time.perf_counter()-t0)/60:.1f} min -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--triples", nargs="*", default=["Tumour_CD8_Fibroblast"],
                    choices=list(triple_defs.TRIPLES) + [],)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    mp.set_start_method("fork", force=True)   # spawn breaks M2S2 on macOS, see REPORT.md
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    for name in args.triples:
        run_triple(name, args.workers, resume=not args.no_resume)


if __name__ == "__main__":
    main()
