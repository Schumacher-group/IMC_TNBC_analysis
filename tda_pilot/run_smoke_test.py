#!/usr/bin/env python
"""Task 4 timing: time ONE (Tumour, CD8) chromatic six-pack on the smoke-test ROI.

This replicates exactly what utils.generate_stats.get_stats_k_chromatic_quotient
does for a 2-label codomain (k=1 -> KChromaticInclusion(filt, 1).sixpack()), but
in-process and single-threaded so the timing and peak memory are clean and include
chalc's native (C++) allocations.

Outputs:
  - per-pair wall-clock (filtration + six-pack), and stats-vectorisation time
  - peak process RSS (resource.getrusage) and tracemalloc Python peak
  - raw six-pack diagrams saved to output/diagrams/<fov>.h5  (+ .pkl of matrices)
  - sanity report: bars per diagram per dimension; finite stats vector + shape
"""

import pickle
import resource
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd
from chalc import chromatic
from chalc.sixpack import KChromaticInclusion

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))
from utils.persistent_stats import get_pers_stats, pers_stats_names  # noqa: E402

FOV = sys.argv[1] if len(sys.argv) > 1 else "Leap066_11"
FILTRATION = "delaunay_cech"  # M2S2 default
MAXDIM = 1  # degrees 0 and 1


def canonicalize_labels(v: np.ndarray) -> np.ndarray:
    label_set = set(v.tolist())
    old_to_new = {label: i for i, label in enumerate(label_set)}
    return np.array([old_to_new[x] for x in v.tolist()], dtype=np.int_)


def maxrss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports kibibytes
    return rss / 1e6 if sys.platform == "darwin" else rss / 1e3


def main() -> None:
    csv = HERE / "data" / f"{FOV}.csv"
    df = pd.read_csv(csv)
    df["label"] = df["label"].astype("category")
    n_tumour = int((df["label"] == "Tumour").sum())
    n_cd8 = int((df["label"] == "CD8").sum())
    print(f"ROI {FOV}: Tumour={n_tumour}  CD8={n_cd8}  total points={len(df)}")
    print(f"filtration={FILTRATION}  max_diagram_dimension={MAXDIM}  workers=1")

    points = df[["x", "y"]].to_numpy().transpose()
    colours = canonicalize_labels(df["label"].cat.codes.to_numpy())

    tracemalloc.start()

    # --- the timed unit: one (Tumour, CD8) six-pack -----------------------
    t0 = time.perf_counter()
    filt = getattr(chromatic, FILTRATION)(points, colours)
    t_filt = time.perf_counter()
    dgms = KChromaticInclusion(filt, 1).sixpack()  # k=1 (2-label codomain)
    t_sixpack = time.perf_counter()
    # ----------------------------------------------------------------------

    # vectorise (persistent statistics) exactly as generate_stats does
    barcodes = {name: dgms.get_matrix(name, list(range(MAXDIM + 1))) for name in dgms}
    stats = {}
    bars_report = {}
    for name, per_dim in barcodes.items():
        for dim, bar in enumerate(per_dim):
            bars_report[f"{name}{dim}"] = 0 if bar is None else int(len(bar))
            for sname, sval in zip(pers_stats_names, get_pers_stats(bar), strict=True):
                stats[f"{name}{dim}-{sname}"] = sval
    t_stats = time.perf_counter()

    _, py_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_rss = maxrss_mb()

    filt_s = t_filt - t0
    sixpack_s = t_sixpack - t_filt
    pair_s = t_sixpack - t0           # filtration + six-pack = the projection unit
    stats_s = t_stats - t_sixpack

    # save raw diagrams (native + portable matrices)
    diag_dir = HERE / "output" / "diagrams"
    diag_dir.mkdir(parents=True, exist_ok=True)
    h5_path = diag_dir / f"{FOV}.h5"
    saved_native = False
    try:
        import h5py

        with h5py.File(h5_path, "w") as fh:
            dgms.save(fh)
        saved_native = h5_path.exists()
    except Exception as e:  # noqa: BLE001
        print(f"  (native SixPack.save failed: {e!r}; falling back to pkl only)")
    matrices = {
        name: {dim: (None if bar is None else np.asarray(bar)) for dim, bar in enumerate(per_dim)}
        for name, per_dim in barcodes.items()
    }
    pkl_path = diag_dir / f"{FOV}_matrices.pkl"
    with pkl_path.open("wb") as fh:
        pickle.dump({"fov": FOV, "filtration": FILTRATION, "matrices": matrices}, fh)

    # ---- report ----------------------------------------------------------
    print("\n=== TIMING (single Tumour-CD8 pair, in-process, 1 thread) ===")
    print(f"  filtration ({FILTRATION}) : {filt_s:8.3f} s")
    print(f"  six-pack persistence      : {sixpack_s:8.3f} s")
    print(f"  --> PAIR COMPUTE (proj. unit): {pair_s:8.3f} s")
    print(f"  stats vectorisation       : {stats_s:8.3f} s")
    print("\n=== MEMORY ===")
    print(f"  peak process RSS          : {peak_rss:8.1f} MB")
    print(f"  tracemalloc Python peak   : {py_peak / 1e6:8.1f} MB")

    print("\n=== SIX-PACK BARS per diagram/dimension ===")
    for name in ("dom", "cod", "ker", "cok", "im", "rel"):
        row = "  ".join(f"{name}{d}={bars_report.get(f'{name}{d}', 0)}" for d in range(MAXDIM + 1))
        print(f"  {row}")

    vec = np.array(list(stats.values()), dtype=float)
    print("\n=== STATS VECTOR ===")
    print(f"  length: {len(vec)}  (expected 6 diagrams x {MAXDIM + 1} dims x {len(pers_stats_names)} stats = {6 * (MAXDIM + 1) * len(pers_stats_names)})")
    print(f"  all finite: {np.isfinite(vec).all()}")

    ker_nonempty = (bars_report.get("ker0", 0) + bars_report.get("ker1", 0)) > 0
    im_nonempty = (bars_report.get("im0", 0) + bars_report.get("im1", 0)) > 0
    cok_present = "cok0" in bars_report and "cok1" in bars_report
    deg1 = bars_report.get("dom1", 0) > 0
    print("\n=== SANITY ===")
    print(f"  kernel non-empty:   {ker_nonempty}")
    print(f"  image non-empty:    {im_nonempty}")
    print(f"  cokernel present:   {cok_present} (dim-0 cokernel is always empty by construction)")
    print(f"  degree-1 features present (loops): {deg1}")
    print(f"  raw diagrams saved: native={saved_native} ({h5_path.name}), matrices={pkl_path.name}")


if __name__ == "__main__":
    main()
