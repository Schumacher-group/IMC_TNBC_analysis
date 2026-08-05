#!/usr/bin/env python
"""Per-ROI label-permutation null: arrangement, conditioned on composition AND geometry.

The question
------------
Every Day-1 adjustment (OLS on cell counts, minority-fraction subsetting) tried to strip
confounds *across* ROIs. That cannot work cleanly here: cell number, class balance, cell
density, imaged area and spatial dispersion are all collinear with each other and with
response, and the relationships are nonlinear (IMBALANCE_CONFOUND.md says so explicitly).

The permutation null sidesteps all of it. For one ROI, keep the cell POSITIONS exactly as
observed and shuffle only the tumour/CD8 LABELS across those positions. Every confound is
held fixed by construction -- same points, same counts of each species, same density, same
hull area, same dispersion, same overall point pattern. The only thing that changes is
which cell is which. Comparing the observed statistic to its own ROI's null distribution
therefore isolates *spatial arrangement* and nothing else.

Each ROI yields z = (observed - null mean) / null sd, plus the observed value's rank in the
null. z is the ROI's arrangement score, on a common scale across ROIs regardless of size,
density or composition. `day2_permutation_test.py` then does the responder contrast on z.

Method
------
Point clouds, filtration and six-pack are built exactly as in `run_full_pilot.py`
(pair_defs -> chromatic.delaunay_cech -> KChromaticInclusion(filt, 1).sixpack()), and
statistics come from M2S2's own `get_stats_from_barcodes_dict` unchanged. The permutation
acts on the canonicalised colour vector only, so species counts are preserved exactly.

Note the filtration is NOT invariant under relabelling -- the chromatic Delaunay-Cech
complex is built on the colour-lifted point set -- so a fresh six-pack is required for every
permutation. That is the whole cost: 1 + B six-packs per ROI.

Scope: pre-treatment ROIs only. Pooled pre+post is confounded by treatment stage (post
resections are 170 NR vs 7 R; see PRETREATMENT_CONFOUND.md), and baseline is the setting
where a predictive claim would live.

Output: output/perm_null_<definition>.parquet, one row per ROI, columns
    <stat>__obs, <stat>__null_mean, <stat>__null_sd, <stat>__z, <stat>__null_rank
for each of the retained statistics.
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
from chalc.sixpack import KChromaticInclusion, SubChromaticInclusion  # noqa: E402
from utils.generate_stats import canonicalize_labels, get_stats_from_barcodes_dict  # noqa: E402

ROIS_DIR = HERE / "data" / "rois"
OUT_DIR = HERE / "output"
FILTRATION = "delaunay_cech"  # as in run_full_pilot.py
MAXDIM = 1
MIN_COUNT = 3

# Pre-specified retained statistics, fixed BEFORE looking at any permutation output.
# Lengths (avg/med/p90) are the scale statistics motivating this run; num_bars is the
# count term that total persistence conflates them with; entropy summarises bar spread.
DIAGRAMS = ("ker", "im", "cok", "dom", "cod")
STATS = ("avg_length", "med_length", "p90_length", "num_bars", "entropy")
RETAINED = tuple(f"{dg}{dim}-{st}" for dg in DIAGRAMS for dim in (0, 1) for st in STATS)

# Primary endpoint, declared in the plan before this script was run.
PRIMARY = "ker1-avg_length"
# Invariant under relabelling (dim-0 bars of the full complex = one per point), so its
# null sd must be exactly 0. A non-zero value means the permutation is not doing what we
# think it is.
INVARIANT_CHECK = "cod0-num_bars"

# Which inclusion K -> L to take the six-pack of.
#
#   symmetric  KChromaticInclusion(filt, 1): domain = monochromatic simplices of BOTH
#              colours. Colour-symmetric, so with tumour outnumbering CD8 ~2:1 the
#              kernel mixes CD8-derived and tumour-derived loops and is tumour-dominated.
#              This is what run_full_pilot.py and the whole Day-1 pilot used.
#
#   other_only SubChromaticInclusion(filt, [[OTHER_COLOUR]]): domain = the non-tumour
#              species alone. DIRECTIONAL, and the statistic an immune-EXCLUSION
#              hypothesis actually names: a tumour nest ringed by CD8 puts an annulus in
#              the CD8-only complex, and that loop dies as soon as the tumour cells
#              filling the nest are added -- i.e. it appears as a long degree-1 KERNEL
#              bar. Exclusion in non-responders therefore predicts NR > R on
#              ker1 length statistics, i.e. Cliff's delta < 0. Note every trend measured
#              under `symmetric` runs the opposite way (delta > 0), so this is a genuine
#              test rather than a rerun of the same thing.
INCLUSIONS = ("symmetric", "other_only")

# Labels are the definition name (e.g. "CD8_primary") and "Tumour"; pandas categorical
# codes are assigned in sorted order, and every definition name sorts before "Tumour",
# so the non-tumour species is always colour 0. Asserted at runtime in _process_roi.
OTHER_COLOUR = 0

_SPECIES: set[str] = set()
_DEFNAME = ""
_NPERM = 0
_SEED = 0
_INCLUSION = "symmetric"


def _init(species, defname, nperm, seed, inclusion) -> None:
    global _SPECIES, _DEFNAME, _NPERM, _SEED, _INCLUSION
    _SPECIES, _DEFNAME, _NPERM, _SEED = species, defname, nperm, seed
    _INCLUSION = inclusion


def _stats_for(points: np.ndarray, colours: np.ndarray) -> dict:
    """One six-pack -> the retained statistics. Same code path as run_full_pilot.py."""
    filt = getattr(chromatic, FILTRATION)(points, colours)
    if _INCLUSION == "other_only":
        dgms = SubChromaticInclusion(filt, [[OTHER_COLOUR]]).sixpack()
    else:
        dgms = KChromaticInclusion(filt, 1).sixpack()
    barcodes = {name: dgms.get_matrix(name, list(range(MAXDIM + 1))) for name in dgms}
    allstats = get_stats_from_barcodes_dict(barcodes)
    return {k: float(allstats[k]) for k in RETAINED if k in allstats}


def _roi_seed(fov: str) -> int:
    """Deterministic per-ROI seed, so a rerun reproduces byte-identical output."""
    return (_SEED + int.from_bytes(fov.encode(), "little")) % (2**32)


def _process_roi(fov: str) -> dict:
    """Observed + null distribution for one ROI. Never raises."""
    try:
        df = pd.read_csv(ROIS_DIR / f"{fov}.csv")
        is_t = df["celltype"].isin(pair_defs.TUMOUR)
        is_s = df["celltype"].isin(_SPECIES)
        df = df[is_t | is_s].copy()
        df["label"] = np.where(df["celltype"].isin(pair_defs.TUMOUR), "Tumour", _DEFNAME)
        n_tumour = int((df["label"] == "Tumour").sum())
        n_other = int((df["label"] == _DEFNAME).sum())

        base = {"roi_id": fov, "n_tumour": n_tumour, "n_other": n_other,
                "n_perm": _NPERM, "status": "ok", "error": ""}
        if n_tumour < MIN_COUNT or n_other < MIN_COUNT:
            base["status"] = "degenerate"
            return base

        points = df[["x", "y"]].to_numpy().transpose()
        cats = df["label"].astype("category")
        colours = canonicalize_labels(cats.cat.codes.to_numpy())
        # The directional inclusion is meaningless if OTHER_COLOUR is not the species we
        # think it is, and this is silent if it ever breaks -- so check it, not assume it.
        if list(cats.cat.categories)[OTHER_COLOUR] == "Tumour":
            msg = (f"{fov}: colour {OTHER_COLOUR} is Tumour, not {_DEFNAME}; "
                   "OTHER_COLOUR assumption violated")
            raise AssertionError(msg)

        obs = _stats_for(points, colours)

        rng = np.random.default_rng(_roi_seed(fov))
        null = {k: np.empty(_NPERM) for k in obs}
        for b in range(_NPERM):
            # Permuting the colour vector preserves both species counts exactly and
            # leaves every point position untouched.
            perm = _stats_for(points, rng.permutation(colours))
            for k in null:
                null[k][b] = perm.get(k, np.nan)

        for k, o in obs.items():
            v = null[k]
            sd = float(np.nanstd(v, ddof=1))
            mean = float(np.nanmean(v))
            base[f"{k}__obs"] = o
            base[f"{k}__null_mean"] = mean
            base[f"{k}__null_sd"] = sd
            base[f"{k}__z"] = (o - mean) / sd if sd > 0 else 0.0
            # fraction of null draws below the observation; 0.5 => indistinguishable
            base[f"{k}__null_rank"] = float(np.nanmean(v < o))
        return base
    except Exception:  # noqa: BLE001 - isolate per-ROI failures, as run_full_pilot.py does
        return {"roi_id": fov, "status": "failed", "error": traceback.format_exc(limit=3)}


def pretreatment_rois() -> list[str]:
    """Cohort ROIs that are pre-treatment and carry a responder label."""
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Sample_Type_(pre/post treatment)": "sample_type"})
    meta = meta[meta["sample_type"].astype(str).str.lower() == "pre"]
    meta = meta[meta["Response"].isin(["Responder", "Non-Responder"])]
    have = {p.stem for p in ROIS_DIR.glob("*.csv")}
    return sorted(r for r in meta["roi_id"]
                  if cohort.is_cohort_member(r) and r in have)


def run_definition(defname: str, workers: int, nperm: int, seed: int, resume: bool,
                   inclusion: str) -> None:
    species = pair_defs.DEFINITIONS[defname]
    # The symmetric run keeps its original filename so existing outputs stay valid.
    tag = "" if inclusion == "symmetric" else f"_{inclusion}"
    out = OUT_DIR / f"perm_null_{defname}{tag}.parquet"
    fail_csv = OUT_DIR / f"perm_null_{defname}{tag}_failures.csv"

    rois = pretreatment_rois()
    prev_rows: list[dict] = []
    done: set[str] = set()
    if resume and out.exists():
        prev = pd.read_parquet(out)
        # only resume rows computed with the same permutation budget and seed
        if int(prev["n_perm"].iloc[0]) == nperm:
            prev_rows = prev.to_dict("records")
            done = set(prev["roi_id"])
    todo = [r for r in rois if r not in done]
    # Cost scales ~n^1.5, and the largest ROI is 120x the smallest. Dispatch big ROIs
    # first so the tail of the run does not leave workers idle.
    sizes = pd.read_csv(HERE / "per_roi_counts.csv").set_index("fov")["n_cells_total"]
    todo.sort(key=lambda r: -float(sizes.get(r, 0)))

    t0 = time.perf_counter()
    print(f"[{defname}] {len(rois)} pre-treatment ROIs, {len(todo)} to compute "
          f"({len(done)} resumed), B={nperm}, workers={workers}, inclusion={inclusion}")
    rows: list[dict] = []
    if todo:
        with mp.Pool(workers, initializer=_init,
                     initargs=(species, defname, nperm, seed, inclusion)) as pool:
            for i, row in enumerate(pool.imap_unordered(_process_roi, todo, chunksize=1), 1):
                rows.append(row)
                if i % 25 == 0 or i == len(todo):
                    el = time.perf_counter() - t0
                    eta = el / i * (len(todo) - i) / 60
                    print(f"  [{defname}] {i}/{len(todo)}  ({el/60:.1f} min elapsed, "
                          f"~{eta:.0f} min left)")

    res = pd.DataFrame(prev_rows + rows)
    failed = res[res["status"] == "failed"]
    res[res["status"] != "failed"].to_parquet(out, index=False)
    if len(failed):
        failed[["roi_id", "error"]].to_csv(fail_csv, index=False)
    elif fail_csv.exists():
        fail_csv.unlink()

    ok = res[res["status"] == "ok"]
    print(f"\n[{defname}] ok={len(ok)} degenerate={int((res.status=='degenerate').sum())} "
          f"failed={len(failed)} in {(time.perf_counter()-t0)/60:.1f} min")
    if len(ok):
        inv_sd = ok[f"{INVARIANT_CHECK}__null_sd"].abs().max()
        print(f"  sanity: max null sd of {INVARIANT_CHECK} = {inv_sd:.3e} (must be 0)")
        zcols = [c for c in ok.columns if c.endswith("__z")]
        print(f"  {len(zcols)} z columns, all finite: {bool(np.isfinite(ok[zcols]).all().all())}")
        print(f"  primary endpoint {PRIMARY}: median z = "
              f"{ok[f'{PRIMARY}__z'].median():+.3f}")
    print("wrote", out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--n-perm", type=int, default=99)
    ap.add_argument("--seed", type=int, default=20260805)
    ap.add_argument("--definitions", nargs="*", default=["CD8_primary"])
    ap.add_argument("--inclusion", choices=INCLUSIONS, default="symmetric",
                    help="'symmetric' reproduces the Day-1/Day-2 six-pack; 'other_only' is "
                         "the directional CD8-vs-tumour inclusion an exclusion hypothesis "
                         "calls for (see INCLUSIONS)")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    mp.set_start_method("fork", force=True)  # see REPORT.md: spawn breaks M2S2 on macOS
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for defname in args.definitions:
        run_definition(defname, args.workers, args.n_perm, args.seed,
                       resume=not args.no_resume, inclusion=args.inclusion)


if __name__ == "__main__":
    main()
