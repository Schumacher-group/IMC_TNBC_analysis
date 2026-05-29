#!/usr/bin/env python
"""Scan the full cell table once to compute per-ROI counts of the tumour aggregate
and the three candidate CD8 definitions, joined with responder status.

Outputs:
  - tda_pilot/per_roi_counts.csv  (one row per fov/ROI)
  - printed summary: how many ROIs clear count thresholds under each CD8 definition,
    and candidate representative ROIs for the smoke test (responder, well populated).

No qc_pass filter (per user: use all cells).
"""

from pathlib import Path

import pandas as pd

import cohort

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CELL_TABLE = ROOT / "CellTable_CleanCohort" / "updated_cell_table.csv"
META = ROOT / "CellTable_CleanCohort" / "CleanCohort_Metadata.csv"

TUMOUR = {"Cancer cell", "B7H4 Cancer cell"}
CD8_PRIMARY = {"CD8 T cell", "Memory CD8 T cell"}
CD8_ALONE = {"CD8 T cell"}
CD8_PLUS_NK = {"CD8 T cell", "Memory CD8 T cell", "NK/CD8"}

USECOLS = ["fov", "cell_meta_cluster"]
CHUNK = 500_000


def main() -> None:
    # per-fov, per-celltype counts accumulated across chunks
    counts = None
    n = 0
    for chunk in pd.read_csv(CELL_TABLE, usecols=USECOLS, chunksize=CHUNK):
        n += len(chunk)
        g = chunk.groupby(["fov", "cell_meta_cluster"]).size()
        counts = g if counts is None else counts.add(g, fill_value=0)
    print(f"scanned rows: {n:,}")

    wide = counts.unstack(fill_value=0)  # rows=fov, cols=cell_meta_cluster
    wide.columns = [str(c) for c in wide.columns]

    def col_sum(cols: set[str]) -> pd.Series:
        present = [c for c in cols if c in wide.columns]
        return wide[present].sum(axis=1) if present else pd.Series(0, index=wide.index)

    out = pd.DataFrame(index=wide.index)
    out.index.name = "fov"
    out["n_cells_total"] = wide.sum(axis=1)
    out["tumour"] = col_sum(TUMOUR)
    out["cd8_primary"] = col_sum(CD8_PRIMARY)
    out["cd8_alone"] = col_sum(CD8_ALONE)
    out["cd8_plus_nk"] = col_sum(CD8_PLUS_NK)

    # join responder status via LEAP_ID (fov = 'Leap001_10' -> 'LEAP001')
    out = out.reset_index()
    out["LEAP_ID"] = out["fov"].str.split("_").str[0].str.upper()
    meta = pd.read_csv(META)
    meta["LEAP_ID"] = meta["LEAP_ID"].str.upper()
    out = out.merge(
        meta[["LEAP_ID", "Patient_ID", "Response", "Biopsy_(Core/Resection)",
              "Sample_Type_(pre/post treatment)"]],
        on="LEAP_ID", how="left",
    )
    # Revision-cohort membership (N = 808). Keep the full 829-row table; downstream
    # analyses filter on this boolean (or use cohort.predictive_cohort_rois).
    out["cohort_member"] = out["fov"].apply(cohort.is_cohort_member)
    out = out.sort_values("fov").reset_index(drop=True)
    out.to_csv(HERE / "per_roi_counts.csv", index=False)

    members = out[out["cohort_member"]]
    print(f"\ntotal ROIs (fov): {len(out)}  [raw cell table]")
    print(f"revision-cohort ROIs (cohort_member): {len(members)}  [expected 808]")
    print(f"revision-cohort patients (Patient_ID): {members['Patient_ID'].nunique()}  [expected 63]")
    print(f"excluded ROIs: {(~out['cohort_member']).sum()}")
    print("Response value counts (per-ROI, full table):")
    print(out["Response"].value_counts(dropna=False).to_dict())

    print("\n=== ROIs clearing thresholds (tumour AND CD8 of given definition) ===")
    for thr in (20, 50):
        print(f"-- threshold >= {thr} cells of EACH --")
        for name in ("cd8_primary", "cd8_alone", "cd8_plus_nk"):
            ok = (out["tumour"] >= thr) & (out[name] >= thr)
            print(f"   tumour & {name:12s}: {ok.sum():4d} / {len(out)} ROIs")

    print("\n=== ROIs dropping below 20 cells of the aggregated CD8 species, by CD8 def ===")
    for name in ("cd8_primary", "cd8_alone", "cd8_plus_nk"):
        below = (out[name] < 20).sum()
        print(f"   {name:12s}: {below:4d} ROIs below 20 ({100*below/len(out):.1f}%)")

    # candidate smoke-test ROIs: responder, well populated, some infiltration
    resp = out[out["Response"].astype(str).str.contains("Responder", case=False, na=False)
               & ~out["Response"].astype(str).str.contains("Non", case=False, na=False)]
    cand = resp[(resp["tumour"] >= 50) & (resp["cd8_primary"] >= 50)].copy()
    cand["min_pair"] = cand[["tumour", "cd8_primary"]].min(axis=1)
    cand = cand.sort_values("min_pair", ascending=False)
    print("\n=== Candidate RESPONDER ROIs (tumour>=50 & cd8_primary>=50), top 15 by min(pair) ===")
    cols = ["fov", "LEAP_ID", "Response", "Sample_Type_(pre/post treatment)",
            "n_cells_total", "tumour", "cd8_primary", "cd8_alone"]
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(cand[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
