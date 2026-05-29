#!/usr/bin/env python
"""Task 3 adapter: extract one ROI's (tumour, CD8) point cloud in M2S2 format.

- Tumour aggregate = {Cancer cell, B7H4 Cancer cell} -> label "Tumour"
- CD8 (primary)    = {CD8 T cell, Memory CD8 T cell} -> label "CD8"
- Coordinates: centroid-1 -> x, centroid-0 -> y, treated as micrometres
  (Hyperion IMC = 1 px = 1 um). No qc_pass filter (per user: use all cells).

Writes:
  tda_pilot/data/<fov>.csv         (columns: x, y, label)
  tda_pilot/data/<fov>_overlay.png (sanity plot)
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import cohort

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CELL_TABLE = ROOT / "CellTable_CleanCohort" / "updated_cell_table.csv"

TUMOUR = {"Cancer cell", "B7H4 Cancer cell"}
CD8 = {"CD8 T cell", "Memory CD8 T cell"}
USECOLS = ["fov", "cell_meta_cluster", "centroid-0", "centroid-1"]
CHUNK = 500_000


def extract(fov: str) -> pd.DataFrame:
    keep = TUMOUR | CD8
    parts = []
    for chunk in pd.read_csv(CELL_TABLE, usecols=USECOLS, chunksize=CHUNK):
        sub = chunk[(chunk["fov"] == fov) & (chunk["cell_meta_cluster"].isin(keep))]
        if len(sub):
            parts.append(sub)
    if not parts:
        msg = f"No tumour/CD8 cells found for fov={fov!r}"
        raise SystemExit(msg)
    df = pd.concat(parts, ignore_index=True)
    label = df["cell_meta_cluster"].map(
        lambda c: "Tumour" if c in TUMOUR else "CD8",
    )
    return pd.DataFrame(
        {"x": df["centroid-1"].astype(float), "y": df["centroid-0"].astype(float), "label": label}
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fov", default="Leap066_11")
    args = ap.parse_args()

    if not cohort.is_cohort_member(args.fov):
        print(
            f"WARNING: {args.fov} is NOT in the revision cohort (N=808) "
            "(excluded patient LEAP149/LEAP150 or a collaborator-flagged ROI). "
            "Proceeding anyway, but it should not feed the pilot."
        )

    out_dir = HERE / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = extract(args.fov)
    counts = df["label"].value_counts().to_dict()
    print(f"fov: {args.fov}")
    print(f"  Tumour cells: {counts.get('Tumour', 0)}")
    print(f"  CD8 cells:    {counts.get('CD8', 0)}")
    print(f"  total points: {len(df)}")
    print(f"  x range: {df['x'].min():.1f} - {df['x'].max():.1f} um")
    print(f"  y range: {df['y'].min():.1f} - {df['y'].max():.1f} um")

    csv_path = out_dir / f"{args.fov}.csv"
    df.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path}")

    fig, ax = plt.subplots(figsize=(7, 7))
    for lab, colour in (("Tumour", "tab:red"), ("CD8", "tab:blue")):
        d = df[df["label"] == lab]
        ax.scatter(d["x"], d["y"], s=6, c=colour, alpha=0.6, label=f"{lab} (n={len(d)})")
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image coordinates
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"{args.fov}: tumour vs CD8")
    ax.legend()
    png_path = out_dir / f"{args.fov}_overlay.png"
    fig.savefig(png_path, dpi=130, bbox_inches="tight")
    print(f"  wrote {png_path}")


if __name__ == "__main__":
    main()
