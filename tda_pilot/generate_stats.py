#!/usr/bin/env python
"""Task 4 wrapper: persistent statistics for the TNBC smoke-test ROI.

Adapted from M2S2_demo/colorectal_cancer_dataset_analysis/generate_stats.py.
Only the override dictionary + filename->metadata map are changed:
  - dataset-dir = tda_pilot/data (single ROI prepared by adapter.py)
  - output-dir  = tda_pilot/output/stats
  - labels-include = Tumour CD8   (the pair only; with only 2 labels present,
    min/max-num-labels 1/2 yields exactly the singletons + the (Tumour,CD8) pair,
    no triples -> "k=2" single gluing map)
  - max-diagram-dimension 1 (degrees 0 and 1), num-workers 1 (clean timing),
    no-resume (force recompute), persistent-statistics vectorisation (default).

macOS/py3.13 note: force the 'fork' start method (M2S2's PoolWithLogger uses a
local-closure pool initialiser that the default 'spawn' cannot pickle). This is
configuration in our launcher, not a change to M2S2 core.
"""

import multiprocessing as mp
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

import pandas as pd

import cohort

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))

from utils import generate_stats, make_cmdline_parser, make_filename_index_map  # noqa: E402

META = HERE.parent / "CellTable_CleanCohort" / "CleanCohort_Metadata.csv"


def _response_map() -> dict[str, str]:
    meta = pd.read_csv(META)
    meta["LEAP_ID"] = meta["LEAP_ID"].str.upper()
    return dict(zip(meta["LEAP_ID"], meta["Response"].astype(str)))


if __name__ == "__main__":
    mp.set_start_method("fork", force=True)
    resp_map = _response_map()

    args_override = {
        "Required arguments": {
            "--dataset-dir": {"default": str(HERE / "data"), "required": False},
            "--output-dir": {"default": str(HERE / "output" / "stats"), "required": False},
        },
        "Dataset-specific options": {
            "--x-column": {"default": "x"},
            "--y-column": {"default": "y"},
            "--label-column": {"default": "label"},
            "--labels-include": {"choices": ["Tumour", "CD8"]},
            "--labels-exclude": {"choices": ["Tumour", "CD8"]},
        },
        "Logging options": {
            "--logfile-dir": {"default": str(HERE / "output" / "logs")},
        },
    }
    Path(args_override["Required arguments"]["--output-dir"]["default"]).mkdir(
        parents=True, exist_ok=True
    )
    Path(args_override["Logging options"]["--logfile-dir"]["default"]).mkdir(
        parents=True, exist_ok=True
    )

    parser = make_cmdline_parser(args_override)
    cmdline_args = parser.parse_args()

    # Restrict to the revision cohort (N=808): if the user did not pass an explicit
    # --files-list, process only the cohort-member CSVs present in the dataset dir.
    if not cmdline_args.files_list:
        dataset_dir = Path(cmdline_args.dataset_dir)
        cohort_files = sorted(
            p.name for p in dataset_dir.glob("*.csv") if cohort.is_cohort_member(p.stem)
        )
        skipped = [p.name for p in dataset_dir.glob("*.csv") if not cohort.is_cohort_member(p.stem)]
        if skipped:
            print(f"Skipping {len(skipped)} non-cohort ROI(s): {skipped}")
        cmdline_args.files_list = cohort_files

    def _mapping_func(filename: str) -> dict[str, Any]:
        # filename like "Leap066_11" -> patient "LEAP066", acquisition "11"
        tokens = filename.split("_")
        leap = tokens[0].upper()
        acq = tokens[-1]
        return {
            "patient_id": leap,
            "acquisition": acq,
            "response": resp_map.get(leap, "unknown"),
            "filtration_algorithm": cmdline_args.filtration_algorithm,
        }

    filename_index_map = make_filename_index_map(_mapping_func)
    generate_stats(cmdline_args, filename_index_map)
