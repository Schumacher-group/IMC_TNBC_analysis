#!/usr/bin/env python
"""Regression test for check 5a's mask-directory identification.

Run with no arguments; it builds its own fixtures in a temp directory.

    python test_mask_identification.py

Why this exists: identifying which segmentation run produced the cell table is
the load-bearing check of stage 1 -- stage 2 extracts uncorrected intensities
using those masks, and a wrong choice would still run, just against the wrong
cells. The first workstation run reported "no candidate matches" when in fact
one matched exactly, because two candidate directories on that machine share
the basename `deepcell_output` and the summary grouped on it. This reproduces
that layout.

The decoy is deliberately hard: same FOV names, same object count, 100% label
coverage. It differs only by a 3 px centroid shift, so name- or count-based
matching accepts it and only the centroid comparison rejects it.
"""

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent


def load_reconcile(path=HERE / "01_reconcile.py"):
    spec = importlib.util.spec_from_file_location("reconcile", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_fixtures(root):
    """Two mask dirs sharing a basename, plus a cell table derived from one."""
    right = root / "sept2024_release" / "deepcell_output"   # made the cell table
    wrong = root / "deepcell_output"                        # a different run
    fovs = ["Leap001_1", "Leap001_2"]
    rows = []

    for mask_dir, offset in [(right, 0), (wrong, 3)]:
        mask_dir.mkdir(parents=True, exist_ok=True)
        for fov in fovs:
            mask = np.zeros((64, 64), dtype=np.uint16)
            obj = 1
            for r in range(4, 60, 12):
                for c in range(4, 60, 12):
                    rr, cc = r + offset, c + offset
                    if rr + 6 < 64 and cc + 6 < 64:
                        mask[rr:rr + 6, cc:cc + 6] = obj
                        obj += 1
            tifffile.imwrite(mask_dir / f"{fov}_whole_cell.tiff", mask)

            if mask_dir is right:
                labels = np.unique(mask)[1:]
                coms = ndi.center_of_mass(np.ones_like(mask, bool), mask, labels)
                for lab, (y, x) in zip(labels, coms):
                    rows.append({"fov": fov, "label": int(lab),
                                 "centroid-0": y, "centroid-1": x})

    table = root / "cell_table.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    return right, wrong, fovs, table


def main():
    recon = load_reconcile()
    with tempfile.TemporaryDirectory() as tmp:
        right, wrong, fovs, table = build_fixtures(Path(tmp))

        recon.MASK_DIR_CANDIDATES = [right, wrong]
        recon.CELL_TABLE = table
        rep = recon.Report()
        chosen, mask_fovs = recon.identify_mask_dir(rep, fovs, tifffile)

        assert chosen == right, f"expected {right}, got {chosen}"
        assert not rep.failures, rep.failures
        assert mask_fovs == set(fovs), mask_fovs

    print("\nPASS: correct mask directory identified; same-basename decoy rejected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
