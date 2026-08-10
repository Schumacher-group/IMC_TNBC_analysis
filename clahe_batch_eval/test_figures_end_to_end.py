#!/usr/bin/env python
"""End-to-end test for 04_figures.py against synthetic data.

    python test_figures_end_to_end.py

The uncorrected tables carry a per-batch intensity offset that the corrected
ones do not, so between-batch dispersion MUST fall after correction for every
marker. That makes the sign of the dispersion summary checkable rather than
merely plausible. Also builds miniature image trees so the example-region
panels are exercised, including the composite and the cell type map.
"""

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

HERE = Path(__file__).resolve().parent
MARKERS = ["Pan-keratin", "Collage-Type_I", "Ki-67", "PD-L2", "Alpha-SMA", "DNA1"]


def load_module():
    spec = importlib.util.spec_from_file_location("figures", HERE / "04_figures.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build(tmp, rng, n_cells=3000):
    import anndata as ad

    n_batches, n_types = 3, 4
    patient = rng.integers(0, 9, n_cells)
    batch = patient % n_batches
    cell_type = rng.integers(0, n_types, n_cells)
    fov = np.array([f"Leap{p:03d}_{c % 2}" for p, c in zip(patient, cell_type)])

    base = (rng.gamma(2.0, 0.4, (n_cells, len(MARKERS)))
            + cell_type[:, None] * 0.3)
    offset = rng.uniform(0.8, 2.5, (n_batches, len(MARKERS)))
    with_batch = base * offset[batch]

    obs = pd.DataFrame({
        "acquisition_ID": fov,
        "Pixie": pd.Categorical([f"type{t}" for t in cell_type]),
        "Patient_ID": patient,
        "Stain_Batch": batch,
        "Response": "Responder",
        "is_pre_treatment": True,
        "in_revision_cohort": True,
    }, index=[f"c{i}" for i in range(n_cells)])
    spatial = rng.random((n_cells, 2)) * 60

    for name, X in [("corrected", base), ("uncorrected", with_batch)]:
        a = ad.AnnData(X=X.astype(np.float32), obs=obs.copy(),
                       var=pd.DataFrame(index=MARKERS),
                       obsm={"spatial": spatial.copy()})
        a.write(tmp / f"adata_{name}.h5ad")

    proc, raw = tmp / "processed", tmp / "non_processed"
    for fov_name in np.unique(fov):
        for d in (proc, raw):
            (d / fov_name).mkdir(parents=True, exist_ok=True)
            for m in MARKERS:
                (d / fov_name / f"{m}.tiff")
        for m in MARKERS:
            img = rng.gamma(2.0, 3.0, (64, 64)).astype(np.float32)
            tifffile.imwrite(raw / fov_name / f"{m}.tiff", img)
            tifffile.imwrite(proc / fov_name / f"{m}.tiff",
                             (img * np.linspace(0.5, 1.5, 64)[:, None]).astype(np.float32))
    return proc, raw


def main():
    rng = np.random.default_rng(0)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        proc, raw = build(tmp, rng)

        mod = load_module()
        mod.PROCESSED_DIR = proc
        mod.NON_PROCESSED_DIR = raw

        out = tmp / "out"
        sys.argv = ["04_figures.py", "--in-dir", str(tmp), "--out-dir", str(out),
                    "--cohort", "all", "--n-rois", "3"]
        rc = mod.main()
        assert rc == 0, f"returned {rc}"

        disp = pd.read_csv(out / "marker_dispersion_all.csv")
        assert len(disp) == len(MARKERS), disp
        print(disp[["marker", "mean_pairwise_emd_uncorrected",
                    "mean_pairwise_emd_corrected", "emd_reduction"]].round(4)
              .to_string(index=False))
        # Ground truth: only the uncorrected copy has a per-batch offset.
        assert (disp["emd_reduction"] > 0).all(), disp
        assert (disp["sd_reduction"] > 0).all(), disp
        print("\nOK: between-batch dispersion falls for every marker, both summaries")

        # The example ROIs must come from DIFFERENT batches -- the figure's whole
        # point is that architecture survives across the batch range.
        import re
        picked = re.findall(r"- `(\S+)` -- batch (\d+)", report_text := (
            out / "figures_report_all.md").read_text())
        assert len(picked) == 3, picked
        assert len({b for _, b in picked}) == 3, f"repeated batch: {picked}"
        print(f"OK: example ROIs drawn from 3 distinct batches: {picked}")

        assert (out / "marker_distributions_all.png").exists()
        assert (out / "example_regions_all.png").exists()
        print("OK: distribution and example-region figures written")

        report = report_text
        assert "cannot bias the conclusion" in report
        assert "36 of" not in report
        print("OK: report written with the all-marker dispersion summary")

    print("\nPASS: stage 4 end-to-end")
    return 0


if __name__ == "__main__":
    sys.exit(main())
