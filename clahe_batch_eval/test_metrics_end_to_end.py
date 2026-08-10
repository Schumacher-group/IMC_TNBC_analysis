#!/usr/bin/env python
"""End-to-end test for 03_metrics.py against synthetic paired AnnData.

    python test_metrics_end_to_end.py

Builds two small stage-2-shaped tables where the ground truth is known by
construction: the "uncorrected" one has a large batch offset added to every
marker, the "corrected" one does not. Both share cell types, patients, batches
and coordinates. So the batch metrics MUST favour the corrected table, and the
cell types are recoverable from both.

This checks the metric directions are not accidentally inverted -- an easy
mistake given ASW-batch and PCR are reported flipped -- and that the script
runs end to end: scanpy embedding, Leiden, silhouettes, the patient-level
bootstrap and the figure.
"""

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def load_module():
    spec = importlib.util.spec_from_file_location("metrics", HERE / "03_metrics.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_paired(tmp, rng, n_cells=4000, n_markers=12):
    """Same cells twice; the uncorrected copy carries a strong batch offset."""
    import anndata as ad
    from sklearn import preprocessing

    n_types, n_batches, n_patients = 4, 4, 12
    cell_type = rng.integers(0, n_types, n_cells)
    patient = rng.integers(0, n_patients, n_cells)
    batch = patient % n_batches            # patients nested in batch, as in the real data

    centres = rng.normal(0, 4.0, (n_types, n_markers))
    base = centres[cell_type] + rng.normal(0, 1.0, (n_cells, n_markers))

    batch_offset = rng.normal(0, 6.0, (n_batches, n_markers))
    with_batch = base + batch_offset[batch]

    obs = pd.DataFrame({
        "acquisition_ID": [f"Leap{p:03d}_{i % 3}" for i, p in enumerate(patient)],
        "label": np.arange(n_cells),
        "Pixie": pd.Categorical([f"type{t}" for t in cell_type]),
        "Patient_ID": patient,
        "Stain_Batch": batch,
        "Response": np.where(patient % 2 == 0, "Responder", "Non-Responder"),
        "is_pre_treatment": True,
        "in_revision_cohort": True,
    }, index=[f"cell{i}" for i in range(n_cells)])

    paths = {}
    for name, X in [("corrected", base), ("uncorrected", with_batch)]:
        a = ad.AnnData(X=X.astype(np.float32), obs=obs.copy(),
                       var=pd.DataFrame(index=[f"M{j}" for j in range(n_markers)]),
                       obsm={"spatial": rng.random((n_cells, 2))})
        a.layers["raw"] = a.X.copy()
        a.layers["scaled"] = preprocessing.StandardScaler().fit_transform(a.X)
        p = tmp / f"adata_{name}.h5ad"
        a.write(p)
        paths[name] = p
    return paths


def main():
    rng = np.random.default_rng(0)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        build_paired(tmp, rng)

        mod = load_module()
        out = tmp / "out"
        sys.argv = ["03_metrics.py", "--in-dir", str(tmp), "--out-dir", str(out),
                    "--cohort", "all", "--n-graph-cells", "4000",
                    "--n-asw-cells", "1500", "--n-bootstrap", "50"]
        rc = mod.main()
        assert rc == 0, f"returned {rc}"

        table = pd.read_csv(out / "metrics_all.csv")
        print("\n" + table[["metric", "corrected", "uncorrected", "difference",
                            "ci_low", "ci_high"]].round(4).to_string(index=False))

        assert len(table) == 7, table
        # The interval must bracket the point estimate. A percentile interval
        # did not, for ARI and NMI, because patient resampling duplicates whole
        # patients and biases those statistics downward.
        finite = table[table["ci_low"].notna()]
        assert (finite["ci_low"] <= finite["difference"]).all(), finite
        assert (finite["difference"] <= finite["ci_high"]).all(), finite
        assert table["corrected"].notna().all(), table
        assert table["difference"].notna().all(), table

        batch = table[table["family"] == "batch"].set_index("metric")
        # Ground truth: the uncorrected copy has a real batch offset, so every
        # batch metric must favour corrected. Catches an inverted sign.
        for metric, row in batch.iterrows():
            assert row["difference"] > 0, f"{metric} should favour corrected: {row.to_dict()}"
        print("\nOK: all 4 batch metrics favour the table without the batch offset")

        bio = table[table["family"] == "biological"].set_index("metric")
        assert (bio["corrected"] > 0).all(), bio
        print("OK: biological conservation metrics computed and positive")

        assert (out / "metric_panel_all.png").exists()
        assert (out / "embedding_all.json").exists()
        report = (out / "metrics_report_all.md").read_text()
        assert "Pixie labels were derived" in report, "caveat missing from report"
        assert "Patient identity is not used" in report, "patient note missing"
        print("OK: figure, params and report written; caveats present")

        params = json.loads((out / "embedding_all.json").read_text())
        assert "clustering" in params, params
        # Both tables must share one clustering method, decided before either ran.
        assert len(set(params["leiden_clusters"])) == 2, params
        print(f"OK: single clustering method for both tables ({params['clustering']})")

    print("\nPASS: stage 3 end-to-end")
    return 0


if __name__ == "__main__":
    sys.exit(main())
