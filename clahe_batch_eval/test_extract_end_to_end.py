#!/usr/bin/env python
"""End-to-end test for 02_extract.py against synthetic data.

    python test_extract_end_to_end.py

Builds a miniature version of the workstation layout -- processed/ and
non_processed/ image trees, Mesmer masks, a cell table whose marker columns
are genuine mean intensities over those masks, and a metadata CSV -- then runs
the real extraction and checks the invariants the paired design depends on:

  * the verification step passes when the extraction matches the table, and
    FAILS when the table is perturbed (otherwise it would rubber-stamp anything)
  * both outputs hold the same cells, in the same order, with the same labels,
    Pixie annotations and coordinates
  * the corrected table reproduces the published intensities
  * the uncorrected table really does differ (it came from different pixels)
  * Carboplatin, present only in processed/, is dropped from both
  * the cache makes a second run reuse rather than re-extract
"""

import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
SHARED = ["DNA1", "Pan-keratin", "CD8a", "Collage-Type_I"]
DERIVED = ["Carboplatin"]
FOVS = ["Leap001_1", "Leap001_2", "Leap002_1"]


def load_module():
    spec = importlib.util.spec_from_file_location("extract", HERE / "02_extract.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_fixtures(root, rng):
    proc, raw, masks = root / "processed", root / "non_processed", root / "masks"
    for d in (proc, raw, masks):
        d.mkdir(parents=True, exist_ok=True)

    rows = []
    for fov in FOVS:
        mask = np.zeros((64, 64), dtype=np.uint16)
        obj = 1
        for r in range(2, 60, 10):
            for c in range(2, 60, 10):
                mask[r:r + 7, c:c + 7] = obj
                obj += 1
        tifffile.imwrite(masks / f"{fov}_whole_cell.tiff", mask)

        (proc / fov).mkdir(exist_ok=True)
        (raw / fov).mkdir(exist_ok=True)
        labels = np.unique(mask)[1:]
        record = {"fov": fov, "label": labels}

        for ch in SHARED:
            raw_img = rng.gamma(2.0, 3.0, size=(64, 64)).astype(np.float32)
            # A deliberately position-dependent transform, as CLAHE would be.
            yy = np.linspace(0.5, 1.5, 64)[:, None]
            proc_img = (raw_img * yy).astype(np.float32)
            tifffile.imwrite(raw / fov / f"{ch}.tiff", raw_img)
            tifffile.imwrite(proc / fov / f"{ch}.tiff", proc_img)
            record[ch] = ndi.mean(proc_img, labels=mask, index=labels)

        for ch in DERIVED:  # exists only in processed/
            img = rng.gamma(2.0, 3.0, size=(64, 64)).astype(np.float32)
            tifffile.imwrite(proc / fov / f"{ch}.tiff", img)
            record[ch] = ndi.mean(img, labels=mask, index=labels)

        area = ndi.sum(np.ones_like(mask, np.float32), labels=mask, index=labels)
        coms = ndi.center_of_mass(np.ones_like(mask, bool), mask, labels)
        record["cell_size"] = area
        record["area"] = area
        record["centroid-0"] = [c[0] for c in coms]
        record["centroid-1"] = [c[1] for c in coms]
        record["cell_meta_cluster"] = rng.choice(["Cancer cell", "CD8 T cell"], len(labels))
        rows.append(pd.DataFrame(record))

    table = pd.concat(rows, ignore_index=True)
    # Column order must match the real table: cell_size, markers..., label, ...
    ordered = (["cell_size"] + SHARED + DERIVED
               + ["label", "area", "centroid-0", "centroid-1", "fov", "cell_meta_cluster"])
    table = table[ordered]
    table_path = root / "cell_table.csv"
    table.to_csv(table_path, index=True)

    meta = pd.DataFrame({
        "LEAP_ID": ["LEAP001", "LEAP002"],
        "Sample_Type_(pre/post treatment)": ["pre", "post"],
        "Response": ["Responder", "Non-Responder"],
        "Patient_ID": [1, 2],
        "Stain_Batch": [1, 2],
    })
    meta_path = root / "metadata.csv"
    meta.to_csv(meta_path, index=False)
    return proc, raw, masks, table_path, meta_path


def run(mod, out_dir, table_path, meta_path, masks, extra=()):
    sys.argv = ["02_extract.py", "--out-dir", str(out_dir),
                "--cell-table", str(table_path), "--metadata", str(meta_path),
                "--mask-dir", str(masks), *extra]
    return mod.main()


def main():
    rng = np.random.default_rng(0)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        proc, raw, masks, table_path, meta_path = build_fixtures(root, rng)

        mod = load_module()
        mod.PROCESSED_DIR = proc
        mod.NON_PROCESSED_DIR = raw

        # --- the verification must reject a table it cannot reproduce ---
        bad = pd.read_csv(table_path, index_col=0)
        bad[SHARED[0]] = bad[SHARED[0]] * 1.5 + 3.0
        bad_path = root / "cell_table_corrupted.csv"
        bad.to_csv(bad_path, index=True)
        rc = run(mod, root / "out_bad", bad_path, meta_path, masks, ["--verify-only"])
        assert rc == 1, "verification passed on a table it should not reproduce"
        print("OK: verification rejects a mismatched cell table")

        # --- and accept the real one ---
        rc = run(mod, root / "out", table_path, meta_path, masks, ["--verify-only"])
        assert rc == 0, "verification failed on a table it should reproduce"
        print("OK: verification accepts the correct cell table")

        # --- full extraction ---
        rc = run(mod, root / "out", table_path, meta_path, masks)
        assert rc == 0, f"extraction returned {rc}"

        import anndata as ad
        cor = ad.read_h5ad(root / "out" / "adata_corrected.h5ad")
        unc = ad.read_h5ad(root / "out" / "adata_uncorrected.h5ad")

        assert list(cor.var_names) == SHARED_SORTED, cor.var_names
        assert "Carboplatin" not in cor.var_names and "Carboplatin" not in unc.var_names
        print(f"OK: {len(cor.var_names)} shared channels, Carboplatin dropped from both")

        assert list(cor.obs_names) == list(unc.obs_names)
        assert (cor.obs["Pixie"].astype(str).values == unc.obs["Pixie"].astype(str).values).all()
        assert np.array_equal(cor.obsm["spatial"], unc.obsm["spatial"])
        assert cor.n_obs == len(pd.read_csv(table_path))
        print(f"OK: {cor.n_obs} cells paired identically across both tables")

        published = pd.read_csv(table_path, index_col=0)
        want = published.sort_values("fov", kind="stable")[SHARED_SORTED].to_numpy()
        assert np.allclose(cor.layers["raw"], want, rtol=1e-5), "corrected != published"
        print("OK: corrected table reproduces the published intensities")

        assert not np.allclose(cor.layers["raw"], unc.layers["raw"]), \
            "uncorrected is identical to corrected -- wrong directory read?"
        print("OK: uncorrected table genuinely differs")

        assert cor.obs["is_pre_treatment"].sum() > 0
        assert (~cor.obs["is_pre_treatment"]).sum() > 0
        print("OK: cohort flags recorded, nothing filtered out")

        # --- cache reuse ---
        report = (root / "out" / "extract_report.md").read_text()
        assert "newly extracted: 3" in report, report
        rc = run(mod, root / "out", table_path, meta_path, masks, ["--skip-verify"])
        report = (root / "out" / "extract_report.md").read_text()
        assert "reused from cache: 3" in report, report
        print("OK: second run reuses the cache instead of re-extracting")

    print("\nPASS: stage 2 end-to-end")
    return 0


SHARED_SORTED = sorted(SHARED)

if __name__ == "__main__":
    sys.exit(main())
