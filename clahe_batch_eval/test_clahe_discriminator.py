#!/usr/bin/env python
"""Regression test for check 1's global-vs-local discriminator.

    python test_clahe_discriminator.py

Why this exists: the first workstation run reported Spearman rho ~0.99 and a
+2.8 bit entropy gain, and those were read as confirming CLAHE. They do not.
The processed data is clipped (max exactly 1.0), and clipping ties pixels,
which lowers rho on its own; separately, stretching a skewed range onto [0,1]
raises binned entropy on its own. A plain global quantile normalisation
reproduces both numbers without a trace of CLAHE.

This builds an image with strong large-scale intensity structure -- the case
where local and global transforms actually diverge -- and puts both candidate
transforms through locality_stats(), asserting it separates them. The global
case uses exactly the transform this pipeline applies elsewhere:
clip(raw / quantile(raw, 0.95), 0, 1).

Two earlier attempts at this statistic failed here and were discarded: binning
raw values gave a global rescale a spurious non-zero score, because proc varies
within a bin whenever raw does; and comparing pixels of exactly equal raw value
collapsed, because CLAHE maps each tile's minimum to zero and these images are
about half zeros, so the only repeated value carries no signal. The isotonic
residual has neither failure mode -- it is exactly zero for any global monotone
map, clipped or not, with no binning and no parametric assumption.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
from skimage.exposure import equalize_adapthist

HERE = Path(__file__).resolve().parent


def load_reconcile(path=HERE / "01_reconcile.py"):
    spec = importlib.util.spec_from_file_location("reconcile", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_raw(size=512, seed=0):
    """Speckled foreground on a strong intensity gradient, like a real IMC field."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:size, 0:size] / size
    background = 4.0 * yy + 1.0 * xx           # large-scale shading
    speckle = rng.gamma(shape=0.6, scale=3.0, size=(size, size))
    img = background * speckle
    img[rng.random((size, size)) < 0.45] = 0.0  # IMC images are mostly zeros
    return img


def global_rescale(raw, q=0.95):
    """The innocent explanation: per-channel quantile normalisation with clipping."""
    denom = np.quantile(raw, q)
    return np.clip(raw / denom, 0, 1)


def clahe(raw):
    scaled = raw / raw.max()
    return equalize_adapthist(scaled, kernel_size=64, clip_limit=0.01)


def main():
    recon = load_reconcile()
    raw = synthetic_raw()

    g = recon.locality_stats(raw.ravel(), global_rescale(raw).ravel())
    c = recon.locality_stats(raw.ravel(), clahe(raw).ravel())

    for name, s in [("global quantile rescale", g), ("CLAHE", c)]:
        print(f"{name:>24}: distinct_at_raw_zero={s['n_distinct_at_raw_zero']}, "
              f"monotone_residual={s['monotone_residual']:.4f}, "
              f"frac_saturated={s['frac_saturated']:.4f}")

    # Any global monotone map is fit exactly by isotonic regression.
    assert g["monotone_residual"] < 1e-9, g
    # CLAHE cannot be: local context decides where a given raw value lands.
    assert c["monotone_residual"] > 0.01, c
    assert c["monotone_residual"] > 1000 * max(g["monotone_residual"], 1e-12), (g, c)

    print("\nPASS: locality_stats separates a global rescale from CLAHE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
