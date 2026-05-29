#!/usr/bin/env python
"""End-to-end install verification for M2S2_demo.

Builds a tiny synthetic 2-colour point cloud and pushes it through the exact
`make_cmdline_parser` + `generate_stats` pipeline (chalc chromatic filtration ->
six-pack persistence -> persistent statistics -> parquet write). This exercises
the full install without depending on downloading the authors' datasets.

NOTE: generate_stats() spawns a worker pool, so on macOS/Python 3.13 (spawn
start method) the body MUST be under `if __name__ == "__main__":`. generate_stats
also ends with sys.exit(0); run as a subprocess and check the exit code.
"""

import multiprocessing as mp
import os
import sys
from pathlib import Path

# M2S2's PoolWithLogger passes a local closure as the worker initialiser, which the
# default macOS/py3.13 "spawn" start method cannot pickle. The authors run on Linux
# (fork). Force fork here (in our launcher, not M2S2 core) to match their assumption.
os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
M2S2 = HERE / "M2S2_demo"
sys.path.insert(0, str(M2S2))

from utils import generate_stats, make_cmdline_parser, make_filename_index_map  # noqa: E402


def main() -> None:
    mp.set_start_method("fork", force=True)
    data_dir = HERE / "_install_test" / "data"
    out_dir = HERE / "_install_test" / "stats"
    log_dir = HERE / "_install_test" / "logs"
    for d in (data_dir, out_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Synthetic ROI: two overlapping Gaussian blobs, labels A and B, ~150 cells each.
    rng = np.random.default_rng(0)
    a = rng.normal(loc=(0, 0), scale=30, size=(150, 2))
    b = rng.normal(loc=(40, 40), scale=30, size=(150, 2))
    df = pd.DataFrame(
        {
            "x": np.concatenate([a[:, 0], b[:, 0]]),
            "y": np.concatenate([a[:, 1], b[:, 1]]),
            "label": ["A"] * len(a) + ["B"] * len(b),
        }
    )
    df.to_csv(data_dir / "synthetic.csv", index=False)

    args_override = {
        "Required arguments": {
            "--dataset-dir": {"default": str(data_dir), "required": False},
            "--output-dir": {"default": str(out_dir), "required": False},
        },
        "Logging options": {"--logfile-dir": {"default": str(log_dir)}},
    }
    parser = make_cmdline_parser(args_override)
    cmdline_args = parser.parse_args(
        [
            "--labels-include", "A", "B",
            "--min-num-labels", "1",
            "--max-num-labels", "2",
            "--num-workers", "1",
            "--no-resume",
            "--disable-progressbar",
            "--max-diagram-dimension", "1",
        ]
    )

    index_map = make_filename_index_map(lambda fn: {"roi": fn})
    generate_stats(cmdline_args, index_map)


if __name__ == "__main__":
    main()
