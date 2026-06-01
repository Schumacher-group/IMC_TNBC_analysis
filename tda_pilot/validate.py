#!/usr/bin/env python
"""Phase 3 validation of the full-pilot stats parquets.

For each of the 7 pair definitions, checks:
  1. completeness  : completed + degenerate + failed == 808
  2. finiteness    : all 456 stat columns finite for every 'ok' row
  3. non-degeneracy: ker0/ker1/im0/im1 num_bars > 0 for >= 90% of 'ok' ROIs
Plus a raw-diagram spot check (one responder + one non-responder, CD8_primary):
overlay of domain/image/kernel/cokernel degree-1 diagrams -> output/spotcheck_diagrams.png

Writes output/VALIDATION.md.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pair_defs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "M2S2_demo"))
from chalc.sixpack import SixPack  # noqa: E402

STATS_DIR = HERE / "output" / "stats"
DIAG_DIR = HERE / "output" / "diagrams" / "CD8_primary"
N_COHORT = 808
STAT_COLS_EXPECTED = 456

meta = pd.read_csv(HERE / "per_roi_counts.csv")[["fov", "Patient_ID", "Response"]]
meta = meta.rename(columns={"fov": "roi_id"})


def df_to_md(df: pd.DataFrame) -> str:
    """Minimal markdown table (no tabulate dependency)."""
    cols = list(df.columns)
    head = "| " + " | ".join(map(str, cols)) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = "\n".join("| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False))
    return "\n".join([head, sep, body])


def validate_definition(defname: str) -> dict:
    df = pd.read_parquet(STATS_DIR / f"{defname}.parquet")
    statcols = [c for c in df.columns if "-" in c]
    fails = pd.read_csv(STATS_DIR / f"{defname}_failures.csv") if (STATS_DIR / f"{defname}_failures.csv").exists() else pd.DataFrame()
    n_failed = len(fails)
    ok = df[df["status"] == "ok"]
    degen = df[df["status"] == "degenerate"]

    # 1. completeness
    total_accounted = len(ok) + len(degen) + n_failed
    complete = total_accounted == N_COHORT and df["roi_id"].nunique() == len(df)

    # 2. finiteness (ok rows only; degenerate have NaN stats by design)
    arr = ok[statcols].to_numpy(dtype=float)
    n_nonfinite_rows = int((~np.isfinite(arr).all(axis=1)).sum())
    finite_ok = n_nonfinite_rows == 0 and len(statcols) == STAT_COLS_EXPECTED

    # 3. non-degeneracy of kernel & image in dim 0 and 1
    nd = {}
    for col in ("ker0-num_bars", "ker1-num_bars", "im0-num_bars", "im1-num_bars"):
        frac = float((ok[col] > 0).mean()) if len(ok) else 0.0
        nd[col] = round(frac, 4)
    nondegen_ok = all(v >= 0.90 for v in nd.values())

    return {
        "definition": defname, "rows": len(df), "ok": len(ok),
        "degenerate": len(degen), "failed": n_failed,
        "completeness_808": complete, "stat_cols": len(statcols),
        "finite_all_ok": finite_ok, "nonfinite_rows": n_nonfinite_rows,
        "nondegen_ok": nondegen_ok, **nd,
    }


def spotcheck_plot() -> tuple[str, str]:
    m = meta.copy()
    resp = m[m["Response"] == "Responder"]["roi_id"]
    nonr = m[m["Response"] == "Non-Responder"]["roi_id"]
    r_roi = next(r for r in ("Leap066_11", *resp) if (DIAG_DIR / f"{r}.h5").exists())
    n_roi = next(r for r in nonr if (DIAG_DIR / f"{r}.h5").exists())

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for row, (roi, label) in enumerate(((r_roi, "Responder"), (n_roi, "Non-Responder"))):
        import h5py
        with h5py.File(DIAG_DIR / f"{roi}.h5", "r") as fh:
            dgms = SixPack.from_file(fh)
        for col, name in enumerate(("dom", "im", "ker", "cok")):
            ax = axes[row, col]
            bar = dgms.get_matrix(name, [1])[0]
            bar = np.asarray(bar) if bar is not None else np.zeros((0, 2))
            fin = bar[np.isfinite(bar).all(axis=1)] if len(bar) else bar
            if len(fin):
                ax.scatter(fin[:, 0], fin[:, 1], s=8, alpha=0.5)
                mx = float(np.nanmax(fin)) if len(fin) else 1.0
                ax.plot([0, mx], [0, mx], "k--", lw=0.6)
            ax.set_title(f"{label} {roi}\n{name} dim1 (n={len(fin)})", fontsize=9)
            ax.set_xlabel("birth"); ax.set_ylabel("death")
    fig.suptitle("CD8_primary degree-1 six-pack diagrams (responder vs non-responder)", fontsize=12)
    fig.tight_layout()
    out = HERE / "output" / "spotcheck_diagrams.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    return r_roi, n_roi


def main() -> None:
    rows = [validate_definition(d) for d in pair_defs.DEFINITIONS]
    res = pd.DataFrame(rows)
    r_roi, n_roi = spotcheck_plot()

    lines = ["# Phase 3 validation — full TDA pilot\n",
             f"Cohort: {N_COHORT} ROIs. Stat columns expected: {STAT_COLS_EXPECTED} "
             "(6 diagrams x 2 dims x 38 stats).\n",
             "## Per-definition checks\n"]
    show = ["definition", "ok", "degenerate", "failed", "completeness_808",
            "finite_all_ok", "nondegen_ok", "ker1-num_bars", "im1-num_bars"]
    lines.append(res[show].pipe(df_to_md))
    lines.append("\n### Non-degeneracy detail (fraction of ok ROIs with >0 bars)\n")
    lines.append(res[["definition", "ker0-num_bars", "ker1-num_bars",
                      "im0-num_bars", "im1-num_bars"]].pipe(df_to_md))
    all_complete = bool(res["completeness_808"].all())
    all_finite = bool(res["finite_all_ok"].all())
    all_nondegen = bool(res["nondegen_ok"].all())
    lines.append("\n## Summary\n")
    lines.append(f"- Completeness (sum==808, unique roi_id) for all 7: **{all_complete}**")
    lines.append(f"- Finiteness (456 finite stats per ok row) for all 7: **{all_finite}**")
    lines.append(f"- Non-degeneracy (ker/im dim0&1 >0 in >=90% ok ROIs) for all 7: **{all_nondegen}**")
    lines.append(f"- Degenerate ROIs (a species had <3 cells): "
                 + ", ".join(f"{r.definition}={r.degenerate}" for r in res.itertuples() if r.degenerate) + ".")
    lines.append(f"- Spot-check diagram overlay saved: `output/spotcheck_diagrams.png` "
                 f"(responder {r_roi}, non-responder {n_roi}).")
    (HERE / "output" / "VALIDATION.md").write_text("\n".join(lines) + "\n")

    token = "ALLPASS" if (all_complete and all_finite and all_nondegen) else "CHECK_FAILED"
    print("VALIDATION_TOKEN:", token)
    res.to_csv(STATS_DIR / "validation_table.csv", index=False)


if __name__ == "__main__":
    main()
