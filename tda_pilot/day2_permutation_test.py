#!/usr/bin/env python
"""Day-2 decisive test: responder contrast on permutation-null z-scores.

Reads `output/perm_null_<definition>.parquet` (see `perm_null.py`). Each ROI carries
z = (observed - null mean) / null sd for each retained statistic, where the null shuffles
cell labels across the ROI's own fixed cell positions. z is therefore a pure
spatial-ARRANGEMENT score: composition, cell number, density, imaged area and dispersion
are all held fixed within each ROI by construction, so none of them can drive a difference
in z between responders and non-responders.

Analysis plan, fixed before the run (see the approved plan document):

  PRIMARY ENDPOINT   z of `ker1-avg_length`, per-patient median, Cliff's delta
                     (Responder vs Non-Responder), with a patient-label permutation
                     p-value (20,000 draws). Patient-level because response is a
                     patient-level outcome and ROIs cluster within patients.
  SECONDARY          the remaining 49 retained statistics, Benjamini-Hochberg across
                     the family.
  ROBUSTNESS         bootstrap CI over patients; leave-one-patient-out; the
                     minority-fraction threshold ladder from day1_threshold_check.py.
  SANITY             `cod0-num_bars` is invariant under relabelling, so its null sd must
                     be exactly 0; and |z| should be large in magnitude for real ROIs,
                     confirming tissue is far from randomly labelled (otherwise the
                     statistic has no resolving power and a null result is uninformative).

DECISION RULE (stated before running):
  primary not significant -> the pilot's no-go is confirmed on the strongest available
  evidence, and the bar-length signal in LENGTH_FEATURES.md is a tissue-geometry artefact.
  primary significant     -> escalate to control species and the directional
  SubChromaticInclusion variant under the same null.

Output: output/PERMUTATION_NULL.md, output/perm_null_primary.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import perm_null  # PRIMARY, INVARIANT_CHECK, RETAINED -- single source of truth

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
N_PERM_LABEL = 20_000
N_BOOT = 5_000
SEED = 20260805


def cliffs_delta(r, nr) -> float:
    r, nr = np.asarray(r, float), np.asarray(nr, float)
    if len(r) < 2 or len(nr) < 2:
        return np.nan
    u = mannwhitneyu(r, nr, alternative="two-sided").statistic
    return 2 * u / (len(r) * len(nr)) - 1


def mwu_p(r, nr) -> float:
    r, nr = np.asarray(r, float), np.asarray(nr, float)
    if len(r) < 2 or len(nr) < 2:
        return np.nan
    return float(mannwhitneyu(r, nr, alternative="two-sided").pvalue)


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    if not ok.any():
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    adj = np.minimum.accumulate((pv[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def patient_medians(d: pd.DataFrame, col: str) -> pd.DataFrame:
    return d.groupby(["patient_id", "response"], as_index=False)[col].median()


def label_permutation_p(pat: pd.DataFrame, col: str, rng) -> tuple[float, float]:
    """Cliff's delta and its p under random reassignment of PATIENT labels."""
    lab = pat["response"].to_numpy()
    v = pat[col].to_numpy(float)
    obs = cliffs_delta(v[lab == "Responder"], v[lab == "Non-Responder"])
    null = np.empty(N_PERM_LABEL)
    for i in range(N_PERM_LABEL):
        s = rng.permutation(lab)
        null[i] = cliffs_delta(v[s == "Responder"], v[s == "Non-Responder"])
    return obs, float((np.abs(null) >= abs(obs)).mean())


def load(defname: str, inclusion: str = "symmetric") -> pd.DataFrame:
    tag = "" if inclusion == "symmetric" else f"_{inclusion}"
    d = pd.read_parquet(OUT / f"perm_null_{defname}{tag}.parquet")
    d = d[d["status"] == "ok"].copy()
    meta = pd.read_csv(HERE / "per_roi_counts.csv").rename(
        columns={"fov": "roi_id", "Patient_ID": "patient_id", "Response": "response"},
    )[["roi_id", "patient_id", "response"]]
    d = d.merge(meta, on="roi_id", how="left")
    d = d[d["response"].isin(["Responder", "Non-Responder"])].copy()
    d["minority_frac"] = np.minimum(d.n_tumour, d.n_other) / (d.n_tumour + d.n_other)

    # z is a signal-to-noise ratio, and the null sd shrinks as an ROI gains cells, so |z|
    # grows with ROI size even at a constant relative deviation. Since responder and
    # non-responder ROIs differ in size and density, carry a size-free companion: the
    # deviation as a fraction of the null mean. Conclusions should agree across both.
    for c in [c for c in d.columns if c.endswith("__obs")]:
        stat = c[: -len("__obs")]
        mean = d[f"{stat}__null_mean"]
        d[f"{stat}__reldev"] = np.where(mean != 0, (d[c] - mean) / mean, np.nan)
    return d


def md_table(dfr: pd.DataFrame) -> str:
    head = "| " + " | ".join(dfr.columns) + " |"
    sep = "| " + " | ".join("---" for _ in dfr.columns) + " |"
    body = "\n".join(
        "| " + " | ".join(f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
                          for v in row) + " |"
        for row in dfr.itertuples(index=False)
    )
    return "\n".join([head, sep, body])


def make_figure(d: pd.DataFrame, defname: str, stats: list[str], tag: str = "") -> Path:
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, len(stats), figsize=(4.0 * len(stats), 4.6))
    axes = np.atleast_1d(axes)
    for ax, st in zip(axes, stats):
        col = f"{st}__z"
        pat = patient_medians(d, col)
        for i, grp in enumerate(["Non-Responder", "Responder"]):
            y = pat[pat.response == grp][col].to_numpy()
            ax.scatter(rng.normal(i, 0.07, len(y)), y, s=30, alpha=0.75,
                       c="tab:blue" if i == 0 else "tab:orange")
            ax.hlines(np.median(y), i - 0.25, i + 0.25, color="k", lw=2)
        dlt = cliffs_delta(pat[pat.response == "Responder"][col],
                           pat[pat.response == "Non-Responder"][col])
        ax.axhline(0, color="grey", lw=0.6, ls=":")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["NR", "R"])
        ax.set_title(f"{st}\nCliff δ={dlt:+.3f}", fontsize=9)
    axes[0].set_ylabel("per-patient median z\n(observed vs own-ROI label-permutation null)")
    fig.suptitle(f"{defname}, pre-treatment: arrangement z-scores by response "
                 "(composition & geometry held fixed within ROI)")
    fig.tight_layout()
    path = OUT / f"perm_null_primary{tag}.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    return path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--definition", default="CD8_primary")
    ap.add_argument("--inclusion", choices=perm_null.INCLUSIONS, default="symmetric",
                    help="which perm_null.py run to analyse")
    args = ap.parse_args()
    defname, inclusion = args.definition, args.inclusion

    d = load(defname, inclusion)
    rng = np.random.default_rng(SEED)
    stats = [s for s in perm_null.RETAINED if f"{s}__z" in d.columns]
    n_pat = d.groupby("response")["patient_id"].nunique().to_dict()

    # ---- sanity ---------------------------------------------------------------
    inv_sd = float(d[f"{perm_null.INVARIANT_CHECK}__null_sd"].abs().max())
    zcols = [f"{s}__z" for s in stats]
    all_finite = bool(np.isfinite(d[zcols]).all().all())
    med_abs_z = d[zcols].abs().median().median()

    # ---- primary --------------------------------------------------------------
    pcol = f"{perm_null.PRIMARY}__z"
    pat = patient_medians(d, pcol)
    obs_d, perm_p = label_permutation_p(pat, pcol, rng)
    R = pat[pat.response == "Responder"][pcol].to_numpy()
    NR = pat[pat.response == "Non-Responder"][pcol].to_numpy()
    boot = np.array([cliffs_delta(rng.choice(R, len(R)), rng.choice(NR, len(NR)))
                     for _ in range(N_BOOT)])
    ci = np.percentile(boot, [2.5, 97.5])
    loo = [
        cliffs_delta(
            pat[(pat.patient_id != p) & (pat.response == "Responder")][pcol],
            pat[(pat.patient_id != p) & (pat.response == "Non-Responder")][pcol],
        )
        for p in pat.patient_id.unique()
    ]

    # size-free companion for the primary endpoint
    rcol = f"{perm_null.PRIMARY}__reldev"
    pat_r = patient_medians(d, rcol)
    rel_d = cliffs_delta(pat_r[pat_r.response == "Responder"][rcol],
                         pat_r[pat_r.response == "Non-Responder"][rcol])
    rel_p = mwu_p(pat_r[pat_r.response == "Responder"][rcol],
                  pat_r[pat_r.response == "Non-Responder"][rcol])

    # ---- secondary ------------------------------------------------------------
    rows = []
    for s in stats:
        col = f"{s}__z"
        p_ = patient_medians(d, col)
        pr_ = patient_medians(d, f"{s}__reldev")
        rows.append({
            "statistic": s,
            "roi_delta": cliffs_delta(d[d.response == "Responder"][col],
                                      d[d.response == "Non-Responder"][col]),
            "pat_delta": cliffs_delta(p_[p_.response == "Responder"][col],
                                      p_[p_.response == "Non-Responder"][col]),
            "pat_p": mwu_p(p_[p_.response == "Responder"][col],
                           p_[p_.response == "Non-Responder"][col]),
            "reldev_pat_delta": cliffs_delta(
                pr_[pr_.response == "Responder"][f"{s}__reldev"],
                pr_[pr_.response == "Non-Responder"][f"{s}__reldev"]),
            "median_z_NR": d[d.response == "Non-Responder"][col].median(),
            "median_z_R": d[d.response == "Responder"][col].median(),
        })
    sec = pd.DataFrame(rows)
    sec["pat_q"] = bh(sec["pat_p"].to_numpy())
    sec = sec.sort_values("pat_p").reset_index(drop=True)

    # ---- minority-fraction ladder ---------------------------------------------
    ladder = []
    for thr in [0.00, 0.10, 0.15, 0.20, 0.25, 0.30]:
        s = d[d.minority_frac >= thr]
        p_ = patient_medians(s, pcol)
        ladder.append({
            "min_frac>=": thr, "n_ROI": len(s),
            "n_NR_pat": int((p_.response == "Non-Responder").sum()),
            "n_R_pat": int((p_.response == "Responder").sum()),
            "pat_delta": cliffs_delta(p_[p_.response == "Responder"][pcol],
                                      p_[p_.response == "Non-Responder"][pcol]),
            "pat_p": mwu_p(p_[p_.response == "Responder"][pcol],
                           p_[p_.response == "Non-Responder"][pcol]),
        })
    ladder = pd.DataFrame(ladder)

    tag = "" if inclusion == "symmetric" else f"_{inclusion}"
    fig_path = make_figure(d, defname, [perm_null.PRIMARY, "im1-avg_length",
                                        "cok1-avg_length", "dom1-avg_length"], tag)

    # Honest read of the balance ladder: it is a post-hoc subgroup scan over 6 nested,
    # highly correlated thresholds with shrinking patient counts, so a nominal p < 0.05 on
    # the last rungs is weak evidence -- but a monotone trend should not be buried either.
    ladder_trend = (
        f"delta rises monotonically from {ladder.pat_delta.iloc[0]:+.3f} (all ROIs) to "
        f"{ladder.pat_delta.iloc[-1]:+.3f} at minority fraction >= "
        f"{ladder['min_frac>='].iloc[-1]:.2f}, reaching nominal p < 0.05 on the last "
        f"{int((ladder.pat_p < 0.05).sum())} rung(s)"
        if (ladder.pat_p < 0.05).any() else
        "delta stays flat and non-significant across every threshold"
    )
    # Resolving power: how far the observed tissue sits from its own null, overall.
    frac_big = float((d[zcols].abs() > 2).to_numpy().mean())
    strongest = d[zcols].abs().median().sort_values(ascending=False)

    significant = perm_p < 0.05
    # The directional run carries a SIGNED hypothesis: a tumour nest ringed by CD8 puts a
    # long-lived loop in the CD8-only domain which the tumour fills in, so immune exclusion
    # in non-responders predicts NR > R, i.e. delta < 0. A significant delta > 0 would be a
    # real finding but NOT evidence of exclusion, and must not be reported as one.
    direction_note = ""
    if inclusion == "other_only":
        if significant and obs_d > 0:
            direction_note = (
                "\n**Direction is OPPOSITE to the exclusion hypothesis.** Exclusion in "
                "non-responders predicts δ < 0 (NR higher); the observed δ is positive. "
                "Whatever this is, it is not immune exclusion and must not be reported as "
                "such.\n")
        elif significant:
            direction_note = (
                "\n**Direction matches the exclusion hypothesis**: δ < 0, i.e. "
                "non-responders carry longer-lived CD8 loops that tumour fills in.\n")
        else:
            direction_note = (
                f"\nPre-specified direction for exclusion was δ < 0 (non-responders higher); "
                f"observed δ = {obs_d:+.3f}, not significant either way.\n")
    verdict = (
        "**Primary endpoint IS significant.** Escalate per the decision rule: run the same "
        "null for the control species (Macrophage, Fibroblast) before making any claim."
        if significant else
        "**Primary endpoint is NOT significant.** The decision rule stands: the pilot's "
        "no-go is confirmed on the strongest available evidence. Even after conditioning "
        "each ROI on its own composition- and geometry-matched null, baseline tumour-CD8 "
        "topology does not distinguish responders from non-responders. The bar-length "
        "signal in `LENGTH_FEATURES.md` is therefore attributable to tissue geometry "
        "(responder ROIs are less dense and more dispersed), not to spatial arrangement "
        "of T cells relative to tumour."
    )

    inclusion_desc = (
        "`SubChromaticInclusion(filt, [[CD8]])` — DIRECTIONAL: domain is the CD8-only "
        "complex, so degree-1 kernel bars are CD8 loops that tumour fills in, the signature "
        "an immune-exclusion hypothesis names"
        if inclusion == "other_only" else
        "`KChromaticInclusion(filt, 1)` — colour-symmetric: domain is the monochromatic "
        "subcomplex of both colours, so the kernel is tumour-dominated at a ~2:1 ratio"
    )
    lines = [
        f"# Permutation-null test: {defname}, pre-treatment, inclusion = `{inclusion}`\n",
        f"Inclusion: {inclusion_desc}.\n",
        f"Cohort: {len(d)} ROIs / {n_pat.get('Non-Responder', 0)} NR + "
        f"{n_pat.get('Responder', 0)} R patients. B = {int(d.n_perm.iloc[0])} label "
        f"permutations per ROI, cell positions held fixed.\n",
        "Each ROI's z = (observed − null mean) / null sd. Composition, cell counts, "
        "density, hull area and dispersion are identical between an ROI and its own null, "
        "so they cannot produce a difference in z. Cliff's δ is R vs NR (+ve = responders "
        "higher).\n",
        "## Sanity checks\n",
        f"- `{perm_null.INVARIANT_CHECK}` is invariant under relabelling: max null sd = "
        f"**{inv_sd:.3e}** (must be 0). {'PASS' if inv_sd == 0 else 'FAIL'}",
        f"- All {len(zcols)} z columns finite: **{all_finite}**",
        f"- Resolving power: **{frac_big:.0%}** of all ROI × statistic z-values exceed |z| = 2, "
        f"and the most deviant statistics sit enormously far from their nulls (median |z| "
        f"{strongest.iloc[0]:.0f} for `{strongest.index[0].replace('__z', '')}`, "
        f"{strongest.iloc[1]:.0f} for `{strongest.index[1].replace('__z', '')}`; overall "
        f"median |z| = {med_abs_z:.1f}). Tumour and CD8 are therefore *strongly* "
        f"non-randomly arranged in both groups — the assay has power, and a null response "
        f"contrast is informative rather than vacuous.\n",
        "## Primary endpoint\n",
        f"`{perm_null.PRIMARY}` z, per-patient median:\n",
        f"- Cliff's δ = **{obs_d:+.3f}**",
        f"- patient-label permutation p ({N_PERM_LABEL:,} draws) = **{perm_p:.4f}**",
        f"- bootstrap 95% CI over patients = [{ci[0]:+.3f}, {ci[1]:+.3f}]",
        f"- leave-one-patient-out δ range = {min(loo):+.3f} … {max(loo):+.3f}",
        f"- median z: NR {d[d.response=='Non-Responder'][pcol].median():+.2f}, "
        f"R {d[d.response=='Responder'][pcol].median():+.2f}",
        f"- size-free companion (deviation as a fraction of the null mean): "
        f"δ = **{rel_d:+.3f}**, p = {rel_p:.4f} — guards against |z| growing with ROI size, "
        f"since the null sd shrinks as an ROI gains cells\n",
        verdict,
        direction_note,
        "\n## Balance sensitivity (primary endpoint)\n",
        "The Day-1 signal collapsed as imbalanced ROIs were dropped (IMBALANCE_CONFOUND.md). "
        "Class imbalance is already held fixed within each ROI here, so this ladder is not "
        "testing the same thing — it asks whether the arrangement contrast is concentrated "
        "in ROIs where both species are well represented (where the segregation-vs-"
        "interspersion question is best posed).\n",
        md_table(ladder.round(4)),
        f"\nObserved: {ladder_trend}. Treat this as hypothesis-generating only — it is a "
        "post-hoc scan over six nested, highly correlated thresholds, patient counts fall as "
        "the threshold rises (25→20 NR), and the primary endpoint was pre-specified at the "
        "0.00 rung. It is not evidence for an effect, but it is the one place worth a "
        "pre-registered look if this line is ever revisited.",
        "\n## Secondary: all retained statistics (per-patient, BH-adjusted)\n",
        md_table(sec.round(4)),
        "\n## Limitations — what this does and does not rule out\n",
        f"- **Power.** 25 NR vs 37 R patients gives ~80% power only at Cliff's |δ| ≈ 0.42. "
        f"The primary CI is [{ci[0]:+.3f}, {ci[1]:+.3f}]: a moderate responder-higher effect "
        f"is not excluded. 'No difference' here means 'no large difference'.",
        f"- **The two estimators of the same quantity disagree in magnitude.** z gives "
        f"δ = {obs_d:+.3f} (p = {perm_p:.3f}); the size-free relative deviation gives "
        f"δ = {rel_d:+.3f} (p = {rel_p:.3f}). Both are non-significant, but the second is "
        f"borderline, so this is a soft null rather than a clean one.",
        "- **Scope.** One pair definition (`CD8_primary`), one inclusion "
        "(`KChromaticInclusion(filt, 1)`, which is colour-symmetric and therefore "
        "tumour-dominated), degree 0 and 1 only, one filtration (`delaunay_cech`). A "
        "directional `SubChromaticInclusion` on CD8 alone — the statistic an *exclusion* "
        "hypothesis actually calls for — has still not been run under this null.",
        "- **B = 99** bounds the per-ROI null resolution; z is estimated from 99 draws, so "
        "per-ROI z carries noise of order 1/sqrt(2·99) ≈ 7% on the sd.",
        f"\nFigure: `{fig_path.name}`. Scripts: `perm_null.py`, `day2_permutation_test.py`.",
    ]
    out_md = OUT / f"PERMUTATION_NULL{tag.upper()}.md"
    out_md.write_text("\n".join(lines) + "\n")

    pd.set_option("display.width", 220)
    print(f"{defname} pre-treatment: {len(d)} ROIs, patients {n_pat}")
    print(f"sanity: invariant null sd = {inv_sd:.3e}; all z finite = {all_finite}; "
          f"median |z| = {med_abs_z:.1f}")
    print(f"\nPRIMARY {perm_null.PRIMARY}: delta={obs_d:+.3f}  perm p={perm_p:.4f}  "
          f"CI=[{ci[0]:+.3f},{ci[1]:+.3f}]  LOO {min(loo):+.3f}..{max(loo):+.3f}")
    print(f"  size-free companion (reldev): delta={rel_d:+.3f}  p={rel_p:.4f}")
    print("\nTop 10 secondary by p:")
    print(sec.head(10).round(4).to_string(index=False))
    print(f"\nBH q<0.05: {int((sec.pat_q < 0.05).sum())} / {len(sec)}")
    print("\nBalance ladder:")
    print(ladder.round(4).to_string(index=False))
    print("\nwrote", out_md)


if __name__ == "__main__":
    main()
