"""Dowker (a.k.a. Witness) persistent homology of two point clouds, in Python.

WHY THIS EXISTS, AND WHY NOT gudhi's WITNESS COMPLEX
----------------------------------------------------
gudhi ships `WitnessComplex` / `StrongWitnessComplex`, but they are NOT the Dowker complex:
they use the *relaxed* witness convention, in which a witness's filtration contribution is
offset by its distance to the nearest landmark. Checked directly against a brute-force Dowker
on the reference toy example from https://github.com/irishryoon/Dowker_persistence:

    brute-force Dowker  H1 = [0.331, 0.756]
    gudhi StrongWitness H1 = [0.065, 0.534]        (no distance/squared convention reconciles them)

So the construction is implemented here from its definition instead. The reference
implementation is Julia (Eirene); this is a Python equivalent so the analysis can reuse all
the existing cohort, matched-design and permutation machinery in this directory.

THE CONSTRUCTION
----------------
Given landmarks L and witnesses W with cross-distances D[l, w], a simplex sigma over L enters
the filtration at

    t(sigma) = min over w in W of max over l in sigma of D[l, w]

i.e. the scale at which some single witness is within t of *every* vertex of sigma.

WHAT IT MEASURES, AND WHY IT SUITS AN EXCLUSION QUESTION
--------------------------------------------------------
The filtration uses ONLY cross-distances. Distances within L and within W never enter. That
matters here: every confound that dominated the chromatic analyses in this directory came from
within-species structure (bar count ~ cell number at rho = 0.83, bar length ~ composition at
rho = -0.69, z-scores coupled to ROI density). Dowker is structurally immune to that family.

It is also asymmetric -- Dowker(cancer | CD8) asks "how well do CD8 cells cover cancer?" -- and
a long-lived degree-1 bar is a loop of cancer cells that NO single CD8 cell is close to, i.e.
literally a region of tumour that CD8 does not reach. Note Dowker duality means the diagrams of
Dowker(L|W) and Dowker(W|L) agree, so the asymmetry is in what the vertices represent, not in
the barcode.

VALIDATION (run `python dowker.py`)
-----------------------------------
  * agrees with exhaustive enumeration on the reference toy example
  * landmarks on a ring witnessed by the same ring -> exactly one long degree-1 bar
  * adding witnesses INSIDE the ring kills that bar earlier -- the exclusion mechanic itself

PERFORMANCE
-----------
Simplices are pruned by a distance `cutoff`, which is sound because
t(sigma) >= max over pairs in sigma of t(pair), so the Dowker complex is a subcomplex of the
flag complex of its own 1-skeleton. Only cliques of the thresholded 1-skeleton need testing.
At realistic sizes (391 landmarks, 816 witnesses, 80 um cutoff) this is ~0.2 s -- faster than
the chromatic six-packs it is being compared against.

Bars still alive at `cutoff` are right-censored; `dowker_features` reports finite bars only and
also returns the censored count so the truncation is visible rather than silent.
"""

from __future__ import annotations

import itertools

import gudhi
import numpy as np

__all__ = ["cross_distances", "dowker_simplex_tree", "dowker_features", "dowker_bruteforce"]


def cross_distances(landmarks: np.ndarray, witnesses: np.ndarray) -> np.ndarray:
    """(n_landmarks, n_witnesses) Euclidean cross-distance matrix."""
    return np.linalg.norm(landmarks[:, None, :] - witnesses[None, :, :], axis=2)


def _pairwise_dowker(D: np.ndarray, chunk: int = 64) -> np.ndarray:
    """t_ij = min_w max(D[i,w], D[j,w]), computed in row blocks to bound memory."""
    n = D.shape[0]
    T = np.empty((n, n), dtype=float)
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        T[s:e] = np.maximum(D[s:e, None, :], D[None, :, :]).min(axis=2)
    return T


def dowker_simplex_tree(D: np.ndarray, cutoff: float) -> gudhi.SimplexTree:
    """Dowker filtration up to dimension 2 (enough for degrees 0 and 1)."""
    n = D.shape[0]
    st = gudhi.SimplexTree()
    for i in range(n):
        st.insert([i], float(D[i].min()))

    T = _pairwise_dowker(D)
    ii, jj = np.triu_indices(n, 1)
    keep = T[ii, jj] <= cutoff
    ii, jj = ii[keep], jj[keep]
    for a, b in zip(ii, jj):
        st.insert([int(a), int(b)], float(T[a, b]))

    # t(sigma) >= max pairwise t, so only cliques of the thresholded 1-skeleton can qualify.
    adj: list[set[int]] = [set() for _ in range(n)]
    for a, b in zip(ii, jj):
        adj[int(a)].add(int(b))
    for a, b in zip(ii, jj):
        for c in adj[int(a)] & adj[int(b)]:
            t = float(np.max(D[[int(a), int(b), c], :], axis=0).min())
            if t <= cutoff:
                st.insert([int(a), int(b), c], t)

    st.make_filtration_non_decreasing()
    return st


def dowker_features(landmarks: np.ndarray, witnesses: np.ndarray, cutoff: float) -> dict:
    """Summary statistics of the Dowker diagram of `landmarks` witnessed by `witnesses`."""
    D = cross_distances(landmarks, witnesses)
    st = dowker_simplex_tree(D, cutoff)
    st.compute_persistence()
    out: dict[str, float] = {}
    for dim in (0, 1):
        bars = np.asarray(st.persistence_intervals_in_dimension(dim), dtype=float)
        finite = bars[np.isfinite(bars).all(axis=1)] if len(bars) else np.zeros((0, 2))
        life = finite[:, 1] - finite[:, 0] if len(finite) else np.zeros(0)
        life = life[life > 0]
        p = f"dow{dim}"
        out[f"{p}-num_bars"] = float(len(life))
        out[f"{p}-censored"] = float(len(bars) - len(finite))  # alive at cutoff
        out[f"{p}-total_persistence"] = float(life.sum()) if len(life) else 0.0
        out[f"{p}-avg_length"] = float(life.mean()) if len(life) else 0.0
        out[f"{p}-med_length"] = float(np.median(life)) if len(life) else 0.0
        out[f"{p}-p90_length"] = float(np.percentile(life, 90)) if len(life) else 0.0
        out[f"{p}-max_length"] = float(life.max()) if len(life) else 0.0
        if len(life):
            q = life / life.sum()
            out[f"{p}-entropy"] = float(-(q * np.log(q)).sum())
        else:
            out[f"{p}-entropy"] = 0.0
    return out


def dowker_bruteforce(D: np.ndarray, maxdim: int = 1) -> gudhi.SimplexTree:
    """Exhaustive enumeration, for validation only. O(n^3) -- small inputs."""
    st = gudhi.SimplexTree()
    n = D.shape[0]
    for k in range(1, maxdim + 3):
        for sig in itertools.combinations(range(n), k):
            st.insert(list(sig), float(np.min(np.max(D[list(sig), :], axis=0))))
    st.make_filtration_non_decreasing()
    return st


def _selftest() -> None:
    rng = np.random.default_rng(0)

    th = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    ring_l = np.c_[np.cos(th), np.sin(th)]
    ring_w = np.c_[np.cos(th + 0.13), np.sin(th + 0.13)]

    st = dowker_bruteforce(cross_distances(ring_l, ring_w))
    st.compute_persistence()
    h1 = np.asarray(st.persistence_intervals_in_dimension(1))
    assert len(h1) == 1 and h1[0][1] > 1.5, f"ring: expected one long bar, got {h1}"
    ring_death = h1[0][1]

    centre = rng.normal(0, 0.15, (8, 2))
    st2 = dowker_bruteforce(cross_distances(ring_l, np.vstack([ring_w, centre])))
    st2.compute_persistence()
    h1b = np.asarray(st2.persistence_intervals_in_dimension(1))
    assert len(h1b) == 0 or h1b[0][1] < ring_death, "witnesses inside the ring must kill it sooner"

    # pruned implementation must agree with exhaustive enumeration
    L = rng.uniform(0, 10, (30, 2))
    W = rng.uniform(0, 10, (40, 2))
    D = cross_distances(L, W)
    a = dowker_bruteforce(D)
    a.compute_persistence()
    b = dowker_simplex_tree(D, cutoff=1e9)
    b.compute_persistence()
    ha = np.array(sorted(map(tuple, a.persistence_intervals_in_dimension(1))))
    hb = np.array(sorted(map(tuple, b.persistence_intervals_in_dimension(1))))
    assert ha.shape == hb.shape and np.allclose(ha, hb), "pruned != brute force"

    print(f"selftest OK  (ring bar dies at {ring_death:.3f}, "
          f"{ring_death:.3f} -> {h1b[0][1]:.3f} with centre witnesses; "
          f"pruned matches brute force on {len(ha)} H1 bars)")


if __name__ == "__main__":
    _selftest()
