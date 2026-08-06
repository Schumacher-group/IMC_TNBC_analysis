# Pre-check: are B7H4+ and B7H4− cancer cells spatially separated?

557 pre-treatment ROIs with at least 30 cells of each subtype (of 593 eligible). k = 10 nearest cancer neighbours, expectation from 20 shuffles of B7H4 status among the same cells.

- Median mixing ratio (observed / shuffled) = **0.847**
- Interquartile range: 0.788 – 0.892
- ROIs with ratio < 0.9 (clearly segregated): **79%**
- Median B7H4+ share of cancer cells: 31.1%

**Verdict: **spatially segregated** — the subtypes occupy distinguishable territory, so a paired CD8 comparison is meaningful.**

Why this matters: it is the design check for the proposed 'is CD8 more excluded from B7H4+ tumour?' comparison. A ratio near 1 would mean any CD8 cell near one subtype is near the other by construction, so both statistics would describe the same tissue and a null result would be uninformative.

Method note: this is nearest-neighbour counting only, no persistent homology, so it is independent of every modelling choice elsewhere in this directory.

Figure: `b7h4_precheck.png`. Script: `day2_b7h4_precheck.py`.
