"""Pair definitions for the Day-1 full TDA pilot (revision cohort, N=808).

Tumour aggregate is fixed; each definition pairs it with one non-tumour species.
Labels are exact `cell_meta_cluster` strings confirmed with the user from the
full 19-label inventory.
"""

# Tumour aggregate (both keratin+ cancer phenotypes)
TUMOUR = {"Cancer cell", "B7H4 Cancer cell"}

# pair_definition name -> set of cell_meta_cluster labels for the non-tumour species
DEFINITIONS: dict[str, set[str]] = {
    # focal pair + two CD8 sensitivity definitions
    "CD8_primary":  {"CD8 T cell", "Memory CD8 T cell"},
    "CD8_strict":   {"CD8 T cell"},
    "CD8_extended": {"CD8 T cell", "Memory CD8 T cell", "NK/CD8"},
    # NK/CD8 alone: the activated GranzymeB+ CD8 cluster the manuscript describes as
    # "regions of NK/CD8 cells infiltrating cancer cell groups" -> cleanest mapping.
    "CD8_NK_only":  {"NK/CD8"},
    # four control pairs (same tumour aggregate, different non-tumour species)
    "CD4":          {"Memory CD4 T cell"},
    "Bcell":        {"B cell", "B cell + HLA"},
    "Macrophage":   {"Macrophage", "Macrophage + HLA", "Macrophage M2"},
    "Fibroblast":   {"Fibroblast"},  # named Fibroblast (not CAF): no CAF-specific annotation
}

# Raw six-pack diagrams (.h5) are retained only for the focal definition.
SAVE_DIAGRAMS = {"CD8_primary"}

# Union of every cell_meta_cluster label we ever need (for the one-pass data build).
ALL_LABELS: set[str] = set(TUMOUR)
for _s in DEFINITIONS.values():
    ALL_LABELS |= _s
