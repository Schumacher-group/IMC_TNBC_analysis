"""Three-species definitions for the chromatic TDA triples.

The whole pilot so far has used a multi-species tool in two-species mode: every run called
`KChromaticInclusion(filt, 1)` on exactly two labels. M2S2 defaults to `--max-num-labels 3`
and uses `KChromaticQuotient(filt, k)` with k = n_labels - 1 for combinations of more than two,
so the triple machinery is theirs, unused by us until now.

Why triples matter here: the hypothesis under review -- "collagen prevents cytotoxic T cells
from reaching the tumour" -- is irreducibly a three-body statement. No pairwise statistic can
express "A is separated from B *by* C". Everything we have tested so far can only say whether
A and B are mixed.

SPECIES ABUNDANCE governs what is feasible. From a full scan of the cell table (4.24M cells,
19 labels), as a fraction of all cells:
    CD8 T cell 18.9%, Cancer cell 18.2%, Fibroblast 14.0%, Macrophage 7.6%,
    B7H4 Cancer cell 7.5%, Endothelial cell 7.4%, Memory CD4 7.2%, APC 6.0%,
    B cell 4.7%, NK/CD8 3.2%, Memory CD8 1.5%
    ... then a long tail: Regulatory T cell 0.42%, Macrophage M2 0.23%.
Regulatory T cells (~21 cells/ROI) and M2 macrophages (~12 cells/ROI) are TOO SPARSE to carry
per-ROI topology -- they would sit near `min_count` and their diagrams would be noise -- so
biologically attractive triples involving them are not viable in this cohort.
"""

TUMOUR = {"Cancer cell", "B7H4 Cancer cell"}
CD8 = {"CD8 T cell", "Memory CD8 T cell"}

# Sentinel for a species defined by MARKER INTENSITY rather than by cell-type label.
# Resolved at runtime in run_triples.py: cells not already claimed by an earlier species in
# the list, whose Collagen-Type-I is at or above the per-ROI COLLAGEN_QUANTILE. The
# per-ROI quantile (not a global threshold) absorbs staining and batch variation.
#
# Why this exists: fibroblasts DEPOSIT collagen but are not collagen, and dense acellular
# stroma is collagen-rich and fibroblast-poor -- exactly the tissue a physical barrier would
# consist of. Measured on a 60-ROI sample, the top collagen quartile is only 23.9%
# fibroblast (against 12.2% fibroblast in tissue overall), and only 48.9% of fibroblasts are
# collagen-high, so the two proxies are far from redundant.
#
# SHARED LIMITATION, worth stating in any writeup: both proxies are cell-based, so neither
# can mark truly ACELLULAR collagen. If the barrier is dense matrix containing no cells,
# no cell-based species can represent it. Only genuine fibre/pixel points could, and fibre
# segmentation did not yield usable ones.
COLLAGEN_HIGH = "__COLLAGEN_HIGH__"
# Real collagen, from the finer fibre segmentation: skeleton points thinned to a fixed
# spacing. Unlike the two cell-based proxies this is not blind to acellular matrix, so it
# is the first representation that can actually test the barrier claim.
COLLAGEN_FIBRE = "__COLLAGEN_FIBRE__"
FIBRE_SPACING = 20.0
# Which segmentation variant to draw fibre points from. non-CLAHE is primary: it shows no
# stain-batch dependence (rho=+0.02 vs -0.11) and does not inflate coverage ~60% by
# amplifying faint signal in collagen-poor regions -- i.e. tumour nest interiors, which is
# exactly where a barrier question is decided.
FIBRE_DIR = "fibre_skeletons_nonCLAHE"
COLLAGEN_COLUMN = "Collage-Type_I"
COLLAGEN_QUANTILE = 0.75

# name -> ordered list of (label, species-set). Colour index follows list order.
TRIPLES: dict[str, list[tuple[str, set[str]]]] = {
    # 1. The barrier hypothesis. Fibroblast stands in for collagen: fibre segmentation did not
    #    yield usable points, and fibroblasts are the cells that deposit collagen. Caveat worth
    #    stating in any writeup -- dense acellular stroma is collagen-rich but fibroblast-poor,
    #    so this proxy misses exactly the regions a physical barrier would occupy.
    "Tumour_CD8_Fibroblast": [
        ("Tumour", TUMOUR), ("CD8", CD8), ("Fibroblast", {"Fibroblast"})],

    # 2. Vascular access. T cells enter tissue through vessels, so apparent "exclusion" can be
    #    a delivery problem rather than a barrier. A genuinely different mechanism.
    "Tumour_CD8_Endothelial": [
        ("Tumour", TUMOUR), ("CD8", CD8), ("Endothelial", {"Endothelial cell"})],

    # 3. Myeloid exclusion. M2 macrophages would be the specific target but are far too sparse
    #    (0.23%), so this uses the whole macrophage compartment.
    "Tumour_CD8_Macrophage": [
        ("Tumour", TUMOUR), ("CD8", CD8),
        ("Macrophage", {"Macrophage", "Macrophage + HLA", "Macrophage M2"})],

    # 4. Checkpoint-specific exclusion. Every analysis so far AGGREGATED 'Cancer cell' and
    #    'B7H4 Cancer cell' into one tumour colour, discarding the distinction that plausibly
    #    matters most for T-cell exclusion -- B7-H4 is an immunosuppressive checkpoint ligand,
    #    and the manuscript's own figure shows it as a separate channel. This triple asks
    #    whether CD8 cells avoid B7H4+ tumour specifically rather than tumour in general.
    "Cancer_B7H4_CD8": [
        ("Cancer", {"Cancer cell"}), ("B7H4_Cancer", {"B7H4 Cancer cell"}), ("CD8", CD8)],

    # 5. Antigen presentation. Where APCs sit relative to tumour and CD8 speaks to priming.
    "Tumour_CD8_APC": [
        ("Tumour", TUMOUR), ("CD8", CD8), ("APC", {"Antigen presenting cell"})],

    # 6. The barrier hypothesis again, with a marker-based collagen proxy rather than a
    #    cell-type one. Tumour and CD8 are claimed first, so the collagen species is
    #    "any OTHER cell sitting in collagen-rich tissue" -- which is closer to the
    #    barrier concept than fibroblast identity, and measurably different from it.
    "Tumour_CD8_CollagenHigh": [
        ("Tumour", TUMOUR), ("CD8", CD8), ("CollagenHigh", {COLLAGEN_HIGH})],

    # 7. The barrier hypothesis with REAL collagen rather than a cell proxy.
    "Tumour_CD8_CollagenFibre": [
        ("Tumour", TUMOUR), ("CD8", CD8), ("CollagenFibre", {COLLAGEN_FIBRE})],
}

# Cell-type labels only; the marker sentinel is not a label and must not reach the data build.
ALL_LABELS: set[str] = set()
for _spec in TRIPLES.values():
    for _name, _labels in _spec:
        ALL_LABELS |= {lab for lab in _labels if lab not in (COLLAGEN_HIGH, COLLAGEN_FIBRE)}
