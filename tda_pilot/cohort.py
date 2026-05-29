"""Cohort definition for the TDA pilot (revision cohort, N = 808).

Cohort sizes (see REPORT.md):
  - 829 ROIs : raw cell table (CellTable_CleanCohort/updated_cell_table.csv)
  - 813 ROIs / 63 patients : originally submitted manuscript cohort
      = 829 minus patients LEAP149 (4 ROIs) and LEAP150 (12 ROIs)
  - 808 ROIs / 63 patients : revision cohort (this pilot)
      = 813 minus five collaborator-flagged ROIs (duplicates / broken / palette)
"""

EXCLUDED_PATIENTS = {"LEAP149", "LEAP150"}
EXCLUDED_ROIS = {
    "Leap005_2_1",  # duplicate of Leap005_1
    "Leap005_2_2",  # duplicate of Leap005_2
    "Leap010_7",    # broken; retaken as Leap010_8
    "Leap094_7",    # broken; retaken as Leap094_10
    "Leap095_13",   # palette region, not tissue
}


def is_cohort_member(roi_id: str) -> bool:
    """True if roi_id belongs to the revision cohort (N = 808)."""
    if roi_id in EXCLUDED_ROIS:
        return False
    patient = roi_id.split("_")[0].upper()  # "Leap066_11" -> "LEAP066"
    return patient not in EXCLUDED_PATIENTS


def predictive_cohort_rois(per_roi_counts_df):
    """Return per_roi_counts filtered to the 808-ROI revision cohort.

    Accepts the ROI identifier in either a ``roi_id`` or ``fov`` column
    (this project's per_roi_counts.csv uses ``fov``).
    """
    id_col = "roi_id" if "roi_id" in per_roi_counts_df.columns else "fov"
    return per_roi_counts_df[per_roi_counts_df[id_col].apply(is_cohort_member)]
