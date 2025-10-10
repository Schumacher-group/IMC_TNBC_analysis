# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains a comprehensive pipeline for analyzing Imaging Mass Cytometry (IMC) data from Triple-Negative Breast Cancer (TNBC) samples. The analysis evaluates patient responses to neoadjuvant chemotherapy by examining cellular phenotypes, spatial relationships, and tissue architecture.

**Related Repository**: This repo contains the **analysis code** for the TNBC IMC study. Data preprocessing (MCD extraction, denoising, normalization) was performed using the [IMC_preprocessing](https://github.com/Schumacher-group/IMC_preprocessing) pipeline.

## Environment Setup

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate imc_analysis
```

The environment includes key dependencies: scanpy, squidpy, scikit-image, xgboost, wandb, and various data science libraries.

## Repository Architecture

### Data Processing Overview

**Preprocessing**: MCD extraction, denoising, and normalization were performed using the [IMC_preprocessing](https://github.com/Schumacher-group/IMC_preprocessing) pipeline, which includes:
- MCD file extraction to OME-TIFF format (steinbock)
- Image denoising (IMC Denoise)
- Multi-stage normalization (per-staining-batch and per-patient)
- Special Carboplatin marker handling (CORE vs RESECTION samples)
- CLAHE contrast enhancement

**Analysis Pipeline** (this repository):
1. **Segmentation**: Cell segmentation using DeepCell's Mesmer model
2. **Phenotyping**: Cell type identification based on marker expression
3. **Spatial Analysis**: Neighborhood enrichment, cell-cell contacts, fiber quantification
4. **Statistical Analysis**: Response prediction and comparative analysis

### Analysis Modules

#### Segmentation (`segmentation/`)

Key notebooks:
- `Mesmer_IMC.ipynb`: Cell segmentation using DeepCell's Mesmer model
- `post_segmentation_analysis.ipynb`: Quality control and segmentation validation
- `visualise_segmentation.ipynb`: Visualization of segmentation masks overlaid on images

The segmentation pipeline produces cell masks that are used in all downstream spatial analyses.

#### Phenotyping (`phenotyping/`)

Located in `phenotyping/pixie/`:
- `2_Pixie_Cluster_Pixels.ipynb`: Pixel-level clustering (intermediate step)
- `3_Pixie_Cluster_Cells.ipynb`: Main cell-level clustering for cell type identification
- `assess_clustering.ipynb`: Quality control for clustering evaluation
- `assessing_batch effects.ipynb`: Batch effect assessment with UMAP visualization
- Uses marker expression to assign cell types (T cells, B cells, macrophages, fibroblasts, tumor cells, etc.)
- Clustering methods include Leiden/Louvain on marker expression space

#### Spatial Analysis (`results_analysis/`)

**Non-spatial analysis** (`results_analysis/non-spatial/`):
- `non_spatial.ipynb`: Cell type composition analysis per FOV and per patient
- `non_spatial_enR.ipynb`: Cell type enrichment analysis comparing responders vs. non-responders
- `predict_response.ipynb`: Machine learning models (XGBoost) predicting treatment response

**Spatial analysis** (`results_analysis/spatial.ipynb`):
- Neighborhood enrichment using Squidpy
- Cell-cell contact analysis
- Distance distribution calculations
- Fiber and collagen architecture quantification


### Core Utilities (`imcsegpipe/`)

A custom package for IMC data handling:
- `_imcsegpipe.py`: Core extraction functions (`extract_mcd_file()`, `create_analysis_stacks()`)
- `utils.py`: Hot pixel filtering, channel sorting, OME-XML generation
- Handles both MCD and TXT file formats for acquisition data recovery

## Key Analysis Concepts

### Data Organization

- **LEAP IDs**: Patient identifiers (e.g., LEAP001, LEAP002)
- **Acquisition IDs**: Individual field-of-view (FOV) within a patient sample
- **Sample Types**: CORE (pre-treatment biopsy) vs. RESECTION (post-treatment surgical sample)
- **Response Status**: RESPONDER vs. NON-RESPONDER based on residual cancer burden (RCB)

### Typical Analysis Workflow

1. Extract raw MCD files → OME-TIFF format
2. Preprocess images (normalization, artifact removal)
3. Segment cells using Mesmer
4. Phenotype cells based on marker expression
5. Perform spatial analysis (neighborhood enrichment, cell-cell contacts)
6. Statistical comparison between responders and non-responders
7. Generate figures for publication

### Important File Paths

When working with notebooks, be aware that paths are often relative to project root:
- Raw data: `../IMC_data/`
- Processed images: `../Img_Denoised/processed/`
- Split channels: `../split_channels_nohpf/`
- Metadata: `../IMC_data/ExtraDocs/`

## Working with Notebooks

Most analysis is performed in Jupyter notebooks. Key notebooks by analysis stage:

**Preprocessing visualization**:
- `explorative_analysis/image_contrast_enhancement.ipynb`: Visualize preprocessing effects (raw → denoised → final)

**Segmentation**:
- `segmentation/Mesmer_IMC.ipynb`: Cell segmentation using Mesmer
- `segmentation/post_segmentation_analysis.ipynb`: Segmentation QC
- `segmentation/visualise_segmentation.ipynb`: Segmentation visualization

**Phenotyping**:
- `phenotyping/pixie/2_Pixie_Cluster_Pixels.ipynb`: Pixel clustering
- `phenotyping/pixie/3_Pixie_Cluster_Cells.ipynb`: Cell phenotyping
- `phenotyping/pixie/assess_clustering.ipynb`: Clustering QC
- `phenotyping/pixie/assessing_batch effects.ipynb`: Batch effect visualization

**Cell counts and composition**:
- `results_analysis/non-spatial/non_spatial.ipynb`: Cell type composition per FOV/patient
- `results_analysis/non-spatial/non_spatial_enR.ipynb`: Enrichment analysis

**Spatial patterns**:
- `results_analysis/spatial.ipynb`: Main spatial analysis
- `nbr_enrich_violins.ipynb`: Neighborhood enrichment plots
- `results_analysis/collagen/`: Fiber segmentation and analysis
- `results_analysis/cell_contact/`: Cell-cell contact analysis

**Prediction modeling**:
- `results_analysis/non-spatial/predict_response.ipynb`: ML models for response prediction

**Manuscript figures**:
- `figures/manuscript/Notebooks/Cell Counts.ipynb`: Cell count comparisons
- `figures/manuscript/Notebooks/Neighbourhood_Enrichment.ipynb`: Spatial enrichment analysis
- `figures/manuscript/Notebooks/KmeansNeighbourhoodsAnalysis.ipynb`: K-means neighborhoods
- `figures/manuscript/Notebooks/Collagen_plots.ipynb`: Collagen/fiber analysis
- `figures/manuscript/Notebooks/stat_tests.ipynb`: Statistical testing

## Common Data Structures

### AnnData Objects

Most analyses use scanpy's AnnData format:
- `.obs`: Cell-level metadata (patient, FOV, cell type, coordinates)
- `.var`: Marker/feature metadata
- `.X`: Expression matrix (cells × markers)
- `.obsm`: Embeddings (UMAP, PCA)
- `.uns`: Unstructured metadata (plotting parameters, statistics)

### Spatial Data with Squidpy

Spatial analyses use Squidpy's spatial data structures:
- Coordinate information in `.obsm['spatial']`
- Spatial graphs in `.obsp`
- Neighborhood enrichment results in `.uns['nhood_enrichment']`

## Data Normalization Strategy

The preprocessing pipeline applies multiple normalization steps:

1. **Per-channel hot pixel filtering** (optional, via `hpf` parameter)
2. **Per-staining-batch quantile normalization** (corrects batch effects across staining dates)
3. **Per-patient intensity normalization** (ensures comparable signal across patients)
4. **Special Carboplatin handling** (sets core samples to zero, normalizes resection samples)

This multi-stage approach accounts for both technical (staining batch) and biological (tissue type) variability.

## Testing and Validation

There are no formal unit tests in this repository. Validation is performed through:
- Visual inspection of segmentation quality
- QC plots in exploratory notebooks
- Comparison of marker distributions across batches
- Cross-validation in prediction models

## Development Notes

- The `imcsegpipe` module is a lightly modified version of external tools adapted for this project
- Many notebooks contain hard-coded paths that may need adjustment if directory structure changes
- Large image files and raw MCD files are not tracked in git (referenced via relative paths outside repo)
- Data preprocessing was performed using the [IMC_preprocessing](https://github.com/Schumacher-group/IMC_preprocessing) pipeline prior to the analysis in this repository
