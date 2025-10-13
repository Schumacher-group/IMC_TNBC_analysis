# IMC Analysis Pipeline

This repository provides a comprehensive workflow for analyzing and visualizing Imaging Mass Cytometry (IMC) data from Triple-Negative Breast Cancer (TNBC) samples. The pipeline covers segmentation, phenotyping, and spatial analysis of IMC data.

## Citation
This pipeline was used in the following preprint:

**"Identifying tissue states by spatial protein patterns related to chemotherapy response in triple-negative breast cancer"**
bioRxiv (2025). DOI: [10.1101/2025.10.06.680783](https://www.biorxiv.org/content/10.1101/2025.10.06.680783)

If you use this pipeline in your work, please consider citing the preprint. 

## Downstream analysis
### Data Preprocessing and Segmentation
**Note**: This repository contains **analysis code only**. Data preprocessing and segmentation are performed using the separate [IMC_preprocessing](https://github.com/Schumacher-group/IMC_preprocessing) pipeline.

The preprocessing pipeline handles:
- MCD extraction to OME-TIFF format
- Image denoising (IMC Denoise)
- CLAHE contrast enhancement for image-level batch effect removal
- Cell segmentation using DeepCell's Mesmer model
- Cell table generation with marker quantification

The output of the preprocessing pipeline (processed images and cell tables) serves as input for the phenotyping and analysis in this repository.

### Using the output of this analysis for prediction of therapy response:
Prediction of therapy response based on the analysed data is handled in: https://github.com/Schumacher-group/ML4SpatialAnalysis

## Repository Structure

- **phenotyping/**
  Notebooks for cellular phenotyping analysis:
  - Pixel and cell-level clustering (Pixie workflow)
  - Cell type identification based on marker expression
  - Quality control and batch effect assessment
  - Dimensionality reduction (UMAP) to visualise batch effects (at a coarse level)

- **figures/manuscript/Notebooks/**
  Notebooks for generating publication figures:
  - **Cell Counts.ipynb**: Cell type composition analysis with statistical testing
  - **Neighbourhood_Enrichment.ipynb**: Spatial neighborhood enrichment analysis
  - **KmeansNeighbourhoodsAnalysis.ipynb**: K-means neighborhood clustering
  - **Collagen_plots.ipynb**: Collagen fiber architecture analysis
  - **nbr_enrich_violins.ipynb**: Neighborhood enrichment visualizations

- **segmentation_QC/**
  Quality control notebooks for cell segmentation validation:
  - Post-segmentation analysis and quality control
  - Segmentation visualization tools

- **Please note repo does not include the cell tables used in the analyses here, as this is a large file (>8GB). Users can generate their own cell table or contact the authors of the publication above to get a copy of the processed cell table file.**

- **Environment**
  - **environment.yml**: Conda environment file with all required dependencies

## Installation

Clone the repository and create the Conda environment using the provided `environment.yml`:

```bash
git clone <repo_url>
cd <repo_folder>
conda env create -f environment.yml
conda activate imc_env
```

## Usage

The analysis pipeline consists of two main stages:

1. **Phenotyping**
   - Run notebooks in `phenotyping/pixie/` to identify cell types based on marker expression
   - Key notebooks:
     - `2_Pixie_Cluster_Pixels.ipynb`: Pixel-level clustering
     - `3_Pixie_Cluster_Cells.ipynb`: Cell phenotyping and clustering
     - `assess_clustering.ipynb`: Quality control for clustering
     - `assessing_batch effects.ipynb`: Batch effect visualization

2. **Analysis and Figure Generation**
   - Use notebooks in `figures/manuscript/Notebooks/` for all publication analyses:
     - Cell type composition analysis
     - Spatial neighborhood enrichment
     - K-means neighborhood clustering
     - Collagen fiber architecture quantification

## Contributing

Contributions are welcome! Please fork the repository and open a pull request with your proposed changes. Make sure to test your changes and update the documentation as necessary.
