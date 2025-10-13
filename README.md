# IMC Analysis Pipeline

This repository provides a comprehensive workflow for analyzing and visualizing Imaging Mass Cytometry (IMC) data from Triple-Negative Breast Cancer (TNBC) samples. The pipeline covers segmentation, phenotyping, and spatial analysis of IMC data.

## Citation

If you use this analysis pipeline, please cite:

**Identifying tissue states by spatial protein patterns related to chemotherapy response in triple-negative breast cancer**
bioRxiv (2025). DOI: [10.1101/2025.10.06.680783](https://doi.org/10.1101/2025.10.06.680783)

## Repository Structure

- **segmentation/**
  Scripts and notebooks for cell segmentation:
  - **Mesmer_IMC.ipynb**: Cell segmentation using DeepCell's Mesmer model
  - Post-segmentation analysis and quality control notebooks
  - Segmentation visualization tools

- **phenotyping/**
  Notebooks for cellular phenotyping analysis:
  - Pixel and cell-level clustering (Pixie workflow)
  - Cell type identification based on marker expression
  - Quality control and batch effect assessment
  - Dimensionality reduction (UMAP) and clustering

- **results_analysis/**
  In-depth spatial and non-spatial analysis:
  - Non-spatial analysis: cell type composition, enrichment analysis, response prediction
  - Spatial analysis: neighborhood enrichment, cell-cell contacts, distance distributions
  - Collagen and fiber architecture analysis

- **figures/manuscript/**
  Notebooks for generating publication figures

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

If needed, install any additional pip dependencies:

```bash
pip install -r requirements.txt
```

## Data Preprocessing

**Note**: This repository contains **analysis code only**. Raw IMC data preprocessing (MCD file extraction, image denoising, normalization, and artifact correction) is performed using the separate [IMC_preprocessing](https://github.com/Schumacher-group/IMC_preprocessing) pipeline.

The preprocessing pipeline handles:
- MCD extraction to OME-TIFF format
- Image denoising (IMC Denoise)
- Multi-stage normalization (per-staining-batch and per-patient)
- Special marker handling (e.g., Carboplatin)
- CLAHE contrast enhancement

The output of the preprocessing pipeline serves as input for this analysis repository.

## Usage

The analysis pipeline consists of three main stages:

1. **Segmentation**
   - Launch `segmentation/Mesmer_IMC.ipynb` to perform cell segmentation using DeepCell's Mesmer model
   - Use post-segmentation notebooks for quality control and visualization

2. **Phenotyping**
   - Run notebooks in `phenotyping/pixie/` to identify cell types based on marker expression
   - Perform clustering, batch effect assessment, and quality control

3. **Spatial and Statistical Analysis**
   - Explore `results_analysis/` notebooks for:
     - Cell type composition and enrichment analysis
     - Spatial neighborhood enrichment
     - Cell-cell contact patterns
     - Response prediction models
   - Generate publication figures using notebooks in `figures/manuscript/`

Each notebook contains detailed documentation.

## Contributing

Contributions are welcome! Please fork the repository and open a pull request with your proposed changes. Make sure to test your changes and update the documentation as necessary.

## License

This project is licensed under the terms specified in the LICENSE file.
