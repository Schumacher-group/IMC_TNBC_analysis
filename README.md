# IMC Analysis Pipeline

This repository provides a comprehensive workflow for processing, analyzing, and visualizing Imaging Mass Cytometry (IMC) data. The pipeline covers everything from data extraction and image preprocessing to segmentation, phenotyping, and spatial analysis.

## Repository Structure

- **imcsegpipe/**  
  Contains utilities and core functions for handling IMC data:
  - Extraction functions for MCD and ZIP files.
  - TXT file matching and analysis stack creation.
  - Export functions to convert data for external tools (e.g., Histocat).
  
- **segmentation/**  
  Includes scripts and notebooks for image segmentation and exploration:
  - **Mesmer_IMC.ipynb**: A Jupyter notebook for cell segmentation using Mesmer.
  - Other segmentation utilities and visualization scripts.

- **phenotyping/**  
  Contains notebooks and utilities for downstream cellular phenotyping analysis:
  - Notebooks for data normalization, quality control, dimensionality reduction (e.g., UMAP), and clustering.
  
- **results_analysis/**  
  Encompasses notebooks and scripts for in-depth analysis of segmentation and spatial data:
  - Notebooks for cell contact analysis, distance distribution, and generating statistical reports.
  
- **Additional Scripts**  
  - **preprocessing.py**: Handles image preprocessing tasks like normalization, contrast enhancement, and hot pixel filtering.
  - **reformatting_all_files.py**: Provides functions to reformat images and standardize file metadata.

- **Environment**  
  - **environment.yml**: A Conda environment file listing all required dependencies for reproducible builds.

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

## Usage

The repository is modularized so you can work on individual components:

- **Data Extraction & Preprocessing:**  
  Run `preprocessing.py` to perform image enhancement and filter artifacts before segmentation.

- **Segmentation:**  
  Launch the segmentation notebooks (e.g., `segmentation/Mesmer_IMC.ipynb`) in Jupyter Notebook to perform and refine cell segmentation.

- **Phenotyping:**  
  Explore the notebooks in the `phenotyping/` folder, which provide tools for data normalization, clustering, and visualization.

- **Results Analysis:**  
  Use the scripts and notebooks in `results_analysis/` for spatial measurements, cell contact analysis, and statistical evaluations.

Each module includes its own documentation or help messages. See individual files and notebooks for further usage instructions.

## Contributing

Contributions are welcome! Please fork the repository and open a pull request with your proposed changes. Make sure to test your changes and update the documentation as necessary.

## License

This project is licensed under the terms specified in the LICENSE file.
