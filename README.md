# Project Title

A comprehensive repository containing tools and analyses for imaging mass cytometry data processing, segmentation, and phenotyping.

## Overview

This repository is organised into multiple modules:
- **AnnoSpat**: Tools for spatial annotation and cell classification.
- **imcsegpipe**: Utilities for image segmentation of mass cytometry data.
- **phenotyping**: Scripts and notebooks for downstream phenotyping analysis.
- **results_analysis**: Notebooks and scripts for visualisation and statistical evaluation.
- **Additional Scripts**: Preprocessing, post-processing, and data management utilities.

## Installation

Clone the repository and install the required dependencies using pip:

```bash
git clone <repo_url>
cd <repo_folder>
pip install -r requirements.txt
```

## Usage

Each module includes its own documentation or README file with usage examples. Some starting points:
- **AnnoSpat**: `python AnnoSpat/AnnoSpat_main/run.py generateLabels --help`
- **imcsegpipe**: `python -m imcsegpipe`
- **phenotyping**: See the notebooks in the `phenotyping/` folder

## Contributing

Contributions are welcome. Please fork the repository and open a pull request with your proposed changes.

## License

This project is licensed under the terms specified in the LICENSE file.
