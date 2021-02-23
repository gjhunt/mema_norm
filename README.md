# Normalization of MEMAs

### Analysis files

Most important analysis files are

1. ``process.ipynb`` (or mirrored ``process.R`` and ``process.Rmd``)

This takes data from raw data to RR transformed data to RR and normalized data and saves it in various formats. Raw level 2 ``.tsv`` files are expected to be in the directory ``analysis/raw_data/MCF10A/``

2. ``analysis_plots.ipynb`` (or mirrored ``analysis_plots.R`` and ``analysis_plots.Rmd``)

This takes processed data and prodcues plots. This needs ``MCF10A_15_df.csv`` to be in the directory ``analysis/processed_data/csv/`` but doesn't require any of the lower-level data. 

### R packages

The underlying R package code may be found in the ``r_packages`` folder