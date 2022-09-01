FROM jupyter/datascience-notebook
#FROM johanngb/rep-int

RUN R -e "install.packages('DEoptim',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('nloptr',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('abind',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dplyr',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('data.table',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('tidyr',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('readr',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('ggplot2',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('reshape2',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('cowplot',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('scales',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('ggplotify',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('GGally',repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('viridis',repos='http://cran.us.r-project.org')"


WORKDIR /home/jovyan

# R PACKAGES
COPY --chown=1000 ./r_packages r_packages

RUN R -e "install.packages('r_packages/rrscale_rpkg/rrscale_1.0.3.tar.gz', repos=NULL,type='source')"
RUN R -e "install.packages('r_packages/mema_norm_rpkg/memanorm_0.0.0.9005.tar.gz', repos=NULL,type='source')"


# ANALYSIS
RUN mkdir analysis

# plot_scripts
COPY --chown=1000 analysis/plot_scripts/make_fmat.R analysis/plot_scripts/make_fmat.R
COPY --chown=1000 analysis/plot_scripts/plot_fns.R analysis/plot_scripts/plot_fns.R

# processed_data
COPY --chown=1000 ./analysis/processed_data/MCF10A analysis/processed_data/MCF10A

# processing_scripts
COPY --chown=1000 analysis/processing_scripts/estimate.R analysis/processing_scripts/estimate.R
COPY --chown=1000 analysis/processing_scripts/make_data_matrices.R analysis/processing_scripts/make_data_matrices.R
COPY --chown=1000 analysis/processing_scripts/make_fmat.R analysis/processing_scripts/make_fmat.R
COPY --chown=1000 analysis/processing_scripts/util.R analysis/processing_scripts/util.R

# raw_data
COPY --chown=1000 ./analysis/raw_data/MCF10A analysis/raw_data/MCF10A


# analysis files to include
COPY --chown=1000 analysis/analysis_plots.ipynb analysis/analysis_plots.ipynb
COPY --chown=1000 analysis/analysis_plots.R analysis/analysis_plots.R
COPY --chown=1000 analysis/analysis_plots.Rmd analysis/analysis_plots.Rmd

COPY --chown=1000 analysis/process.ipynb analysis/process.ipynb
COPY --chown=1000 analysis/process.R analysis/process.R
COPY --chown=1000 analysis/process.Rmd analysis/process.Rmd

COPY --chown=1000 analysis/d_sensitivity.ipynb analysis/d_sensitivity.ipynb

# other files
COPY --chown=1000 LICENSE LICENSE
COPY --chown=1000 README.md README.md
