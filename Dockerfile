FROM johanngb/rep-int

RUN R -e "install.packages('DEoptim')"
RUN R -e "install.packages('nloptr')"
RUN R -e "install.packages('abind')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('tidyr')"
RUN R -e "install.packages('readr')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('reshape2')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('scales')"
RUN R -e "install.packages('ggplotify')"
RUN R -e "install.packages('GGally')"


WORKDIR /home/rep

# raw data
COPY ./analysis/raw_data analysis/raw_data

# processed data
COPY ./analysis/processed_data analysis/processed_data

# packages to include
COPY ./r_packages r_packages

RUN R -e "install.packages('r_packages/rrscale_rpkg/rrscale_1.0.3.tar.gz', repos=NULL,type='source')"
RUN R -e "install.packages('r_packages/mema_norm_rpkg/memanorm_0.0.0.9004.tar.gz', repos=NULL,type='source')"

# analysis files to include
COPY analysis/analysis_plots.R analysis/analysis_plots.R
COPY analysis/analysis_plots.Rmd analysis/analysis_plots.Rmd
COPY analysis/analysis_plots.ipynb analysis/analysis_plots.ipynb

COPY analysis/plot_scripts/make_fmat.R analysis/plot_scripts/make_fmat.R
COPY analysis/plot_scripts/plot_fns.R analysis/plot_scripts/plot_fns.R

COPY analysis/process.R analysis/process.R
COPY analysis/process.Rmd analysis/process.Rmd
COPY analysis/process.ipynb analysis/process.ipynb
COPY analysis/processing_scripts/estimate.R analysis/processing_scripts/estimate.R
COPY analysis/processing_scripts/make_data_matrices.R analysis/processing_scripts/make_data_matrices.R
COPY analysis/processing_scripts/make_fmat.R analysis/processing_scripts/make_fmat.R
COPY analysis/processing_scripts/util.R analysis/processing_scripts/util.R

# other files
COPY LICENSE LICENSE
COPY README.md README.md
