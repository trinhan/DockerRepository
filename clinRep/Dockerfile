FROM rocker/tidyverse:4.0.5

# Dockerfile for Clinical Reporting scripts
# For use on MAF files annotated with VEP and oncokb

RUN apt-get update -y && apt-get install -y \
	libbz2-dev \
	liblzma-dev
RUN R -e 'install.packages(c("docopt", "vcfR", "kableExtra", "dplyr", "gplots", "DT", "yaml", "knitr", "circlize", "plotly"))'
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
RUN R -e 'BiocManager::install("GSEABase")'
RUN R -e 'install.packages(c("DT", "matrixStats"))'
RUN R -e 'BiocManager::install("karyoploteR")'
RUN R -e 'BiocManager::install("biomaRt")'


RUN mkdir template
COPY *.R /opt/
COPY *.Rmd /template/
COPY PathwayList.csv /opt/

COPY Dockerfile /opt/
