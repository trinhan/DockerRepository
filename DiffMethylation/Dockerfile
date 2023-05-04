FROM rocker/tidyverse:4.0.3

# Dockerfile for R functions needed for DifferentialMethylation project

RUN apt-get update -y && apt-get install -y \
	libbz2-dev \
	liblzma-dev
RUN R -e 'install.packages(c("naturalsort", "dplyr", "tibble", "cowplot", "knitr", "stringr", "reshape2", "ggplot2", "picante", "RColorBrewer", "pheatmap"))'
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("minfi", "GEOquery", "GenomicRanges", "TCGAbiolinks", "SummarizedExperiment", "rtracklayer", "IlluminaHumanMethylationEPICmanifest","IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))'

RUN mkdir /opt/function

COPY functions/* /opt/function
COPY *.Rmd /opt
COPY *.R /opt

COPY Dockerfile /opt/

RUN mkdir /Rlibs
RUN export R_LIBS=/Rlibs

COPY .Rprofile .
