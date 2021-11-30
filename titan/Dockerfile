FROM rocker/tidyverse:4.0.5

# TitanCNA

RUN apt-get update -y && apt-get install -y \
	libbz2-dev \
	liblzma-dev
RUN R -e 'install.packages(c("doMC", "optparse"))'
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'

RUN R -e 'BiocManager::install("Rhtslib")'
RUN R -e 'BiocManager::install("TitanCNA")'
##RUN R -r 'BiocManager::install("SNPchip")'
RUN R -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'

RUN mkdir -p TitanCNA TitanCNA/R TitanCNA/scripts/

COPY titanCNA.R TitanCNA/scripts/titanCNA.R
COPY selectSolution.R TitanCNA/scripts/selectSolution.R
COPY R/* TitanCNA/R/*

WORKDIR .