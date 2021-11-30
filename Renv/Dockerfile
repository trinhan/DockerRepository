FROM r-base
WORKDIR /app    
RUN Rscript -e "install.packages(c('docopt', 'vcfR', 'matrixStats'))"
COPY MergeVCFRecall.R /app/
COPY Dockerfile /opt/
