FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev git
RUN apt-get install -y r-base

RUN mkdir -p /oncokb
RUN mkdir -p /data
RUN git clone https://github.com/oncokb/oncokb-annotator.git /oncokb
WORKDIR /oncokb
RUN pip3 install -r requirements/common.txt -r requirements/pip3.txt

## INSTALL R packages
RUN Rscript -e "install.packages(c('docopt', 'matrixStats'))"
RUN Rscript -e "install.packages('vcfR')"

COPY vepVCF2maf4Oncokb.R /opt/
COPY annotateProteins.R /opt/

COPY Dockerfile /opt/
