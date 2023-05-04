FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.8 python3-pip python3-setuptools python3-dev git
RUN apt-get install -y r-base

RUN mkdir -p /oncokb
RUN mkdir -p /data
RUN git clone -b v3.2.4 https://github.com/oncokb/oncokb-annotator.git /oncokb
WORKDIR /oncokb
RUN pip3 install wheel setuptools pip --upgrade
RUN pip3 install -r requirements/common.txt -r requirements/pip3.txt

## INSTALL R packages
RUN Rscript -e "install.packages(c('docopt', 'matrixStats'))"
RUN Rscript -e "install.packages('vcfR')"

COPY vepVCF2maf.R /opt/

COPY Dockerfile /opt/
