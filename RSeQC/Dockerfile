####################################################
#RNA-seq Tools
#Dockerfile to build a container with RSeQC package and python
#Ubuntu:latest
####################################################
#Build the image based on Ubuntu
FROM ubuntu:18.04

#Install all the required packages in python
#Scripts are in: /usr/local/bin

RUN apt-get update && apt-get upgrade -y && apt-get clean
RUN apt-get install -y curl wget python3.7 python3.7-dev python3.7-distutils samtools
RUN apt-get update && apt-get install --yes default-jdk default-jre libncurses5-dev libbz2-dev liblzma-dev libsqlite3-dev sqlite3 pigz unzip bzip2 g++ make ncurses-dev 
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1
RUN update-alternatives --set python /usr/bin/python3.7
RUN curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py --force-reinstall && \
    rm get-pip.py

#RUN apt-get update -y && apt-get install -y wget git unzip bzip2 g++ make ncurses-dev 
#RUN apt-get update && apt-get install --yes default-jdk default-jre libncurses5-dev libbz2-dev liblzma-dev libsqlite3-dev sqlite3 pigz
#RUN apt-get update \#
#	&& apt-get install -y gcc \
#	&& apt-get install -y python3.8-venv \
#	&& apt-get install -y zlib1g-dev libssl-dev \
#	&& apt-get install -y samtools \
#	&& apt-get install -y python3-pip \
#	&& rm -rf /var/lib/apt/lists/*
#RUN python3.8 -m venv /venv
#ENV PATH=/venv/bin:$PATH

# install hts-lib  
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 
RUN tar -vxjf htslib-1.9.tar.bz2 && cd htslib-1.9 && make && make lib-static && make install 

RUN pip3 install numpy
RUN pip3 install RSeQC==4.0.0 
RUN pip3 install HTSeq

# copy the dockerfile
ADD Dockerfile /opt

#Change working directory

WORKDIR /