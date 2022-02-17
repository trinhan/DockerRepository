FROM ubuntu:18.04

LABEL \
  version="3.1" \
  description="Docker image to run AnnotSV"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && apt-get install -y \
 build-essential \
  git \
  curl \
  g++ \
  make \
  python \
  tar \
  tcl \
  tcllib \
  unzip \
  wget \
  zlib1g-dev \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
  libperl-dev \
libssl-dev

ENV BEDTOOLS_INSTALL_DIR=/opt/bedtools2
ENV BEDTOOLS_VERSION=2.28.0

WORKDIR /tmp
RUN wget https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VERSION/bedtools-$BEDTOOLS_VERSION.tar.gz && \
  tar -zxf bedtools-$BEDTOOLS_VERSION.tar.gz && \
  rm -f bedtools-$BEDTOOLS_VERSION.tar.gz

WORKDIR /tmp/bedtools2
RUN make && \
  mkdir --parents $BEDTOOLS_INSTALL_DIR && \
  mv ./* $BEDTOOLS_INSTALL_DIR

WORKDIR /
RUN ln -s $BEDTOOLS_INSTALL_DIR/bin/* /usr/bin/ && \
  rm -rf /tmp/bedtools2

ENV ANNOTSV_VERSION=3.1
ENV ANNOTSV_COMMIT=eb519723c8445187281762ef3c059c11630ad7d0
ENV ANNOTSV=/opt/AnnotSV_$ANNOTSV_VERSION

WORKDIR /opt
RUN wget https://github.com/lgmgeo/AnnotSV/archive/${ANNOTSV_COMMIT}.zip && \
  unzip ${ANNOTSV_COMMIT}.zip && \
  mv AnnotSV-${ANNOTSV_COMMIT} ${ANNOTSV} && \
  rm ${ANNOTSV_COMMIT}.zip && \
  cd ${ANNOTSV} && \
  make PREFIX=. install 
## && don't install this at the time being
##  make PREFIX=. install-human-annotation


ARG htsversion=1.13

WORKDIR /opt
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig && \
    curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && make install) && \
    git clone --depth 1 git://github.com/samtools/htslib-plugins && \
    (cd htslib-plugins && make PLUGINS='hfile_cip.so hfile_mmap.so' install)


ENV PATH="${ANNOTSV}/bin:${PATH}"
WORKDIR /