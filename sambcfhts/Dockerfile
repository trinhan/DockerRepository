FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
	build-essential \
	curl \
	git \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
	libperl-dev \
	libssl-dev \
	zlib1g-dev \
	openjdk-8-jdk \
	python3 \
	python3-pip \
	python3-setuptools \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /tmp
ARG htsversion=1.13
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig && \
    curl -L https://github.com/samtools/samtools/releases/download/${htsversion}/samtools-${htsversion}.tar.bz2 | tar xj && \
    (cd samtools-${htsversion} && ./configure --with-htslib=system && make install) && \
    curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && make install) && \
    git clone --depth 1 git://github.com/samtools/htslib-plugins && \
    (cd htslib-plugins && make PLUGINS='hfile_cip.so hfile_mmap.so' install)

ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

RUN pip3 install -U pip setuptools 
RUN pip3 install scikit-allel pandas

ADD picard.jar /tmp

ADD M1_txt2vcf.py /usr/local/bin/
ADD filter_indels.py/ /usr/local/bin/
ADD filter_passed_mutations.py /usr/local/bin/
ADD strelka_allelic_count_indel.py /usr/local/bin/
ADD strelka_allelic_count_snv.py /usr/local/bin/
ADD merge_callstats.py /usr/local/bin/
ADD vcf2mafbed.py /usr/local/bin/
ADD Dockerfile /opt/