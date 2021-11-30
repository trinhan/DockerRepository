FROM mcr.microsoft.com/dotnet/core/sdk:2.1 AS build

#get source
#RUN git clone https://github.com/Illumina/Pisces.git

WORKDIR app

ADD Pisces/binaries/5.2.10.49/Pisces_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/CreateGenomeSizeFile_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/Hygea_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/Psara_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/Scylla_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/Stitcher_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/VariantQualityRecalibration_5.2.10.49.tar.gz /app
ADD Pisces/binaries/5.2.10.49/VennVcf_5.2.10.49.tar.gz /app

RUN chmod -R a+xrw /app