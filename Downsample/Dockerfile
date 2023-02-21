FROM rocker/tidyverse:4.0.3

# R container needed for downsampling analysis

RUN apt-get update -y && apt-get install -y \
	libbz2-dev \
	liblzma-dev

RUN R -e 'install.packages("R.utils")'
RUN R -e 'BiocManager::install("GenomicAlignments")'

# Create a directory to mount extra R libraries incase they are needed later on. E.g. for HPC use
# If using singularity use singularity exec --no-home -B Rlibs:/Rlibs

RUN mkdir /Rlibs
RUN export R_LIBS=/Rlibs
COPY .Rprofile .

# Copy the required scripts and this Dockerfile to /opt/

COPY *.R /opt
COPY Dockerfile /opt/


