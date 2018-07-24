FROM bioconductor/release_base
MAINTAINER Mark Fernandes<mark.fernandes@cruk.cam.ac.uk>
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update 
RUN apt-get install --fix-missing -y git
RUN mkdir -p /home/participant/
RUN git clone --recursive https://github.com/bioinformatics-core-shared-training/cruk-autumn-school-2017.git /home/participant/Course_Materials


WORKDIR /home//participant/Course_Materials/
