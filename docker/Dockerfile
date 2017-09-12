FROM bioconductor/release_base
MAINTAINER Mark Dunning<mark.dunning@cruk.cam.ac.uk>
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update 
RUN apt-get install --fix-missing -y git
RUN mkdir -p /home/participant/
RUN git clone --recursive https://github.com/bioinformatics-core-shared-training/cruk-autumn-school-2017.git /home/participant/Course_Materials



# Install fastqc
WORKDIR /tmp
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip -P /tmp
RUN unzip fastqc_v0.11.3.zip
RUN sudo chmod 755 FastQC/fastqc
RUN ln -s $(pwd)/FastQC/fastqc /usr/bin/fastqc

# Install various alignment tools and bedtools
RUN apt-get install -y bowtie2 samtools bedtools bwa

## installing latest version of SRA toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz
RUN gunzip sratoolkit.2.8.1-3-ubuntu64.tar.gz
RUN tar xvf sratoolkit.2.8.1-3-ubuntu64.tar
RUN ln -s /tmp/sratoolkit.2.8.1-3-ubuntu64/bin/* /usr/bin/

## Install cutadapt and macs2
RUN apt-get install -y python-dev
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN sudo python get-pip.py
RUN sudo pip install cython
RUN sudo pip install --user --upgrade cutadapt
RUN rm get-pip.py
RUN pip install numpy
RUN pip install MACS2

RUN chmod +x ~/.local/bin/cutadapt
RUN ln -s ~/.local/bin/cutadapt /usr/bin/cutadapt

# Install samstat
RUN wget https://sourceforge.net/projects/samstat/files/latest/samstat-1.5.1.tar.gz
RUN tar -zxvf samstat-1.5.1.tar.gz
WORKDIR samstat-1.5.1
RUN ./configure
RUN make
RUN make check
RUN make install

# Install meme-chip
WORKDIR /tmp
RUN wget http://meme-suite.org/meme-software/4.12.0/meme_4.12.0.tar.gz
RUN tar zxf meme_4.12.0.tar.gz
WORKDIR meme_4.12.0
RUN ./configure  --prefix=/usr/ --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt
RUN make
#RUN make test
RUN make install

# Get data for the bulk RNA-seq materials

WORKDIR /tmp
RUN apt-get install unzip
RUN wget https://ndownloader.figshare.com/articles/3219673?private_link=f5d63d8c265a05618137 -O fastq.zip
RUN unzip fastq.zip -d /home/participant/Course_Materials/RNASeq/data/
RUN rm fastq.zip
RUN wget https://ndownloader.figshare.com/articles/3219685?private_link=1d788fd384d33e913a2a -O raw.zip
RUN unzip raw.zip -d /home/participant/Course_Materials/RNASeq/data/
RUN rm raw.zip

# Install R packages for different sections
COPY install_rna_packages.R /home/participant/Course_Materials/
COPY install_single_cell.R /home/participant/Course_Materials/
COPY install_chip.R /home/participant/Course_Materials/

# Add a couple of packages required by intro section
# Commented-out for now so I can test the command line tools first
#

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('genefilter')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('breastCancerVDX')"

RUN R -f /home/participant/Course_Materials/install_rna_packages.R
RUN R -f /home/participant/Course_Materials/install_single_cell.R
RUN R -f /home/participant/Course_Materials/install_chip.R

RUN rm -rf /tmp
RUN chown rstudio /home/participant/Course_Materials/
WORKDIR /home//participant/Course_Materials/
