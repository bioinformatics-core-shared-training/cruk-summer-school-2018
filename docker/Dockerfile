FROM bioconductor/release_base
MAINTAINER Mark Dunning<mark.dunning@cruk.cam.ac.uk>
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update 
RUN apt-get install --fix-missing -y git
###Get repository of the course. Install data and R packages
#RUN apt-get install -y sra-toolkit
RUN mkdir -p /home/participant/
RUN git clone --recursive https://github.com/bioinformatics-core-shared-training/cruk-autumn-school-2017.git /home/participant/Course_Materials
RUN R -f /home/participant/Course_Materials/RNASeq/install_bioc_packages.R
RUN R -f /home/participant/Course_Materials/install_single_cell.R
WORKDIR /tmp
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip -P /tmp
RUN unzip fastqc_v0.11.3.zip
RUN sudo chmod 755 FastQC/fastqc
RUN ln -s $(pwd)/FastQC/fastqc /usr/bin/fastqc
RUN apt-get install -y bowtie2 samtools bedtools
## installing latest version of SRA toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz
RUN gunzip sratoolkit.2.8.1-3-ubuntu64.tar.gz
RUN tar xvf sratoolkit.2.8.1-3-ubuntu64.tar
RUN ln -s /tmp/sratoolkit.2.8.1-3-ubuntu64/bin/* /usr/bin/
## Install cutadapt
RUN apt-get install -y python-dev
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN sudo python get-pip.py
RUN sudo pip install cython
RUN sudo pip install --user --upgrade cutadapt
RUN rm get-pip.py
RUN wget https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2
RUN pip install MACS2
#RUN apt-get install -y bzip2 tabix bedtools build-essential git-core cmake zlib1g-dev libncurses-dev
RUN bzip2 -d htslib-1.3.1.tar.bz2
RUN tar xvf htslib-1.3.1.tar
WORKDIR htslib-1.3.1
RUN make install 
RUN chmod +x ~/.local/bin/cutadapt
RUN ln -s ~/.local/bin/cutadapt /usr/bin/cutadapt
RUN wget https://sourceforge.net/projects/samstat/files/latest/samstat-1.5.1.tar.gz
RUN tar -zxcf samstat.tgz 
WORKDIR samstat
RUN make clean
RUN make
RUN sudo cp samstat /usr/local/bin/ 
WORKDIR /tmp
RUN wget http://meme-suite.org/meme-software/4.12.0/meme_4.12.0.tar.gz
RUN tar zxf meme_4.12.0.tar.gz
WORKDIR meme_4.12.0
RUN ./configure --prefix=$HOME/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt
RUN make
RUN make test
RUN make install
WORKDIR /tmp
RUN R -f /home/participant/Course_Materials/install_chip.R
RUN apt-get install unzip
RUN wget https://ndownloader.figshare.com/articles/3219673?private_link=f5d63d8c265a05618137 -O fastq.zip
RUN unzip fastq.zip -d /home/participant/Course_Materials/RNASeq/data/
RUN rm fastq.zip
RUN wget https://ndownloader.figshare.com/articles/3219685?private_link=1d788fd384d33e913a2a -O raw.zip
RUN unzip raw.zip -d /home/participant/Course_Materials/RNASeq/data/
RUN rm raw.zip
RUN chown rstudio /home/participant/Course_Materials/
WORKDIR /home//participant/Course_Materials/
