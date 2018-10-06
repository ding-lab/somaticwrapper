FROM ubuntu:xenial

LABEL author="Matt Wyczalkowski" \
      maintainer="m.wyczalkowski@wustl.edu"

RUN apt-get update && apt-get install -y \
    autoconf \
    bedtools \
    build-essential \
    bwa \
    cmake \
    cpanminus \
    curl \
    default-jre \
    gawk \
    git \
    libbz2-dev \ 
    liblzma-dev \
    libncurses-dev \
    python \
    python-dev \
    python-pip \
    sshpass \
    subversion \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    && apt-get clean

# Installing HTSLIB
RUN cd /usr/local \
    && git clone https://github.com/samtools/htslib.git \
    && cd /usr/local/htslib \
    && make \
    && make lib-static \
    && make install

ENV LD_LIBRARY_PATH=/usr/local/htslib

# Install samtools
RUN cd /usr/local/ \
    && git clone https://github.com/samtools/samtools.git \
    && cd samtools \
    && make \
    && make install
    
# Various Python packages
RUN pip install \
    pyvcf \
    pysam \
    numpy

# Perl packags
RUN cpanm XML::XPath

# Notes on Strelka2 installation
# there are issues with installing strelka in /usr/local, because strelka1 already exists
# for now, removing /usr/local/strelka and /usr/local/strelka_workflow-1.0.15
# cd /usr/local
# git clone https://github.com/Illumina/strelka.git  strelka2
# mkdir strelka2_build; cd strelka2_build
# /usr/local/strelka2/configure --prefix=/usr/local/strelka2
# make -C /usr/local/strelka2_build
# make install

# Note that SomaticWrapper v1.1 uses strelka 1.0.15.  We download that too
RUN cd /usr/local \
    && wget https://sites.google.com/site/strelkasomaticvariantcaller/home/download/strelka_workflow-1.0.15.tar.gz \
    && tar -zxf strelka_workflow-1.0.15.tar.gz \
    && cd /usr/local/strelka_workflow-1.0.15 \
    && ./configure --prefix=/usr/local/strelka \
    && make \
    && make install \
    && ln -s /usr/local/samtools /usr/local/strelka_workflow-1.0.15/src/perl/opt \
    && rm /usr/local/strelka_workflow-1.0.15.tar.gz


# Install varscan 2.3.8 to match SomaticWrapper
# Note that newer versions available here: https://github.com/dkoboldt/varscan
RUN cd /usr/local \
    && wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.8.jar \
    && ln -s VarScan.v2.3.8.jar VarScan.jar

# install bam-readcount, as per https://gist.github.com/ckandoth/87ba44948cb747916f8d
RUN cd /usr/local \
    && wget https://github.com/genome/bam-readcount/archive/master.tar.gz \
    && tar -zxf master.tar.gz \
    && rm master.tar.gz \
    && mkdir bam-readcount-master/build \
    && cd bam-readcount-master/build \
    && cmake ../ \
    && make \
    && cp bin/bam-readcount /usr/local/bin

# install SnpEff, per http://snpeff.sourceforge.net/download.html
RUN cd /usr/local \
    && wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip snpEff_latest_core.zip \
    && rm snpEff_latest_core.zip

# Install GATK and picard
RUN cd /usr/local \
    && wget https://software.broadinstitute.org/gatk/download/auth?package=GATK -O GATK.tar.bz2 \
    && bunzip2 GATK.tar.bz2 \
    && tar -xf GATK.tar \
    && ln -s GenomeAnalysisTK* GenomeAnalysisTK \
    && wget https://github.com/broadinstitute/picard/releases/download/2.12.1/picard.jar

# DEVELOPMENT 
# Install mysql so that can use db-based VEP annotation 
# NOTE: this can be removed if using VEP cache on e.g. production servers
# based on
# https://stackoverflow.com/questions/38356219/mysql-in-docker-frozen-at-root-password-config
#    && apt-get install -y mysql-server libmysqlclient-dev libnss-sss \
RUN apt-get update \
    && { \
        echo debconf debconf/frontend select Noninteractive; \
        echo mysql-community-server mysql-community-server/data-dir select ''; \
        echo mysql-community-server mysql-community-server/root-pass password ''; \
        echo mysql-community-server mysql-community-server/re-root-pass password ''; \
        echo mysql-community-server mysql-community-server/remove-test-db select true; \
    } | debconf-set-selections \
    && apt-get install -y mysql-server libmysqlclient-dev \
    && apt-get clean
RUN cpanm DBI DBD::mysql 

# install VEP - we do not install cache into the image, it will be downloaded during setup
# See here for discussion of caches: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
# Cache dir will be /data/D_Ensembl

RUN cd /usr/local \
    && git clone https://github.com/Ensembl/ensembl-vep.git  \
    && cd /usr/local/ensembl-vep \
    && perl INSTALL.pl /data/D_Ensembl -a a

# Install Strelka2 
RUN cd /usr/local \
    && git clone https://github.com/Illumina/strelka.git strelka2 \
    && mkdir strelka2_build \
    && cd strelka2_build \
    && /usr/local/strelka2/configure --prefix=/usr/local/strelka2 \
    && make -C /usr/local/strelka2_build \
    && make install

# Install Pindel.  Using ding-lab branch
RUN cd /usr/local/ \
    && git clone --recursive https://github.com/ding-lab/pindel.git \
    && cd /usr/local/pindel/ \
    && ./INSTALL /usr/local/htslib

# Install vcf2maf: https://github.com/mskcc/vcf2maf
RUN cd /usr/local \
    && export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4` \
    && curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL \
    && tar -zxf mskcc-vcf2maf.tar.gz \
    && ln -s mskcc-vcf2maf-* mskcc-vcf2maf 

# Note that this creates directory /usr/local/mskcc-vcf2maf-decbf60.
# TODO: ln -s /usr/local/mskcc-vcf2maf-decbf60 /usr/local/mskcc-vcf2maf

COPY vimrc /root/.vimrc

# change the hash below to force clone
# get it with git rev-parse HEAD, though it doesn't matter what it is as long as it changes
ENV latest_commit boo102
RUN cd /usr/local \
    && git clone -b cwl-dev https://github.com/ding-lab/somaticwrapper 

WORKDIR /usr/local/somaticwrapper

ENV PERL5LIB=/usr/local/somaticwrapper
ENV LC_ALL="C"


# following invocation will create common data dir in ./
# docker run -v $PWD/data:/data -it breakpoint_surveyor
#VOLUME ["/data"]

CMD ["/bin/bash"]
