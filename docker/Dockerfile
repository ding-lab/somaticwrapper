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
    
# Install Pindel
RUN cd /usr/local/ \
    && git clone --recursive https://github.com/genome/pindel.git  \
    && cd /usr/local/pindel/ \
    && ./INSTALL /usr/local/htslib

# Various Python packages
RUN pip install \
    pyvcf \
    pysam \
    numpy

# Perl packags
RUN cpanm XML::XPath

# Strelka v 2.7.1
#  -> incorporate the newest Strelka at some point 
#RUN cd /usr/local \
#    && wget https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2 \
#    && bzip2 -d strelka-2.7.1.centos5_x86_64.tar.bz2 \
#    && tar -xf strelka-2.7.1.centos5_x86_64.tar \
#    && rm strelka-2.7.1.centos5_x86_64.tar

# Note that SomaticWrapper v1.1 uses strelka 1.0.14.  We download that too
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

COPY vimrc /root/.vimrc

RUN cd /usr/local \
    &&  git clone -b docker https://github.com/ding-lab/somaticwrapper 

WORKDIR /usr/local/somaticwrapper


# following invocation will create common data dir in ./
# docker run -v $PWD/data:/data -it breakpoint_surveyor
VOLUME ["/data"]

CMD ["/bin/bash"]
