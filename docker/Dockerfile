FROM mwyczalkowski/breakpoint_surveyor:latest

# Based on /home/mwyczalk_test/data/BPS.TCGA_Virus.Lite/docker/Dockerfile

LABEL author="Matt Wyczalkowski" \
      maintainer="m.wyczalkowski@wustl.edu"

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

RUN apt-get update && apt-get install -y \
    default-jre \
    libncurses-dev \
    sshpass \
    vim \
    && apt-get clean

RUN cd /usr/local \
    && git clone -b docker https://mwyczalkowski:Elohim,598@github.com/ding-lab/somaticwrapper 


RUN apt-get update && apt-get install -y \
    cmake \
    && apt-get clean

# Install varscan 2.3.8 to match SomaticWrapper
# Note that newer versions available here: https://github.com/dkoboldt/varscan
RUN cd /usr/local \
    && wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.8.jar \
    && ln -s VarScan.v2.3.8.jar VarScan.jar

# install bam-readcount, as per https://gist.github.com/ckandoth/87ba44948cb747916f8d
RUN cd /usr/local \
    wget https://github.com/genome/bam-readcount/archive/master.tar.gz \
    tar -zxf master.tar.gz \
    rm master.tar.gz \
    mkdir bam-readcount-master/build \
    cd bam-readcount-master/build \
    cmake ../ \
    make \
    cp bin/bam-readcount /usr/local/bin

# install SnpEff, per http://snpeff.sourceforge.net/download.html
RUN cd /usr/local \
    wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    unzip snpEff_latest_core.zip \
    rm snpEff_latest_core.zip

# install VEP - this downloads a cache, which should probably be installed as a runtime step (like reference)
# Rather than installed in docker imate.  See here for discussion of caches: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
# So the plan is to 1) install VEP in docker, with no cache, and 2) when starting first container, download cache then
# TODO: (I did CTRL-C to move on)  This will take an hour or so prob
#  - use CTRL-C to cancel if you do not wish to convert this cache now (you may run convert_cache.pl later)

#cd /usr/local
#cpanm DBI
#git clone https://github.com/Ensembl/ensembl-vep.git
#cd ensembl-vep
#perl INSTALL.pl
#  - y, y
#  - homo_sapiens_vep_90_GRCh37.tar.gz homo_sapiens_vep_90_GRCh37.tar.gzhomo_sapiens_vep_90_GRCh37.tar.gz
#  - this is typically 59 60
# cache files stored in /root/.vep  -> this should probably be set to /data/D_Ensembl
# See also here for how to install cache separately in a docker environment: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#docker

WORKDIR /usr/local/somaticwrapper

COPY vimrc /root/.vimrc

# following invocation will create common data dir in ./
# docker run -v $PWD/data:/data -it breakpoint_surveyor
VOLUME ["/data"]

CMD ["/bin/bash"]
