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
    && ln -s /usr/local/samtools /usr/local/strelka_workflow-1.0.15/src/perl/opt
#    && rm strelka_workflow-1.0.15.tar.gz

RUN apt-get update && apt-get install -y \
    libncurses-dev \
    vim \
    default-jre \
    && apt-get clean

RUN cd /usr/local \
    && git clone -b docker https://mwyczalkowski:Elohim,598@github.com/ding-lab/somaticwrapper 

# required for running bsub_varscan script - apparently not
#RUN apt-get update && apt-get install -y \
#    bc \
#    && apt-get clean

# Install varscan 2.3.8 to match SomaticWrapper
# Note that newer versions available here: https://github.com/dkoboldt/varscan
RUN cd /usr/local \
    && wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.8.jar


WORKDIR /usr/local/somaticwrapper

# following invocation will create common data dir in ./
# docker run -v $PWD/data:/data -it breakpoint_surveyor
VOLUME ["/data"]

CMD ["/bin/bash"]
