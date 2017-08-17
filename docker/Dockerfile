FROM mwyczalkowski/breakpoint_surveyor:latest

# Based on /home/mwyczalk_test/data/BPS.TCGA_Virus.Lite/docker/Dockerfile

LABEL author="Matt Wyczalkowski" \
      maintainer="m.wyczalkowski@wustl.edu"
# Strelka
RUN cd /usr/local \
    && wget https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2 \
    && bzip2 -d strelka-2.7.1.centos5_x86_64.tar.bz2 \
    && tar -xf strelka-2.7.1.centos5_x86_64.tar \
    && rm strelka-2.7.1.centos5_x86_64.tar

RUN apt-get update && apt-get install -y \
    libncurses-dev \
    vim \
    default-jre \
    && apt-get clean

RUN cd /usr/local \
    && git clone -b docker https://mwyczalkowski:Elohim,598@github.com/ding-lab/somaticwrapper 


WORKDIR /usr/local/somaticwrapper

# following invocation will create common data dir in ./
# docker run -v $PWD/data:/data -it breakpoint_surveyor
VOLUME ["/data"]

CMD ["/bin/bash"]
