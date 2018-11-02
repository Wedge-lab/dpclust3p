FROM ubuntu:16.04

USER root

RUN apt-get update && apt-get -y install r-base libxml2 libxml2-dev libcurl4-gnutls-dev libssl-dev curl

RUN R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("optparse","VariantAnnotation","GenomicRanges","Rsamtools","IRanges","S4Vectors","ggplot2","reshape2"))'

RUN mkdir -p /opt/dpclust3p
COPY . /opt/dpclust3p/
RUN R -q -e 'install.packages("/opt/dpclust3p", repos=NULL, type="source")'

RUN mkdir /tmp/downloads

RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/samtools/htslib/archive/1.7.tar.gz && \
    mkdir /tmp/downloads/htslib && \
    tar -C /tmp/downloads/htslib --strip-components 1 -zxf tmp.tar.gz && \
    make -C /tmp/downloads/htslib && \
    rm -f /tmp/downloads/tmp.tar.gz

ENV HTSLIB /tmp/downloads/htslib

RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/cancerit/alleleCount/archive/v4.0.0.tar.gz && \
    mkdir /tmp/downloads/alleleCount && \
    tar -C /tmp/downloads/alleleCount --strip-components 1 -zxf tmp.tar.gz && \
    cd /tmp/downloads/alleleCount/c && \
    mkdir bin && \
    make && \
    cp /tmp/downloads/alleleCount/c/bin/alleleCounter /usr/local/bin/. && \
    cd /tmp/downloads && \
    rm -rf /tmp/downloads/alleleCount /tmp/downloads/tmp.tar.gz

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]