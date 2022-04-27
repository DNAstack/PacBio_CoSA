# TAG dnastack/cosa:c0a2fa8
FROM ubuntu:bionic

MAINTAINER Heather Ward <heather@dnastack.com>

ENV COSA_GIT_HASH c0a2fa8b538f226195d61f91b4fbfdaf41d9a99e
ENV SAMTOOLS_VERSION 1.12
ENV MINIMAP2_VERSION 2.18
ENV BAMTOOLS_GIT_HASH de86007b944bba1e5305fcf8cf24e84262b5fcad

RUN apt-get -qq update && \
	apt-get -qq install \
	git \
	python3.8-dev \
	python3-pip \
	python3-setuptools \
	tabix \
	wget \
	make \
	gcc \
	build-essential \
	unzip \
	libz-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	cmake \
	libjsoncpp-dev

# CoSA
RUN python3.8 -m pip install pyvcf Cython

RUN cd /opt && git clone https://github.com/PacificBiosciences/CoSA.git
RUN cd /opt/CoSA && \
	git checkout ${COSA_GIT_HASH} && \
	python3.8 setup.py build && \
	python3.8 setup.py install

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt/ && \
	rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} && \
	./configure --prefix=/usr/local/ && \
	make && \
	make install

RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 && \
	tar -jxvf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 --directory /opt/ && \
	rm minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2

ENV PATH $PATH:/opt/minimap2-${MINIMAP2_VERSION}_x64-linux/

RUN cd /opt && git clone https://github.com/pezmaster31/bamtools.git
RUN cd /opt/bamtools && \
	git checkout ${BAMTOOLS_GIT_HASH} && \
	mkdir build && \
	cd build && \
	cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
	make && \
	make install
