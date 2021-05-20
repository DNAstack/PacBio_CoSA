# TAG dnastack/pacbio_tools:1
FROM ubuntu:bionic

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		curl

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
	rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH $PATH:/opt/miniconda3/bin/

RUN conda config --add channels bioconda
RUN conda install pbccs lima pbaa

