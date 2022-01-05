# TAG dnastack/pangolin:a8ee773
# Based on https://github.com/cov-lineages/pangolin/blob/master/Dockerfile
FROM continuumio/miniconda3:4.9.2-alpine

MAINTAINER Heather Ward <heather@dnastack.com>

ENV PANGOLIN_GIT_HASH a8ee7732d9c2c1bc1422c0900cff0e80aa051cea

# Install git for pangolin
RUN apk update && \
    apk add git bash

RUN cd /opt && git clone https://github.com/cov-lineages/pangolin.git
WORKDIR /opt/pangolin
RUN git checkout ${PANGOLIN_GIT_HASH}

# Python 3.8.5 already installed along with recent version of pip
# so remove Python and pip deps from environment.yml before installation
RUN sed -i "$(grep -n 'python>=3.7' ./environment.yml | cut -f1 -d:)d" ./environment.yml && \
    sed -i "$(grep -n pip= ./environment.yml | cut -f1 -d:)d" ./environment.yml

# Install the conda environment
RUN conda env create --quiet -f ./environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/pangolin/bin:$PATH

# Install Pangolin
RUN pip install . && rm -rf /root/.cache/pip
RUN pangolin --version &> /pangolin-version.txt

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name pangolin > /pangolin.yml
