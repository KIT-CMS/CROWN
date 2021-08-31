FROM mschnepf/slc7-condocker:latest

LABEL maintainer="SEBASTIAN BROMMER  <sebastian.brommer@cern.ch>"

# install conda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/anaconda/ \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && chmod ugo+w /opt/anaconda

ARG conda_env=KingMaker
ENV PATH="/opt/anaconda/bin:${PATH}"
ARG PATH="/opt/anaconda/bin:${PATH}"

# Initialize conda in bash config fiiles:
RUN conda init --system bash

# get environment from KingMaker repository
RUN wget https://raw.githubusercontent.com/harrypuuter/KingMaker/main/environment.yml \
    && conda env create -f environment.yml