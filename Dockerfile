# conda setup
FROM continuumio/miniconda3:4.7.10
RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

# creating conda environment
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# installs focus using bioconda
RUN conda install focus

# Entrypoint
CMD /bin/bash
