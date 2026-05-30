FROM --platform=linux/amd64 ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies
RUN apt-get -y update \
 && apt-get -y install --no-install-recommends \
      build-essential curl wget git unzip \
      python3 python3-pip python3-venv \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install miniforge (provides conda + mamba)
RUN curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
      -o /tmp/miniforge.sh \
 && bash /tmp/miniforge.sh -b -p /opt/miniforge3 \
 && rm /tmp/miniforge.sh

ENV PATH=/opt/miniforge3/bin:$PATH

# Install Jellyfish via bioconda
RUN mamba install -y -c bioconda -c conda-forge jellyfish \
 && conda clean -af -y

# Install focus
RUN pip install --no-cache-dir metagenomics-focus

# Extract database
RUN python -c "import focus, zipfile, pathlib; \
      pkg = pathlib.Path(focus.__file__).parent; \
      z = pkg / 'db.zip'; \
      zipfile.ZipFile(z).extractall(pkg) if z.exists() else None"

RUN mkdir -p /data
WORKDIR /data

ENTRYPOINT ["focus"]
CMD ["--help"]
