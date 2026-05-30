FROM continuumio/miniconda3:24.1.2-0

# Install Jellyfish via bioconda, then focus via pip
RUN conda config --add channels defaults \
 && conda config --add channels bioconda \
 && conda config --add channels conda-forge \
 && conda install -y jellyfish \
 && conda clean -af -y

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir -e .

# Extract database if zip is present
RUN if [ -f focus/db.zip ]; then \
      unzip focus/db.zip -d focus/ && rm focus/db.zip; \
    fi

ENTRYPOINT ["focus"]
CMD ["--help"]
