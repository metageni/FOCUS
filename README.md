![](logo/focus_small_logo.png "Logo")

#### FOCUS: An Agile Profiler for Metagenomic Data
* [Installation](#installation)
  * [dependencies](#dependencies)
  * [pip3](#pip3)
  * [bioconda](#bioconda)
  * [git](#git)
* [Usage](#usage)
* [Output](#output)
* [Citing](#citing)

## Installation
### Dependencies
  - [Python 3.6](http://www.python.org/download)
  - [Setuptools 36.0.1](https://setuptools.readthedocs.io/en/latest/)
  - [Jellyfish 2.2.6](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6). if using macOS, use [bioconda](https://anaconda.org/bioconda/jellyfish)
  - [Numpy 1.12.1](https://github.com/numpy/numpy)
  - [SciPy 0.19.0](https://github.com/scipy/scipy)
  - unzip/curl

### pip3
	# pip3 also install numpy and scipy
	pip3 install metagenomics-focus

### Bioconda
You can now easily install FOCUS using [conda](https://conda.io) via the
[Bioconda](https://bioconda.github.io/) channel. It is as easy as:

    # bioconda should handle all the dependencies
    conda create -n focus -c bioconda focus
	source activate focus

This will create a conda environment called `focus` (as specified by the
`-n` argument), and install FOCUS along with all its dependencies. The second
line activates the newly created `focus` conda environment.

### Git

	# these steps should install Numpy and Scipy

	# clone focus
	git clone git@github.com:metageni/FOCUS.git

	# unzip database and move it to folder
	unzip FOCUS/focus_app/db.zip && mv db FOCUS/focus_app/

	# install focus
	cd FOCUS && python setup.py install

## Usage
    focus [-h] [-v] -q QUERY -o OUTPUT_DIRECTORY [-k KMER_SIZE]
             [-b ALTERNATE_DIRECTORY] [-p OUTPUT_PREFIX] [-t THREADS]
             [--list_output] [-l LOG]

    FOCUS: An Agile Profiler for Metagenomic Data

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -q QUERY, --query QUERY
                            Path to directory with FAST(A/Q) files
      -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                            Path to output files
      -k KMER_SIZE, --kmer_size KMER_SIZE
                            K-mer size (6 or 7) (Default: 6)
      -b ALTERNATE_DIRECTORY, --alternate_directory ALTERNATE_DIRECTORY
                            Alternate directory for your databases
      -p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output prefix (Default: output)
      -t THREADS, --threads THREADS
                            Number Threads used in the k-mer counting (Default: 4)
      --list_output         Output results as a list
      -l LOG, --log LOG     Path to log file (Default: STDOUT).

    example > focus -q samples_directory


## Output
FOCUS generates a tabular output per taxonomic level (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`, and `Strain`) and one with all levels which can be used as [STAMP](http://kiwi.cs.dal.ca/Software/STAMP)'s input for statistical analysis.

## Citing
FOCUS was written by Genivaldo G. Z. Silva. Feel free to [contact me](mailto:genivaldo.gueiros@gmail.com)

If you use FOCUS, please cite it:

    Silva, G. G. Z., D. A. Cuevas, B. E. Dutilh, and R. A. Edwards, 2014:
    FOCUS: An alignment-free model to identify organisms in metagenomes
    using non-negative least squares. PeerJ.
