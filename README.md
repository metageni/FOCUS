![](logo/focus_small_logo.png "Logo")

#### FOCUS: An Agile Profiler for Metagenomic Data
* [Installation](#installation)
    * [Dependencies](#dependencies)
    * [Bioconda](#bioconda)
    * [Pip](#pip)
    * [Git](#git)
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
- unzip


### Bioconda
You can now easily install FOCUS using [conda](https://conda.io) via the
[Bioconda](https://bioconda.github.io/) channel. It is as easy as:

    conda create -n focus -c bioconda focus
	source activate focus

This will create a conda environment called `focus` (as specified by the
`-n` argument), and install FOCUS along with all its dependencies. The second
line activates the newly created `focus` conda environment.

### Pip
	pip install metagenomics-focus

### Git

	# clone focus
	git clone git@github.com:metageni/FOCUS.git

	# unzip database and move it to folder
	unzip FOCUS/focus_app/db.zip && mv db FOCUS/focus_app/

	# install focus
	cd FOCUS && python setup.py install

## Usage
	focus -q INPUT_DIR -k -o OUTPUT_DIR

		-q Path to directory with FASTA/FASTQ file(s)

		-o Path to output directory

		-k k-mer size (6 or 7 avaliable) (Default: 6)

    	-p Ouput prefix (Default: output)


## Output
FOCUS generates a tabular output per taxonomic level (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`, and `Strain`) and one with all levels which can be used as [STAMP](http://kiwi.cs.dal.ca/Software/STAMP)'s input for statistical analysis.

## Citing
FOCUS was written by Genivaldo G. Z. Silva. Feel free to [contact me](mailto:genivaldo.gueiros@gmail.com)

If you use FOCUS, please cite it:

    Silva, G. G. Z., D. A. Cuevas, B. E. Dutilh, and R. A. Edwards, 2014:
    FOCUS: An alignment-free model to identify organisms in metagenomes
    using non-negative least squares. PeerJ.
