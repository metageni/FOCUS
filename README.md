![](logo/focus_small_logo.png "FOCUS")

# FOCUS: An Agile Profiler for Metagenomic Data

* [Installation](#installation)
* [Usage](#usage)
* [Output](#output)
* [Docker](#docker)
* [Changelog](#changelog)
* [Citing](#citing)

---

## Binning vs Profiling

FOCUS is a **profiling** tool — it estimates the relative abundance of organisms in a metagenome without assembling or binning reads.
See this [blog post](https://onestopdataanalysis.com/metagenome-profile/) for the distinction.

---

## Installation

### Dependencies

| Dependency | Version |
|---|---|
| Python | ≥ 3.8 |
| [Jellyfish](https://github.com/gmarcais/Jellyfish) | 2.x |
| NumPy | ≥ 1.12.1 |
| SciPy | ≥ 0.19.0 |
| Click | ≥ 8.0 |

### Bioconda (recommended)

```bash
conda install -c bioconda -c conda-forge metagenomics-focus
```

This installs FOCUS and Jellyfish automatically.

### pip

```bash
pip install metagenomics-focus

# also install Jellyfish separately (macOS)
brew install jellyfish
# or via bioconda
conda install -c bioconda jellyfish
```

### From source

```bash
git clone https://github.com/your-org/focus.git
cd focus

# extract the database (one-time)
cd focus_app && unzip db.zip && cd ..

# install with pip
pip install -e .

# or with uv (developer tool — not required for end users)
uv pip install -e .
```

> **Note:** `uv` is a fast developer tool for managing Python environments. It is **not** a runtime dependency and is not required to install or use FOCUS. The package installs with standard `pip` or `conda`.

---

## Usage

```
focus [OPTIONS]

  FOCUS: An Agile Profiler for Metagenomic Data.

  Example: focus -q samples/ -o results/

Options:
  -v, --version                   Show version and exit.
  -q, --query TEXT                FAST(A/Q) file or directory (repeatable).  [required]
  -o, --output_directory TEXT     Directory for output CSV files.  [required]
  -k, --kmer_size [6|7]           K-mer size.  [default: 6]
  -b, --alternate_directory TEXT  Alternate directory containing the db/ folder.
  -p, --output_prefix TEXT        Prefix for output files.  [default: output]
  -t, --threads TEXT              Jellyfish thread count.  [default: 4]
  --list_output                   Return results as a list of lists.
  -l, --log TEXT                  Log file path (default: STDOUT).
  -h, --help                      Show this message and exit.
```

### Examples

```bash
# single file
focus -q sample.fastq -o results/

# multiple files
focus -q sample1.fasta -q sample2.fastq -o results/

# directory of files, k-mer size 7, 8 threads
focus -q samples/ -o results/ -k 7 -t 8

# log to file
focus -q samples/ -o results/ -l run.log
```

---

## Output

One CSV per taxonomic level plus an all-levels file:

| File | Contents |
|---|---|
| `output_All_levels.csv` | All 8 taxonomy columns + abundance per sample |
| `output_Kingdom_tabular.csv` | Kingdom-level aggregated abundances |
| … | Class, Order, Family, Genus, Species |
| `output_Strain_tabular.csv` | Strain-level abundances |

Abundances are expressed as **percentages** (0–100).
The all-levels file is compatible with [STAMP](http://kiwi.cs.dal.ca/Software/STAMP).

---

## Docker

See [`Docker/README.md`](Docker/README.md) for build and Singularity instructions.

```bash
docker build -t focus .
docker run --rm \
  -v $(pwd)/samples:/data/samples \
  -v $(pwd)/results:/data/results \
  focus -q /data/samples -o /data/results
```

---

## Development

```bash
# install with dev dependencies (pip)
pip install -e ".[dev]"

# or with uv (faster)
uv pip install -e ".[dev]"

# run tests
pytest

# run tests with coverage
pytest --cov=focus_app --cov-report=term-missing
```

---

## Changelog

See [`CHANGELOG`](CHANGELOG).

---

## Citing

If you use FOCUS, please cite:

> Silva, G. G. Z., D. A. Cuevas, B. E. Dutilh, and R. A. Edwards, 2014:
> FOCUS: An alignment-free model to identify organisms in metagenomes
> using non-negative least squares. *PeerJ*.
