#!/usr/bin/env python3
"""Unit tests for focus_app.focus core logic.

Covers: normalise, is_wanted_file, load_database, run_nnls, aggregate_level,
        refine_results, write_results, and count_kmers (when Jellyfish available).
"""

import csv
import random
import tempfile

import pytest

from pathlib import Path
from shutil import which

import numpy as np

from focus_app.focus import (
    normalise,
    is_wanted_file,
    load_database,
    run_nnls,
    aggregate_level,
    refine_results,
    write_results,
    count_kmers,
)

DATA = Path(__file__).parent / "data"
DB_SMALL = DATA / "k6_small_sample.txt"
DB_EMPTY_ROW = DATA / "k6_small_sample_empty_record.txt"
FASTA = DATA / "mock_sample.fasta"
FASTA_EMPTY = DATA / "mock_sample_empty.fasta"
FASTQ = DATA / "mock_sample.fastq"


# ---------------------------------------------------------------------------
# normalise
# ---------------------------------------------------------------------------

def test_normalise_basic():
    """Equal values should produce uniform distribution."""
    result = normalise([1, 1])
    np.testing.assert_array_almost_equal(result, [0.5, 0.5])


def test_normalise_unequal():
    """Unequal values should sum to 1."""
    result = normalise([2, 2, 2, 2])
    np.testing.assert_array_almost_equal(result, [0.25, 0.25, 0.25, 0.25])


def test_normalise_zero_raises():
    """All-zero input must raise RuntimeWarning."""
    with pytest.raises(RuntimeWarning):
        normalise([0, 0, 0])


# ---------------------------------------------------------------------------
# is_wanted_file
# ---------------------------------------------------------------------------

def test_is_wanted_file_filters_and_sorts():
    """Only .fasta/.fastq/.fna files should be returned, sorted."""
    paths = ["b.fastq", "a.fasta", "n.fna", "x.png", "y.txt"]
    result = is_wanted_file(paths)
    assert [p.name for p in result] == ["a.fasta", "b.fastq", "n.fna"]


def test_is_wanted_file_case_insensitive():
    """Extension matching must be case-insensitive."""
    paths = ["X.FASTQ", "Y.FASTA", "Z.FNA"]
    result = is_wanted_file(paths)
    assert len(result) == 3


def test_is_wanted_file_empty():
    """Non-sequence files should return empty list."""
    assert is_wanted_file(["a.png", "b.csv"]) == []


# ---------------------------------------------------------------------------
# load_database
# ---------------------------------------------------------------------------

def test_load_database_shape():
    """Matrix shape should be (n_kmers, n_organisms)."""
    matrix, organisms, kmer_order = load_database(DB_SMALL)
    assert len(kmer_order) == 3
    assert len(organisms) == 2
    assert matrix.shape == (3, 2)


def test_load_database_kmer_order():
    """K-mer order must match the header of the database file."""
    _, _, kmer_order = load_database(DB_SMALL)
    assert kmer_order == ["GAACGC", "GAACGA", "CACCCA"]


def test_load_database_organisms():
    """Organism labels must be tab-joined taxonomy strings."""
    _, organisms, _ = load_database(DB_SMALL)
    assert organisms[0].startswith("Bacteria\tSpirochaetes")
    assert organisms[1].startswith("Bacteria\tFirmicutes")


def test_load_database_normalised():
    """Each organism column should sum to 1 (normalised)."""
    matrix, _, _ = load_database(DB_SMALL)
    for col in matrix.T:
        np.testing.assert_almost_equal(col.sum(), 1.0)


def test_load_database_empty_row_raises():
    """A row with all-zero k-mer counts must raise RuntimeWarning."""
    with pytest.raises(RuntimeWarning):
        load_database(DB_EMPTY_ROW)


# ---------------------------------------------------------------------------
# run_nnls
# ---------------------------------------------------------------------------

def test_run_nnls_sums_to_one():
    """NNLS output must be normalised (sum ≈ 1)."""
    matrix, _, _ = load_database(DB_SMALL)
    random.seed(42)
    query = normalise([random.randint(1000, 100000) for _ in range(3)])
    result = run_nnls(matrix, query)
    np.testing.assert_almost_equal(result.sum(), 1.0)


def test_run_nnls_non_negative():
    """All NNLS abundances must be ≥ 0."""
    matrix, _, _ = load_database(DB_SMALL)
    random.seed(7)
    query = normalise([random.randint(1000, 100000) for _ in range(3)])
    result = run_nnls(matrix, query)
    assert (result >= 0).all()


def test_run_nnls_known_seed():
    """Regression: known seed must produce stable output."""
    matrix, _, _ = load_database(DB_SMALL)
    random.seed(1128)
    query = normalise([random.randint(10000, 200000) for _ in range(3)])
    result = run_nnls(matrix, query)
    np.testing.assert_almost_equal(result[0], 0.11743935706399153, decimal=6)
    np.testing.assert_almost_equal(result[1], 0.88256064293600844, decimal=6)


# ---------------------------------------------------------------------------
# aggregate_level
# ---------------------------------------------------------------------------

def test_aggregate_level_kingdom():
    """All organisms are Bacteria — kingdom sum should equal 1."""
    matrix, organisms, _ = load_database(DB_SMALL)
    random.seed(500)
    query = normalise([random.randint(10000, 200000) for _ in range(3)])
    abundances = run_nnls(matrix, query)
    results = {org: float(abundances[i]) for i, org in enumerate(organisms)}
    kingdom = aggregate_level(results, 0)
    np.testing.assert_almost_equal(kingdom["Bacteria"], 1.0)


def test_aggregate_level_phylum():
    """Phylum-level aggregation should split into two distinct phyla."""
    matrix, organisms, _ = load_database(DB_SMALL)
    random.seed(500)
    query = normalise([random.randint(10000, 200000) for _ in range(3)])
    abundances = run_nnls(matrix, query)
    results = {org: float(abundances[i]) for i, org in enumerate(organisms)}
    phylum = aggregate_level(results, 1)
    assert set(phylum.keys()) == {"Spirochaetes", "Firmicutes"}
    np.testing.assert_almost_equal(sum(phylum.values()), 1.0)


# ---------------------------------------------------------------------------
# write_results
# ---------------------------------------------------------------------------

def test_write_results_creates_file():
    """write_results must create a CSV with correct header and non-zero rows."""
    matrix, organisms, _ = load_database(DB_SMALL)
    random.seed(1)
    query = normalise([random.randint(10000, 200000) for _ in range(3)])
    abundances = run_nnls(matrix, query)
    results = {org: [float(abundances[i])] for i, org in enumerate(organisms)}

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as tmp:
        tmp_path = Path(tmp.name)

    write_results(results, tmp_path, ["sample.fasta"], ["Kingdom", "Phylum", "Class",
                                                         "Order", "Family", "Genus",
                                                         "Species", "Strain"])
    with open(tmp_path) as fh:
        rows = list(csv.reader(fh))

    assert rows[0][-1] == "sample.fasta"
    assert len(rows) > 1  # at least one data row
    tmp_path.unlink()


# ---------------------------------------------------------------------------
# refine_results
# ---------------------------------------------------------------------------

def test_refine_results_structure():
    """refine_results must return header + data rows with strain-level names."""
    matrix, organisms, _ = load_database(DB_SMALL)
    random.seed(3)
    query = normalise([random.randint(10000, 200000) for _ in range(3)])
    abundances = run_nnls(matrix, query)
    results = {org: [float(abundances[i])] for i, org in enumerate(organisms)}
    levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
    rows = refine_results(results, ["sample.fasta"], levels)
    assert rows[0] == ["Strain", "sample.fasta"]
    assert len(rows) >= 2


# ---------------------------------------------------------------------------
# count_kmers (requires Jellyfish)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not which("jellyfish"), reason="Jellyfish not installed")
def test_count_kmers_fasta():
    """count_kmers on mock_sample.fasta must return expected counts."""
    kmer_order = ["AAAAAA", "AAAAAT", "TTTTTT"]
    result = count_kmers(FASTA, "6", "1", kmer_order)
    assert result == [19, 3, 0]


@pytest.mark.skipif(not which("jellyfish"), reason="Jellyfish not installed")
def test_count_kmers_fastq():
    """count_kmers on mock_sample.fastq must return a list of ints."""
    kmer_order = ["AAAAAA", "AAAAAT", "TTTTTT"]
    result = count_kmers(FASTQ, "6", "1", kmer_order)
    assert isinstance(result, list)
    assert all(isinstance(v, int) for v in result)


@pytest.mark.skipif(not which("jellyfish"), reason="Jellyfish not installed")
def test_count_kmers_empty_raises():
    """count_kmers on an empty file must raise an Exception."""
    with pytest.raises(Exception):
        count_kmers(FASTA_EMPTY, "6", "1", ["AAAAAA"])
