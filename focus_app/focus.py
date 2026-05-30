#!/usr/bin/env python3
"""Core logic for focus_app — an alignment-free metagenomics profiler.

Uses k-mer counting (via Jellyfish) and non-negative least squares (NNLS)
to estimate the relative abundance of organisms in a metagenome.
"""

import os
import csv
import random
import logging

from pathlib import Path
from collections import defaultdict

import numpy as np

from scipy.optimize import nnls

TAXONOMY_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
VALID_EXTENSIONS = {".fna", ".fasta", ".fastq"}

logger = logging.getLogger(__name__)


def normalise(counts):
    """Normalise a count vector to relative abundances summing to 1.

    Args:
        counts (array-like): Raw count values.

    Returns:
        numpy.ndarray: Normalised array.

    Raises:
        RuntimeWarning: If all counts are zero.
    """
    arr = np.array(counts, dtype=float)
    total = arr.sum()
    if total == 0:
        raise RuntimeWarning("All values in input are 0.")
    return arr / total


def is_wanted_file(paths):
    """Filter a list of paths to only accepted sequence file extensions.

    Args:
        paths (list): List of str or Path objects.

    Returns:
        list[Path]: Sorted list of Path objects with valid extensions.
    """
    result = [Path(p) for p in paths if Path(p).suffix.lower() in VALID_EXTENSIONS]
    result.sort()
    return result


def load_database(database_path):
    """Load a FOCUS k-mer database file into a numpy matrix.

    The file is tab-separated: 8 taxonomy columns followed by k-mer count columns.
    Each organism row is normalised before storage.

    Args:
        database_path (str or Path): Path to the database file.

    Returns:
        tuple:
            numpy.ndarray: Transposed matrix (kmers × organisms).
            list[str]: Organism labels (tab-joined taxonomy strings).
            list[str]: Ordered k-mer strings.
    """
    rows = {}
    with open(database_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        kmer_order = next(reader)[8:]
        for row in reader:
            key = "\t".join(row[:8])
            rows[key] = normalise(np.array(row[8:], dtype=int))

    organisms = list(rows.keys())
    matrix = np.array([rows[o] for o in organisms])
    return matrix.T, organisms, kmer_order


def count_kmers(query_file, kmer_size, threads, kmer_order):
    """Count k-mers in a FASTA/FASTQ file using Jellyfish.

    Args:
        query_file (str or Path): Input sequence file.
        kmer_size (str): K-mer length (e.g. "6" or "7").
        threads (str): Number of threads for Jellyfish.
        kmer_order (list[str]): Ordered list of k-mers to report counts for.

    Returns:
        list[int]: Count for each k-mer in kmer_order.

    Raises:
        Exception: If Jellyfish fails or produces no output.
    """
    suffix = str(random.random())
    count_file = Path("kmer_counting_{}".format(suffix))
    dump_file = Path("kmer_dump_{}".format(suffix))

    os.system("jellyfish count -m {} -o {} -s 100M -t {} -C {}".format(
        kmer_size, count_file, threads, query_file))

    if not count_file.exists():
        raise Exception(
            "Jellyfish failed to count k-mers. Ensure Jellyfish 2.x is installed.")

    os.system("jellyfish dump {} -c > {}".format(count_file, dump_file))
    os.remove(count_file)

    if not dump_file.exists():
        raise Exception("Jellyfish dump step failed.")

    if dump_file.stat().st_size == 0:
        os.remove(dump_file)
        raise Exception("{} produced no k-mer counts — file may be invalid.".format(query_file))

    counts = defaultdict(int)
    with open(dump_file) as fh:
        for kmer, count in csv.reader(fh, delimiter=" "):
            counts[kmer] = int(count)
    os.remove(dump_file)

    return [counts[k] for k in kmer_order]


def run_nnls(database_matrix, query_counts):
    """Solve NNLS to estimate organism abundances from k-mer counts.

    Args:
        database_matrix (numpy.ndarray): k-mers × organisms matrix.
        query_counts (numpy.ndarray): Normalised k-mer counts for the query.

    Returns:
        numpy.ndarray: Normalised abundance per organism.
    """
    raw, _ = nnls(database_matrix, query_counts)
    return normalise(raw)


def aggregate_level(results, position):
    """Aggregate per-organism abundances to a given taxonomy level.

    Args:
        results (dict): Mapping of tab-joined taxonomy string → abundance (float or array).
        position (int): Column index (0–7) corresponding to the taxonomy level.

    Returns:
        dict: Mapping of taxon name → summed abundance.
    """
    grouped = defaultdict(list)
    for taxa_key, abundance in results.items():
        taxon = taxa_key.split("\t")[position]
        grouped[taxon].append(abundance)
    return {t: np.sum(grouped[t], axis=0) for t in grouped}


def write_results(results, output_path, query_files, header):
    """Write profiling results to a CSV file.

    Args:
        results (dict): Mapping of taxonomy string → list of abundances per query.
        output_path (str or Path): Destination CSV file.
        query_files (list): Query file paths (used as column headers).
        header (list[str]): Taxonomy level column names.
    """
    file_names = [Path(f).name for f in query_files]
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header + file_names)
        for taxa, abundances in results.items():
            if np.sum(abundances) > 0:
                row = taxa.split("\t") + [a * 100 for a in np.atleast_1d(abundances)]
                writer.writerow(row)


def refine_results(results, query_files, taxonomy_levels):
    """Return results as a list of lists (for programmatic use / --list_output).

    Args:
        results (dict): Mapping of taxonomy string → list of abundances per query.
        query_files (list): Query file paths.
        taxonomy_levels (list[str]): Full taxonomy level names.

    Returns:
        list[list]: Header row followed by data rows (strain-level only, non-zero).
    """
    file_names = [Path(f).name for f in query_files]
    rows = [[taxonomy_levels[-1]] + file_names]
    for taxa, abundances in results.items():
        if sum(abundances) > 0:
            rows.append([taxa.split("\t")[-1]] + [a * 100 for a in abundances])
    return rows


def profile(query_files, database_path, kmer_size="6", threads="4"):
    """Run the full FOCUS profiling pipeline on a list of query files.

    Args:
        query_files (list[Path]): Validated input sequence files.
        database_path (Path): Path to the k-mer database file.
        kmer_size (str): K-mer size ("6" or "7").
        threads (str): Jellyfish thread count.

    Returns:
        tuple:
            dict: Mapping of organism taxonomy string → list of abundances per query.
            list[str]: Organism labels.
    """
    logger.info("Loading reference database from %s", database_path)
    db_matrix, organisms, kmer_order = load_database(database_path)
    logger.info("Database loaded: %d reference genomes", len(organisms))

    results = {taxa: [0.0] * len(query_files) for taxa in organisms}

    for idx, query in enumerate(query_files):
        logger.info("Processing (%d/%d): %s", idx + 1, len(query_files), query)
        raw_counts = count_kmers(query, kmer_size, threads, kmer_order)
        query_vec = normalise(raw_counts)
        abundances = run_nnls(db_matrix, query_vec)
        for pos, org in enumerate(organisms):
            results[org][idx] = float(abundances[pos])

    return results, organisms
