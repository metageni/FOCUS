# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import os
import random

from pathlib import Path
from shutil import which
from collections import defaultdict

from focus_app import version

from numpy import array
from numpy import sum as numpy_sum
from scipy.optimize import nnls

LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'


def normalise(raw_counts):
    """Normalise raw counts into proportions.

    Args:
        raw_counts (numpy.ndarray): Array with raw count.

    Returns:
        numpy.ndarray: Normalised data.

    """
    sum_values = numpy_sum(raw_counts)

    if not sum_values:
        raise RuntimeWarning('All values in input are 0.')
    return raw_counts / sum_values


def is_wanted_file(queries):
    """Remove files from query files that not have extension .fasta/.fastq/.fna

    Args:
        queries (list): List with query names.

    Returns:
        list: Sorted list with only .fasta/.fastq/.fna files.

    """
    queries = [query for query in queries if query.split(".")[-1].lower() in ["fna", "fasta", "fastq"]]
    queries.sort()

    return queries


def load_database(database_path):
    """Load database.

    Args:
        database_path (PosixPath): Path to database.

    Returns:
        numpy.ndarray: Matrix with loaded database.
        list: List of organisms in the database.
        list: K-mer database order.

    """
    database_results = {}
    with open(database_path) as database_file:
        database_reader = csv.reader(database_file, delimiter='\t')
        kmer_order = next(database_reader, None)[8:]
        for row in database_reader:
            # if k-mers are all 0, normalise will raise an expection
            database_results["\t".join(row[:8])] = normalise(array(row[8:], dtype='i'))

    organisms = list(database_results.keys())
    database_results = array([database_results[organism] for organism in organisms])

    return database_results.T, organisms, kmer_order


def count_kmers(query_file, kmer_size, threads, kmer_order):
    """Count k-mers on FAST(A/Q) file.

    Args:
        query_file (PosixPath): Query in FAST(A/Q) file.
        kmer_size (str): K-mer size.
        threads (str): Number of threads to use in the k-mer counting.
        kmer_order (list): List with k-mers database order.

    Returns:
        numpy.ndarray: K-mer counts.

    """
    suffix = str(random.random())
    output_count = Path("kmer_counting_{}".format(suffix))
    output_dump = Path("kmer_dump_{}".format(suffix))

    # count and dump kmers counts
    os.system("jellyfish count -m {} -o {} -s 100M -t {} -C --disk {}".format(kmer_size, output_count, threads,
                                                                              query_file))
    os.system("jellyfish dump {} -c > {}".format(output_count, output_dump))

    # delete binary counts
    if output_count.exists():
        os.remove(output_count)
    else:
        raise Exception('Jellyfish failed to count the k-mers. Make sure you have installed its required version')

    # checks if k-mer counting correctly was dumped
    if output_dump.exists():
        # not empty file
        if output_dump.stat().st_size:
            counts = defaultdict(int)
            with open(output_dump) as counts_file:
                counts_reader = csv.reader(counts_file, delimiter=' ')
                for kmer, count in counts_reader:
                    counts[kmer] = int(count)

            # delete dump file
            os.remove(output_dump)

            return [counts[kmer_temp] for kmer_temp in kmer_order]

        else:
            os.remove(output_dump)
            raise Exception('{} has no k-mers count. Probably not valid file'.format(query_file))

    else:
        raise Exception('Something went wrong when trying to dump the k-mer counting.')


def refine_results(results, query_files, taxonomy_level):
    """Result result removing rows with 0 counts.

     Args:
         results (dict): Profile for all the metagenomes.
         query_files (list): List with list of files profiled.
         taxonomy_level (list): Taxonomy level(s).

    Returns:
        list of list: List with refined results.

     """
    refined_results = [[taxonomy_level[-1:] + [Path(targeted_file).parts[-1] for targeted_file in query_files]]]
    for taxa in results:
        if sum(results[taxa]) > 0:
            refined_results.append(taxa.split("\t")[-1:] + [abundance * 100 for abundance in results[taxa]])

    return refined_results


def write_results(results, output_directory, query_files, taxonomy_level):
    """Write FOCUS results.

     Args:
         results (dict): Profile for all the metagenomes.
         output_directory (PosixPath): Path to output file.
         query_files (list): List with list of files profiled.
         taxonomy_level (list): Taxonomy level(s).

     """
    with open(output_directory, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
        writer.writerow(taxonomy_level + [Path(targeted_file).parts[-1] for targeted_file in query_files])

        for taxa in results:
            if sum(results[taxa]) > 0:
                writer.writerow(taxa.split("\t") + [abundance * 100 for abundance in results[taxa]])


def aggregate_level(results, position):
    """Aggregate abundance of metagenomes by taxonomic level.

    Args:
        results (dict): Path to results.
        position (int): Position of level in the results.

    Returns:
        dict: Aggregated result targeting chosen level.

    """
    level_results = defaultdict(list)

    for all_taxa in results:
        taxa = all_taxa.split("\t")[position]
        abundance = results[all_taxa]
        level_results[taxa].append(abundance)

    return {temp_taxa: numpy_sum(level_results[temp_taxa], axis=0) for temp_taxa in level_results}


def run_nnls(database_matrix, query_count):
    """Run Non-negative least squares (NNLS) algorithm.

    Args:
        database_matrix (ndarray): Matrix with count for organisms in the database.
        query_count (ndarray): Metagenome k-mer count.

    Returns:
        numpy.ndarray: Abundances of each organism.

    """
    return normalise(nnls(database_matrix, query_count)[0])


def parse_args():
    """Parse args entered by the user.

    Returns:
        argparse.Namespace: Parsed arguments.

    """
    parser = argparse.ArgumentParser(description="FOCUS: An Agile Profiler for Metagenomic Data",
                                     epilog="example > focus -q samples_directory")
    parser.add_argument('-v', '--version', action='version', version='FOCUS {}'.format(version))
    parser.add_argument("-q", "--query", help="Path to directory with FAST(A/Q) files", required=True)
    parser.add_argument("-o", "--output_directory", help="Path to output files", required=True)
    parser.add_argument("-k", "--kmer_size", help="K-mer size (6 or 7) (Default: 6)", default="6")
    parser.add_argument("-b", "--alternate_directory", help="Alternate directory for your databases", default="")
    parser.add_argument("-p", "--output_prefix", help="Output prefix (Default: output)", default="output")
    parser.add_argument("-t", "--threads", help="Number Threads used in the k-mer counting (Default: 4)", default="4")
    parser.add_argument('--list_output', help='Output results as a list',
                        action='store_true', required=False)
    parser.add_argument('-l', '--log', help='Path to log file (Default: STDOUT).', required=False)

    return parser.parse_args()


def main(args=False):

    if not args:
        args = parse_args()

    # parameters and other variables
    query = Path(args.query)
    prefix = args.output_prefix
    output_directory = Path(args.output_directory)
    kmer_size = args.kmer_size
    work_directory = Path(args.alternate_directory) if args.alternate_directory else Path(__file__).parents[0]
    database_path = Path(work_directory, "db/k" + kmer_size)
    threads = args.threads
    jellyfish_path = which("jellyfish")
    jellyfish_version = os.popen("jellyfish count --version").read().split(".")[0] if jellyfish_path else None

    if args.log:
        logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO, filename=args.log)
    else:
        logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)

    logger = logging.getLogger(__name__)

    logger.info("FOCUS: An Agile Profiler for Metagenomic Data")

    # check if output_directory is exists - if not, creates it
    if not output_directory.exists():
        Path(output_directory).mkdir(parents=True, mode=511)
        logger.info("OUTPUT: {} does not exist - just created it :)".format(output_directory))

    # check if at least one of the queries is valid
    if not query.is_dir():
        logger.critical("QUERY: {} is not a directory".format(query))

    # check if database exists
    if not database_path.exists():
        logger.critical("DATABASE: {} does not exist. Did you extract db.zip?".format(database_path))
        compressed_db = Path(work_directory, "db.zip")
        uncompress_path = Path(work_directory)

        # check if unzip is installed
        if not which("unzip"):
            logger.critical("Install unzip !!!")
        # try to uncompress database for you
        elif compressed_db.exists():
            logger.info("DATABASE: Uncompressing Database for you :)")
            os.system("unzip {} -d {}".format(str(compressed_db), uncompress_path))
        else:
            logger.critical("{} was not found".format(compressed_db))

    # check if at least one of the queries is valid
    if is_wanted_file(os.listdir(query)) == []:
        logger.critical("QUERY: {} does not have any Fasta/Fna/Fastq file".format(query))

    # check if k-mer counter is installed
    elif not jellyfish_path:
        logger.critical("K-MER COUNTER: Jellyfish is not installed. Please install 2.XX")

    # check jellyfish installed is the correct version
    elif jellyfish_version != "2":
        logger.critical("K-MER COUNTER: Jellyfish needs to be version 2.XX. You have version {}".
                        format(jellyfish_version))

    # check if query is exists
    elif not query.exists():
        logger.critical("QUERY: {} does not exist".format(query))

    # check if work directory exists
    elif work_directory != work_directory or not work_directory.exists():
        logger.critical("WORK_DIRECTORY: {} does not exist".format(work_directory))

    # check k-mer size
    elif kmer_size not in ["6", "7"]:
        logger.critical("K-MER SIZE: {} is not a valid k-mer size for this program - "
                        "please choose 6 or 7".format(kmer_size))

    else:
        logger.info("1) Loading Reference DB")
        database_path = Path(work_directory, "db/k" + kmer_size)
        database_matrix, organisms, kmer_order = load_database(database_path)

        logger.info("2) Reference DB was loaded with {} reference genomes".format(len(organisms)))
        # get fasta/fastq files
        query_files = is_wanted_file([temp_query for temp_query in os.listdir(query)])

        results = {taxa: [0] * len(query_files) for taxa in organisms}

        for counter, temp_query in enumerate(query_files):
            logger.info("3.{}) Working on: {}".format(counter + 1, temp_query))

            logger.info("   Counting k-mers")
            # count k-mers
            query_count = normalise(count_kmers(Path(query, temp_query), kmer_size, threads, kmer_order))
            # find the best set of organisms that reconstruct the user metagenome using NNLS
            logger.info("   Running FOCUS")
            organisms_abundance = run_nnls(database_matrix, query_count)

            # store results
            query_index = query_files.index(temp_query)
            for pos in range(len(organisms)):
                results[organisms[pos]][query_index] = organisms_abundance[pos]

        taxomy_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]

        logger.info('5) Writing Results to {}'.format(output_directory))
        # All taxonomy levels in one output
        output_file = Path(output_directory, prefix + "_All_levels.xls")
        write_results(results, output_file, query_files, taxomy_levels)

        # write output for each taxonomy level
        for pos, level in enumerate(taxomy_levels):
            logger.info('  5.{}) Working on {}'.format(pos + 1, level))
            level_result = aggregate_level(results, pos)
            output_file = Path(output_directory, prefix + "_" + level + "_tabular.xls")
            write_results(level_result, output_file, query_files, [level])

        if args.list_output:
            logger.info('6) Creating list of lists with output')
            return refine_results(results, query_files, taxomy_levels)

    logger.info('Done'.format(output_directory))


if __name__ == "__main__":
    main()
