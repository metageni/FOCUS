# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import os
import random
import sys

from pathlib import Path

from numpy import array
from numpy import sum as numpy_sum
from scipy.optimize import nnls


LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'
logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)
LOGGER = logging.getLogger(__name__)

WORK_DIRECTORY = 'focus_app'


def normalise(raw_counts):
    """Normalise raw counts into proportions.

    Args:
        raw_counts (numpy.ndarray): array with raw count

    Returns:
        numpy.ndarray: Normalised data

    """
    return raw_counts / (numpy_sum(raw_counts) * 1.)


def load_database(database_path):
    """Load database.

    Args:
        database_path (str): Path to database

    Returns:
        numpy.ndarray: matrix with loaded database
        list: list of organisms in the database

    """
    database_results = {}
    with open(database_path) as database_file:
        database_reader = csv.reader(database_file, delimiter='\t')
        next(database_reader, None)
        for row in database_reader:
            if '0' in row[9] and numpy_sum(array(row[8:], dtype='i')) == 0:
                sys.stderr.write("There are no kmers found for " + "\t".join(row[:8]))
                continue
            database_results["\t".join(row[:8])] = normalise(array(row[8:], dtype='i'))

    organisms = list(database_results.keys())
    database_results = array([database_results[organism] for organism in organisms])

    return database_results.T, organisms


def count_kmers(kmer_size, method="jellyfish"):
    """Count k-mers on FAST(A/Q) file.

    Args:
        kmer_size (str): k-mer size
        method (str): software to count k-mers

    Returns:
        numpy.ndarray: raw count of k-mers

    """
    pass


def count_kmers_fake(kmer_size, method="jellyfish"):
    """Count k-mers on FAST(A/Q) file.

    Args:
        kmer_size (str): k-mer size
        method (str): software to count k-mers

    Returns:
        numpy.ndarray: raw count of k-mers

    """
    return normalise([random.randint(10000, 200000) for _ in range(2080)])


def write_results(results, output_directory, query_files):
    """Write FOCUS results.

     Args:
         results (dict): profile for all the metagenomes
         output_directory (str): Path to output directory
         query_files (list): List with list of files profiled

     """
    with open("{}/STAMP_tabular.xls".format(output_directory), 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        writer.writerow(["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"] +
                        [targeted_file for targeted_file in query_files])

        for taxa in results:
            if sum(results[taxa]) > 0:
                writer.writerow(taxa.split("\t") + [abundance * 100 for abundance in results[taxa]])


def main():
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-q", "--query", help="Path to FAST(A/Q) file or directory with these files", required=True)
    #parser.add_argument("-o", "--output_directory",  help="Path to output files", required=True)
    #parser.add_argument("-k", "--kmer_size",  help="K-mer size (6 or 7)", default="6")
    #parser.add_argument("-d", "--work_directory",  help="Work directory", default="focus_app")
    #parser.add_argument("-b", "--alternate_directory",  help="Alternate directory for your databases", default="")

    #args = parser.parse_args()

    #query = Path(args.query)
    #output_directory = Path(args.output_directory)
    #kmer_size = args.kmer_size
    #WORK_DIRECTORY = Path(args.alternate_directory) if args.alternate_directory else Path(args.work_directory)

    query = Path("../queries/")
    output_directory = Path("../output/")
    kmer_size = "6"
    WORK_DIRECTORY = Path("")

    # check if query is exists
    if not query.exists():
        LOGGER.critical("QUERY: {} does not exist".format(query))

    # check if output_directory is exists
    elif not output_directory.exists():
        LOGGER.critical("OUTPUT: {} does not exist".format(output_directory))

    # check if work directory exists
    elif WORK_DIRECTORY != WORK_DIRECTORY or not output_directory.exists():
        LOGGER.critical("WORK_DIRECTORY: {} does not exist".format(WORK_DIRECTORY))

    # check k-mer size
    elif kmer_size not in ["6", "7"]:
        LOGGER.critical("K-MER SIZE: {} is not a valid k-mer size for this program - please choose 6 or 7".format(kmer_size))

    else:
        LOGGER.info("FOCUS: An Agile Profiler for Metagenomic Data")
        LOGGER.info("   1) Loading Reference DB")
        database_path = Path(WORK_DIRECTORY, "db/k" + kmer_size)

        if not database_path.exists ():
            sys.stderr.write ("ERROR: Database {} was not found. Did you extract db.zip?".format(database_path))
            sys.exit()

        database_matrix, organisms = load_database(database_path)

        LOGGER.info("   2) Reference DB was loaded with {} reference genomes".format(len(organisms)))

        # get fasta/fastq files
        query_files = [query] if query.is_file() else [temp_query for temp_query in os.listdir(query)]
        query_files.sort()

        results = {taxa:[0] * len(query_files) for taxa in organisms}

        counter = 1
        for temp_query in query_files:
            LOGGER.info("   3.{}) Working on: {}".format(counter, temp_query))

            LOGGER.info("      Counting k-mers")

            #query_count = count_kmers(kmer_size)
            query_count = count_kmers_fake(kmer_size)
            # find the best set of organisms that reconstruct the user metagenome using NNLS
            LOGGER.info("      Running FOCUS")
            organisms_abundance = normalise(nnls(database_matrix, query_count)[0])

            # store results
            query_index = query_files.index(temp_query)
            for pos in range(len(organisms)):
                results[organisms[pos]][query_index] = organisms_abundance[pos]

            counter += 1

        LOGGER.info('   5) Writing Results to {}'.format(output_directory))
        write_results(results, output_directory, query_files)

    LOGGER.info('Done')


if __name__ == "__main__":
    main()