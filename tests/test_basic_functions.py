# -*- coding: utf-8 -*-

import random

from focus_app.focus import normalise, load_database, is_wanted_file, count_kmers, which, run_nnls, aggregate_level

import pytest

from numpy import array


def test_normalise():
    assert list(normalise(array([1, 1,]))) == [0.5, 0.5]
    assert list(normalise(array([2, 2, 2, 2]))) == [0.25, 0.25, 0.25, 0.25]

    # bad input
    with pytest.raises (RuntimeWarning):
        normalise([0, 0, 0])


def test_is_wanted_file():
    assert is_wanted_file(["a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna', 'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png", "a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna', 'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png" , "queries/"]) == []


def test_count_kmers():
    if which("jellyfish"):
        query_file = "tests/data/mock_sample.fasta"
        threads = "1"

        kmer_size = "6"
        kmer_order = ["AAAAAA", "AAAAAT", "TTTTTT"]
        assert count_kmers(query_file, kmer_size, threads, kmer_order) == [19, 3, 0]

        kmer_size = "7"
        kmer_order = ["AAAAAAA", "AAAGAAT", "ATTTTTT"]
        assert count_kmers(query_file, kmer_size, threads, kmer_order) == [17, 0, 0]

        # testing for empty input
        query_file = "tests/data/mock_sample_empty.fasta"
        kmer_size = "6"
        kmer_order = ["AAAAAA", "AAAAAT", "TTTTTT"]
        with pytest.raises(Exception):
            count_kmers(query_file, kmer_size, threads, kmer_order)


def test_load_database():
    # simple database
    database_matrix, organisms, kmer_order = load_database("tests/data/k6_small_sample.txt")
    assert kmer_order == ['GAACGC', 'GAACGA', 'CACCCA']
    assert organisms == ['Bacteria\tSpirochaetes\tSpirochaetia\tSpirochaetales\tSpirochaetaceae\tSphaerochaeta\tSphaerochaeta_pleomorpha\tSphaerochaeta_pleomorpha_Grapes_uid82365',
                  'Bacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\tStreptococcus\tStreptococcus_suis\tStreptococcus_suis_98HAH33_uid58665']
    assert [list(xx) for xx in database_matrix] == [[0.26407041877834087, 0.27930682976554538],
                                                    [0.38383568951720459,0.38379204892966362],
                                                    [0.35209389170445454, 0.33690112130479105]]

    # bad database with all keys 0 in row = raise error
    with pytest.raises (RuntimeWarning):
        load_database("tests/data/k6_small_sample_empty_record.txt")


def test_run_nnls():
    database_matrix, organisms, kmer_order = load_database("tests/data/k6_small_sample.txt")

    # fake metagenomic count
    # example 1
    random.seed(1128)
    fake_query_count = normalise([random.randint(10000, 200000) for _ in range (3)])
    assert list(run_nnls(database_matrix, fake_query_count)) == [0.11743935706399153, 0.88256064293600844]

    # example 2
    random.seed(2)
    fake_query_count = normalise([random.randint(10000, 200000) for _ in range (3)])
    assert list(run_nnls(database_matrix, fake_query_count)) == [0.79079139795692821, 0.20920860204307179]

    # example 3
    random.seed(500)
    fake_query_count = normalise([random.randint(10000, 200000) for _ in range (3)])
    assert list(run_nnls(database_matrix, fake_query_count)) == [0.40650798820847034, 0.59349201179152955]


def test_aggregate_level():

    true_answer = [
        {'Bacteria': [0.54337724854388625, 0.76762773844686816,
                      0.68899501300924559]},
        {'Spirochaetes': [0.26407041877834087, 0.38383568951720459,
                          0.35209389170445454], 'Firmicutes': [0.27930682976554538,
                                                               0.38379204892966362, 0.33690112130479105]},
        {'Spirochaetia': [0.26407041877834087, 0.38383568951720459,
                          0.35209389170445454], 'Bacilli': [0.27930682976554538,
                                                            0.38379204892966362, 0.33690112130479105]},
        {'Spirochaetales': [0.26407041877834087, 0.38383568951720459,
                            0.35209389170445454], 'Lactobacillales': [0.27930682976554538,
                                                                      0.38379204892966362, 0.33690112130479105]},
        {'Spirochaetaceae': [0.26407041877834087, 0.38383568951720459,
                             0.35209389170445454], 'Streptococcaceae': [0.27930682976554538,
                                                                        0.38379204892966362, 0.33690112130479105]},
        {'Sphaerochaeta': [0.26407041877834087, 0.38383568951720459,
                           0.35209389170445454], 'Streptococcus': [0.27930682976554538,
                                                                   0.38379204892966362, 0.33690112130479105]},
        {'Sphaerochaeta_pleomorpha': [0.26407041877834087,
                                      0.38383568951720459, 0.35209389170445454],
         'Streptococcus_suis': [0.27930682976554538, 0.38379204892966362,
                                0.33690112130479105]},
        {'Sphaerochaeta_pleomorpha_Grapes_uid82365': [0.26407041877834087,
                                                      0.38383568951720459, 0.35209389170445454],
         'Streptococcus_suis_98HAH33_uid58665': [0.27930682976554538,
                                                 0.38379204892966362, 0.33690112130479105]},
    ]
    database_matrix, organisms, kmer_order = load_database ("tests/data/k6_small_sample.txt")
    results = {organisms[pos]: profile for pos, profile in enumerate(database_matrix.T)}

    producted_result = []

    # runs aggregate function in all the levels
    for pos in range(8):
        aggregate_results = aggregate_level(results, pos)

        # convert results to list rather the numpy.array
        aggregate_results_list = {level:list(aggregate_results[level]) for level in aggregate_results}

        producted_result.append(aggregate_results_list)

    assert producted_result == true_answer


def test_write_results():
    pass
