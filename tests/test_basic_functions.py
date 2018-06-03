# -*- coding: utf-8 -*-

from numpy import array

from focus_app.focus import normalise, load_database, is_wanted_file, count_kmers


def test_normalise():
    assert list(normalise(array([1, 1,]))) == [0.5, 0.5]
    assert list(normalise(array([2, 2, 2, 2]))) == [0.25, 0.25, 0.25, 0.25]


def test_load_database():
    database_path = "tests/data/k6_subset.txt"
    database_matrix, organisms, kmers_order = load_database(database_path)

    assert list(organisms) == ['Bacteria\tSpirochaetes\tSpirochaetia\tSpirochaetales\tSpirochaetaceae\tSphaerochaeta\tSphaerochaeta_pleomorpha\tSphaerochaeta_pleomorpha_Grapes_uid82365',
                               'Bacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\tStreptococcus\tStreptococcus_suis\tStreptococcus_suis_98HAH33_uid58665',
                               'Bacteria\tProteobacteria\tGammaproteobacteria\tAlteromonadales\tAlteromonadaceae\tAlteromonas\tAlteromonas_macleodii\tAlteromonas_macleodii__Deep_ecotype__uid58251'
                               ]
    # look into small subset of the matrix
    assert list(database_matrix[0])[0] == 0.00027570089293671023
    assert list(database_matrix[1])[0] == 0.00040074099488477375
    assert list(database_matrix[2])[0] == 0.00036760119058228027


def test_is_wanted_file():
    assert is_wanted_file(["a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna', 'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png", "a.fasta", "b.fastq", "x.FASTq", "y.FASTA", "n.fna"]) == ['a.fasta', 'b.fastq', 'n.fna', 'x.FASTq', 'y.FASTA']
    assert is_wanted_file(["f.png" , "queries/"]) == []


def test_count_kmers():
    pass


def test_write_results():
    pass
