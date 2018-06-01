# Test data set

This is a test data set provided by Rob. The run in the fasta file comes from the
[IBD Multi 'Omics project](https://ibdmdb.org/) and you can find the original on the SRA as project
[SRX3105022](https://www.ncbi.nlm.nih.gov/sra/SRX3105022).
I chose this because it is a human gut data set and should have lots of well characterized bacteria.

To run the test use:

```
python focus.py test/MSM5LLDU.fasta
```

Here is the output table for species level that you should get, and the full output file is [also included](MSM5LLDU.fasta__output.txt)


### Species Level
| Rank | Predicted Organism |   Estimated Abundance (%)
|---|-----------------------|---------------------------------------------------
| 1 |    Bacteroides_helcogenes |  19.1211636317
| 2 |       Bacteroides_xylanisolvens |       15.5031066261
| 3 |       Bacteroides_salanitronis |        11.3255480363
| 4 |       Bacteroides_thetaiotaomicron |    8.36145300543
| 5 |       Bacteroides_vulgatus |    6.59710893868
| 6 |       Prevotella_denticola |    4.42671781981
| 7 |       Bacteroides_fragilis |    4.27881060106
| 8 |       Eubacterium_siraeum |     3.42323728348
| 9 |       Parabacteroides_distasonis |      3.00810957398
| 10 |      Alistipes_finegoldii |    2.51926524262
| 11 |      Chitinophaga_pinensis |   2.04409920634
| 12 |      Mahella_australiensis |   1.91722582299
| 13 |      Cenarchaeum_symbiosum |   1.89285520872
| 14 |      Prevotella_intermedia |   1.69865846183
| 15 |      Faecalibacterium_prausnitzii |    1.66073709429
| 16 |      Thermocrinis_albus |      1.48830645646
| 17 |      Desulfurococcus_fermentans |      1.30994756082
| 18 |      Treponema_brennaborense |         1.25593593996
| 19 |      Anaplasma_centrale |      1.12299034792
| 20 |      Spirochaeta_africana |    1.06189465851
| - |       Others (abundance < 1%) |         5.98


# Testing the stamp output

The directory test_stamp has several fasta files from the same series as described above ([SRX3105022](https://www.ncbi.nlm.nih.gov/sra/SRX3105022)). You can use this dataset to output tab-separated files that you can open in Libre Office, Excel, or STAMP.

To run that test use:

```
python focus.py -q test/test_stamp
```


