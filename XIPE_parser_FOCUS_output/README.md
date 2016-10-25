xipe_comparison.py 0.1
----------------------------------------------
This program uses FOCUS output (http://edwards.sdsu.edu/FOCUS) into XIPE (http://edwards.sdsu.edu/cgi-bin/xipe.cgi)
For a a non-parametric statistical analysis of the distribution of samples to determine which samples are statistically
significantly different.
The script parses FOCUS output and compare the each pair of samples and reports which targeted level (Genus for example)
were statistically different in a confidence level


xipe_comparison.py: Uses XIPE for a non-parametric statistical comparision of the samples in the FOCUS ouput
	-q FOCUS output (*__STAMP_tabular.spf
	
	-c Minimum Confidence Level (Default: 95)
	
	-l Comparison Level [kingdom, phylum, class, order, family, genus, species, strain, all](Default: genus)
	
	-o Output Name (Default: my_xipe_output.xls)
     
	 example> python xipe_comparison.py -q input__STAMP_tabular.spf.xls -c 95 -l genus -o xipe__genus___sharks_ouput