#FOCUS
An Alignment-Free Model To Identify Organisms In Metagenomes Using Non-Negative Least Squares 
version: 0.3
-----

(c)            Silva, G. G. Z., D. A. Cuevas, B. E. Dutilh, and R. A. Edwards, 2014: FOCUS: an alignment-free model
			   to identify organisms in metagenomes using non-negative least squares. PeerJ, 2, e425,doi:10.7717/peerj.425.
website: 	   http://edwards.sdsu.edu/FOCUS


(1) USAGE
-----
python focus.py -q query_sequence.fna [-k] [-m] [-s]

	-q Specify input file. 
		Required: Input files should be in FASTA/FASTQ format. 

	-k Specify k-mer frequency used on the profile (default: 7)
		6,7 and 8 frequencies are available. 
	
	-m minimum relative abundance to show in the results (default: 1%)
	
	-s output in the STAMP format for multiple files [0=False (default) and 1=True]
	
	-l Split STAMP output in different levels (default: all; options: kingdom, phylum, class, order, family, genus, or species)

	
(2) ADDING DATA INTO THE DATABASE
-----
In order to insert data into the FOCUS database, you have to run focus only with -d parameter
	For example, python focus.py -d myGenomicData

	The genomic files has to be in the following format having 9 tabulated(\t) fields:
	 1. File location: You fasta/fastq file location
	 2-9.For your organisms: Kingdom, phylum, class, order, family, genus, species, and strain
	 
(3) OUTPUT

Single file:
---------------------------------------------------------------------------------------------------------------
- FOCUS prints the output after it finishes the running, and it also write a tabular file (YOURQUERY__output.txt)

Multiple files:
---------------------------------------------------------------------------------------------------------------
- FOCUS also run the program and generate a tabular output comparing all the samples. This output can be used
as STAMP's input for statistical analysis.

- Split STAMP output in different levels by setting -l to the wanted level 
(default: all; options: kingdom, phylum, class, order, family, genus, or species)
-----
For further info please visit the FOCUS website at http://edwards.sdsu.edu/FOCUS


(4) DEPENDENCIES
------------
Python >= 2.7.X and 3.0: http://www.python.org/download
Jellyfish: http://www.cbcb.umd.edu/software/jellyfish
Numpy: http://sourceforge.net/projects/numpy/files/NumPy
SciPy: http://sourceforge.net/projects/scipy

COPYRIGHT AND LICENSE
---------------------
Copyright (C) 2013-2015  Genivaldo Gueiros Z. Silva

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
