- RBNS pipeline built Peter Freese (freese.peter@gmail.com)

- Basic functionality includes taking in a .fastq file and:

	- Splitting the reads into different concentrations based on the sequencing barcode
	- Calculating QC metrics (reads at expected length, library complexity)
	- Computing kmer counts and frequencies for specified k's (typically k=4-7), including:

		- Simple enrichments (R=ratio of frequency in each PD concentration/input)
		- Streaming kmer assignment (SKA) library fractions and B values (see Lambert et al. 2014, Mol. Cell )





- Optional additional functionalities include:
	- Producing motif sequence logos from alignment of weighted enriched kmers (requires the weblogo command line program)
	- Producing C+G matched files at each PD concentration that can be fed into RNAfold for secondary structure analysis
	- Producing .bed files of enriched kmer occurrences in the hg19 transcriptome that can be loaded into the UCSC Genome Browser (requires download of additional large files of precomputed kmer occurrences in the transcriptome)


- Package is developed for Python 2.7* on a Linux high performance computing cluster

	- Dependencies:
		- scipy
		- numpy
		- matplotlib
		- simplejson

	- Optional dependencies for additional functionalities:
		- WebLogo3 command line executable on your $PATH (tested with WebLogo 3.3: https://pypi.python.org/pypi/weblogo/3.3)
		- RNAfold command line executable (tested with RNAfold 2.1.6: https://www.tbi.univie.ac.at/RNA/index.html)




- To run with provided RBFOX3 test data:


