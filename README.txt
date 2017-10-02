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


