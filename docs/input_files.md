![Logo](../img/SHAPEware-logo.png)

***
# Input files for SHAPEware

### SHAPEware-MaP (analyze_shape-map.py)

The recommended use for SHAPEware-MaP is to first create a directory that contains the input files in standard locations.
The analysis results will be placed in the same directory under specific file and folder paths. You can override this 
behavior to some extent by giving additional arguments. Run "./analyze_shape-map.py -h" for more details. 

For now, let's assume you created a directory called **my_shape_data/**

#### File: my_shape_data/raw_data/

This is the default directory where the FASTQ files should be placed. Currently, SHAPEware-MaP expects paired end reads,
but this requirement will be relaxed in the near future. In general, SHAPEware is designed to analyze batches of experiments that are either
biologically related or were generated during the same sequencing run, but you can also analyze single datasets.

#### File: my_shape_data/input_sheet.csv (or .xlsx file, if given as argument)

This is a spreadsheet that describes how sequence read files (FASTQ files) relate to reference sequence and to the three types 
of SHAPE experiments that are commonly performed on each RNA (one sample with SHAPE reagent, one without, and one
treated with SHAPE reagent but chemically denatured). In the example below, these are columns E-G.

![sheet_example](../img/input_sheet_example.png)

Sample names follow the [Illumina FASTQ naming convention](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm),
although there is some built-in flexibility. However, in general you should expect the values of the three sample columns in the spreadsheet to match the "SampleName" portion of the file.

In other words, if your first sample is named "SHAPE-example", the FASTQ files are expected to adhere to a similar naming scheme as these:
- raw_data/SHAPE-example_S1_L001_R1_001.fastq.gz
- raw_data/SHAPE-example_S1_L001_R2_001.fastq.gz

** Important: ** Sample names must only include alphanumeric characters and dashes. This is because Snakemake relies
on wildcards in filenames. You can also enter the sample numbers (S1, S2, S3) into the spreadsheet instead if this is an issue.

It is important that the name provided in column C ("Reference sequence") exactly matches the reference sequence names provided 
in the reference FASTA file (see below).

#### File: my_shape_data/ref_seqs.fa

A simple FASTA file that contains the sequences of the RNAs assayed in each SHAPE-MaP experiment. They can be provided as either DNA or RNA.

For the TPP riboswitch example, this is the file:

	>TPP
	GGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTC

Each row in the spreadsheet must be matched (by name/header) to one of the reference sequences described in the FASTA file.
Otherwise, this will create an error.

#### File: my_shape_data/ref_masks.fa

A simple FASTA file that contains the sequences of the RNAs assayed in each SHAPE-MaP experiment, with positions that should
be masked (i.e., positions for which SHAPE values will not be used for folding and will be highlighted as masked in the plots)
replaced with an X.

For the TPP riboswitch example, this is the file:

	>TPP
	XXXXXXXXGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTC

Note that the first few bases have been replaced by "X" characters. This is because the experiment was done with a primer
that anneals to that sequence, reducing the apparent mutation rate and leading to non-representative SHAPE values. 
The same approach can be used to indicate bases in 
a sequence for which the SHAPE-MaP data is not reliable due to a high mutation rate or other issues with reverse transcription,
PCR, or sequencing. 

As before, each row in the spreadsheet must be matched (by name/header) to one of the reference sequences described in the FASTA file.
Otherwise, this will create an error.

#### _Optional file:_ my_shape_data/ref_structures.dot

A FASTA-like file that contains the known 2D structure of the corresponding RNA in dot-bracket notation. If this is present,
the software will compute statistics that describe the agreement between the predicted and observed folds.

This is an example file for TPP:

	>TPP
	GGACUCGGGGUGCCCUUCUGCGUGAAGGCUGAGAAAUACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC
	(((((((((.((((.((((...)))))))(...).)..)))).(..(((((..((((......)))).).))))))))))

As before, sequence names must match the FASTA files and the names provided in the spreadsheet.

### Running analyze_shape-map.py

Once these files are in place, SHAPEware-MaP can be run in place with the following commands (with paths adjusted for your system):
	
	source activate shapeware
	shapeware/analyze_shape-map.py my_shape_data/