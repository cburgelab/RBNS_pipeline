![Logo](../img/RBNS_logo.png)

***
# Input files for the RBNS_pipeline:

### FASTQ file

The FASTQ file must contain reads from at least two libraries: and input library, and one or more pulldown libraries, each of which will be compared to the input library. Each library must have its own sequencing barcode (for example, GGCTAC and TCGGCA as in the two reads below).

	@HWI-ST:7:1101:6120:4119#GGCTAC/1
	CCGTAACTTACTGAACAGCATGGAATTCTCGGGTGCCAAG
	+HWI-ST:7:1101:6120:4119#GGCTAC/1
	BBCDBDFFHHHHHJJJJIJJJJJHJIJJJJJJJHHHIGJJ
	@HWI-ST:7:1101:2531:2734#TCGGCA/1
	TTTACATTTAAATTGGCCTCTGGAATTCTCGGGTGCCAAG
	+HWI-ST:7:1101:2531:2734#TCGGCA/1
	CCCFFFFFHHHHHIIJIJJJJJJJJJJJJJJJJIJJJJJJ


### The settings.json file

This file contains all of the settings for the RBNS experiment, and it allows you to specify what should be calculated in a particular run of RBNS_main.py (e.g., 'naive' enrichment calculations on which k's, SKA library fractions, enrichments by read position, RBNS logo generation, RBNS folding, etc.).

An example minimal settings file is provided as test_data/settings.RBFOX3.json. See [here](settings_information.md) for a more detailed description of the required and optional arguments within the settings.json file.


## Running RBNS_main.py

Once the .settings file & FASTQ are in place, the RBNS_pipeline can be run in place with the following commands:
	
	source activate rnbs_pipeline (should have already been done upon installation)
	python src/RBNS_main.py test_data/settings.RBFOX3.json (or wherever your future settings.json file is)
