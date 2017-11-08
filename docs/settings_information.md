![Logo](../img/RBNS_logo.png)

***
# The settings file is a .json file that contains the required information for an RBNS experiments (e.g., file and output paths, protein name and concentrations, FASTQ barcodes), and it also tells RBNS_main.py which functionalities to perform.

### REQUIRED information for each RBNS experiment

The recommended use for SHAPEware-MaP is to first create a directory that contains the input files in standard locations.
The analysis results will be placed in the same directory under specific file and folder paths. You can override this 
behavior to some extent by giving additional arguments. Run "./analyze_shape-map.py -h" for more details. 

For now, let's assume you created a directory called **my_shape_data/**

#### REQUIRED parameter specifications

Specifying each of the following in the .json file is required; if left out, an error will be immediately returned.


#### OPTIONAL parameter specificiations

Specifying any of the following in the .json file is optional; if left out, they will simply not be computed.


### Running RBNS_main.py

Once these files are in place, the RBNS_pipeline can be run in place with the following commands (with paths adjusted for your system):
	
	source activate rbns_pipeline
	python src/RBNS_main.py test_data/settings.RBFOX3.json
