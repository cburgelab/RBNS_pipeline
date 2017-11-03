
![Logo](img/SHAPEware-logo.png)

***


Welcome to the repository for SHAPEware! 

SHAPEware is a set of bioinformatics tools to analyze data from high-throughput sequencing experiments that chemically
probe RNA to investigate its structure. In developing SHAPEware, our goal is to create a unified package for all types
of SHAPE data that provides highly accurate and customizable analysis capabilities.

This is a quick start guide. For more background information, visit the [wiki](https://bitbucket.org/arrakistx/shapeware/wiki/Home).

## Installation

##### Requirements

SHAPEware is designed to run on Linux and MacOS. In addition, it requires the following software to be pre-installed on your computing environment:

- Python
- R
- The [Miniconda](https://conda.io/miniconda.html) package manager
- The "git" command line tool

If you need help installing any of these tools, see the [detailed documentation](docs/installation.md). When installing dependencies, make sure you
agree with the corresponding licenses of various software tools.

Note that if you run SHAPEware on MacOS, you will need to manually install the [RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) package if you would like to
obtain pseudoknot predictions from the ShapeKnots program.

## Download and install SHAPEware

The easiest way to get the SHAPEware software is to clone our repository. This will ensure you always have access to the
latest version. 

	git clone https://bitbucket.org/arrakistx/shapeware.git

After cloning our repository, run the included installation script: 

	cd shapeware
	./install.sh

This will use the Conda package manager to ensure that 
you have all of the key dependencies. This step will create a stable environment in which to run analysis jobs. 
This approach will keep your results reproducible and will not affect other software that you have installed on your system.
While you can also use SHAPEware without this step by installing all necessary packages manually, it is not recommended.

The installation script will also prompt you to download example files (~40Mb). We recommend you download these files as they 
will allow you to properly test your SHAPEware installation.

If no errors were reported, test the new environment.
	
	source activate shapeware

## Test SHAPEware with example data

In this public version, SHAPEware is able to analyze SHAPE-MaP (Smola et al., 2015) data. You can find example input files in the example_data/
directory within the repository. These were derived from experiments that assayed the TPP riboswitch (Haller et al, 2013) in the presence or 
absence of its cognate ligand. 

If you didn't already download the example files with the installation script, run:

	./example_data/download_files.sh 

Now, try running the SHAPEware-MaP on the example data. By default, this will use all available processors, but this 
can be configured with the --num_cores option.

	source activate shapeware
	./analyze_shape-map.py example_data/

Once the script has finished running, you can find the output from the software in the example_data/results/ folder.

## General use

SHAPEware is designed to analyze batches of samples, such as multiplexed experiments coming from the same sequencing run 
or replicate experiments for the same RNA in the presence and absence of ligands. 

### SHAPEware-MaP

The inputs to SHAPEware-MaP are described in more detail [here](docs/input_files.md). Here is a quick summary:

- A spreadsheet/CSV file describing the samples (input_sheet.csv) 
- A FASTA file specifying the reference sequences corresponding to the constructs used in SHAPE-MaP (ref_seqs.fa)
- A FASTA file specifying sequence masks (with 'X' characters in positions that should be ignored in downstream analysis) (ref_masks.fa)
- (Optional): a file containing reference 2D structures in dot-bracket notation (ref_structures.dot)

You can find examples of all of these files in the example_data/ folder. It is probably easiest to just take a look at these files first. You can 
run SHAPEware-MaP on this example and reproduce the figures below.

#### Example output: TPP

The output includes several types of files that provide useful visualizations. The most direct output is a plot of
the inferred SHAPE reactivities at each position. For a more detailed description of SHAPE reactivity values, see the paper by 
Smola et al., referenced below. 

SHAPEware analyzes ambiguous and unambiguous mutations to identify positions where SHAPE values may be 
over- or underestimated. These uncertainties are reflected by the thin vertical bars in the following type of plot. 
We are actively investigating error models for SHAPE reactivities that consider read depth, mutation ambiguity and
inter-sample variability. We will post updates here and on the Wiki page as we update these models.

![TPP_shape_bar](img/TPP_shape_values.png)

You can also see the SHAPE reactivity values superimposed on various folding predictions (some made using the SHAPE data as constraints, 
some made while blinded to any SHAPE data).

Presently, SHAPEware integrates with the RNAfold program from the [ViennaRNA suite](https://www.tbi.univie.ac.at/RNA/) and the ShapeKnots program from
David Matthews' [RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) package. You can find both types of folds in the output directory.
This example was generated by providing SHAPE values from our example dataset to computationally fold the TPP riboswitch in the absence of TPP ligand. 

![TPP_example](img/TPP_example.png)

Plots are also generated to compare the SHAPE reactivities under the presence or absence of its cognate ligand (or other differential conditions)
and the level of inter-sample variability.

![SHAPE_signal](img/SHAPE_signal_example.png)

You can find the derived SHAPE reactivity values and other useful plots and information in the output/ and plots/ 
subdirectories that are created during the analysis process. The files with the prefix output/agg_reactivity_values are the most direct
place to find the numeric values.

## License

SHAPEware is developed at [Arrakis Therapeutics](http://arrakistx.com/) and released under a GPL v3 license.

## Contact

For any questions or comments about SHAPEware, contact Luis Barrera (lbarrera [at] arrakistx {dot} com)

## References

- Haller A, et al. **Folding and ligand recognition of the TPP riboswitch aptamer at single-molecule resolution** _PNAS_. 2013 Mar 12; 110(11): 4188-4193. doi:  10.1073/pnas.1218062110

- Smola MJ, et al. **Selective 2'-hydroxyl acylation analyzed by primer extension and mutational profiling (SHAPE-MaP) for direct, versatile and accurate RNA structure analysis.** _Nat Protoc_. 2015 Nov;10(11):1643-69. doi: 10.1038/nprot.2015.103. Epub 2015 Oct 1.

- Eddy SR, et al. **Computational Analysis of Conserved RNA Secondary Structure in Transcriptomes and Genome** _Ann. Rev. Biophys_. 2015, 43, 433-456.
