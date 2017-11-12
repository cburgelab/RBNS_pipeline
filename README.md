![Logo](img/RBNS_logo.png)

***


Welcome to the repository for RNA Bind-N-Seq Analysis! 

The RBNS pipeline is a set of bioinformatics tools to analyze data from high-throughput sequencing experiments of protein-bound RNAs. The current version includes read splitting, calculation of kmer frequencies and enrichments, QC metrics, production of motif sequence logos, and RNA secondary structure analysis. Incorporation of functions to compute presence of bipartite motifs & flanking nucleotide context preferences are forthcoming.


## Installation

##### Requirements

The RBNS pipeline is designed to run on a Linux computing cluster, optionally with jobs parallelized by submitting them to a PBS/Torque queue. In addition, it requires the following software to be pre-installed on your computing environment:

- Python (tested on version 2.7.11)
- The [Miniconda](https://conda.io/miniconda.html) or [anaconda](https//docs.anaconda.com/) package manager (if you download this now, be sure to 'source ~/.bashrc' after so that 'conda' is on your $PATH).
- The [Weblogo](http://weblogo.threeplusone.com/manual.html) program (if sequence motif logos are to be produced.)
- The [forgi](https://viennarna.github.io/forgi/) library (if RNA secondary secondary structure analysis is performed).
- The [RNAfold](https://www.tbi.univie.ac.at/RNA/) program (if RNA secondary secondary structure analysis is performed).

If you need help installing any of these tools, see the [detailed documentation](docs/installation.md). When installing dependencies, make sure you agree with the corresponding licenses of various software tools.


## Download and install the RBNS_pipeline

The easiest way to get the RBNS pipeline software is to clone this repository. This will ensure you always have access to the latest version. 

	https://pfreese@bitbucket.org/pfreese/rbns_pipeline.git

After cloning the repository, run the included installation script: 

	cd rbns_pipeline
	./install.sh

This will use the Conda package manager to ensure that you have all of the key dependencies. This step will create a stable environment in which to run analysis jobs. 
This approach will keep your results reproducible and will not affect other software that you have installed on your system.
While you can also use the RBNS pipeline without this step by installing all necessary packages manually, it is not recommended.


## Test the RBNS pipeline with example data

In this version, the RBNS_pipeline is able to analyze RBNS (Dominguez et al., 2017) data. You can find example input files in the test_data/ directory within the repository. These were derived from experiments that assayed the RBFOX3 protein.

Once the script has finished running, you can find the output from the pipeline in the results_dir given in the settings.RBFOX3.json file.


### RBNS_pipeline

The inputs to the RBNS_pipeline are described in more detail [here](docs/input_files.md). Here is a quick summary:

- A settings .json file describing the experiment, different libraries assayed, and what counts & optional additional functionalities are to be performed by the pipeline.
- A FASTQ file containing the multiplexed sequencing reads from the different libraries to be split & analyzed.

You can find examples of all of these files in the test_data/ folder. It is probably easiest to just take a look at these files first. You can run the RBNS_pipeline on this example and reproduce the logo below.

#### Example output: RBFOX3

The complete set of output files is described in [here](docs/output_information.md). Briefly, output should include:

- Split read files for each library, including files containing QC stats about the number of reads in each library & library complexity.
- Enrichment tables of kmers
- Pickled intermediate files of kmer counts & frequencies.

Additional output files depending on functionalities requested include:

- SKA (Streaming kmer Assignment) library fraction and nt. frequencies & kmer enrichments at different positions of the random region
- Sequence motif logos as shown below
- RNA secondary structure analyses of the top enriched kmers

An example of an output logo for RBFOX3, derived from computing the significantly enriched 5mers through an interative procedure at a Z-score > 3 cutoff, is:

![RBFOX3_logo](img/RBFOX3_5mer_seqlogos.png)


## License

RBNS_pipeline is developed by Peter Freese and released under a GPL v3 license.

## Contact

For any questions or comments about the RBNS pipeline, contact Peter Freese (pfreese [at] mit {dot} edu).

## References

- Lambert, et al. **RNA Bind-n-Seq: quantitative assessment of the sequence and structural binding specificity of RNA binding proteins** _Mol Cell_. 2014 Jun 5;54(5):887-900. doi:  [10.1016/j.molcel.2014.04.016](https://www.ncbi.nlm.nih.gov/pubmed/24837674)
- Lambert, et al. **RNA Bind-n-Seq: Measuring the Binding Affinity Landscape of RNA-Binding Proteins** _Methods Enzymol_. 2015 558:465-93. doi:  [10.1016/bs.mie.2015.02.007](https://www.ncbi.nlm.nih.gov/pubmed/26068750)
- Dominguez et al. **Sequence, Structure and Context Preferences of Human RNA Binding Proteins** _bioRxiv_. 2017 Oct 12. doi:  [10.1101/201996](https://www.biorxiv.org/content/early/2017/10/12/201996)

