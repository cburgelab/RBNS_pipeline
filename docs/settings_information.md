![Logo](../img/RBNS_logo.png)

***
## The settings file is a .json file that contains the required information for an RBNS experiment (e.g., file and output paths, protein name and concentrations, FASTQ barcodes), and it also tells RBNS_main.py which functionalities to perform.

### REQUIRED information for each RBNS experiment

- fastq: the file that contains the contains the sequencing reads. Preferably an absolute path to the file, although if just a basename is given (as is the case for the test_data), it will assume the .fastq is: 'test_data/$FASTQ')
- read_len: An int of the random region read length (e.g., 20 or 40). Note that the reads in the fastq file can be longer than this (e.g., if the random reads are of length 20 but 40bp was sequenced, with the last 20 bases of each being the 3' adapter) - in this case, only the first read_len positions of each will be used.
- protein_name: the name of the protein analyzed. Highly recommended that this be alpha-numeric  only (no underscores, dots, etc. as these may cause an error in name parsing within some scripts)
- barcodes: the 'barcodes', and 'concentrations' parameters must be lists of the same length and in the same order. The barcodes correspond to the sequencing barcode reported in the 1st line of each of the 4 read lines in the .fastq (for example, if the 1st line is like '+HWI-ST:7:1101:6316:2864#AGTTCC/1', the barcode corresponding to this concentration would be 'AGTTCC'. Note that the get_barcode() function in RBNS_utils.py extracts the barcode according to 'line.split(begin_barcode_symb)[-1].split(end_barcode_symb)[0]', where begin_barcode_symb is '#' and end_barcode_symb is '/' by default, though these can be changed via optional parameter specifications, or the get_barcode() function can be changed to suit your needs).
- concentrations: A list of integers of the protein concentrations (typically in nM) corresponding to the same order of the barcodes. The 'input' concentration should be 0. These must be unique (with the possible exception of two 0's if one corresponds to input and the other to a 0 nM library) as file names will commonly be of the form $PROTEINNAME_$CONCENTRATION, with 'input' replacing the 0 $CONCENTRATION for that library). If you happen to have two replicates of the same concentration, rename one of them for this purpose (e.g., 200 and 201).
- input_barcode: This must be the same as one of the entries in the 'barcodes' list, with the corresponding 'concentration' being 0.
- results_dir: The main directory where all files will be output.
- scratch_dir: For some I/O intensive functions (such as motif logo generation and RNAfolding), many temporary files will be written to and read from temporary directories within here. It will be made and deleted upon successful completion of the pipeline run.


### REQUIRED parameter specifications

Specifying each of the following in the .json file is required; if left out, an error will be immediately returned.

- ks_to_test_naive: By default, [5]. A list of k sizes to compute enrichments for; note that 5 should always be in this list as kmers are 'sorted' by decreasing 5mer enrichment values.


### OPTIONAL functionalities & parameter specificiations

Specifying any of the following in the .json file is optional; if left out, they will simply not be computed or the defaults will be used.

### Functionalities
### ============

- nt_freqs_by_position: If True, it will calculate and plot the nucleotide frequencies at each position in each of the libraries, so you can see if your input library is biased.
- ks_to_test_by_position: If True, it will calculate the frequency of each kmer at each of the read_len - k + 1 positions in the random region to see if kmers are evenly or unevenly distributed from 5' -> 3'.
- stream_count and ks_to_test_stream: If stream_count = True, Streaming Kmer Assignment (SKA) value library fractions will be calculated and reported. See https://bitbucket.org/marjens/cska for a faster version implemented in Cython.
- naive_max_once_count and ks_to_test_naive_max_once: By default False, but if naive_max_once_count = True, it will compute kmer frequencies and enrichments only counting up to 1 occurrence of each kmer within each read (e.g., UUUUU for a read containing UUUUUUUUUU would only get one count rather than 6 as it would with normal naive enrichments.
- ks_to_test_logos and Z_scores_for_logos: Lists which specify the k's for which to perform logo generation, with the Z_scores_for_logos being the lower bound for Z-score enriched kmers to be included for inclusion into logo alignments. ks_to_test_logos should be a list of integer(s) (k=5 and/or 6 recommended), with Z_scores_for_logos being a list of integers and/or floats (can mix & match, e.g., [2, 2.5, 3]). Logos will be made for all pairwise combinations of k & Z_score (logos used in Dominguez et al. 2017 were for k=5, Z_score = 3).
- fold_each_reads_F: If True, will perform RNAfold'ing on the libraries (default = False) by folding each reads file, creating a file of 5 million input sequences as well as subsets of each pulldown library (up to 5 million reads each) that have a matched C+G content distribution within the random region.

### Advanced parameters
### ==================

- rna_5p_adapter and rna_3p_adapter: By default, "GGGGAGTTCTACAGTCCGACGATC" as "TGGAATTCTCGGGTGTCAAGG" as these were used in all Lambert et al. 2014 & Dominguez et al. 2017 experiments. All RNAfolding will fold the reads as rna_5p_adapter + random region + rna_3p_adapter, but only analyze the positions out_f the returned structure corresponding to the random region.
- start_of_adapter_seq: If the read_len is significantly shorter than the .fastq read length, the beginning of the 3' sequencing adapter can be included as  start_of_adapter_seq, and only reads that contain it at the expected position (i.e., no indels in the read) will be included. For example, if reads of length 40 were sequenced for random 20mer oligos with the 3' sequencing adapter being 'TGGAATTCTCGGGTGTCAAGG', setting start_of_adapter_seq = 'TGGAAT' would search for this sequence from the 3' end of reads, and only reads that include this sequence at positions 20-26 would be analyzed).
- temp: the temperature at which the RBNS experiment was performed, in degrees Celsius (default = 4). This parameter is only used as a setting for RNAfold (e.g., the X for 'RNAfold -p --temp=X')
- num_reads_per_folding_block: by default, 1000000. RNA folding of each library will be split up into jobs, each containing num_reads_per_folding_block reads per job. If you are submitting your jobs to a cluster queue, you may want to decrease this parameter depending on the length each job is allowed to run for. A folding job with 1000000 may take anywhere from single-digit to ~20 hours depending on your cluster, so if a job length is limited to 12 hours, you may want to set this to 500000 to be safe, for example.
- mismatches_allowed_in_barcode: by default, 0. Greater than 0 will allow barcode mismatches within a Hamming distance of mismatches_allowed_in_barcode to one of the specified barcodes to count toward the library.
- num_reads_for_logos: If logo generation is requested, this parameter tells how many reads OR what proprotion of reads to use for each of the two sets on which logos will be made (with the final output logo only including kmers appearing in both logos). This argument can be either an integer (the number of reads, e.g., 1000000) or a proportion (e.g., 0.1). By default, this argument is 0.5 (the two logos will each be performed on half of the reads), but for testing purposes, you may want to initially set this low (e.g., 50000) to make sure everything works before running it on a larger read set.
- fold_all_or_mostenrichedconc_only: If RNAfolding is requested, this argument should either be 'all' or 'most_enriched_only' (default: 'most_enriched_only'). If 'all', all protein libraries will be folded, while if 'most_enriched_only', only the input library and the pulldown library with the highest enrichment will be folded.
- begin_barcode_symb and end_barcode_symb: '#' and '/', respectively (see 'barcodes' in the 'REQUIRED' information above).
- kmers_to_ignore: If your enrichments are low and kmers you know to be spurious are interspersed with the top kmers, you can enter a list of short sequences with this argument and any kmers containing any of them will not be included in enrichment tables or logo generation (e.g., if you are confident your motif is C-rich but some low A-rich kmers appear interspersed, kmers_to_ignore = ['AAA'] would cause any kmers containing an A3 to be ignored).
- ks_to_fold: By default, [6] (i.e., will analyze 6mers if fold_each_reads_F is True).



