![Logo](../img/RBNS_logo.png)

***
# The output directories & files of an RBNS_pipeline run are (assuming you ran the pipeline on the RBFOX3 test_data, which contains the 20 nM and 1300 nM pulldown libraries as well as the input library):

## split_reads/
- RBFOX3_barcode_log.txt contains the number of reads assigned to each barcode, as well as the observed insert lengths if the 'start_of_adapter_seq' argument was specified in the .settings file. It also includes a summary of the 'bad' barcodes (i.e., those that did not match the ones specified).
- 'RBFOX3_1300.reads', 'RBFOX3_1300.reads', and 'RBFOX3_1300.reads' contain just the random portion of each read; RBFOX3_1300.fastq.gz, etc. contain the 4 lines for each read that were in the original FASTQ file.
- 'RBFOX3_input.lib_complexity', etc. contain the library complexity of each library (how many sequences were unique, how many showed up 2x, 3x, etc. in the library)
- 'RBFOX3_1300.wrong_insert_len.reads': if 'start_of_adapter_seq' was specified, these files contain reads with the library barcode but with an incorrect insert length (indel in the read; these reads were not used).

## unique/

- This directory contains the same files as the 'split_reads' directory, except that each read library contains only 1 occurrence of each read (i.e., even if occurs at more than 1x in 'RBFOX3_input.lib_complexity', it only occurs once in these split read files). Our own analysis shows enrichments, etc. from these files are largely the same as those on the original split_read files, but if you have very low library complexity, you may want to use these read files instead.

## 

## errors/

- If all goes well, there will be nothing in here, but if you launch jobs and something goes wrong, you may want to investigate here.

