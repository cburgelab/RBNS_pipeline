#!/usr/bin/env python
import os
import datetime
import pprint
import subprocess

import RBNS_utils




def get_sig_enriched_kmers_from_txt_R_F(
        txt_R_F,
        most_enriched_lib_conc = None,
        num_std_for_sig = 2):
    """
    INPUT:
        - txt_R_F: a file of enrichments from the pipeline (e.g.
        /net/uorf/data/nobackup/pfreese/RBNS_results/Fox_1_7_14/tables/RBFox2_enrichment_R.6mers.txt
        - most_enriched_lib_conc: the column header in the first row of the txt_F (e.g. "80");
            can pass this in manually. If nothing is passed in, the concentration
            with the highest enrichment will be chosen

    RETURNS:
        - return_D =
            {"sig_enriched_kmers_L": sig_enriched_kmers_L,
            "sig_enrichments_by_kmer_D": sig_enrichments_by_kmer_D}
    """
    #### get the dictionary of enrichments
    enriches_by_kmer_D = return_D_of_enrichments_from_txt_F(
            txt_R_F,
            most_enriched_lib_conc = most_enriched_lib_conc)

    #### a list of enrichments
    enrichments_L = enriches_by_kmer_D.values()
    mean, std = RBNS_utils.mean_std( enrichments_L )
    sig_threshold = mean + (num_std_for_sig * std)

    #### a list of the sig. enriched kmers and R values
    sig_enriched_kmer_R_tuples_L = []
    #### a dictionary of the enrichments, containing ONLY the sig. enriched
    ####    kmers
    sig_enrichments_by_kmer_D = {}
    for kmer in enriches_by_kmer_D:
        R = enriches_by_kmer_D[kmer]
        if (R >= sig_threshold):
            sig_enriched_kmer_R_tuples_L.append( (kmer, R) )
            sig_enrichments_by_kmer_D[kmer] = enriches_by_kmer_D[kmer]
    sig_enriched_kmer_R_tuples_L.sort( key = lambda x: -1*x[1] )
    sig_enriched_kmers_L = [tupl[0] for tupl in sig_enriched_kmer_R_tuples_L]

    return_D = {"sig_enriched_kmers_L": sig_enriched_kmers_L,
            "sig_enrichments_by_kmer_D": sig_enrichments_by_kmer_D}

    return return_D



def get_least_enriched_kmers_from_txt_R_F(
        txt_R_F,
        most_enriched_lib_conc = None,
        num_to_return = 25 ):
    """
    INPUT:
        - txt_R_F: a file of enrichments from the pipeline (e.g.
        /net/uorf/data/nobackup/pfreese/RBNS_results/Fox_1_7_14/tables/RBFox2_enrichment_R.6mers.txt
        - most_enriched_lib_conc: the column header in the first row of the txt_F (e.g. "80");
            can pass this in manually. If nothing is passed in, the concentration
            with the highest enrichment will be chosen

    RETURNS:
        - return_D =
            {"sig_enriched_kmers_L": sig_enriched_kmers_L,
            "sig_enrichments_by_kmer_D": sig_enrichments_by_kmer_D}
    """
    #### get the dictionary of enrichments
    enriches_by_kmer_D = return_D_of_enrichments_from_txt_F(
            txt_R_F,
            most_enriched_lib_conc = most_enriched_lib_conc)

    kmers_enriches_T_L = [(kmer, enriches_by_kmer_D[kmer]) for kmer in\
            enriches_by_kmer_D ]
    kmers_enriches_T_L.sort( key = lambda x: x[1] )

    kmers_to_return_L = [x[0] for x in kmers_enriches_T_L[:num_to_return]]

    return_D = {"least_enriched_kmers_L": kmers_to_return_L}

    return return_D





def return_D_of_enrichments_from_txt_F(
        txt_F,
        most_enriched_lib_conc = None):
    """
    INPUT:
        - txt_F: a file of enrichments from the pipeline (e.g.
        /net/uorf/data/nobackup/pfreese/RBNS_results/Fox_1_7_14/tables/enrich_naive.6.txt
        - most_enriched_lib_conc: the column header in the first row of the txt_F (e.g. "80");
            can pass this in manually. If nothing is passed in, the concentration
            with the highest enrichment will be chosen
    RETURNS;
        - kmer_to_enrichments_D, a dictionary with
        kmer_to_enrichments_D["20."] = 3.6
    """
    print txt_F
    #### get the enrichments for each of the concentrations
    enriches_by_conc_D = get_D_from_txt_table( txt_F )

    #### IF a most_enriched_lib_conc was passed in, use that concentration
    if (most_enriched_lib_conc != None):
        most_enriched_D = enriches_by_conc_D[most_enriched_lib_conc]
    #### Otherwise, calculate which conc has the highest enrichment
    else:
        most_enriched_lib = ""
        highest_enrichment_all_libs = 0.0

        # go through all libraries
        for lib_key in enriches_by_conc_D.keys():
            highest_enrichment_this_lib = 0.0
            # go through all 4^k kmers for that library
            for kmer in enriches_by_conc_D[lib_key].keys():
                enrich = enriches_by_conc_D[lib_key][kmer]
                if (enrich > highest_enrichment_this_lib):
                    highest_enrichment_this_lib = enrich
            # see if that library's highest enrichment is higher than all previous
            # libraries
            if (highest_enrichment_this_lib > highest_enrichment_all_libs):
                most_enriched_lib = lib_key
                highest_enrichment_all_libs = highest_enrichment_this_lib

        most_enriched_D = enriches_by_conc_D[most_enriched_lib]

    return most_enriched_D




def get_D_from_txt_table(
        txt_table_F ):
    """
    - Returns a dictionary of
        enriches_by_conc_D,
        where, e.g., enriches_by_conc_D["1300"]["AAAAA"] = 2.1123
    """
    header_line = True
    enriches_by_conc_D = {}
    with open( txt_table_F ) as f:
        for line in f:
            stripped_L = line.strip().split("\t")
            if header_line:
                concs_L = stripped_L[1:]
                #### all except the first column are library concentrations
                ####    (the first column is the protein name)
                for conc in concs_L:
                    enriches_by_conc_D[conc] = {}
                header_line = False
            #### for all subsequent lines
            else:
                kmer = stripped_L[0]
                for lib_num, val in enumerate(stripped_L[1:]):
                    conc_this_val = concs_L[lib_num]
                    try:
                        enriches_by_conc_D[conc_this_val][kmer] = float( val )
                    except KeyError:
                        enriches_by_conc_D[conc_this_val] = {kmer: float( val )}

    k = len( kmer )
    # make sure that each library has the proper number of values (note -
    ### this is now skipped because if there are kmers_to_ignore, it will not
    ### match
    for conc in enriches_by_conc_D.keys():
        if (len(enriches_by_conc_D[conc]) == 0):
            del enriches_by_conc_D[conc]
        #else:
        #    assert( len(enriches_by_conc_D[conc]) == (4**k) )
    # make sure there are an adequate number of libraries in the enriches_by_conc_D
    # (with maximally one library (input) having been removed)
    #assert((len(enriches_by_conc_D) == len(concs_L)) or\
    #        (len(enriches_by_conc_D) == len(concs_L)-1))
    return enriches_by_conc_D




def make_temp_reads_F(
        orig_reads_F,
        target_reads_DIR,
        read_length_to_use = "full_length",
        num_reads_to_use = "all",
        target_reads_basename = None,
        start_frac = 0.0 ):
    """
    - Given an orig_reads_F, makes a new reads file containing only reads that
        DON'T contain any N's

    - INPUT:
        - orig_reads_F: a file (e.g. .reads) of the reads
        - target_reads_DIR: where the output .reads file will be written
        - rd_length_to_use: the length of reads that will be included in the
            output reads file
            - "full_length": it will be the full read
        - num_reads_to_use:
            - "all": all reads
            - an int: will use up to that many reads
            - a float (from 0.0 to 1.0): that proportion of the total reads in
                orig_reads_F
        - target_reads_basename: the output reads basename; if not passed in,
            it will be the orig_reads_F basename with a time_stamp appended
        - start_frac: how far through the orig_reads_F to start getting reads

    - RETURNS:
        return_D = {"out_reads_F": out_reads_F,
                "num_reads_in_out_reads_F": rd_num}
    """
    #### Make the out_reads_F
    if (target_reads_basename == None):
        reads_basename = os.path.basename( orig_reads_F ).rsplit( ".", 1 )[0] +\
            ".{}.reads".format(datetime.datetime.now().strftime("%Hh_%Mm_%Ss" ))
    else:
        if (target_reads_basename[-6:] == ".reads"):
            reads_basename = target_reads_basename
        else:
            reads_basename = target_reads_basename + ".reads"

    out_reads_F = os.path.join( target_reads_DIR, reads_basename )
    os.system( "mkdir -p {}".format( target_reads_DIR ) )
    out_reads_f = open( out_reads_F, "w" )

    #### Get the num_reads_to_use
    total_lines_in_file = RBNS_utils.return_num_lines_in_F( orig_reads_F )
    if (num_reads_to_use == "all"):
        reads_to_use = 10000000000
    elif (type( num_reads_to_use ) is int ):
        reads_to_use = num_reads_to_use
    elif (type( num_reads_to_use ) is float):
        assert( num_reads_to_use <= 1.0 and num_reads_to_use > 0.0 )
        #### use the helper function in file_helpers.py
        reads_to_use = int(num_reads_to_use * total_lines_in_file)


    #### Get the line_to_start_at from the start_frac passed in
    if (start_frac == 0.0):
        line_to_start_at = 0
    else:
        line_to_start_at = int(total_lines_in_file * start_frac)

    line_lower = line_to_start_at
    line_upper = line_to_start_at + reads_to_use

    #### Get the read_length_to_use if it's "all"
    if (read_length_to_use == "full_length"):
        rd_length_to_use = get_readlength( orig_reads_F )
    else:
        assert( type(read_length_to_use) is int )
        rd_length_to_use = read_length_to_use

    #### Now populate the out_reads_F
    with open( orig_reads_F ) as f:
        reads_written = 0
        this_read = -1
        for line in f:
            this_read += 1
            if (this_read >= line_to_start_at):
                ln = line.strip()
                #### Only use reads that have no N's
                if ((len(ln) >= rd_length_to_use) and (ln.find("N") == -1)):
                    out_reads_f.write(ln[:rd_length_to_use] + "\n")
                    reads_written += 1
                if (reads_written >= reads_to_use):
                    break
    out_reads_f.close()

    return_D = {"out_reads_F": out_reads_F,
            "num_reads_in_out_reads_F": reads_written}
    return return_D





def get_readlength(
            reads_F ):
    """
    - Gets the read length from the first line
    """
    with open( reads_F ) as f:
        for line in f:
            read_len = len( line.strip() )
            break
    return read_len




def get_kmer_freqs_from_reads_F(
        reads_F,
        k,
        vals_sum_to = "sumto1"):
    """
    - Returns the kmer counts & freqs in reads_F

    - INPUTs:
        - vals_sum_to
            "sumto1": all 4^k entries sum to 1
            "sumto4^k": all 4^k entries sum to 4^k
    """
    counts_by_kmer_D = {}

    for kmer in RBNS_utils.yield_kmers( k ):
        counts_by_kmer_D[kmer] = 0
    with open( reads_F ) as f:
        for line in f:
            read = line.strip()
            for start_pos in range(len(read) - k + 1):
                kmer = read[start_pos:(start_pos+k)]
                # only include it if it doesn't have an N
                try:
                    counts_by_kmer_D[kmer] += 1
                except KeyError: pass
    return_D = {"counts_by_kmer_D": counts_by_kmer_D}

    if (vals_sum_to == "none"):
        return counts_by_kmer_D

    #### Normalize using the helper function in dict_helpers
    elif (vals_sum_to == "sumto1"):
        freqs_by_kmer_D = RBNS_utils.normalize_D( counts_by_kmer_D )
    elif ( vals_sum_to == "sumto4^k" ):
        freqs_by_kmer_D = RBNS_utils.normalize_D(
                counts_by_kmer_D,
                vals_sum_to = "sumto4^k")
    else:
        print "{0} IS NOT A VALID vals_sum_to ARGUMENT. REPLACE AND TRY AGAIN\n".format(
                vals_sum_to )

    return_D["freqs_by_kmer_D"] = freqs_by_kmer_D
    return return_D



def return_frequency_and_number_of_reads_kmer_in_reads_F(
        reads_F,
        kmer ):
    """
    - For a reads_F, makes a new out_reads_F in the same directory
        in which each occurrence of the kmer is replaced with "X"s
    - Called by functions in ~/RBNS/RBNS_motifs.py

    - RETURNS:
            return_D = {"out_reads_F": out_reads_F,
                    "tot_num_reads": tot_num_reads,
                    "num_reads_w_kmer": num_reads_w_kmer,
                    "freq_reads_w_kmer": freq_reads_w_kmer,
                    "tot_num_kmer_occurs" : tot_num_kmer_occurs,
                    "counts_by_kmer_D": counts_by_kmer_D}
                    "freqs_by_kmer_D": freqs_by_kmer_D}
    """
    k = len( kmer )
    read_len = get_readlength( reads_F )

    orig_reads_DIR = os.path.dirname( reads_F )
    orig_reads_basename = os.path.basename( reads_F )
    out_basename = orig_reads_basename.rsplit(".", 1)[0] +\
            "_{}.reads".format( kmer )
    #### If the file name is over 100 characters, shorten it
    if (len( out_basename ) >= 100):
        out_basename = "{}.reads".format( kmer )
    out_reads_F = os.path.join( orig_reads_DIR, out_basename )

    #### The number of reads and number of times a kmer was found
    tot_num_reads = 0
    num_reads_w_kmer = 0
    tot_num_kmer_occurs = 0

    #### A dictionary of kmer frequencies for the reads written out
    counts_by_kmer_D = {}
    for this_kmer in RBNS_utils.yield_kmers( k ):
        counts_by_kmer_D[this_kmer] = 0

    reads_f = open( reads_F )
    out_reads_f = open( out_reads_F, "w" )

    reads_to_write_out_L = []

    for line in reads_f:

        tot_num_reads += 1
        if ( len( reads_to_write_out_L ) == 10000 ):
            for read in reads_to_write_out_L:
                out_reads_f.write( read + "\n" )
            reads_to_write_out_L = []
        read = line.strip()

        cont = True
        found_any = False

        while (cont == True):
            kmer_pos = read.find( kmer )
            if (kmer_pos == -1):

                if (found_any == True):
                    num_reads_w_kmer += 1
                for start_pos in range( read_len - k + 1 ):
                    this_kmer = read[start_pos:(start_pos + k)]
                    try:
                        counts_by_kmer_D[this_kmer] += 1
                    except KeyError:
                        pass
                reads_to_write_out_L.append( read )
                #out_reads_f.write( read + "\n" )
                cont = False

            #### If an occurrence of this kmer was found, replace it with X's
            ####    and write out the read
            else:
                found_any = True
                tot_num_kmer_occurs += 1
                read = read[:kmer_pos] + "X"*k + read[(kmer_pos+k):]

    for read in reads_to_write_out_L:
        out_reads_f.write( read + "\n" )

    reads_f.close()
    out_reads_f.close()


    #### Normalize the counts_by_kmer_D into freqs using the helper function
    ####    in dict_helpers.py
    freqs_by_kmer_D = RBNS_utils.normalize_D( counts_by_kmer_D )
    freq_reads_w_kmer = float( num_reads_w_kmer ) / tot_num_reads

    return_D = {"out_reads_F": out_reads_F,
            "tot_num_reads": tot_num_reads,
            "num_reads_w_kmer": num_reads_w_kmer,
            "freq_reads_w_kmer": freq_reads_w_kmer,
            "tot_num_kmer_occurs" : tot_num_kmer_occurs,
            "counts_by_kmer_D": counts_by_kmer_D,
            "freqs_by_kmer_D": freqs_by_kmer_D}
    #### Remove the old reads_F
    #os.system( "rm {}".format( reads_F ) )
    return return_D





def get_num_reads_by_barcode_and_conc_D(
        barcode_log_txt_F ):
    """
    {'barcode_to_numreads_D': {'ACTGAT': 26787411,
                               'CGTACG': 21423352,
                               'GAGTGG': 27350419,
                               'GGTAGC': 16016161,
                               'GTTTCG': 15968650,
                               'TCGGCA': 11921053},
        'conc_to_numreads_D': {'5 nM': 26787411,
                                '20 nM': 16016161,
                                '80 nM': 27350419,
                                '320 nM': 21423352,
                                '1300 nM': 15968650,
                                'Input': 11921053}}
    """
    barcode_to_numreads_D = {}
    conc_to_numreads_D = {}

    with open( barcode_log_txt_F ) as f:
        next( f )
        for line in f:
            if ( line.find( 'total_reads' ) != -1 ):
                return_D = {
                        'barcode_to_numreads_D': barcode_to_numreads_D,
                        'conc_to_numreads_D': conc_to_numreads_D }
                return return_D
            try:
                barcode, conc, num_reads = line.strip().split("\t")[:3]
                try:
                    conc_descrip = "{} nM".format( int( float( conc ) ) )
                except:
                    conc_descrip = 'Input'
                num_reads = int( num_reads.replace(",","") )

                barcode_to_numreads_D[barcode] = num_reads
                conc_to_numreads_D[conc_descrip] = num_reads
            except:
                pass






