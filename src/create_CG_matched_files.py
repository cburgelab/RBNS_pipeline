#!/usr/bin/env python
import os, sys
import time, gzip, copy, pprint

import RBNS_utils

##### For RNA folding applications, these functions create sets of Pulldown reads
#####   files that match the distribution of C+G content within input reads



def make_out_Fs_of_PD_reads_that_match_input_CplusG_content(
        input_reads_w_str_F,
        PD_reads_w_str_Fs_L,
        fiveP_adapt_len,
        threeP_adapt_len,
        min_perc_reads_tokeep_CGbin = 4. ):
    """
    - Given reads files that contains the RNA folding output (4 lines per read)
        for the input and various pulldown libraries (PD_reads_w_str_Fs_L),
        along with the 5' and 3' adapter sequencing lengths (i.e., the read
        in input_reads_w_str_F should be
        random_len + len( fiveP_adapt_len ) + len( threeP_adapt_len ),

        calculates the 'bins' of total C+G content to use in the input,
        considering only bins that have at least min_perc_reads_tokeep_CGbin
        of the total reads

    - The returned dictionary is like:

    {'CGbin_normedfreqofreads_T_L': [(4, 0.011543667774174992),
                                    (5, 0.03146770413989078),
                                    (6, 0.06686856820248833),
                                    (7, 0.11379754117391035),
                                    (8, 0.15852543985619572),
                                    (9, 0.18160974445675013),
                                    (10, 0.16983774528197618),
                                    (11, 0.12929578756854745),
                                    (12, 0.08037548189544161),
                                    (13, 0.040380307164332135),
                                    (14, 0.01629801248629254)],
        'numreads_by_numCG_D': {0: 16,
                                1: 329,
                                2: 2825,
                                3: 15092,
                                4: 57129,
                                5: 155732,
                                6: 330929,
                                7: 563178,
                                8: 784534,
                                9: 898777,
                                10: 840518,
                                11: 639878,
                                12: 397774,
                                13: 199840,
                                14: 80658,
                                15: 25256,
                                16: 6209,
                                17: 1170,
                                18: 136,
                                19: 19,
                                20: 1}}
    """
    orig_DIR = os.path.dirname( input_reads_w_str_F )
    split_reads_DIR = orig_DIR.rsplit( "/", 1 )[0]
    out_DIR = os.path.join( split_reads_DIR, "fld_CG_match" )
    logs_DIR = os.path.join( out_DIR, "logs" )
    os.system( "mkdir -p {}".format( logs_DIR ) )

    ################### < SEE IF ALL FILES ALREADY EXIST > ####################
    all_exist = True
    out_input_F = os.path.join( out_DIR, os.path.basename( input_reads_w_str_F ) )
    all_out_Fs_L = [out_input_F]
    if not os.path.exists( out_input_F ) or ( os.stat( out_input_F ).st_size < 100000 ):
        all_exist = False
    for PD_reads_F in PD_reads_w_str_Fs_L:
        out_PD_F = os.path.join( out_DIR, os.path.basename( PD_reads_F ) )
        all_out_Fs_L.append( out_PD_F )
        if not os.path.exists( out_PD_F ) or ( os.stat( out_PD_F ).st_size < 100000 ):
            all_exist = False
    if all_exist:
        print "INPUT and all PD already exist:\n\n{0}\n\t{1}\n".format( out_DIR,
            "\n\t".join( [os.path.basename(x) for x in all_out_Fs_L] ) )
        return
    ################### </ SEE IF ALL FILES ALREADY EXIST > ###################

    making_F = os.path.join( out_DIR, "making.txt" )
    if os.path.exists( making_F ):
        return
    with open( making_F, 'w' ) as f:
        pass

    ##### Get the distribution of C+G content in the input reads, for that
    # input_returned_D = {
    #        "numreads_by_numCG_D": numreads_by_numCG_D,
    #       'CGbin_normedfreqofreads_T_L': CGbin_normedfreqofreads_T_L }
    input_returned_D = get_num_reads_by_CplusG_content(
            input_reads_w_str_F,
            fiveP_adapt_len,
            threeP_adapt_len,
            min_perc_reads_for_CGbin = min_perc_reads_tokeep_CGbin )

    #### The total number of input reads for each of the bins
    input_numreads_to_keep_by_bin_D = {}
    input_out_F = os.path.join( logs_DIR, "input.CGbins_used.txt" )
    with open( input_out_F, "w" ) as f:
        pass
    for num_CG, perc in input_returned_D['CGbin_normedfreqofreads_T_L']:
        ##### Write out the frequency of reads in each bin
        with open( input_out_F, "a" ) as f:
            f.write( "CG_bin: {0}\t{1:.1f}%\n".format( num_CG, 100. * perc ) )

        input_numreads_to_keep_by_bin_D[num_CG] = input_returned_D['numreads_by_numCG_D'][num_CG]

    ########################### < GET THE INPUT READS > #######################
    write_output_F_from_inputreads_and_numreadstokeepbyCG_D(
            input_reads_w_str_F,
            out_input_F,
            input_numreads_to_keep_by_bin_D,
            fiveP_adapt_len,
            threeP_adapt_len )
    ########################## </ GET THE INPUT READS > #######################

    ########## < GET EACH OF THE PD READS TO MATCH INPUT C/G> #################
    for PD_reads_F in PD_reads_w_str_Fs_L:
        out_PD_F = os.path.join( out_DIR, os.path.basename( PD_reads_F ) )
        write_output_F_from_inputreads_props_max_given_PD_readstokeepbyCG_D(
                PD_reads_F,
                out_PD_F,
                input_numreads_to_keep_by_bin_D,
                fiveP_adapt_len,
                threeP_adapt_len )
    ########## </ GET EACH OF THE PD READS TO MATCH INPUT C/G> ################

    os.system( "rm -rf {}".format( making_F ) )








def write_output_F_from_inputreads_and_numreadstokeepbyCG_D(
        in_reads_w_str_F,
        out_reads_w_str_F,
        numreads_to_keep_by_bin_D,
        fiveP_adapt_len,
        threeP_adapt_len ):
    """
    - Given a pulldown folded reads file (in_reads_w_str_F), will write out
        a new file (out_reads_w_str_F) containing a subset of reads according
        to numreads_to_keep_by_bin_D, which dictates how many reads containing
        each number of C+G bases should be included
    """
    #### If the out_reads_w_str_F already exist, PASS
    if ( os.path.exists( out_reads_w_str_F ) and\
            os.stat( out_reads_w_str_F ).st_size > 100000 ):
        return

    #### Make a copy of the numreads_to_keep_by_bin_D to run down
    copied_numreads_to_keep_by_bin_D = copy.copy( numreads_to_keep_by_bin_D )

    out_DIR = os.path.dirname( out_reads_w_str_F )
    os.system( "mkdir -p {}".format( out_DIR ) )

    log_DIR = os.path.join( out_DIR, "logs" )

    out_log_F = os.path.join( log_DIR,
            os.path.basename( out_reads_w_str_F ).split(".gz")[0] + ".log.txt" )

    out_f = gzip.open( out_reads_w_str_F, 'wb' )

    #### Get the read length
    for four_lines_L in RBNS_utils.iterNlines(
            in_reads_w_str_F, 4, strip_newlines = True ):
        read_len = len( four_lines_L[0] ) - fiveP_adapt_len - threeP_adapt_len
        break
    assert( read_len in [20, 40] )

    #### The total number of reads to write out
    num_reads_to_write = sum( copied_numreads_to_keep_by_bin_D.values() )

    #### Get the read length
    for four_lines_L in RBNS_utils.iterNlines(
            in_reads_w_str_F, 4, strip_newlines = True ):

        rand_read = four_lines_L[0][fiveP_adapt_len:(fiveP_adapt_len+read_len)]
        num_CG = len( [x for x in rand_read if x in ["C", "G"]] )
        if ( num_CG in copied_numreads_to_keep_by_bin_D ):
            if ( copied_numreads_to_keep_by_bin_D[num_CG] > 0 ):
                out_f.write( "\n".join( four_lines_L ) + "\n" )
                copied_numreads_to_keep_by_bin_D[num_CG] -= 1
                num_reads_to_write -= 1

        if ( num_reads_to_write == 0 ):
            break


    out_f.close()
    with open( out_log_F, 'w' ) as f:
        CG_bins_L = copied_numreads_to_keep_by_bin_D.keys()
        CG_bins_L.sort()
        f.write( "REMAINING READS:\n" )
        for CG_bin in CG_bins_L:
            f.write( "CG {0}:\t{1}\n".format(
                CG_bin, copied_numreads_to_keep_by_bin_D[CG_bin] ) )

        f.write( "\n\n reads used:\n" )
        for CG_bin in CG_bins_L:
            f.write( "CG_num_reads {0}:\t{1}\n".format(
                CG_bin, numreads_to_keep_by_bin_D[CG_bin] ) )



def get_num_reads_by_CplusG_content(
        reads_w_str_F,
        fiveP_adapt_len,
        threeP_adapt_len,
        min_perc_reads_for_CGbin = 4. ):
    """
    - Give a file of reads, calculates the number of reads in each C+G bin
        (counting the number of C+G's in the random region) and determines
        which bins will be used to match PD reads to
    - Only C+G bins that have at least min_perc_reads_for_CGbin % of reads will
        be used (i.e., don't want to use very lowly populated C+G bins)
    """
    ##['GAGTTCTACAGTCCGACGATCTGAACCGAACATATTCTACGTGGAATTCTCGGGTGCCAAGG',
    ## '0.819 0.818 0.867 0.867 0.866 0.900 0.857 0.821 0.983 0.182 0.930 0.953 0.176 0.076 0.834 0.848 0.017 0.957 0.196 0.868 0.948 0.202 0.200 0.009 0.012 0.089 0.102 0.968 0.947 0.093 0.787 0.029 0.042 0.120 0.910 0.950 0.232 0.119 0.777 0.781 0.822 0.890 0.988 0.998 0.965 0.961 0.989 0.898 0.851 0.226 0.185 0.874 0.877 0.096 0.074 0.077 0.756 0.746 0.005 0.091 0.103 0.056',
    ##  '(((((((((.((..((.(.((......)).)...))..)))))))))))..((...))....',
    ##  'sssssssssissiissisisshhhhhhssisiiissiisssssssssssmmsshhhssffff']
    numreads_by_numCG_D = {}
    start_time = time.time()

    #### Get the read length
    for four_lines_L in RBNS_utils.iterNlines( reads_w_str_F, 4, strip_newlines = True ):
        read_len = len( four_lines_L[0] ) - fiveP_adapt_len - threeP_adapt_len
        break
    for i in range( read_len + 1 ):
        numreads_by_numCG_D[i] = 0

    #### Go through all of the reads
    for four_lines_L in RBNS_utils.iterNlines( reads_w_str_F, 4, strip_newlines = True ):
        rand_read = four_lines_L[0][fiveP_adapt_len:(fiveP_adapt_len+read_len)]
        num_CG = len( [x for x in rand_read if x in ["C", "G"]] )
        numreads_by_numCG_D[num_CG] += 1

    pprint.pprint( numreads_by_numCG_D )
    end_time = time.time()

    #### Get the number of reads in each CG bin
    total_num_reads = float( sum( numreads_by_numCG_D.values() ) )
    bin_percreads_T_L = [(num_CG, numreads_by_numCG_D[num_CG]*100./total_num_reads)\
            for num_CG in numreads_by_numCG_D ]

    #### Prune for those CG-bins that have at least x%
    bin_percreads_T_L = [tupl for tupl in bin_percreads_T_L if (tupl[1] >= min_perc_reads_for_CGbin)]
    #### Renormalize so the sum of the pruned bins adds up to 1.
    total_after_pruned = sum( [tupl[1] for tupl in bin_percreads_T_L] )
    CGbin_normedfreqofreads_T_L = [(tupl[0], tupl[1] / total_after_pruned)\
            for tupl in bin_percreads_T_L]
    CGbin_normedfreqofreads_T_L.sort( key = lambda x: x[0] )

    return_D = {
            "numreads_by_numCG_D": numreads_by_numCG_D,
            'CGbin_normedfreqofreads_T_L': CGbin_normedfreqofreads_T_L }

    return return_D






def write_output_F_from_inputreads_props_max_given_PD_readstokeepbyCG_D(
        in_reads_w_str_F,
        out_reads_w_str_F,
        numreads_to_keep_by_bin_D,
        fiveP_adapt_len,
        threeP_adapt_len ):
    """
    - Given a set of input reads (in_reads_w_str_F) and the number of reads
        to keep for each C+G bin, gets the required number of reads in each bin
        and writes them out to out_reads_w_str_F
    """
    #### If the out_reads_w_str_F already exist, RETURN
    if ( os.path.exists( out_reads_w_str_F ) and\
            os.stat( out_reads_w_str_F ).st_size > 100000 ):
        return

    #### Make a copy of the numreads_to_keep_by_bin_D to run down
    cp_numreads_to_keep_by_bin_D = copy.copy( numreads_to_keep_by_bin_D )

    for four_lines_L in RBNS_utils.iterNlines(
            in_reads_w_str_F, 4, strip_newlines = True ):
        read_len = len( four_lines_L[0] ) - fiveP_adapt_len - threeP_adapt_len
        break
    assert( read_len in [20, 40] )

    #### FOr the in_reads_w_str_F, GET the actual number of CG reads in each bin
    ####    to see if there are enough for each; if not, we will need to
    ####    downsample each of the CG bins
    num_reads_in_F_by_CG_bin_D = {}
    for CG_bin in cp_numreads_to_keep_by_bin_D:
        num_reads_in_F_by_CG_bin_D[CG_bin] = 0
    for four_lines_L in RBNS_utils.iterNlines(
            in_reads_w_str_F, 4, strip_newlines = True ):

        rand_read = four_lines_L[0][fiveP_adapt_len:(fiveP_adapt_len+read_len)]
        num_CG = len( [x for x in rand_read if x in ["C", "G"]] )
        try:
            num_reads_in_F_by_CG_bin_D[num_CG] += 1
        except KeyError:
            pass

    ##### Go through and get the lowest CG_bin factor - that is, which needed C+G
    #####   bin is least populated and will be the factor multiplied by the
    #####   values of numreads_to_keep_by_bin_D to determine how many in each
    ####    C+G bin will __actually__ be able to be written to out_reads_w_str_F
    min_CG_factor = 1.
    any_limiting_CG_bin = 'No limiting CG bin'
    for CG_bin, num_target_reads in cp_numreads_to_keep_by_bin_D.iteritems():

        num_in_reads = num_reads_in_F_by_CG_bin_D[CG_bin]
        prop_of_reads_needed_in_inF = float( num_in_reads ) / num_target_reads
        if ( prop_of_reads_needed_in_inF < min_CG_factor ):
            any_limiting_CG_bin = 'CG bin {0}: desired {1:,}, have only {2:,} -> min_CG_factor = {3:.2f}'.format(
                    CG_bin, num_target_reads, num_in_reads, prop_of_reads_needed_in_inF )
        min_CG_factor = min( min_CG_factor, prop_of_reads_needed_in_inF )

    #### Now go through and downsample the copied_numreads_to_keep_by_bin_D if
    ####    necessary (i.e., if the min_CG_factor is less than 1
    copied_numreads_to_keep_by_bin_D = {}
    for CG_bin, orig_num_target_reads in cp_numreads_to_keep_by_bin_D.iteritems():
        copied_numreads_to_keep_by_bin_D[CG_bin] = int( min_CG_factor * orig_num_target_reads )

    out_DIR = os.path.dirname( out_reads_w_str_F )
    os.system( "mkdir -p {}".format( out_DIR ) )

    log_DIR = os.path.join( out_DIR, "logs" )

    out_log_F = os.path.join( log_DIR,
            os.path.basename( out_reads_w_str_F ).split(".gz")[0] + ".log.txt" )

    out_f = gzip.open( out_reads_w_str_F, 'wb' )

    #### Get the read length
    for four_lines_L in RBNS_utils.iterNlines(
            in_reads_w_str_F, 4, strip_newlines = True ):
        read_len = len( four_lines_L[0] ) - fiveP_adapt_len - threeP_adapt_len
        break

    #### The total number of reads to write out
    num_reads_to_write = sum( copied_numreads_to_keep_by_bin_D.values() )

    #### Get the read length
    for four_lines_L in RBNS_utils.iterNlines(
            in_reads_w_str_F, 4, strip_newlines = True ):

        rand_read = four_lines_L[0][fiveP_adapt_len:(fiveP_adapt_len+read_len)]
        num_CG = len( [x for x in rand_read if x in ["C", "G"]] )
        if ( num_CG in copied_numreads_to_keep_by_bin_D ):
            if ( copied_numreads_to_keep_by_bin_D[num_CG] > 0 ):
                out_f.write( "\n".join( four_lines_L ) + "\n" )
                copied_numreads_to_keep_by_bin_D[num_CG] -= 1
                num_reads_to_write -= 1

        if ( num_reads_to_write == 0 ):
            break

    out_f.close()
    with open( out_log_F, 'w' ) as f:
        CG_bins_L = copied_numreads_to_keep_by_bin_D.keys()
        CG_bins_L.sort()
        f.write( any_limiting_CG_bin )

        f.write( "\n\nREMAINING READS:\n" )
        for CG_bin in CG_bins_L:
            f.write( "CG {0}:\t{1}\n".format( CG_bin, copied_numreads_to_keep_by_bin_D[CG_bin] ) )








