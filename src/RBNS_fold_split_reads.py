#!/usr/bin/env python
import os, sys, time, gzip
import subprocess, pprint, glob
import random
import shutil
import socket
import inspect
import cPickle as pickle
import numpy as np

import RBNS_utils
import RBNS_cluster_utils

import forgi.graph.bulge_graph as cgb





def submit_get_Ppaired_DotBracket_andletters_for_reads_F_for_block(
        in_reads_F,
        temp,
        block_idx,
        starting_scratch_DIR ,
        fiveP_adapter,
        threeP_adapter ):
    """
    - For an in_reads_F split_reads file,
        submits a job to run the
        get_Ppaired_DotBracket_andletters_for_reads_F_for_block() function
           below
    """
    #### Get this file path
    filename = inspect.getframeinfo( inspect.currentframe() ).filename
    this_script_path = os.path.abspath( filename )

    output_DIR = os.path.dirname( in_reads_F )
    errors_outputs_DIR = os.path.join( output_DIR, "errors_outputs" )
    RBNS_utils.make_dir( errors_outputs_DIR )

    command = ('hostname ; python %(this_script_path)s '
            'get_Ppaired_DotBracket_andletters_for_reads_F_for_block '
            '%(in_reads_F)s '
            '%(temp)s '
            '%(block_idx)s '
            '%(fiveP_adapter)s '
            '%(threeP_adapter)s '
            '%(starting_scratch_DIR)s ' % locals())

    job_name = "{0}_block{1}_get_Ppaired_DotBracket_andletters_for_reads_F_for_block".format(
            os.path.basename( in_reads_F ).split(".reads")[0], block_idx )
    RBNS_cluster_utils.launch(
            command,
            out_file = os.path.join( errors_outputs_DIR,
                "{}.log".format( job_name ) ),
            jobname = job_name,
            error_dir = errors_outputs_DIR,
            q = 'long' )








def get_Ppaired_DotBracket_andletters_for_reads_F_for_block(
        in_reads_F,
        temp,
        block_idx,
        starting_scratch_DIR,
        fiveP_adapter,
        threeP_adapter,
        num_reads_per_block = 1000000,
        num_reads_to_get_status_after = 1000,
        copy_back_every_x_seconds = 10000000 ):
    """
    - For an in_reads_F in the /split_reads directory, folds a block of
        num_reads_per_block reads, where this is for the block_idx, at the
        reported temp

    - starting_scratch_DIR is a directory in which temporary I/O files will go
    - fiveP_adapter & threeP_adapter will be added to each of the reads in
        in_reads_F

    - Makes an output file in /split_reads/fld/by_block
    """
    temp = int( temp )
    block_idx = int( block_idx )

    start_basename = os.path.basename( in_reads_F ).split(".reads")[0] +\
            ".block_{}".format( block_idx )

    scratch_DIR = os.path.join( starting_scratch_DIR,
            "{0}_{1}".format( start_basename, random.randint( 0, 100000 ) ) )
    os.system( "mkdir -p {}".format( scratch_DIR ) )

    #### If the out_DIR is different than the one simply passed in
    out_DIR = os.path.join( os.path.dirname( in_reads_F ), "fld", "by_block" )
    out_logs_DIR = os.path.join( out_DIR, "logs" )
    os.system( "mkdir -p {}".format( out_logs_DIR ) )

    prev_copies_DIR = os.path.join( os.path.join( out_DIR, "prev_copy" ) )
    os.system( "mkdir -p {}".format( prev_copies_DIR ) )

    #### copy over the reads to the scratch space
    scratch_in_reads_F = os.path.join( scratch_DIR,
            os.path.basename( in_reads_F ) )
    print "copying to {}".format( scratch_in_reads_F )

    lower_this_block = block_idx * num_reads_per_block
    upper_this_block = ( block_idx + 1 ) * num_reads_per_block

    RBNS_utils.copy_lines_lower_through_upper_to_another_F_without_N(
            in_reads_F,
            scratch_in_reads_F,
            lower_this_block,
            upper_this_block )
    #shutil.copyfile( in_reads_F, scratch_in_reads_F )
    print "\tDONE copying {}".format( os.path.basename( in_reads_F ) )

    #### Get the read length
    with open( scratch_in_reads_F ) as f:
        for line in f:
            read_len = len( line.strip() )
            break

    scratch_reads_F = os.path.join( scratch_DIR,
            "{}.w_struc.reads.gz".format( start_basename ) )
    out_reads_F = os.path.join( out_DIR,
            "{}.w_struc.reads.gz".format( start_basename ) )
    log_F = os.path.join( out_logs_DIR,
            "{}.w_struc.log.txt".format( start_basename ) )

    with open( log_F, "w" ) as log_f:
        log_f.write( "Starting at: {}\n".format(
            RBNS_utils.return_nice_datetime_str_for_filename() ) )

        log_f.write( "hostname is: {}\n".format( socket.gethostname() ) )

        log_f.write( "\n\n" + "="*80 + "\n" + "="*80 + "\n\n" )

        log_f.write( "in_reads_F: {}\n".format( in_reads_F ) )
        log_f.write( "temp: {}\n".format( temp ) )
        log_f.write( "fiveP_adapter: {}\n".format( fiveP_adapter ) )
        log_f.write( "threeP_adapter: {}\n".format( threeP_adapter ) )

        log_f.write( "\n\n" + "="*80 + "\n" + "="*80 + "\n\n" )

    num_reads_covered = 0
    start_time = time.time()
    last_copy_back_time = time.time()

    copy_from_file = os.path.join( scratch_DIR,
        "{}.w_struc.reads.gz".format( start_basename ) )
    copy_to_file = os.path.join( out_DIR,
        "{}.w_struc.reads.gz".format( start_basename ) )
    with gzip.open( scratch_reads_F, "wb" ) as out_f:
        pass

    with open( scratch_in_reads_F ) as f:

        reads_this_set_L = []

        num_reads_this_set = 0
        for line in f:

            rand_RNA_seq = line.strip()
            reads_this_set_L.append( rand_RNA_seq )
            num_reads_this_set += 1
            num_reads_covered += 1

            if ( num_reads_this_set == num_reads_to_get_status_after ):

                avg_bp_probs_Ls_L = get_avg_bp_probs_w_adapters(
                    reads_this_set_L,
                    read_len,
                    scratch_DIR,
                    temp,
                    fiveP_adapter,
                    threeP_adapter )

                # Ds_by_read_D = {
                #           "dot_struc": ".....(((((.(((((.((((((((((..(....)..))))))).))).))))))))..)).",
                #           "dotbracket_str": "fffffsssssisssssissssssssssiishhhhsiisssssssisssissssssssiissf"}
                Ds_by_read_D = get_MFE_DotBracket_and_symbol_over_reads(
                    reads_this_set_L,
                    scratch_DIR,
                    temp,
                    fiveP_adapter,
                    threeP_adapter )

                #### Now go through each of the seqs
                with gzip.open( scratch_reads_F, "ab" ) as out_f:
                    for read_idx, rand_read in enumerate( reads_this_set_L ):

                        read_w_adapt = fiveP_adapter + rand_read + threeP_adapter

                        avg_bp_probs_L = avg_bp_probs_Ls_L[read_idx]
                        out_Ppaired_str = " ".join( ["{0:.3f}".format( Ppaired )\
                                for Ppaired in avg_bp_probs_L] )

                        try:
                            dotstruc_D = Ds_by_read_D[read_w_adapt]

                            this_line = "\n".join( [
                                read_w_adapt,
                                out_Ppaired_str,
                                dotstruc_D['dot_struc'],
                                dotstruc_D['dotbracket_str']] ) + "\n"
                            out_f.write( this_line )
                        except KeyError:
                            pass

                reads_this_set_L = []
                num_reads_this_set = 0

                if ( num_reads_covered in [5000, 10000, 50000] + [i*100000 for i in range( 1, 10 )] ):

                    curr_time = time.time()

                    secs_per_million_reads = (curr_time - start_time) * 1000000. / num_reads_covered
                    readable_string_secs_per_mil = RBNS_utils.return_human_readable_time_from_secs(
                            secs_per_million_reads )
                    with open( log_F, "a" ) as log_f:
                        curr_pprint_str = RBNS_utils.return_nice_datetime_str_for_filename()
                        curr_str = "{0:,} at {1} for {2}: {3} sec per million rds\n".format(
                            num_reads_covered,
                            curr_pprint_str,
                            start_basename,
                            readable_string_secs_per_mil )
                        print curr_str
                        log_f.write( curr_str )

                #### If it's been more than copy_back_every_x_seconds
                ####    since the last copy, do it now
                #if ( ( time.time() - last_copy_back_time ) >= copy_back_every_x_seconds ):

                #    with open( log_F, "a" ) as log_f:
                #        curr_pprint_str = RBNS_utils.return_nice_datetime_str_for_filename()
                #        curr_str = "\tCopying {0} to {1} w/ {2} reads at {3}\n".format(
                #                scratch_reads_F,
                #                out_reads_F,
                #                num_reads_covered,
                #                curr_pprint_str )
                #        print curr_str
                #        log_f.write( curr_str )

                #    try:
                #        shutil.copy( out_reads_F, prev_copies_DIR )
                #    except IOError:
                #        pass

                #    shutil.copy( scratch_reads_F, out_reads_F )
                #    last_copy_back_time = time.time()


    shutil.copyfile( scratch_reads_F, out_reads_F )
    with open( log_F, "a" ) as log_f:
        total_sec = time.time() - start_time
        curr_pprint_str = RBNS_utils.return_nice_datetime_str_for_filename()
        curr_str = "FINISHED SUCCESSFULLY! Copied {0:,} reads to {1} at {2}\n\tTook {3:,} sec.".format(
                num_reads_covered,
                out_reads_F,
                curr_pprint_str,
                int( total_sec ) )
        print curr_str
        log_f.write( curr_str )

    os.system( "rm -rf {}".format( scratch_DIR ) )
    making_F = out_reads_F + ".making"
    if os.path.exists( making_F ):
        os.system( "rm {}".format( making_F ) )












def get_suboptimal_sampled_DotBracket_reads_F(
        in_struct_gz_F,
        temp,
        starting_scratch_DIR,
        num_subopt_to_get = 20,
        num_reads_to_get_status_after = 1000,
        copy_back_every_x_seconds = 10000000 ):
    """
    - For an in_struct_gz_F that already has 4 lines per read (e.g.,
            RBFOX3_input.w_struc.reads.gz
            RBFOX3_20.w_struc.reads.gz, etc.),
        will get each read in in_struct_gz_F and sample num_subopt_to_get
        suboptimal structures from the thermodynamic ensemble

    - Makes an output file with num_subopt_to_get + 1 lines per read
        (the read, followed by the num_subopt_to_get structures)
    """
    temp = int( temp )

    start_basename = os.path.basename( in_struct_gz_F ).split(".reads")[0]
    out_DIR = os.path.join( os.path.dirname( in_struct_gz_F ),
            "subopt_DB_{0}reads".format( num_subopt_to_get ) )
    out_logs_DIR = os.path.join( out_DIR, "logs" )
    os.system( "mkdir -p {}".format( out_logs_DIR ) )

    out_reads_F = os.path.join( out_DIR,
            "{}.subopt_DB.gz".format( start_basename ) )
    if os.path.exists( out_reads_F ):
        print "\n\t{} already exists - skipping".format( out_reads_F )
        return

    scratch_DIR = os.path.join( starting_scratch_DIR, "{0}_{1}".format(
            start_basename, random.randint( 0, 100000 ) ) )
    os.system( "mkdir -p {}".format( scratch_DIR ) )

    #### copy over the reads to the scratch space
    scratch_in_reads_F = os.path.join( scratch_DIR,
            os.path.basename( in_struct_gz_F ) )
    print "\n\tCopying {0} to {1}".format( in_struct_gz_F, scratch_in_reads_F )

    shutil.copyfile( in_struct_gz_F, scratch_in_reads_F )
    print "\t\tDONE"

    out_F_to_append_results_to = os.path.join( scratch_DIR,
            "{}.subopt_DB.gz".format( start_basename ) )
    log_F = os.path.join( out_logs_DIR,
            "{}.subopt_DB.log.txt".format( start_basename ) )

    with open( log_F, "w" ) as log_f:
        log_f.write( "Starting at: {}\n".format(
            RBNS_utils.return_nice_datetime_str_for_filename() ) )

        log_f.write( "hostname is: {}\n".format( socket.gethostname() ) )

        log_f.write( "\n\n" + "="*80 + "\n" + "="*80 + "\n\n" )

        log_f.write( "in_struct_gz_F: {}\n".format( in_struct_gz_F ) )
        log_f.write( "temp: {}\n".format( temp ) )

        log_f.write( "\n\n" + "="*80 + "\n" + "="*80 + "\n\n" )

    num_reads_covered = 0
    start_time = time.time()
    last_copy_back_time = time.time()

    reads_this_set_L = []
    num_reads_this_set = 0
    for lines_L in RBNS_utils.iterNlines(
            scratch_in_reads_F, 4, strip_newlines = True ):

        rand_RNA_seq = lines_L[0]
        reads_this_set_L.append( rand_RNA_seq )
        num_reads_this_set += 1
        num_reads_covered += 1

        if ( num_reads_this_set == num_reads_to_get_status_after ):

            get_subopt_folding_of_reads(
                reads_this_set_L,
                scratch_DIR,
                temp,
                out_F_to_append_results_to,
                num_to_return_for_each_read = num_subopt_to_get )

            reads_this_set_L = []
            num_reads_this_set = 0

    get_subopt_folding_of_reads(
        reads_this_set_L,
        scratch_DIR,
        temp,
        out_F_to_append_results_to,
        num_to_return_for_each_read = num_subopt_to_get )

    shutil.copyfile( out_F_to_append_results_to, out_reads_F )
    with open( log_F, "a" ) as log_f:
        total_sec = time.time() - start_time
        curr_pprint_str = RBNS_utils.return_nice_datetime_str_for_filename()
        curr_str = "\nFINISHED SUCCESSFULLY! Copied {0:,} reads to {1} at {2}\n\tTook {3:,} seconds\n".format(
                num_reads_covered,
                out_reads_F,
                curr_pprint_str,
                int( total_sec ) )
        print curr_str
        log_f.write( curr_str )

    os.system( "rm -rf {}".format( scratch_DIR ) )
    making_F = os.path.join( os.path.dirname( out_reads_F ),
        "{0}.subopt_DB.gz.making".format(
            os.path.basename( in_struct_gz_F ).split(".reads")[0] ) )
    if os.path.exists( making_F ):
        os.system( "rm {}".format( making_F ) )








def calc_Ppaired_over_top_enriched_kmers_and_flanking(
        reads_w_struct_F,
        k,
        fiveP_adapter,
        threeP_adapter,
        random_read_len,
        num_bins = 5 ):
    """
    - For an input reads_struct_F like:
            RBFOX3_input.w_struc.reads.gz
            RBFOX3_20.w_struc.reads.gz,

        gets all occurrences of each of the kmers and
        calculates the Ppaired over each position of the motif & 10 bases
        flanking it upstream & downstream, as well as the average Ppaired
        over each motif occurrence (i.e., which of the num_bins Ppaired bins
        it should go into for later calculating the R by Ppaired bin)

    - Pickles an output dictionary in out_Ds_DIR for later loading & analysis
    """
    assert( num_bins in [5, 10] )
    starting_basename = os.path.basename( reads_w_struct_F ).split('.w_st')[0]

    fiveP_len = len( fiveP_adapter )
    threeP_len = len( threeP_adapter )

    out_Ds_DIR = os.path.join( os.path.dirname( reads_w_struct_F ),
            'Ppaired_Ds', str( k ) )
    RBNS_utils.make_dir( out_Ds_DIR )

    out_D_F = os.path.join( out_Ds_DIR, "{0}.D.pkl".format( starting_basename ) )
    if os.path.exists( out_D_F ):
        return

    ##### Make sure that the adapter lengths & random read length match up
    for lines_L in RBNS_utils.iterNlines( reads_w_struct_F, 4, strip_newlines = True ):

        read_w_adapter = lines_L[0]
        calculated_random_read_len = len( read_w_adapter ) - fiveP_len - threeP_len
        assert( calculated_random_read_len == random_read_len )
        break
    random_idx_L = range( random_read_len )

    upper_index_of_random = random_read_len + fiveP_len
    num_kmers_each_read = random_read_len - k + 1

    D = {"num_reads": 0,
        "Ppair_and_count_by_kmer_idx_D": {},
        "counts_by_kmer_binidx_D": {} }

    for kmer in RBNS_utils.yield_kmers( k ):
        D["counts_by_kmer_binidx_D"][kmer] = {}
        for i in range( num_bins ):
            D["counts_by_kmer_binidx_D"][kmer][i] = 0
        D["Ppair_and_count_by_kmer_idx_D"][kmer] = {}
        for idx in range( -10, k + 10 ):
            D["Ppair_and_count_by_kmer_idx_D"][kmer][idx] = {'counts': 0, 'Ppaired_sum': 0.}

    for lines_L in RBNS_utils.iterNlines( reads_w_struct_F, 4, strip_newlines = True ):

        read_w_adapter = lines_L[0]
        random_seq = read_w_adapter[fiveP_len:upper_index_of_random]
        seq_L = [x for x in random_seq]

        Ppaired_L = lines_L[1].split(" ")
        pruned_Ppaired_L = [float(x) for x in Ppaired_L[fiveP_len:upper_index_of_random]]

        seq_Ppaired_T_L = zip( seq_L, pruned_Ppaired_L )

        D["num_reads"] += 1

        for start_idx in range( num_kmers_each_read ):

            #### Get the kmers in the read
            kmer = random_seq[start_idx:(start_idx + k)]
            Ppaired_kmer_L = pruned_Ppaired_L[start_idx:(start_idx + k)]

            mean_Ppaired = np.mean( Ppaired_kmer_L )
            if ( num_bins == 5 ):
                bin_idx = get_bin_of_5_from_mean_Ppaired( mean_Ppaired )
            elif ( num_bins == 10 ):
                bin_idx = get_bin_of_10_from_mean_Ppaired( mean_Ppaired )

            D["counts_by_kmer_binidx_D"][kmer][bin_idx] += 1

            ##### Go through and get all of the Ppaired flanking
            for rel_idx in range( -10, 10 + k ):

                this_idx = start_idx + rel_idx
                if this_idx in random_idx_L:
                    D["Ppair_and_count_by_kmer_idx_D"][kmer][rel_idx]['counts'] += 1
                    D["Ppair_and_count_by_kmer_idx_D"][kmer][rel_idx]['Ppaired_sum'] +=\
                            pruned_Ppaired_L[this_idx]

    ##### Pickle to out_D_F
    RBNS_utils.pkl_with_formatfile( D, out_D_F )







def plot_RBNS_Ppaired_ratio_w_sig(
        readswstruct_startingbasename_myannot_L,
        read_len,
        k,
        effective_R_D,
        kmers_to_do = "top_10" ):
    """
    - Makes a plot of the Ppaired ratios for the desired set of kmers

    - fld_CG_match_DIR is the directory that contains the Ppaired_Ds directory:
        /net/eofe-data010/data001/burgelab/nevermind/data/nm/pfreese/RBFOX3_test/split_reads/fld_CG_match

    """
    import random
    import RBNS_plots

    assert( kmers_to_do in ["top_10"] )

    fld_CG_match_DIR = os.path.dirname( readswstruct_startingbasename_myannot_L[0][0] )
    RBP = readswstruct_startingbasename_myannot_L[0][1].split('_')[0]

    #### Get the list of all kmers, which to motif_num matches
    all_kmers_L = RBNS_utils.return_all_kmers_L( k )
    kmer_to_motifidx_D = {}
    motifidx_to_kmer_D = {}
    for idx, kmer in enumerate( all_kmers_L ):
        kmer_to_motifidx_D[kmer] = idx
        motifidx_to_kmer_D[idx] = kmer

    out_DIR = fld_CG_match_DIR.split("/split_reads")[0]
    os.system( "mkdir -p {}".format( out_DIR ) )

    out_DIR_this_RBP_k = os.path.join( out_DIR, "{0}mer_plots".format( k ) )

    out_Ds_DIR = os.path.join( out_DIR_this_RBP_k, 'Ds' )
    os.system( "mkdir -p {}".format( out_Ds_DIR ) )

    tables_DIR = os.path.join( out_DIR_this_RBP_k, 'tables' )
    os.system( "mkdir -p {}".format( tables_DIR ) )

    out_F_start = os.path.join( out_DIR_this_RBP_k, RBP )

    #### Get the most enriched concentration
    most_enriched_conc_str = ""
    for T in readswstruct_startingbasename_myannot_L:
        if ( T[2] == 'Most enriched' ):
            most_enriched_conc_str = "{0} nM".format( T[1].split("_")[-1] )

    Ds_DIR = os.path.join( fld_CG_match_DIR, "Ppaired_Ds/{0}".format( k ) )

    most_enriched_lib_annotation = ""
    #### annots_L will be like:
    ####    ['input', '5_nM', '20_nM', '80_nM', '320_nM', '1300_nM']
    annots_L = []
    ####    D_by_annot_D will have keys like 'input', '5_nM', etc. and values:
    ## {'AAAAAA': {-10: {'Ppaired_sum': 2.061,
    ##                    'counts': 4},
    ##              -9: {'Ppaired_sum': 1.027,
    ##                    'counts': 5},
    D_by_annot_D = {}
    for reads_w_struct_F, starting_basename, my_annot in\
        readswstruct_startingbasename_myannot_L:

        if ( my_annot == "Input" ):
            lib_annot = 'input'
        else:
            lib_annot = starting_basename.split("_")[-1] + "_nM"
        if ( my_annot == "Most enriched" ):
            most_enriched_lib_annotation = lib_annot

        annots_L.append( lib_annot )

        D_F = os.path.join( Ds_DIR, "{0}.D.pkl".format( starting_basename ) )
        D_by_annot_D[lib_annot] = pickle.load( open( D_F ) )['Ppair_and_count_by_kmer_idx_D']

    #### Get the desired kmers (e.g., the top 5 )
    if ( kmers_to_do == "top_10" ):
        kmer_R_T_L = [(kmer, effective_R_D[kmer]) for kmer in effective_R_D]
        kmer_R_T_L.sort( key = lambda x: -1 * x[1] )
        top_kmers_L = [tupl[0] for tupl in kmer_R_T_L[:1000]]

    #top_kmers_L = RBNS_exp.return_top_X_kmers( k, int( 4 ** k ) )
    for kmer_idx, kmer_to_plot in enumerate( top_kmers_L ):

        R = effective_R_D[kmer_to_plot]
        title = r"{0}, {1} (\#{2}: $R={3:.2f}$)".format(
                RBP, kmer_to_plot.replace("T","U"),
                kmer_idx + 1, R )
        print title

        ##### Make an output .txt table of the Ppaireds of this motif and
        ####    the Ppaired ratio for the upstream & downstream flanking
        ####    positions
        out_txt_F = os.path.join( tables_DIR,
                "{0}.{1}.Ppaired_w_flanking_ratios.txt".format(
                    RBP, kmer_to_plot.replace("T","U") ) )
        out_txt_f = open( out_txt_F, 'w' )
        out_txt_f.write( "{0} {1} Ppaired, {2}mer plus 10 flanking".format(
            RBP, kmer_to_plot.replace("T","U"), k ) )
        for idx in range( -10, 10 + k ):
            out_txt_f.write( "\t{}".format( idx ) )

        Ppaired_L_by_libannot_D = {}
        to_plot_avg_fold_probs_by_motif_and_pos_D = {}
        for lib_annot, D in D_by_annot_D.iteritems():

            to_plot_avg_fold_probs_by_motif_and_pos_D[lib_annot] = {}
            Ppaired_L_by_libannot_D[lib_annot] = []
            #### The original motif_num of this kmer
            motif_num = kmer_to_motifidx_D[kmer_to_plot]

            #### Since we're only plotting one motif, give it index 0
            to_plot_avg_fold_probs_by_motif_and_pos_D[lib_annot][0] = []
            Ppaired_D = D[kmer_to_plot]
            #pprint.pprint( Ppaired_D )
            for idx in range( k ):
                Ppaired = Ppaired_D[idx]['Ppaired_sum'] / Ppaired_D[idx]['counts']
                to_plot_avg_fold_probs_by_motif_and_pos_D[lib_annot][0].append( Ppaired )
                Ppaired_L_by_libannot_D[lib_annot].append( Ppaired )

            ##### Go through the 10 flanking upstream & downstream positions
            out_txt_f.write( "\n{}".format( lib_annot ) )
            for idx in range( -10, 10 + k ):
                Ppaired = Ppaired_D[idx]['Ppaired_sum'] / Ppaired_D[idx]['counts']
                out_txt_f.write( "\t{0:.3f}".format( Ppaired ) )

        print "\nSaving to: {}".format( out_F_start )
        returned_D = RBNS_plots.plot_Ppaired_ratio_of_motif(
                to_plot_avg_fold_probs_by_motif_and_pos_D,
                annots_L,
                [kmer_to_plot],
                most_enriched_lib_annotation,
                out_F_start,
                read_len,
                title = title,
                #plot_signif = True,
                plot_signif = False,
                #skip_if_already_exists = True )
                skip_if_already_exists = False )
        ratios_by_motif_pos_D = returned_D['ratios_by_motif_pos_D']

        out_ratio_D_F = os.path.join( out_Ds_DIR,
                "{}.ratios_by_pos_D.pkl".format( kmer_to_plot ) )
        try:
            D = ratios_by_motif_pos_D[kmer_to_plot]
            RBNS_utils.pkl_with_formatfile(
                    D,
                    out_ratio_D_F,
                    num_to_include_in_format = "all" )
        except KeyError:
            pass

        out_D_F = os.path.join( out_Ds_DIR,
                "{}.Ppaired_L_by_libannot_D.pkl".format( kmer_to_plot ) )
        RBNS_utils.pkl_with_formatfile(
                Ppaired_L_by_libannot_D,
                out_D_F,
                num_to_include_in_format = "all" )

        if ( 'sigB_by_motif_pos_D' in returned_D ):
            sigB_by_motif_pos_D = returned_D['sigB_by_motif_pos_D']
            out_sigB_D_F = os.path.join( out_Ds_DIR,
                    "{}.sigB_by_pos_D.pkl".format( kmer_to_plot ) )
            try:
                D = sigB_by_motif_pos_D[kmer_to_plot]
                RBNS_utils.pkl_with_formatfile(
                        D,
                        out_sigB_D_F,
                        num_to_include_in_format = "all" )
            except KeyError:
                pass

        out_txt_f.close()




def plot_R_by_Ppaired_bin_w_sig(
        readswstruct_startingbasename_myannot_L,
        read_len,
        k,
        effective_R_D,
        kmers_to_do = "top_10" ):
    """
    - Makes a plot of the Ppaired ratios for the desired set of kmers

    - fld_CG_match_DIR is the directory that contains the Ppaired_Ds directory:
        /net/eofe-data010/data001/burgelab/nevermind/data/nm/pfreese/RBFOX3_test/split_reads/fld_CG_match

    """
    import random
    import RBNS_plots


    assert( kmers_to_do in ["top_10"] )

    fld_CG_match_DIR = os.path.dirname( readswstruct_startingbasename_myannot_L[0][0] )
    RBP = readswstruct_startingbasename_myannot_L[0][1].split('_')[0]

    #### Get the list of all kmers, which to motif_num matches
    all_kmers_L = RBNS_utils.return_all_kmers_L( k )
    kmer_to_motifidx_D = {}
    motifidx_to_kmer_D = {}
    for idx, kmer in enumerate( all_kmers_L ):
        kmer_to_motifidx_D[kmer] = idx
        motifidx_to_kmer_D[idx] = kmer

    out_DIR = fld_CG_match_DIR.split("/split_reads")[0]
    os.system( "mkdir -p {}".format( out_DIR ) )

    out_DIR_this_RBP_k = os.path.join( out_DIR, "{0}mer_plots".format( k ) )

    out_Ds_DIR = os.path.join( out_DIR_this_RBP_k, 'Ds' )
    os.system( "mkdir -p {}".format( out_Ds_DIR ) )

    tables_DIR = os.path.join( out_DIR_this_RBP_k, 'tables' )
    os.system( "mkdir -p {}".format( tables_DIR ) )

    out_F_start = os.path.join( out_DIR_this_RBP_k, RBP )

    #### Get the most enriched concentration
    most_enriched_conc_str = ""
    for T in readswstruct_startingbasename_myannot_L:
        if ( T[2] == 'Most enriched' ):
            most_enriched_conc_str = "{0} nM".format( T[1].split("_")[-1] )

    Ds_DIR = os.path.join( fld_CG_match_DIR, "Ppaired_Ds/{0}".format( k ) )

    most_enriched_lib_annotation = ""
    #### annots_L will be like:
    ####    ['input', '5_nM', '20_nM', '80_nM', '320_nM', '1300_nM']
    annots_L = []
    ####    D_by_annot_D will have keys like 'input', '5_nM', etc. and values:
    ## {'AAAAAA': {-10: {'Ppaired_sum': 2.061,
    ##                    'counts': 4},
    ##              -9: {'Ppaired_sum': 1.027,
    ##                    'counts': 5},
    D_by_annot_D = {}
    for reads_w_struct_F, starting_basename, my_annot in\
        readswstruct_startingbasename_myannot_L:

        if ( my_annot == "Input" ):
            lib_annot = 'input'
        else:
            lib_annot = starting_basename.split("_")[-1] + "_nM"
        if ( my_annot == "Most enriched" ):
            most_enriched_lib_annotation = lib_annot

        annots_L.append( lib_annot )

        D_F = os.path.join( Ds_DIR, "{0}.D.pkl".format( starting_basename ) )
        D_by_annot_D[lib_annot] = pickle.load( open( D_F ) )['counts_by_kmer_binidx_D']

    #### Get the desired kmers (e.g., the top 5 )
    if ( kmers_to_do == "top_10" ):
        kmer_R_T_L = [(kmer, effective_R_D[kmer]) for kmer in effective_R_D]
        kmer_R_T_L.sort( key = lambda x: -1 * x[1] )
        top_kmers_L = [tupl[0] for tupl in kmer_R_T_L[:1000]]

    #top_kmers_L = RBNS_exp.return_top_X_kmers( k, int( 4 ** k ) )
    for kmer_idx, kmer_to_plot in enumerate( top_kmers_L ):

        R = effective_R_D[kmer_to_plot]
        title = r"{0}, {1} (\#{2}: $R={3:.2f}$)".format(
                RBP, kmer_to_plot.replace("T","U"),
                kmer_idx + 1, R )
        print title

        ##### Make an output .txt table of the Ppaireds of this motif and
        ####    the Ppaired ratio for the upstream & downstream flanking
        ####    positions
        out_txt_F = os.path.join( tables_DIR,
                "{0}.{1}.R_by_Ppaired_bin.txt".format(
                    RBP, kmer_to_plot.replace("T","U") ) )
        out_txt_f = open( out_txt_F, 'w' )
        out_txt_f.write( "{0} {1} R by Ppaired bin".format(
            RBP, kmer_to_plot.replace("T","U") ) )
        out_txt_f.write( "\t0-0.2\t0.2-0.4\t0.4-0.6\t0.6-0.8\t0.8-1.0" )

        ##### First get the INPUT frequency in each of the 5 Ppaired bins
        input_D = D_by_annot_D['input']
        input_freq_by_bin_D = {}
        input_kmer_counts_all_bins = 0.
        input_all_counts_all_bins = 0.
        for bin_idx in range( 5 ):
            total_counts_this_bin = sum( [input_D[all_kmer][bin_idx] for all_kmer in input_D] )
            kmer_counts_this_bin = input_D[kmer_to_plot][bin_idx]
            kmer_freq = float( kmer_counts_this_bin ) / total_counts_this_bin
            input_freq_by_bin_D[bin_idx] = kmer_freq
            input_kmer_counts_all_bins += kmer_counts_this_bin
            input_all_counts_all_bins += total_counts_this_bin

        enrichments_by_kmer_conc_bin_D = {kmer_to_plot: {}}
        for lib_annot, D in D_by_annot_D.iteritems():

            if ( lib_annot == 'input' ):
                continue

            kmer_counts_all_bins = 0.
            all_counts_all_bins = 0.

            enrichments_by_kmer_conc_bin_D[kmer_to_plot][lib_annot] = {}
            #### The original motif_num of this kmer
            motif_num = kmer_to_motifidx_D[kmer_to_plot]

            this_D = D_by_annot_D[lib_annot]

            out_txt_f.write( "\n{}".format( lib_annot ) )
            for bin_idx in range( 5 ):
                total_counts_this_bin = sum( [this_D[all_kmer][bin_idx] for all_kmer in this_D] )
                kmer_counts_this_bin = this_D[kmer_to_plot][bin_idx]
                kmer_freq = float( kmer_counts_this_bin ) / total_counts_this_bin
                kmer_counts_all_bins += kmer_counts_this_bin
                all_counts_all_bins += total_counts_this_bin

                try:
                    kmer_R = kmer_freq / input_freq_by_bin_D[bin_idx]
                except ZeroDivisionError:
                    kmer_R = 1.
                out_txt_f.write( "\t{0:.3f}".format( kmer_R ) )

                enrichments_by_kmer_conc_bin_D[kmer_to_plot][lib_annot][bin_idx] = kmer_R

            ##### Get the OVERALL (over all bins) R
            overall_R = (kmer_counts_all_bins/all_counts_all_bins) / (input_kmer_counts_all_bins/input_all_counts_all_bins)
            enrichments_by_kmer_conc_bin_D[kmer_to_plot][lib_annot]['overall'] = kmer_R


        print "\nSaving to: {}".format( out_F_start )
        Ppaired_upper_bins_L = [0.2, 0.4, 0.6, 0.8, 1.]
        returned_D = RBNS_plots.plot_enrichment_by_5_Ppaired_bins(
                enrichments_by_kmer_conc_bin_D,
                annots_L,
                [kmer_to_plot],
                Ppaired_upper_bins_L,
                out_F_start,
                read_len,
                title = title,
                #plot_signif = True,
                plot_signif = False )
        sigB_by_kmer_conc_bin_D = returned_D['sigB_by_kmer_conc_bin_D']

        out_ratio_D_F = os.path.join( out_Ds_DIR,
                "{}.sigB_by_kmer_conc_bin_D.pkl".format( kmer_to_plot ) )
        try:
            D = sigB_by_kmer_conc_bin_D[kmer_to_plot]
            RBNS_utils.pkl_with_formatfile(
                    D,
                    out_ratio_D_F,
                    num_to_include_in_format = "all" )
        except KeyError:
            pass

        out_D_F = os.path.join( out_Ds_DIR,
                "{}.enrichments_by_conc_bin_D.pkl".format( kmer_to_plot ) )
        RBNS_utils.pkl_with_formatfile(
                enrichments_by_kmer_conc_bin_D[kmer_to_plot],
                out_D_F,
                num_to_include_in_format = "all" )

        out_txt_f.close()







###############################################################################
#################################### < UTILS > ################################




def get_avg_bp_prob_with_flankingprimers(
        RNA_seq,
        read_length,
        tmp_DIR,
        fiveP_adapter,
        threeP_adapter,
        temp = 37 ):
    """
    - Calculates the average bp. probability for each base of RNA_seq (which
        does not have the 5' / 3' adapter sequence flanking it)
    - tmp_DIR is an already-made scratch directory for temporary file I/0
    - temp is degrees C which will be passed into the RNAfold command
    - fiveP_adapter & threeP_adapter are the 5' and 3' sequencing adapters that
        will be added to & folded along with the RNA_seq
    """
    # Get the primer lengths
    fiveP_adapter_len = len( fiveP_adapter )
    threeP_adapter_len = len( threeP_adapter )

    # First, make a temporary .fa file of the read
    tmp_read_fasta_F = os.path.join( tmp_DIR, "read.fa" )
    with open( tmp_read_fasta_F, "w" ) as f:
        f.write( ">read\n{0}".format(
            fiveP_adapter + RNA_seq + threeP_adapter ) )
    os.system("cd {0}".format( tmp_DIR ))
    # make the .ps files
    os.chdir( tmp_DIR )
    fold_CMD = "RNAfold -p --temp={0} < {1} > /dev/null".format(
            temp, tmp_read_fasta_F )
    fold = subprocess.Popen(fold_CMD, shell=True)
    stdoutdata, stderrdata = fold.communicate()
    # get pos
    probs_ps_F = os.path.join( tmp_DIR, "read_dp.ps" )
    avg_bp_probs_L = parse_probs_ps_F(
            probs_ps_F,
            fiveP_adapter_len + read_length + threeP_adapter_len )

    avg_bp_probs_randomregion_L =\
        avg_bp_probs_L[fiveP_adapter_len:(-1*threeP_adapter_len)]

    avg_bp_probs_D = {"all_bp_probs_L": avg_bp_probs_L,
            "randomregion_bp_probs_L": avg_bp_probs_randomregion_L}

    return avg_bp_probs_D






def parse_probs_ps_F(
        probs_ps_F,
        read_length ):
    """
    - Parses probs_ps_F, returning a list of length read_length,
        which contains the posterior probability of each base being
        paired
            - NOTE that the base numbers in probs_ps_F are 1 higher than the
                list indices in probs_ps_F, so there must be a change-by-1
    """
    avg_bp_probs_L = [0.] * read_length

    with open( probs_ps_F ) as f:
        start_parsing = False
        for line in f:
            ln = line.strip()
            if not start_parsing and (ln == "%start of base pair probability data"):
                start_parsing = True
            elif start_parsing:
                try:
                    bases_and_probs_L = ln.split(" ")[:4]
                    # it's only a probability if the 4th column is "ubox" (not "lbox")
                    if ( bases_and_probs_L[3] == "ubox" ):
                        # subtract 1 when going from base in the .ps file to the list
                        # index
                        try:
                            prob = (float( bases_and_probs_L[2] ) ** 2)
                            avg_bp_probs_L[ int( bases_and_probs_L[0] ) - 1 ] += prob
                            avg_bp_probs_L[ int( bases_and_probs_L[1] ) - 1 ] += prob
                        except IndexError:
                            pass
                    elif ( bases_and_probs_L[3] == "lbox" ):
                        return avg_bp_probs_L
                except IndexError:
                    continue
            if ( ln == "showpage" ):
                return avg_bp_probs_L
                start_parsing = False
    return avg_bp_probs_L









def get_avg_bp_probs_w_adapters(
        RNA_seqs_L,
        read_length,
        tmp_DIR,
        temperature,
        fiveP_adapter,
        threeP_adapter ):
    """
    - Calculates the average bp. probability for each base for each of the
        strings in RNA_seqs_L

    INPUTS:
        - RNA_seqs_L: a list of sequences (e.g.: ["AACAA","CCAAG"])
        - read_length: the length of each of the RNAs in RNA_seqs_L

    RETURNS:
        - avg_bp_probs_Ls_L: a list of length len(RNA_seqs_L), which contains
            lists of the average bp. probability for each base of that seq
    """

    avg_bp_probs_Ls_L = []
    len_5p = len( fiveP_adapter )
    len_3p = len( threeP_adapter )

    # First, make a temporary .fa file of the read
    tmp_read_fasta_F = os.path.join( tmp_DIR, "read.fa" )
    with open( tmp_read_fasta_F, "w" ) as f:
        for seq_num, RNA_seq in enumerate( RNA_seqs_L ):
            f.write( ">{0}\n{1}\n".format(
                seq_num,
                fiveP_adapter + RNA_seq + threeP_adapter ) )
    os.system("cd {0}".format( tmp_DIR ))
    # make the .ps files
    os.chdir( tmp_DIR )
    fold_CMD = "RNAfold -p --temp={0} < {1} > /dev/null".format(
            temperature, tmp_read_fasta_F )
    fold = subprocess.Popen(fold_CMD, shell=True)
    stdoutdata, stderrdata = fold.communicate()

    for seq_num, RNA_seq in enumerate( RNA_seqs_L ):
        probs_ps_F = os.path.join( tmp_DIR, "{0}_dp.ps".format( seq_num ) )
        avg_bp_probs_L = parse_probs_ps_F( probs_ps_F,
                read_length + len_5p + len_3p )
        avg_bp_probs_Ls_L.append( avg_bp_probs_L )

    return avg_bp_probs_Ls_L









def get_MFE_DotBracket_and_symbol_over_reads(
        random_reads_L,
        scratch_DIR,
        temp,
        fiveP_adapter,
        threeP_adapter ):
    """
    - Returns a dictionary like:
        return_D = {"seq_w_adapters": seq_w_adapters,
                "dot_struc": dot_struc,
                "dotbracket_str": dotbracket_str }
    """
    fa_to_fold_F = os.path.join( scratch_DIR, "rd.fa" )
    with open( fa_to_fold_F, "w" ) as fa_f:

        for read_idx, random_read in enumerate( random_reads_L ):

            seq_w_adapters = fiveP_adapter + random_read + threeP_adapter

            fa_f.write( ">" + str( read_idx ) + "\n" + seq_w_adapters + "\n" )

    #### Now fold all of the reads
    out_MFE_F = os.path.join( scratch_DIR, "rd.folded.fa" )

    # Get the MFE structure for each of the reads
    os.chdir( scratch_DIR )
    lines_L = subprocess.check_output( "RNAfold -T {0} < {1}".format( temp, fa_to_fold_F ), stderr=subprocess.STDOUT, shell = True ).split("\n")

    Ds_by_read_D = {}
    for idx in range( len( random_reads_L ) ):
        #### lower = 3 * idx; then 3*idx + 1; 3*idx + 2
        read_w_adapters = lines_L[(3*idx)+1].replace("U","T")
        dot_struc = lines_L[(3*idx)+2].split( " " )[0]
        dotbracket_str = get_elementstring_from_DotBracket( dot_struc ).replace( "t", "f" )
        Ds_by_read_D[read_w_adapters] = {
                "dot_struc": dot_struc,
                "dotbracket_str": dotbracket_str }

    return Ds_by_read_D





def get_elementstring_from_DotBracket( DotBracket_str ):
    """
    INPUT:
        DotBracket_string =
        (((((((((...((((((.........))))))........((((((.......))))))..)))))))))
    RETURNS:
        sssssssssmmmsssssshhhhhhhhhssssssmmmmmmmmsssssshhhhhhhssssssmmsssssssss
    """
    bg = cgb.BulgeGraph()
    bg.from_dotbracket( DotBracket_str )
    return bg.to_element_string()








def get_subopt_folding_of_reads(
        reads_L,
        scratch_DIR,
        temp,
        out_F_to_append_results_to,
        num_to_return_for_each_read = 20 ):
    """
    - Given a list of reads (with adapters) reads_L, will get
        num_to_return_for_each_read suboptimal DotBracket structures sampled
        with probabilities equal to their Boltzmann weights
    """
    tmp_read_fasta_F = os.path.join( scratch_DIR, "reads.fa" )
    out_F = os.path.join( scratch_DIR, "reads.out.txt" )

    read_by_readwindex_D = {}
    with open( tmp_read_fasta_F, "w" ) as f:
        for idx, read in enumerate( reads_L ):
            read_w_index = "read{0}".format( idx )
            read_by_readwindex_D[read_w_index] = read
            f.write( ">read{0}\n{1}\n".format( idx, read ) )

    os.chdir( scratch_DIR )
    fold_CMD = "RNAsubopt --temp={0} --stochBT={1} < {2} > {3}".format(
        temp, num_to_return_for_each_read, tmp_read_fasta_F, out_F )
    # make the .ps files
    fold = subprocess.Popen( fold_CMD, shell = True )
    stdoutdata, stderrdata = fold.communicate()

    #### Now go through all of the reads
    this_read = ""
    out_f_to_append_results_to = gzip.open( out_F_to_append_results_to, 'ab' )
    for lines_L in RBNS_utils.iterNlines( out_F,
            num_to_return_for_each_read + 2,
            strip_newlines = True ):

        this_read = lines_L[1]
        out_f_to_append_results_to.write( ">" + this_read + "\n" )

        #### Now go through all of the num_to_return_for_each_read DotBracket
        ####    structures
        for DB_str in lines_L[2:]:
            element_string = get_elementstring_from_DotBracket( DB_str )
            out_f_to_append_results_to.write(
                    element_string + "\n" )

    out_f_to_append_results_to.close()



def get_bin_of_5_from_mean_Ppaired( mean_Ppaired ):
    if ( mean_Ppaired <= 0.2 ):
        return 0
    elif ( mean_Ppaired <= 0.4 ):
        return 1
    elif ( mean_Ppaired <= 0.6 ):
        return 2
    elif ( mean_Ppaired <= 0.8 ):
        return 3
    else:
        return 4



def get_bin_of_10_from_mean_Ppaired( mean_Ppaired ):
    if ( mean_Ppaired <= 0.1 ):
        return 0
    elif ( mean_Ppaired <= 0.2 ):
        return 1
    elif ( mean_Ppaired <= 0.3 ):
        return 2
    elif ( mean_Ppaired <= 0.4 ):
        return 3
    elif ( mean_Ppaired <= 0.5 ):
        return 4
    elif ( mean_Ppaired <= 0.6 ):
        return 5
    elif ( mean_Ppaired <= 0.7 ):
        return 6
    elif ( mean_Ppaired <= 0.8 ):
        return 7
    elif ( mean_Ppaired <= 0.9 ):
        return 8
    else:
        return 9

################################# </ UTILS > ##################################
###############################################################################







if __name__ == '__main__':

    fxn = sys.argv[1].strip('-')
    args = ['"' + arg + '"' for arg in sys.argv[2:]]
    python_command = fxn + '(' + ','.join( args ) + ')'

    print python_command
    eval( python_command )







