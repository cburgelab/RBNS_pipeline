#!/usr/bin/env python
import os, sys, time, gzip
import subprocess, pprint, glob
import random
import shutil
import socket
import signal
import inspect

os.sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__)) ) )

import RBNS_utils
import RBNS_cluster_utils

import forgi.graph.bulge_graph as cgb

import cPickle as pickle







def submit_get_Ppaired_DotBracket_andletters_for_reads_F_for_block(
        in_reads_F,
        temp,
        block_idx ):
    """
    - For an in_reads_F split_reads file,
        submits a job to run the
        get_Ppaired_DotBracket_andletters_for_reads_F_for_block() function
           below
    10/3/16
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
            '%(block_idx)s ' % locals())

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
        num_reads_per_block = 1000000,
        #z#fiveP_adapt = "GAGTTCTACAGTCCGACGATC",
        fiveP_adapt = "GGGGAGTTCTACAGTCCGACGATC",
        #threeP_adapt = "TGGAATTCTCGGGTGCCAAGG",
        threeP_adapt = "TGGAATTCTCGGGTGTCAAGG",
        num_reads_to_get_status_after = 1000,
        #num_reads_to_get_status_after = 1000000,
        #copy_back_every_x_seconds = 7500 ):
        copy_back_every_x_seconds = 10000000 ):
    """
    - For an in_reads_F like:
        /net/utr/data/atf/pfreese/RBNS_results/igf2bp1/split_reads/IGF2BP1_input.reads

        for a block_idx, makes a final out_F like:
            /net/utr/data/atf/pfreese/RBNS_results/igf2bp1/split_reads/w_struc/by_block/
            IGF2BP1_input.block_0.w_struc.reads.gz OR
            IGF2BP1_320.block_2.w_struc.reads.gz,
        while also copying back temporary files after every
            copy_back_every_x_seconds; previous versions go in the
            /prev_copy directories

    1/3/17
    """
    temp = int( temp )
    block_idx = int( block_idx )

    start_basename = os.path.basename( in_reads_F ).split(".reads")[0] +\
            ".block_{}".format( block_idx )

    scratch_DIR = "/scratch/pfreese/{0}_{1}".format(
            start_basename, random.randint( 0, 100000 ) )
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
        log_f.write( "fiveP_adapt: {}\n".format( fiveP_adapt ) )
        log_f.write( "threeP_adapt: {}\n".format( threeP_adapt ) )

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

            ############## < if this read is not to be done > #########

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
                    fiveP_adapt,
                    threeP_adapt )

                # Ds_by_read_D = {
                #           "dot_struc": ".....(((((.(((((.((((((((((..(....)..))))))).))).))))))))..)).",
                #           "dotbracket_str": "fffffsssssisssssissssssssssiishhhhsiisssssssisssissssssssiissf"}
                Ds_by_read_D = get_MFE_DotBracket_and_symbol_over_reads(
                    reads_this_set_L,
                    scratch_DIR,
                    temp,
                    fiveP_adapt,
                    threeP_adapt )

                #### Now go through each of the seqs
                with gzip.open( scratch_reads_F, "ab" ) as out_f:
                    for read_idx, rand_read in enumerate( reads_this_set_L ):

                        read_w_adapt = fiveP_adapt + rand_read + threeP_adapt

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

            #if ( num_reads_covered == 100000 ):
            #    break

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










def submit_get_subopt_DotBracket_reads_F_for_block(
        in_struct_gz_F,
        temp,
        block_idx ):
    """
    - For an in_reads_F split_reads file,
        submits a job to run the
        get_Ppaired_DotBracket_andletters_for_reads_F_for_block() function
           below
    10/3/16
    """
    #### Get this file path
    filename = inspect.getframeinfo( inspect.currentframe() ).filename
    this_script_path = os.path.abspath( filename )

    output_DIR = os.path.dirname( in_struct_gz_F )
    errors_outputs_DIR = os.path.join( output_DIR, "errors_outputs" )
    RBNS_utils.make_dir( errors_outputs_DIR )

    command = ('hostname ; python %(this_script_path)s '
            'get_subopt_DotBracket_reads_F_for_block '
            '%(in_struct_gz_F)s '
            '%(temp)s '
            '%(block_idx)s ' % locals())

    job_name = "{0}_block{1}_get_subopt_DotBracket_reads_F_for_block".format(
            os.path.basename( in_struct_gz_F ).split(".")[0], block_idx )
    RBNS_cluster_utils.launch(
            command,
            out_file = os.path.join( errors_outputs_DIR,
                "{}.log".format( job_name ) ),
            jobname = job_name,
            error_dir = errors_outputs_DIR,
            q = 'long' )






def get_subopt_DotBracket_reads_F_for_block(
        in_struct_gz_F,
        temp,
        block_idx,
        num_subopt_to_get = 20,
        num_reads_per_block = 500000,
        num_reads_to_get_status_after = 1000,
        #num_reads_to_get_status_after = 1000000,
        #copy_back_every_x_seconds = 7500 ):
        copy_back_every_x_seconds = 10000000 ):
    """
    - For an in_struct_gz_F:
        /net/utr/data/atf/pfreese/RBNS_results/igf2bp1/split_reads/fld_CG_match

    5/17/17
    """
    temp = int( temp )
    block_idx = int( block_idx )

    start_basename = os.path.basename( in_struct_gz_F ).split(".reads")[0] +\
            ".block_{}".format( block_idx )

    scratch_DIR = "/scratch/pfreese/{0}_{1}".format(
            start_basename, random.randint( 0, 100000 ) )
    os.system( "mkdir -p {}".format( scratch_DIR ) )

    #### If the out_DIR is different than the one simply passed in
    #z#out_DIR = os.path.join( os.path.dirname( in_reads_F ), "w_str", "by_block" )
    #out_DIR = os.path.join( os.path.dirname( in_reads_F ), "str", "by_block" )
    out_DIR = os.path.join( os.path.dirname( in_struct_gz_F ),
            "subopt_DB_{0}reads".format( num_subopt_to_get ) )
    out_logs_DIR = os.path.join( out_DIR, "logs" )
    os.system( "mkdir -p {}".format( out_logs_DIR ) )

    #prev_copies_DIR = os.path.join( os.path.join( out_DIR, "prev_copy" ) )
    #os.system( "mkdir -p {}".format( prev_copies_DIR ) )

    #### copy over the reads to the scratch space
    scratch_in_reads_F = os.path.join( scratch_DIR,
            os.path.basename( in_struct_gz_F ) )
    print "copying to {}".format( scratch_in_reads_F )

    lower_this_block = block_idx * ( num_reads_per_block * 4 )
    upper_this_block = ( block_idx + 1 ) * ( num_reads_per_block * 4 )

    RBNS_utils.copy_gzip_lines_lower_through_upper_to_another_F(
            in_struct_gz_F,
            scratch_in_reads_F,
            lower_this_block,
            upper_this_block )
    #shutil.copyfile( in_reads_F, scratch_in_reads_F )
    print "\tDONE copying {}".format( os.path.basename( in_struct_gz_F ) )

    out_reads_F = os.path.join( out_DIR,
            "{}.subopt_DB.gz".format( start_basename ) )
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

    #copy_from_file = os.path.join( scratch_DIR,
    #    "{}.w_struc.reads.gz".format( start_basename ) )
    #out_DIR = os.path.join( os.path.dirname( in_reads_F ), "w_str" )
    #z#out_DIR = os.path.join( os.path.dirname( in_reads_F ), "w_str" )
    #copy_to_file = os.path.join( out_DIR,
    #    "{}.w_struc.reads.gz".format( start_basename ) )
    #with gzip.open( scratch_reads_F, "wb" ) as out_f:
    #    pass

    #abs_read_idx = 0

    reads_this_set_L = []
    num_reads_this_set = 0
    for lines_L in RBNS_utils.iterNlines_strip( scratch_in_reads_F, 4 ):

        ############## < if this read is not to be done > #########
        #if ( abs_read_idx < lower_this_block ):
        #    abs_read_idx += 1
        #    continue
        #if ( abs_read_idx >= upper_this_block ):
        #    break
        #abs_read_idx += 1

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
        #if ( num_reads_covered in [5000, 10000, 50000] + [i*100000 for i in range( 1, 10 )] ):

        #    curr_time = time.time()

        #    secs_per_million_reads = (curr_time - start_time) * 1000000. / num_reads_covered
        #    readable_string_secs_per_mil = RBNS_utils.return_human_readable_time_from_secs(
        #            secs_per_million_reads )
        #    with open( log_F, "a" ) as log_f:
        #        curr_pprint_str = RBNS_utils.return_nice_datetime_str_for_filename()
        #        curr_str = "{0:,} at {1} for {2}: {3} sec per million rds\n".format(
        #            num_reads_covered,
        #            curr_pprint_str,
        #            start_basename,
        #            readable_string_secs_per_mil )
        #        print curr_str
        #        log_f.write( curr_str )

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

            #if ( num_reads_covered == 100000 ):
            #    break

    shutil.copyfile( out_F_to_append_results_to, out_reads_F )
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
    making_F = os.path.join( os.path.dirname( out_reads_F ),
        "{0}.block_{1}.subopt_DB.gz.making".format(
            os.path.basename( in_struct_gz_F ).split(".reads")[0], block_idx ) )
    if os.path.exists( making_F ):
        os.system( "rm {}".format( making_F ) )










###############################################################################
###############################################################################
#################################### < UTILS > ################################




def get_avg_bp_prob_with_flankingprimers(
        RNA_seq,
        read_length,
        tmp_DIR,
        fivePprimer,
        threePprimer,
        temp = 37 ):
    """
    - Calculates the average bp. probability for each base of RNA_seq
    """
    # Get the primer lengths
    fivePprimer_len = len( fivePprimer )
    threePprimer_len = len( threePprimer )

    # First, make a temporary .fa file of the read
    tmp_read_fasta_F = os.path.join( tmp_DIR, "read.fa" )
    with open( tmp_read_fasta_F, "w" ) as f:
        f.write( ">read\n{0}".format(
            fivePprimer + RNA_seq + threePprimer ) )
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
            fivePprimer_len + read_length + threePprimer_len )

    avg_bp_probs_randomregion_L =\
        avg_bp_probs_L[fivePprimer_len:(-1*threePprimer_len)]

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
    # go through the file
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
                    #print "INDEX ERROR FOR {0} in {1}".format(
                    #        ln, probs_ps_F )
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
        - read_length: the length of the seqs to keep (e.g. 39 of the 40
            bases of each read)
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
        #avg_bp_probs_wo_adapters_L = avg_bp_probs_L[len_5p:(-1*len_3p)]
        #avg_bp_probs_Ls_L.append( avg_bp_probs_wo_adapters_L )
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
    #rand_str = str( random.random() )

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
        DotBracket_string = '(((((((((...((((((.........))))))........((((((.......))))))..)))))))))'
    RETURNS:
        sssssssssmmmsssssshhhhhhhhhssssssmmmmmmmmsssssshhhhhhhssssssmmsssssssss

    10/3/16
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
    5/12/17
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
    #os.system("cd {0}".format( scratch_DIR ))
    # make the .ps files
    fold = subprocess.Popen( fold_CMD, shell = True )
    stdoutdata, stderrdata = fold.communicate()

    #### Now go through all of the reads
    this_read = ""
    out_f_to_append_results_to = gzip.open( out_F_to_append_results_to, 'ab' )
    with open( out_F ) as f:
        for line in f:
            if ( line[0] == '>' ):
                read_w_idx = line.strip()[1:]
                this_read = read_by_readwindex_D[read_w_idx]
                out_f_to_append_results_to.write( ">" + this_read + "\n" )
            else:
                DB_str = line.strip()
                element_string = get_elementstring_from_DotBracket( DB_str )
                out_f_to_append_results_to.write(
                        element_string + "\n" )

    out_f_to_append_results_to.close()

    #print out_F










################################# </ UTILS > ##################################
###############################################################################
###############################################################################




def test():

    scratch_DIR = "/scratch/pfreese/test_subopt"
    os.system( "mkdir -p {}".format( scratch_DIR ) )
    out_F = os.path.join( scratch_DIR, "reads.output.txt" )
    reads_L = ["ATGTAATCCCCACTATATTCAACTCAAAGGTATAGCTTTA",
            "GTGACGGAACCATGAAATACGTCTAATTTTTTGCACGATG"]
    temp = 37
    #get_subopt_folding_of_reads(
    #    reads_L,
    #    scratch_DIR,
    #    temp,
    #    out_F,
    #    num_to_return_for_each_read = 20 )
    get_subopt_DotBracket_reads_F_for_block(
            "/net/utr/data/atf/pfreese/RBNS_results/igf2bp1/split_reads/fld_CG_match/IGF2BP1_0.w_struc.reads.gz",
            4,
            0,
            num_subopt_to_get = 20,
            num_reads_per_block = 10000 )



if __name__ == '__main__':

    fxn = sys.argv[1].strip('-')
    args = ['"' + arg + '"' for arg in sys.argv[2:]]
    python_command = fxn + '(' + ','.join( args ) + ')'

    print python_command
    eval( python_command )







