#!/usr/bin/env python
import os, sys
import pprint, tempfile, glob
import glob
import datetime, time
import shutil
import argparse
import math, random
import cPickle as pickle
from math import log


import RBNS_Config_helpers
from RBNS_cluster_utils import launch
import RBNS_logo_helpers
import file_IO
import RBNS_plots
import RBNS_utils




def return_Config_D(
        config_F ):
    """
    - Returns a settings_D for an input Config_F
    - uses the helper function in RBNS_Config_helpers as a template to get config_D
    """
    #### Parses the config_F to include the following keys in the returned
    ####    config_D, using the function in RBNS_Config_helpers
    int_settings_L = ["read_len_to_use",
                        "starting_k",
                        "ending_k" ]
    float_settings_L = ["Zscore_kmers_to_keep"]
    list_settings_L = ["kmers_to_ignore"]
    boolean_settings_L = ["remove_all_previous_DIRs"]
    int_or_float_settings_L = ["num_input_pulldown_rds_to_use"]

    #### calls the helper function in RBNS_Config_helpers.py
    config_D = RBNS_Config_helpers.get_settings_template(
        config_F,
        int_settings_L,
        float_settings_L,
        int_or_float_settings_L,
        list_settings_L,
        boolean_settings_L )

    return config_D







class RBNS_motifs:
    """
    - Performs the motif logo pipeline
    """
    def __init__( self, config_F ):
        self.config_F = config_F
        self.process_config_F( config_F )


    def process_config_F( self, config_F ):

        #### Get a time_stamp
        self.time_stamp = datetime.datetime.now().strftime("%m_%d_%Y__%Hh:%Mm:%Ss")

        config_D = return_Config_D( config_F )

        ############# < GET THE VARIABLES FROM THE config_D > #################
        ########################## < REQUIRED ARGS > ##########################
        self.input_reads_F = config_D["input_reads_F"]
        self.pulldown_reads_F = config_D["pulldown_reads_F"]
        self.protein_name = config_D["protein_name"]
        self.starting_k = config_D["starting_k"]

        self.input_orig_read_len = file_IO.get_readlength(
                self.input_reads_F )
        self.input_pulldown_read_len = file_IO.get_readlength(
                self.pulldown_reads_F )
        assert( self.input_orig_read_len == self.input_pulldown_read_len )
        ########################## </ REQUIRED ARGS > #########################

        ########################## < OPTIONAL ARGS > ##########################
        #### The lowest k to look for motifs ( by default, =4 )
        try:
            ending_k = config_D["ending_k"]
        except KeyError:
            ending_k = 4
        self.ending_k = ending_k

        #### The Z-score criterion for when to stop searching for motifs
        try:
            Zscore_kmers_to_keep = config_D["Zscore_kmers_to_keep"]
        except KeyError:
            Zscore_kmers_to_keep = 2.
        self.Zscore_kmers_to_keep = Zscore_kmers_to_keep


        #### The output DIR
        try:
            out_DIR = config_D["output_DIR"]
        except KeyError:
            out_DIR = os.path.join(
                os.path.dirname(os.path.dirname(self.pulldown_reads_F)), "logos" )

        out_ks_DIR = os.path.join( out_DIR,
                "k_{0}_to_{1}_Zscoretokeep_{2:.1f}".format(
                    self.starting_k, self.ending_k,
                    Zscore_kmers_to_keep ))
        self.out_DIR = out_ks_DIR
        #### If all the previously made directories for this range of
        ####    k's should be deleted
        try:
            remove_all_previous_DIRs = config_D["remove_all_previous_DIRs"]
        except KeyError:
            remove_all_previous_DIRs = True
        if (remove_all_previous_DIRs == True):
            try:
                os.system( "rm -rf {}".format( out_ks_DIR ) )
                time.sleep( 1 )
            except:
                pass
        #### Now make the out_DIR
        try:
            os.makedirs( self.out_DIR )
        except OSError:
            pass

        #### The directory with the logos made from the kmers consistent
        ####    between the two halves
        self.consistent_DIR = os.path.join( self.out_DIR, "consistent" )
        try:
            os.makedirs( self.consistent_DIR )
        except OSError:
            pass

        #### The directory with the logos made from the kmers that are in
        ####    both halves
        self.in_both_DIR = os.path.join( self.out_DIR, "in_both" )
        try:
            os.makedirs( self.in_both_DIR )
        except OSError:
            pass


        #### If a trimmed read length should be used
        try:
            read_len_to_use = config_D["read_len_to_use"]
        except KeyError:
            read_len_to_use = self.input_orig_read_len
        self.read_len_to_use = read_len_to_use


        #### The number of input & pulldown reads to use
        try:
            num_input_pulldown_rds_to_use = config_D["num_input_pulldown_rds_to_use"]
        except KeyError:
            num_input_pulldown_rds_to_use = 0.5
        self.num_input_pulldown_rds_to_use = num_input_pulldown_rds_to_use


        #### The protein name to put on plots
        try:
            protein_name_for_plotting = config_D["protein_name_for_plotting"]
        except KeyError:
            protein_name_for_plotting = self.protein_name
        self.protein_name_for_plotting = protein_name_for_plotting

        #### Any kmers which should be ignored
        try:
            kmers_to_ignore_L = config_D["kmers_to_ignore"]
        except KeyError:
            kmers_to_ignore_L = []
        self.kmers_to_ignore_L = kmers_to_ignore_L

        #### A title to override
        try:
            barplot_title = config_D["barplot_title"]
        except KeyError:
            barplot_title = "{0}: {1}mers".format(
                self.protein_name_for_plotting,
                self.starting_k )
        self.barplot_title = barplot_title
        ########################## </ OPTIONAL ARGS > #########################
        ############### </ GET THE VARIABLES FROM THE config_D > ##############

        #### Make a scratch directory to use, and sub_DIRs within it:
        ####    /reads_Fs
        ####    /logos
        scratch_DIR = config_D['scratch_DIR']
        self.temp_DIR = os.path.join( scratch_DIR,
                "{0}_{1}".format( self.time_stamp, random.random()) )

        self.temp_reads_DIR = os.path.join( self.temp_DIR, "reads_Fs" )
        os.makedirs( self.temp_reads_DIR )

        #self.temp_logos_DIR_1 = os.path.join( self.temp_DIR, "logos_1", "indiv" )
        self.temp_logos_DIR_1 = os.path.join( self.temp_DIR, "logos_1" )
        os.makedirs( self.temp_logos_DIR_1 )
        #self.temp_logos_DIR_1_to_move = os.path.dirname( self.temp_logos_DIR_1 )

        #self.temp_logos_DIR_2 = os.path.join( self.temp_DIR, "logos_2", "indiv" )
        self.temp_logos_DIR_2 = os.path.join( self.temp_DIR, "logos_2" )
        os.makedirs( self.temp_logos_DIR_2 )
        #self.temp_logos_DIR_2_to_move = os.path.dirname( self.temp_logos_DIR_2 )


        #######################################################################
        ########## < Record all of the parameters in the log_F > ##############

        log_F_1 = os.path.join( self.out_DIR, "log_1.txt" )
        self.log_F_1 = log_F_1
        print "\nlog_F_1 is: {}\n\n".format( self.log_F_1 )
        self.log_f_1 = open( log_F_1, "w" )

        log_F_2 = os.path.join( self.out_DIR, "log_2.txt" )
        self.log_F_2 = log_F_2
        print "\nlog_F_2 is: {}\n\n".format( self.log_F_2 )
        self.log_f_2 = open( log_F_2, "w" )

        ########## </ Record all of the parameters in the log_F > #############
        #######################################################################


        #######################################################################
        ### < Files that will contain the kmers aligned to make the logos > ###

        self.kmers_to_make_logos_Fs_by_half_D = {}
        for half in [1, 2]:
            log_F = os.path.join( self.out_DIR, "{0}_{1}mers_for_logo.{2}.txt".\
                format( config_D["protein_name"], config_D["starting_k"], half ) )
            self.kmers_to_make_logos_Fs_by_half_D[half] = log_F
        #### Also make one for the "consistent" kmers
        log_F = os.path.join( self.out_DIR, "{0}_{1}mers_for_logo.txt".\
            format( config_D["protein_name"], config_D["starting_k"] ) )
        self.kmers_to_make_logos_Fs_by_half_D["consistent"] = log_F
        #### Also make one for the kmers in both
        log_F = os.path.join( self.out_DIR, "{0}_{1}mers_for_logo.in_both.txt".\
            format( config_D["protein_name"], config_D["starting_k"] ) )
        self.kmers_to_make_logos_Fs_by_half_D["in_both"] = log_F

        ## </ Files that will contain the kmers aligned to make the logos > ###
        #######################################################################

        os.system( "rm -rf {}".format( self.temp_DIR ) )








def execute_get_logos(
        config_F,
        verbose = False):
    """
    - Performs the generation of RBNS motif logos

    REQUIRED in the config_F:
        - input_reads_F
        - pulldown_reads_F
        - protein_name: a unique identifier of this experiment (will be in file
            names, so it should not have any spaces)
        - starting_k: the largest k to start with; will go from this k to
            ending_k

    OPTIONAL ARGUMENTS, with default values:
        - out_DIR: the directory above that containing the pulldown_reads_F,
            plus the sub_DIR "logos"
        - read_len_to_use: the read length of input_reads_F / pulldown_reads_F
        - protein_name_for_plotting: a name for use in plots (can include
            spaces)
        - ending_k: the lowest length kmers that will be looked for: 4
        - num_input_pulldown_rds_to_use: 0.5
        - Zscore_kmers_to_keep: 2.
    """

    settings = RBNS_motifs( config_F )

    ############## < GET A LOGO FOR THE FIRST HALF OF THE READS > #############
    settings.log_f_1.write( "config_F: {}\n\n".format( config_F ) )

    settings.log_f_1.write( "="*31 + " < SETTINGS > " + "="*31 + "\n\n" )
    attributes = vars( settings )
    #### Write all non-dictionary attributes to the log_f (kmer freq.
    ####    dictionaries are too large to write out)
    settings.log_f_1.write('\n\n'.join("%s: %s" % item for item in attributes.items()\
            if type(item[1]) != dict))
    settings.log_f_1.write( "\n\n" + "="*30 + " < / SETTINGS > " +\
            "="*30 + "\n\n" )

    #### Call the main function that gets the logos
    returned_1_D = get_logos( settings, 1 )

    #### Once finished, copy the logos_DIR from the scratch space to the
    ####    out_DIR
    out_thisrep_DIR = os.path.join( settings.out_DIR,
            os.path.basename( settings.temp_logos_DIR_1 ) )
    if os.path.exists( out_thisrep_DIR ):
        os.system( "rm -rf {}".format( out_thisrep_DIR ) )
    shutil.move( settings.temp_logos_DIR_1, settings.out_DIR )

    #### Close the log_F
    os.system( "rm -rf {}".format( settings.temp_reads_DIR ) )
    settings.log_f_1.close()
    ############# </ GET A LOGO FOR THE FIRST HALF OF THE READS > #############

    ############# < GET A LOGO FOR THE SECOND HALF OF THE READS > #############
    os.makedirs( settings.temp_reads_DIR )
    settings.log_f_2.write( "config_F: {}\n\n".format( config_F ) )

    settings.log_f_2.write( "="*31 + " < SETTINGS > " + "="*31 + "\n\n" )
    attributes = vars( settings )
    #### Write all non-dictionary attributes to the log_f (kmer freq.
    ####    dictionaries are too large to write out)
    settings.log_f_2.write('\n\n'.join("%s: %s" % item for item in attributes.items()\
            if type(item[1]) != dict))
    settings.log_f_2.write( "\n\n" + "="*30 + " < / SETTINGS > " +\
            "="*30 + "\n\n" )

    #### Call the main function that gets the logos
    returned_2_D = get_logos( settings, 2 )

    #### Once finished, copy the logos_DIR from the scratch space to the
    ####    out_DIR
    out_thisrep_DIR = os.path.join( settings.out_DIR,
            os.path.basename( settings.temp_logos_DIR_2 ) )
    if os.path.exists( out_thisrep_DIR ):
        os.system( "rm -rf {}".format( out_thisrep_DIR ) )
    shutil.move( settings.temp_logos_DIR_2, settings.out_DIR )

    #### Close the log_F
    settings.log_f_2.close()
    ############# </ GET A LOGO FOR THE SECOND HALF OF THE READS > ############

    os.system( "rm -rf {}".format( settings.temp_reads_DIR ) )
    print "\n\nFINISHED! See files in: {}".format( settings.out_DIR )

    ########## < kmers in a consistent ordering between the 2 halves > ########
    if os.path.exists( settings.consistent_DIR ):
        os.system( "rm -rf {}".format( settings.consistent_DIR ) )
    get_one_set_logos_based_on_consistent_kmer_orderings(
            returned_1_D["kmer_offset_R_T_Ls_by_classnum_D"],
            returned_2_D["kmer_offset_R_T_Ls_by_classnum_D"],
            settings )
    ######### </ kmers in a consistent ordering between the 2 halves > ########

    ################# < kmers that are in both of the 2 halves > ##############
    kmers_F_1 = settings.kmers_to_make_logos_Fs_by_half_D[1]
    kmers_F_2 = settings.kmers_to_make_logos_Fs_by_half_D[2]
    if os.path.exists( settings.in_both_DIR ):
        os.system( "rm -rf {}".format( settings.in_both_DIR ) )
    get_one_set_logos_based_on_kmers_in_both(
            kmers_F_1,
            kmers_F_2,
            settings )
    ################# < kmers that are in both of the 2 halves > ##############

    #### Get the .PWM files for each logo
    logos_1_DIR = os.path.join( settings.out_DIR, "logos_1" )
    make_PWM_files_for_all_aligned_kmers_Fs_in_DIR( logos_1_DIR )

    logos_2_DIR = os.path.join( settings.out_DIR, "logos_2" )
    make_PWM_files_for_all_aligned_kmers_Fs_in_DIR( logos_2_DIR )

    make_PWM_files_for_all_aligned_kmers_Fs_in_DIR( settings.consistent_DIR )
    make_PWM_files_for_all_aligned_kmers_Fs_in_DIR( settings.in_both_DIR )








def get_logos( settings,
        half ):
    """
    - Takes in a settings instance from the RBNS_motifs class and gets logos
        for 'half' 1 or 2
    """
    #### Remove any reads already in the temp_reads_DIR
    os.system( "rm {}".format( os.path.join( settings.temp_reads_DIR, "*.reads" ) ) )
    assert( half in [1, 2] )
    if (half == 1):
        settings.log_f = settings.log_f_1
        settings.temp_logos_DIR = settings.temp_logos_DIR_1
        start_frac = 0.
    else:
        settings.log_f = settings.log_f_2
        settings.temp_logos_DIR = settings.temp_logos_DIR_2
        start_frac = 0.5

    #### Copy over the original input & pulldown reads using the helper
    ####    function in file_IO
    ########################### < INPUT > #####################################
    returned_input_reads_D = file_IO.make_temp_reads_F(
            settings.input_reads_F,
            settings.temp_reads_DIR,
            read_length_to_use = settings.read_len_to_use,
            num_reads_to_use = settings.num_input_pulldown_rds_to_use,
            start_frac = start_frac,
            target_reads_basename = "input.reads" )

    #### add the returned .reads file and the number of reads in it to settings
    settings.input_reads_F_created = returned_input_reads_D["out_reads_F"]
    settings.num_input_reads_F_created = returned_input_reads_D[
            "num_reads_in_out_reads_F"]

    ########################## < PULLDOWN > ###################################
    returned_pulldown_reads_D = file_IO.make_temp_reads_F(
            settings.pulldown_reads_F,
            settings.temp_reads_DIR,
            read_length_to_use = settings.read_len_to_use,
            num_reads_to_use = settings.num_input_pulldown_rds_to_use,
            start_frac = start_frac,
            target_reads_basename = "pulldown.reads" )

    #### add the returned .reads file and the number of reads in it to settings
    settings.pulldown_reads_F_created = returned_pulldown_reads_D["out_reads_F"]
    settings.num_pulldown_reads_F_created = returned_pulldown_reads_D[
            "num_reads_in_out_reads_F"]

    #### A list of all kmers which have significant R's but are not included in
    ####    the logo since they contain one of the kmers_to_ignore
    kmers_ignored_in_logo_L = []

    #################### < ORIG. FREQS & ENRICHMENTS > ########################
    settings.orig_input_freqs_Ds_by_k = {}
    settings.orig_pulldown_freqs_Ds_by_k = {}

    settings.orig_enrich_Ds_by_k = {}
    settings.orig_kmer_enrich_tuples_Ls_by_k = {}

    settings.orig_R_mean_std_by_k = {}

    for k in range( settings.ending_k, settings.starting_k + 1 ):

        ############################### < INPUT > #############################
        ####    Use the helper function in file_IO to get the kmer_freqs_Ds
        input_kmer_freqs_D = file_IO.get_kmer_freqs_from_reads_F(
                settings.input_reads_F_created,
                k )["freqs_by_kmer_D"]
        settings.orig_input_freqs_Ds_by_k[k] = input_kmer_freqs_D

        ############################## < PULLDOWN > ###########################
        pulldown_kmer_freqs_D = file_IO.get_kmer_freqs_from_reads_F(
                settings.pulldown_reads_F_created,
                k )["freqs_by_kmer_D"]
        settings.orig_pulldown_freqs_Ds_by_k[k] = pulldown_kmer_freqs_D

        ################## < ENRICHMENTS and EXCESS FREQS > ###################
        enrichments_D = {}
        kmer_enrich_tuples_L = []
        for kmer in RBNS_utils.yield_kmers( k ):
            try:
                R = pulldown_kmer_freqs_D[kmer] / input_kmer_freqs_D[kmer]
            except ZeroDivisionError:
                R = 1.
            enrichments_D[kmer] = R
            if ( kmer not in settings.kmers_to_ignore_L ):
                kmer_enrich_tuples_L.append( (kmer, R) )

        settings.orig_enrich_Ds_by_k[k] = enrichments_D

        kmer_enrich_tuples_L.sort( key = lambda x: -1*x[1] )

        settings.orig_kmer_enrich_tuples_Ls_by_k[k] = kmer_enrich_tuples_L

        #### Get the mean & st. deviation of the Rs, and add these to
        ####    settings.orig_R_mean_std_by_k
        mean_R, std_R = RBNS_utils.mean_std(
                [ x[1] for x in kmer_enrich_tuples_L ] )
        settings.orig_R_mean_std_by_k[k] =\
            {"mean": mean_R,
            "std": std_R,
            "thresh_Zscore2": mean_R + 2*std_R,
            "thresh_Zscore3": mean_R + 3*std_R,
            "thresh_R":\
                mean_R + (settings.Zscore_kmers_to_keep * std_R)}

    #### Write out the settings.orig_R_mean_std_by_k to the log_F
    settings.log_f.write( "\n\n\norig_R_mean_std_by_k:\n" )
    pprint.pprint( settings.orig_R_mean_std_by_k, settings.log_f )
    settings.log_f.write( "\n\n\n" )

    #################### </ ORIG. FREQS & ENRICHMENTS > #######################


    ################ < GO THROUGH EACH k, GET KMER FREQUENCIES > ##############

    #### A dictionary that will be used to make each of the motifs, and
    ####    eventually used to make the "consistent" logos
    kmer_offset_R_T_Ls_by_classnum_D = {}

    #### The "current" pulldown and input files, which will be updated with the
    ####    created reads files that have "X"s for kmers removed
    current_pulldown_reads_F = settings.pulldown_reads_F_created
    current_input_reads_F = settings.input_reads_F_created

    #### A list of lines of the aligned kmers to write at the bottom of the
    ####    log_F
    aligned_lines_end_of_log_L = []

    #### A list of all kmers (including all previous k's) that have been
    ####    used
    all_kmers_used_L = []

    #### The maximum number of kmers to include in the logo
    max_num_kmers = int( 4 ** (k - 2 ) )

    num_kmers_covered = 0
    for k in range( settings.starting_k , settings.ending_k - 1, -1 ):

        print "\n\t\t\t\t\tk = {}".format( k )
        settings.log_f.write( "\n\n" + "="*34 + " < k = {} > ".format( k ) +\
                "="*34 + "\n\n" )

        aligned_lines_end_of_log_L.append(
                "\n\nk={}\n".format( k ) )

        cont_this_k = True
        if ( num_kmers_covered >= max_num_kmers ):
            cont_this_k = False
        kmer_num_this_k = 0
        class_num = 0

        ####    The excess_freq's above this threshhold will be added to
        ####    the logos
        thresh_R = settings.orig_R_mean_std_by_k[k]["thresh_R"]

        ####    Get the current kmer freq. counts
        current_pulldown_kmer_freqs_D = file_IO.get_kmer_freqs_from_reads_F(
                current_pulldown_reads_F,
                k )["freqs_by_kmer_D"]
        current_input_kmer_freqs_D = file_IO.get_kmer_freqs_from_reads_F(
                current_input_reads_F,
                k )["freqs_by_kmer_D"]

        kmer_R_tuples_L = []
        for kmer, pulldown_freq in current_pulldown_kmer_freqs_D.iteritems():
            try:
                R = pulldown_freq / current_input_kmer_freqs_D[kmer]
            except ZeroDivisionError:
                R = 1.
            if kmer not in settings.kmers_to_ignore_L:
                kmer_R_tuples_L.append( (kmer, R) )
        kmer_R_tuples_L.sort( key = lambda x: -1*x[1] )
        current_most_enriched_kmer = kmer_R_tuples_L[0][0]
        settings.log_f.write( "\nTop 10 enriched kmers for the first run of k={}:\n".format(
            k) )
        pprint.pprint( kmer_R_tuples_L[:10], settings.log_f )
        settings.log_f.write( "\n\n" )

        while (cont_this_k == True):

            ##### Get the kmer to count & replace; if this is the first kmer
            #####   for this k, use the highest originall enriched kmer;
            #####   otherwise, use the most enriched kmer calculated at the
            #####   end of the previous round using all of the X'ed out kmers
            kmer_to_count = current_most_enriched_kmer
            orig_enrich =\
                settings.orig_enrich_Ds_by_k[k][kmer_to_count]
            settings.log_f.write( "\n\nCONSIDERING: {0} (orig R: {1:.3f})\n".format(
                kmer_to_count, orig_enrich ) )

            #### If the R exceeds the threshold cutoff calculated
            ####    from the original reads file, use this kmer and add it to
            ####    a motif class;
            pulldown_freq_kmer = current_pulldown_kmer_freqs_D[kmer_to_count]
            input_freq_kmer = current_input_kmer_freqs_D[kmer_to_count]
            try:
                kmer_R = pulldown_freq_kmer / input_freq_kmer
            except ZeroDivisionError:
                kmer_R = 1.
            if (kmer_R >= thresh_R):

                #### Even if this kmer has a high R, make sure that it doesn't
                ####    contain any of the kmers_to_ignore_L
                no_kmers_to_ignore = True
                for kmer_to_ignore in settings.kmers_to_ignore_L:
                    if (kmer_to_count.find( kmer_to_ignore ) != -1):
                        no_kmers_to_ignore = False
                        settings.log_f.write( "IGNORING {0} since it contains the kmer_to_ignore {1}\n\n".format(
                            kmer_to_count, kmer_to_ignore ))
                        kmers_ignored_in_logo_L.append( kmer_to_count )
                        print "\nIGNORING {0} (contains {1})\n".format(
                                kmer_to_count,
                                kmer_to_ignore )
                        break

                no_previous_kmers = True
                for prev_kmer in all_kmers_used_L:
                    if (prev_kmer.find( kmer_to_count ) != -1):
                        no_previous_kmers = False

                ##### < IF THIS kmer doesn't contain any kmers_to_ignore > ####
                if (no_kmers_to_ignore == True) and (no_previous_kmers == True):

                    all_kmers_used_L.append( kmer_to_count )

                    settings.log_f.write( "\t- ACCEPTED! Adding {0} since the R {1:.4g} is greater than the thresh_R (={2:.4g})\n\n".format(
                            kmer_to_count, kmer_R, thresh_R ))
                    #### the excess R is kmer_R - 1
                    excess_R = kmer_R - 1
                    # THE IMPORTANT FUNCTION THAT ADDS THE KMER TO EXISTING MOTIFS
                    returned_D = RBNS_logo_helpers.add_kmer_to_existing_motif_class_or_start_new(
                            kmer_to_count,
                            excess_R,
                            kmer_offset_R_T_Ls_by_classnum_D )
                    #< THE IMPORTANT FUNCTION THAT ADDS THE KMER TO EXISTING MOTIFS

                    kmer_offset_R_T_Ls_by_classnum_D = returned_D[
                            "kmer_offset_R_T_Ls_by_classnum_D"]
                    motif_added_to = returned_D["motif_class_added_to"]
                    offset_used = returned_D["offset_used"]

                    #### Print out this kmer aligned to its motif class
                    excess_R_str = "{:.3f}".format( excess_R )
                    excess_R_str += (" " * (8 - len( excess_R_str )))
                    spaces_before = (3 * k * motif_added_to) + int(k/2) + offset_used
                    aligned_line = excess_R_str + ( " "*spaces_before ) + kmer_to_count
                    print aligned_line
                    aligned_lines_end_of_log_L.append( "\n" + aligned_line )

                #### < / IF THIS kmer doesn't contain any kmers_to_ignore > ###


                #### Update the input & pulldown reads to the new reads files
                ####    with kmer occurrences X'ed out
                ########## Replace this kmer in the PULLDOWN file #################
                #### - Uses the helper function in file_IO, which returns:
                ####    return_D = {"out_reads_F": out_reads_F,
                ####                "tot_num_reads": tot_num_reads,
                ####                "num_reads_w_kmer": num_reads_w_kmer,
                ####                "freq_reads_w_kmer": freq_reads_w_kmer,
                ####                "tot_num_kmer_occurs" : tot_num_kmer_occurs,
                ####                "counts_by_kmer_D": counts_by_kmer_D,
                ####                "freqs_by_kmer_D": freqs_by_kmer_D}
                pulldown_returned_D = file_IO.return_frequency_and_number_of_reads_kmer_in_reads_F(
                        current_pulldown_reads_F,
                        kmer_to_count )

                ########### < Replace this kmer in the INPUT file > ###############
                input_returned_D = file_IO.return_frequency_and_number_of_reads_kmer_in_reads_F(
                        current_input_reads_F,
                        kmer_to_count )

                current_pulldown_reads_F = pulldown_returned_D["out_reads_F"]
                current_input_reads_F = input_returned_D["out_reads_F"]

                #### Calculate the new enrichments using the pulldown & input
                ####    files for which the kmer's have been X'ed out; the
                ####    kmer used in the next round will be the most enriched
                ####    from this calculation
                current_pulldown_kmer_freqs_D = pulldown_returned_D["freqs_by_kmer_D"]
                current_input_kmer_freqs_D = input_returned_D["freqs_by_kmer_D"]
                kmer_R_tuples_L = []
                for kmer, pulldown_freq in current_pulldown_kmer_freqs_D.iteritems():
                    try:
                        R = pulldown_freq / current_input_kmer_freqs_D[kmer]
                    except ZeroDivisionError:
                        R = 1.
                    kmer_R_tuples_L.append( (kmer, R) )
                kmer_R_tuples_L.sort( key = lambda x: -1*x[1] )
                current_most_enriched_kmer = kmer_R_tuples_L[0][0]
                settings.log_f.write( "\nTop 10 enriched kmers after removing {}:\n".format(
                    kmer_to_count) )
                pprint.pprint( kmer_R_tuples_L[:10], settings.log_f )
                settings.log_f.write( "\n\n" )

                kmer_num_this_k += 1

            #### Otherwise, DON'T use it (use the reads files from before this
            ####    kmer was X'ed out)
            else:
                cont_this_k = False
                settings.log_f.write( "\t- REJECTED! NOT ADDING {0} since the R {1:.4g} is less than the thresh_R (={2:.4g})\n\n".format(
                    kmer_to_count, kmer_R, thresh_R ))

            num_kmers_covered += 1
            if ( num_kmers_covered >= max_num_kmers ):
                cont_this_k = False
                break
        settings.log_f.write( "\n\n" + "="*33 + " < / k = {} > ".format( k ) +\
                "="*33 + "\n\n\n" )

    #### Write out all of the aligned kmers in the log_F
    for aligned_line in aligned_lines_end_of_log_L:
        settings.log_f.write( aligned_line )

    #### If any kmers_ignored_in_logo_L were encountered
    if ( len( kmers_ignored_in_logo_L ) > 0 ):
        settings.log_f.write( "The following kmers had significant R's, but were NOT included since they contain one of the kmers_to_ignore (={}):\n\n".format( settings.kmers_to_ignore_L ) )
        pprint.pprint( kmers_ignored_in_logo_L, settings.log_f )
        settings.log_f.write( "\n\n" )
    elif (len( settings.kmers_to_ignore_L ) > 0):
        settings.log_f.write( "No kmers w/ significant R contained any of the kmers_to_ignore (={}).\n".format( settings.kmers_to_ignore_L ) )

    #### Now make all of the logos using the helper function in
    ####    RBNS_logo_helpers.py
    all_logos_L = []
    all_seq_logos_L = []
    logos_excess_Rs_L = []

    temp_logos_DIR = os.path.join( settings.temp_logos_DIR, "indiv" )
    os.makedirs( temp_logos_DIR )

    #### Call the helper function to make the logos
    returned_D = RBNS_logo_helpers.make_logos_from_kmer_offset_R_T_Ls_by_classnum_D(
            kmer_offset_R_T_Ls_by_classnum_D,
            temp_logos_DIR,
            protein_name = settings.protein_name_for_plotting )

    prob_logos_by_classnum_D = returned_D["prob_logos_by_classnum_D"]
    seq_logos_by_classnum_D = returned_D["seq_logos_by_classnum_D"]
    pruned_pos_to_include_by_classnum_D = returned_D["pruned_pos_to_include_by_classnum_D"]

    #### Go through each motif_class and add the logo & its excess_R
    ####     to all_logos_L & logos_excess_Rs_L, respectively
    motif_classes_L = prob_logos_by_classnum_D.keys()
    motif_classes_L.sort()
    for motif_class in motif_classes_L:

        prob_logo_F = prob_logos_by_classnum_D[motif_class]
        seq_logo_F = seq_logos_by_classnum_D[motif_class]
        excess_R = sum([x[2] for x in kmer_offset_R_T_Ls_by_classnum_D[motif_class]])

        all_logos_L.append( prob_logo_F )
        all_seq_logos_L.append( seq_logo_F )
        logos_excess_Rs_L.append( excess_R )

    #### Now that all logos have been made, plot them all together in one
    ####    stacked bargraph w/ logos using the helper function in plots_bar
    out_combined_plots_F = os.path.join( settings.temp_logos_DIR,
            "{0}_{1}mer_logos_{2}.pdf".format(
                settings.protein_name_for_plotting,
                settings.starting_k,
                half ) )
    returned_D = RBNS_plots.plot_stacked_bargraph_with_logos(
            logos_excess_Rs_L,
            all_logos_L,
            out_combined_plots_F,
            title = settings.barplot_title )
    out_combined_seq_plots_F = os.path.join( settings.temp_logos_DIR,
            "{0}_{1}mer_seqlogos_{2}.pdf".format(
                settings.protein_name_for_plotting,
                settings.starting_k,
                half ) )
    returned_D = RBNS_plots.plot_stacked_bargraph_with_logos(
            logos_excess_Rs_L,
            all_seq_logos_L,
            out_combined_seq_plots_F,
            title = settings.barplot_title )

    #### Make the final file of all the aligned kmers that went into
    ####    each motif class
    out_aligned_kmers_F = settings.kmers_to_make_logos_Fs_by_half_D[half]
    RBNS_logo_helpers.make_final_file_aligned_kmers_weights_per_motif_class(
        kmer_offset_R_T_Ls_by_classnum_D,
        out_aligned_kmers_F )

    #### The dictionary to return
    return_D = { "kmer_offset_R_T_Ls_by_classnum_D":
            kmer_offset_R_T_Ls_by_classnum_D }

    return return_D










def get_one_set_logos_based_on_consistent_kmer_orderings(
        kmer_offset_R_T_Ls_by_classnum_D_1,
        kmer_offset_R_T_Ls_by_classnum_D_2,
        settings ):
    """
    - Using the constituent kmers, offsets, and R values of two sets of logos,
        each derived from half of the data, makes a
        composite logo consisting only of kmers that are consistent in the two
        halve (i.e., kmers that are in the same relative ordering between them)

    - This function is called at the end of execute_get_logos()
    """

    all_logos_L = []
    all_seq_logos_L = []
    logos_excess_Rs_L = []

    consistent_class_idx = 0
    ### A dictionary that will be passed into
    ####    make_aligned_kmers_F_from_kmer_PWMs_by_k_classnum_D() in
    ####    RBNS_logo_helpers.py
    kmer_offset_R_tuples_Ls_by_classnum_D = {}

    num_motif_classes = len( kmer_offset_R_T_Ls_by_classnum_D_1 )

    for motif_class_num in range( num_motif_classes ):

        kmer_offset_R_tuples_L = []
        cont = True
        num_kmers = len( kmer_offset_R_T_Ls_by_classnum_D_1[motif_class_num] )
        for kmer_idx in range( num_kmers ):
            kmer_1, offset_1, excess_R_1 = kmer_offset_R_T_Ls_by_classnum_D_1[motif_class_num][kmer_idx]
            try:
                if (cont == True):
                    kmer_2, offset_2, excess_R_2 = kmer_offset_R_T_Ls_by_classnum_D_2[motif_class_num][kmer_idx]
                    if (kmer_2 == kmer_1):
                        avg_excess_R = (excess_R_1 + excess_R_2 ) / 2.
                        kmer_offset_R_tuples_L.append(
                                (kmer_1, offset_1, avg_excess_R ) )
                    else:
                        cont = False

            except (KeyError, IndexError):
                cont = False
                continue

        if (len( kmer_offset_R_tuples_L ) > 0):
            kmer_offset_R_tuples_Ls_by_classnum_D[consistent_class_idx] =\
                    kmer_offset_R_tuples_L
            consistent_class_idx += 1
            excess_R = sum( [x[2] for x in kmer_offset_R_tuples_L] )
            logos_excess_Rs_L.append( excess_R )

    #### Record the kmer_offset_R_tuples_Ls_by_classnum_D in a file
    ####    that can manually be adjusted if desired
    RBNS_logo_helpers.make_final_file_aligned_kmers_weights_per_motif_class(
            kmer_offset_R_tuples_Ls_by_classnum_D,
            settings.kmers_to_make_logos_Fs_by_half_D["consistent"] )

    returned_D = RBNS_logo_helpers.make_logos_from_kmer_offset_R_T_Ls_by_classnum_D(
            kmer_offset_R_tuples_Ls_by_classnum_D,
            settings.consistent_DIR,
            protein_name = settings.protein_name_for_plotting )
    pruned_pos_to_include_by_classnum_D = returned_D[
            "pruned_pos_to_include_by_classnum_D"]
    prob_logos_by_classnum_D = returned_D["prob_logos_by_classnum_D"]
    seq_logos_by_classnum_D = returned_D["seq_logos_by_classnum_D"]

    #### Pickle the pruned_pos_to_include_by_classnum_D
    out_D_F = settings.kmers_to_make_logos_Fs_by_half_D["consistent"].split(".txt")[0] +\
                    ".pruned_pos_to_include_by_classnum_D.pkl"
    with open( out_D_F, "wb" ) as out_f:
        pickle.dump( pruned_pos_to_include_by_classnum_D, out_f )

    #### Make the kmers_in_logo_F including ONLY the pruned positions
    RBNS_logo_helpers.make_final_file_aligned_kmers_weights_per_motif_class_only_pruned_positions(
            settings.kmers_to_make_logos_Fs_by_half_D["consistent"],
            pruned_pos_to_include_by_classnum_D )

    #### Go through each motif_class and add the logo & its excess_R
    ####     to all_logos_L & logos_excess_Rs_L, respectively
    motif_classes_L = prob_logos_by_classnum_D.keys()
    motif_classes_L.sort()
    for motif_class in motif_classes_L:

        prob_logo_F = prob_logos_by_classnum_D[motif_class]
        seq_logo_F = seq_logos_by_classnum_D[motif_class]

        all_logos_L.append( prob_logo_F )
        all_seq_logos_L.append( seq_logo_F )

    #### Now that all logos have been made, plot them all together in one
    ####    stacked bargraph w/ logos using the helper function in plots_bar
    out_combined_plots_F = os.path.join( settings.consistent_DIR,
            "{0}_{1}mer_logos.pdf".format(
                settings.protein_name_for_plotting,
                settings.starting_k ) )
    returned_D = RBNS_plots.plot_stacked_bargraph_with_logos(
            logos_excess_Rs_L,
            all_logos_L,
            out_combined_plots_F,
            title = settings.barplot_title )
    out_combined_seq_plots_F = os.path.join( settings.consistent_DIR,
            "{0}_{1}mer_seqlogos.pdf".format(
                settings.protein_name_for_plotting,
                settings.starting_k ) )
    returned_D = RBNS_plots.plot_stacked_bargraph_with_logos(
            logos_excess_Rs_L,
            all_seq_logos_L,
            out_combined_seq_plots_F,
            title = settings.barplot_title )






def get_one_set_logos_based_on_kmers_in_both(
        kmers_F_1,
        kmers_F_2,
        settings ):
    """
    - Using the constituent kmers, offsets, and R values of two sets of logos,
        each derived from half of the data, makes a
        composite logo consisting only of kmers that are IN BOTH (in any order)
        in the two halves

    - This function is called at the end of execute_get_logos()
    """

    all_logos_L = []
    all_seq_logos_L = []
    logos_excess_Rs_L = []

    kmer_offset_R_tuples_Ls_by_classnum_D =\
            RBNS_logo_helpers.get_in_both_kmers_F_from_2_indiv_notrequiring_same_logonum(
                    kmers_F_1, kmers_F_2 )
    num_motifs = len( kmer_offset_R_tuples_Ls_by_classnum_D )

    for motif_class in range( num_motifs ):
        kmer_offset_R_tuples_L = kmer_offset_R_tuples_Ls_by_classnum_D[motif_class]
        excess_R = sum( [x[2] for x in kmer_offset_R_tuples_L] )
        logos_excess_Rs_L.append( excess_R )

    #### Record the kmer_offset_R_tuples_Ls_by_classnum_D in a file
    ####    that can manually be adjusted if desired
    RBNS_logo_helpers.make_final_file_aligned_kmers_weights_per_motif_class(
            kmer_offset_R_tuples_Ls_by_classnum_D,
            settings.kmers_to_make_logos_Fs_by_half_D["in_both"] )

    returned_D = RBNS_logo_helpers.make_logos_from_kmer_offset_R_T_Ls_by_classnum_D(
            kmer_offset_R_tuples_Ls_by_classnum_D,
            settings.in_both_DIR,
            protein_name = settings.protein_name_for_plotting )
    pruned_pos_to_include_by_classnum_D = returned_D[
            "pruned_pos_to_include_by_classnum_D"]
    prob_logos_by_classnum_D = returned_D["prob_logos_by_classnum_D"]
    seq_logos_by_classnum_D = returned_D["seq_logos_by_classnum_D"]
    #### Pickle the pruned_pos_to_include_by_classnum_D
    out_D_F = settings.kmers_to_make_logos_Fs_by_half_D["in_both"].split(".txt")[0] +\
                    ".pruned_pos_to_include_by_classnum_D.pkl"
    with open( out_D_F, "wb" ) as out_f:
        pickle.dump( pruned_pos_to_include_by_classnum_D, out_f )

    #### Make the kmers_in_logo_F including ONLY the pruned positions
    RBNS_logo_helpers.make_final_file_aligned_kmers_weights_per_motif_class_only_pruned_positions(
            settings.kmers_to_make_logos_Fs_by_half_D["in_both"],
            pruned_pos_to_include_by_classnum_D )

    #### Go through each motif_class and add the logo & its excess_R
    ####     to all_logos_L & logos_excess_Rs_L, respectively
    motif_classes_L = prob_logos_by_classnum_D.keys()
    motif_classes_L.sort()
    for motif_class in motif_classes_L:

        prob_logo_F = prob_logos_by_classnum_D[motif_class]
        seq_logo_F = seq_logos_by_classnum_D[motif_class]

        all_logos_L.append( prob_logo_F )
        all_seq_logos_L.append( seq_logo_F )

    #### Now that all logos have been made, plot them all together in one
    ####    stacked bargraph w/ logos using the helper function in plots_bar
    out_combined_plots_F = os.path.join( settings.in_both_DIR,
            "{0}_{1}mer_logos.pdf".format(
                settings.protein_name_for_plotting,
                settings.starting_k ) )
    returned_D = RBNS_plots.plot_stacked_bargraph_with_logos(
            logos_excess_Rs_L,
            all_logos_L,
            out_combined_plots_F,
            title = settings.barplot_title )
    out_combined_seq_plots_F = os.path.join( settings.in_both_DIR,
            "{0}_{1}mer_seqlogos.pdf".format(
                settings.protein_name_for_plotting,
                settings.starting_k ) )
    returned_D = RBNS_plots.plot_stacked_bargraph_with_logos(
            logos_excess_Rs_L,
            all_seq_logos_L,
            out_combined_seq_plots_F,
            title = settings.barplot_title )





def make_PWM_file_from_aligned_kmers_F(
        aligned_kmers_F ):
    """
    - Given a (pruned) aligned kmers file, makes a .PWM flat text file
    """
    out_PWM_F = os.path.join( os.path.dirname( aligned_kmers_F ),
            os.path.basename( aligned_kmers_F ).split( "pruned_" )[-1].\
                    split(".aligned")[0] + ".PWM" )

    out_f = open( out_PWM_F, 'w' )

    with open( aligned_kmers_F ) as f:
        for line in f:
            num_pos = len( line.strip() )
            break

    #### Make a PWM_D
    PWM_counts_D = {}
    for i in range( num_pos ):
        PWM_counts_D[i] = {"A":0, "C":0, "G": 0, "U": 0}
    #### Populate the PWM_counts_D
    with open( aligned_kmers_F ) as f:
        for line in f:
            aligned_kmer = line.strip()
            for pos, nt in enumerate( aligned_kmer ):
                try:
                    PWM_counts_D[pos][nt] += 1
                except KeyError:
                    pass
    #### Convert the counts into freqs
    PWM_D = {}
    for i in PWM_counts_D:
        total_counts = sum( PWM_counts_D[i].values() )

        PWM_D[i] = {}
        for nt, count in PWM_counts_D[i].iteritems():
            PWM_D[i][nt] = float( count ) / total_counts

    pos_L = PWM_D.keys()
    pos_L.sort()
    nts_L = PWM_D[pos_L[0]].keys()
    nts_L.sort()
    #### RNA logos must have "U" in the 2nd row header line
    nts_L_header = ["A","C","G","U"]

    out_f.write( 'ID test\n' )
    out_f.write( "\t".join(["PO"] + nts_L_header) + "\n" )
    for position in pos_L:
        out_f.write(str(position))
        for letter in nts_L:
            out_f.write("\t" + str(PWM_D[position][letter]))
        out_f.write("\n")

    out_f.close()





def make_PWM_files_for_all_aligned_kmers_Fs_in_DIR(
        DIR ):
    """
    - Given a DIR, gets all pruned*.aligned_kmers.txt files, and makes a .PWM
        for each using the make_PWM_file_from_aligned_kmers_F() function above
    """

    #### Go through all of the pruned*.aligned_kmers.txt files in the DIR
    aligned_kmers_Fs_L = glob.glob( os.path.join( DIR, "pruned*.aligned_kmers.txt" ) )
    for aligned_kmers_F in aligned_kmers_Fs_L:
        make_PWM_file_from_aligned_kmers_F( aligned_kmers_F )


###############################################################################
################################## < UTILS > ##################################


def is_homopolymer( kmer ):
    """
    - Returns True or False, depending on if kmer is a homopolymer
    """
    nts_S = set( [x for x in kmer] )
    if (len( nts_S ) == 1):
        return True
    else:
        return False



def is_homopolymer_w_1_intervening( kmer ):
    """
    - Returns a dictionary depending on if kmer is a homopolymer with
        exactly 1 intervening nt

    - return_D ALWAYS has the key: "is_homopolymer_w_1_intervening"
    """

    nts_S = set( [x for x in kmer] )
    if (len( nts_S ) == 2):
        num_by_nt_D = {}
        for nt in kmer:
            try:
                num_by_nt_D[nt] += 1
            except KeyError:
                num_by_nt_D[nt] = 1

        nt_numtimes_T_L = [(nt, num_by_nt_D[nt]) for nt in num_by_nt_D]
        nt_numtimes_T_L.sort( key = lambda x: -1*x[1] )
        if (nt_numtimes_T_L[0][1] == len( kmer ) - 1) and\
                (nt_numtimes_T_L[1][1] == 1):
            major_nt = nt_numtimes_T_L[0][0]
            minor_nt = nt_numtimes_T_L[1][0]
            pos_of_minor_nt = kmer.find( minor_nt )
            return_D = {"is_homopolymer_w_1_intervening": True,
                    "major_nt": major_nt,
                    "minor_nt": minor_nt,
                    "pos_of_minor_nt": pos_of_minor_nt}
        else:
            return_D = {"is_homopolymer_w_1_intervening": False}
    else:
        return_D = {"is_homopolymer_w_1_intervening": False}
    return return_D

################################# </ UTILS > ##################################
###############################################################################






def parse_args():
    parser = argparse.ArgumentParser(
            description = 'Generates RNA Bind-n-Seq motif logos.')

    parser.add_argument("config_F")

    parser.add_argument("--num-reads", metavar='N', type = float, nargs = 1,
            default = [0.5],
            help = "Number of reads to use from each of input & pulldown .reads files.")

    #### args.starting_k will be a list of ints like: [[5], [6]]
    parser.add_argument("--starting-k", metavar='k', type = int, nargs = 1,
            action = "append",
            help = "The highest value of k to look for enriched kmers. Can repeat this option to run multiple jobs, each with a different starting k, e.g.:  --starting-k 5 --starting-k 6")

    #### args.ending_k will be a list of ints like: [[4], [3]]
    parser.add_argument("--ending-k", metavar='k', type = int, nargs = 1,
            action = "append",
            help = "The lowest value of k to look for enriched kmers; will search from k = starting-k -> ending-k, inclusive. Can repeat this option to run multiple jobs, each with a different ending k, e.g.:  --ending-k 3 --ending-k 4")

    #### args.Zscore_kmers_to_keep will be a list of floats like: [[2.], [3.]]
    parser.add_argument("--Zscore-kmers-to-keep", metavar='x', type = float, nargs = 1,
            action = "append",
            help = "The Z-score of frequency of reads that contain a kmer to use that kmer in one of the motifs. Default is 2, can increase to 3 or higher for only including more highly enriched kmers. Can repeat this option to run multiple jobs, each with a different Z-score, e.g.:  --Zscore-kmers-to-keep 2 -Zscore-kmers-to-keep 3")

    parser.add_argument("--launch-onto-cluster",
            action = "store_true",
            help = "If the jobs should separately each be launched onto the cluster.")

    args = parser.parse_args()
    return args












def main():
    """
    - Parses arguments and either runs logos jobs sequentially, or launches
        a separate job for each onto the cluster
    """
    args = parse_args()

    ##### Go through each of the args
    #####   Get the config_DIR
    config_DIR = os.path.dirname( args.config_F )
    config_starting_basename = os.path.basename( args.config_F ).rsplit(".", 1)[0]

    config_with_args_DIR = os.path.join( config_DIR, "with_args_added" )
    try:
        os.makedirs( config_with_args_DIR )
    except OSError: pass

    ########## < Provide default args. if they weren't passed in > ############

    ######################## < starting-k > ###################################
    if (args.starting_k is None):
        args.starting_k = [[6]]
        print "\t- Using default starting-k = 6"
    ######################## </ starting-k > ##################################

    ########################## < ending-k > ###################################
    if (args.ending_k is None):
        args.ending_k = [[4]]
        print "\t- Using default ending-k = 4"
    ########################## </ ending-k > ##################################

    ###### < Z-score enrichment cutoff for kmers to include in logos > ########
    if (args.Zscore_kmers_to_keep is None):
        args.Zscore_kmers_to_keep = [[2.]]
        print "\t- Using default Zscore-kmers-to-keep = 2.0"
    ##### < / Z-score enrichment cutoff for kmers to include in logos > #######

    ######### < / Provide default args. if they weren't passed in > ###########


    #### Launch a job for each k / Z-score combination
    for starting_k_L in args.starting_k:
        starting_k = starting_k_L[0]
        basename = config_starting_basename + "_k_{0}".format( starting_k )
        for ending_k_L in args.ending_k:
            ending_k = ending_k_L[0]
            basename += "_to_{0}".format( ending_k )
            if (ending_k > starting_k):
                continue
            for Zscore_kmers_to_keep_L in args.Zscore_kmers_to_keep:
                Zscore_kmers_to_keep = Zscore_kmers_to_keep_L[0]
                basename += "_Zscore_{0:.1f}".format(
                        Zscore_kmers_to_keep )
                #### add the arguments to the config_F
                this_config_F = os.path.join( config_with_args_DIR,
                        basename + ".config_F" )
                config_f = open( this_config_F, "w" )

                for line in open( args.config_F ):
                    if (line.find( "starting_k" ) != -1) or\
                            (line.find( "ending_k" ) != -1) or\
                            (line.find( "Zscore_kmers_to_keep" ) != -1) or\
                            (line.find( "num_" ) != -1):
                        continue
                    config_f.write( line )

                #### Write out the arguments
                config_f.write( "starting_k = {}\n".format( starting_k ) )
                config_f.write( "ending_k = {}\n".format( ending_k ) )
                config_f.write( "Zscore_kmers_to_keep = {}\n".format(
                    Zscore_kmers_to_keep ) )
                config_f.write( "num_input_pulldown_rds_to_use: {}\n".format(
                    args.num_reads[0] ) )

                config_f.close()

                ################# < RUN EXECUTE_GET_LOGOS > ###################
                if ( args.launch_onto_cluster == False ):
                    execute_get_logos(
                        this_config_F )
                else:
                    error_output_DIR = os.path.join( config_with_args_DIR,
                            basename, "errors_outputs" )
                    if not os.path.exists( error_output_DIR ):
                        os.makedirs( error_output_DIR )
                    #output_F = os.path.join( error_output_DIR, "output.txt" )
                    #error_F = os.path.join( error_output_DIR, "error.txt" )
                    #cmd = 'python {0} {1} --num-reads {2} --starting-k {3} --ending-k {4} --Zscore-kmers-to-keep {5} 1> {6} 2> {7}'.format(
                    #        os.path.abspath( sys.argv[0] ),
                    #        args.config_F,
                    #        args.num_reads[0],
                    #        starting_k,
                    #        ending_k,
                    #        Zscore_kmers_to_keep,
                    #        output_F,
                    #        error_F )
                    cmd = 'python {0} {1} --num-reads {2} --starting-k {3} --ending-k {4} --Zscore-kmers-to-keep {5}'.format(
                            os.path.abspath( sys.argv[0] ),
                            args.config_F,
                            args.num_reads[0],
                            starting_k,
                            ending_k,
                            Zscore_kmers_to_keep )
                    print cmd

                    launch(
                        cmd,
                        jobname = basename,
                        time_mins = 690,
                        error_DIR = error_output_DIR )

                ################# < RUN EXECUTE_GET_LOGOS > ###################






if __name__ == '__main__':
    main()





