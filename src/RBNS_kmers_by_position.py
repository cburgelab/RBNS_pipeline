#!/usr/bin/env python
import os, sys
import inspect
import math
import cPickle as pickle
import glob
import pprint


import RBNS_utils
import RBNS_plots





def analyze_freqs_by_position_one_barcodes_ordered_kmers_to_consider(
        ordered_kmers_to_consider_Ls_by_k_D,
        protein,
        main_DIR,
        conc_for_fastq,
        ks_L,
        ordered_kmers_description_fnames = "",
        make_output_Fs = True,
        num_controls = 20,
        max_log2_val_colormap = None):
    """
    - A helper function called by the
    analyze_freqs_by_position_all_barcodes_one_protein_top_enriched_kmers()
    function below

    - Using the previously pickled dictionaries from the
        get_counts_freqs_by_pos_one_F() function above, loads them and
        calculates the KL divergence of
        (Uniform across read || Observed freqs. across read) for each kmer,
        and outputs a .txt table with kmers in descending order of
        KL Divergence

    - INPUTs:
        - counts_by_pos_DIR:
            - directory that has pickled dictionaries, like:
                /net/uorf/data/backup/RBNS_results/srsf8/counts/by_position
            - basename_start:
                - e.g., "80", or "input_library"
            - pprint_lib_name:
                - a "nice" name to use for the title (e.g., "Input Library")
            - ordered_kmers_to_consider_L:
                - the kmers to consider, e.g., the sig. enriched kmers
            - add_to_end_title_str:
                - e.g., ", Perfect Adapter Reads Only"
            - make_output_Fs:
                - If True, makes a .txt out file and a plot;
                - If False, doesn't make .txt/.pdf (this is used the first time
                    around to get the maximum absolute log2 value so that on
                    the second time around when plots are made, all the
                    colorbars can be coordinated togeter)
            - ordered_kmers_description_fnames:
                - a name that will be added to the output PDF & .txt files to
                    distinguish it (e.g., "3std" to denote these kmers as those
                    with Z-score >= 3)
                    - "3.0std"
                    - "2.0std"
                    - "least"
                    - "Adapter 1" or "Adapter 2"
            - max_log2_val_colormap:
                - If passed in, the heatmap colorbar will go from
                    -max_log2_val_colormap to max_log2_val_colormap
    """
    return_D = {}

    frequency_Ds_DIR = os.path.join( main_DIR, "frequency_Ds" )
    RBNS_utils.make_dir( frequency_Ds_DIR )

    if (conc_for_fastq == "input"):
        conc_label = "Input lib."
    else:
        conc_label = "{} nM lib.".format( conc_for_fastq )

    #### go through each of the k's
    for k in ks_L:

        ordered_kmers_to_consider_L = ordered_kmers_to_consider_Ls_by_k_D[k]
        if ( len( ordered_kmers_to_consider_L ) == 0 ):
            return_D[k] = {}
            continue

        #### Load the previously pickled dictionary of kmer frequencies at each
        ####    position
        D_F = os.path.join( main_DIR,
                "frequency_Ds/{0}_{1}.{2}mer.frequencies.by_position.pkl".format(
                    protein, conc_for_fastq, k ) )
        with open( D_F ) as f:
            freqs_by_pos_kmer_D = pickle.load( f )
        num_kmers_per_rd = len( freqs_by_pos_kmer_D.keys() )

        #### A uniform distribution over all positions in the read
        uniform_L = [1./num_kmers_per_rd] * num_kmers_per_rd

        #### Now go through each of the kmers and get the KL divergence of
        ####    KLDiv( uniform || observed freqs. across read )
        kmer_KLDiv_tuples_L = []
        kmer_to_KLDiv_D = {}
        for kmer in RBNS_utils.yield_kmers( k ):
            obs_freqs_L = [freqs_by_pos_kmer_D[x][kmer] for x in\
                    range( num_kmers_per_rd )]
            sum_obs_freqs = sum( obs_freqs_L )
            #### Normalize the obs_freqs_L
            obs_freqs_L = [x/sum_obs_freqs for x in obs_freqs_L]

            #### Get the KL Divergence
            KL = RBNS_utils.KL_divergence(uniform_L, obs_freqs_L)
            kmer_to_KLDiv_D[kmer] = KL
            #### Also get the log2(OBSERVED/UNIFORM) at each position
            observed_over_unif_L = []
            for i in range( num_kmers_per_rd ):
                try:
                    observed_over_unif_L.append(
                        math.log(obs_freqs_L[i] / uniform_L[i], 2) )
                except ValueError:
                    observed_over_unif_L.append( 0. )

            kmer_KLDiv_tuples_L.append( (kmer, KL, observed_over_unif_L) )

        #### Sort the kmers by decreasing
        kmer_KLDiv_tuples_L.sort( key = lambda x: -1*x[1] )

        #### Get the mean KL divergence
        KL_divs = [x[1] for x in kmer_KLDiv_tuples_L]
        mean_KL, std_KL = RBNS_utils.mean_std( KL_divs )

        #### Extract the ordered_kmers_to_consider
        kmer_KLDiv_tuples_desiredkmers_L = []
        for kmer in ordered_kmers_to_consider_L:
            #### Go through kmer_KLDiv_tuples_L and get the tuple for this kmer
            for tupl in kmer_KLDiv_tuples_L:
                if (tupl[0] == kmer):
                    kmer_KLDiv_tuples_desiredkmers_L.append( tupl )
                    break

        #### Go through each of the significant kmers
        kmers_to_KL_Div_D = {}
        kmers_to_log2_Obs_over_Exp_L_D = {}
        #### the maximum absolute log2 value plotted, so that all of the
        ####    libraries can have the same colormap scale
        max_abs_log2_plotted = 0.
        for kmer, KL, observed_over_unif_L in kmer_KLDiv_tuples_desiredkmers_L:
            kmers_to_KL_Div_D[kmer] = KL
            kmers_to_log2_Obs_over_Exp_L_D[kmer] = observed_over_unif_L
            #### update max_abs_log2_plotted
            max_abs_log2_plotted = max( max_abs_log2_plotted,
                    max( observed_over_unif_L ), -1*min( observed_over_unif_L ))

        #### Get "control" kmer distributions that have KL divergence below
        ####    the mean
        ctrl_STD_tuples_L = [x for x in kmer_KLDiv_tuples_L if x[1] <
                mean_KL ]
        control_delta = int(len( ctrl_STD_tuples_L ) / float( num_controls ))
        #### go through and get the 20 evently spaced controls to plot
        control_tuples_L = [ctrl_STD_tuples_L[control_delta*x] for x in
                range(num_controls)]
        #### Go through each of the control kmers
        control_kmers_to_KL_Div_D = {}
        control_kmers_to_log2_Obs_over_Exp_L_D = {}
        for kmer, KL, observed_over_unif_L in control_tuples_L:
            control_kmers_to_KL_Div_D[kmer] = KL
            control_kmers_to_log2_Obs_over_Exp_L_D[kmer] = observed_over_unif_L
            #### update max_abs_log2_plotted
            max_abs_log2_plotted = max( max_abs_log2_plotted,
                    max( observed_over_unif_L ), -1*min( observed_over_unif_L ))

        #### add the max_abs_log2_plotted to the return_D
        return_D[k] = {"max_abs_log2_plotted": max_abs_log2_plotted,
                "kmers_to_KL_Div_D": kmers_to_KL_Div_D,
                "kmers_to_log2_Obs_over_Exp_L_D": kmers_to_log2_Obs_over_Exp_L_D,
                "control_kmers_to_KL_Div_D": control_kmers_to_KL_Div_D,
                "control_kmers_to_log2_Obs_over_Exp_L_D": control_kmers_to_log2_Obs_over_Exp_L_D}


        if make_output_Fs:

            #### Make the out_F
            out_DIR = os.path.join( main_DIR, "tables/by_position" )
            RBNS_utils.make_dir( out_DIR )

            #### < Make a table of KL div by decreasing R of sig. R kmers > ###
            out_basename = "{0}mers.{1}sig_R_{2}.KL_div_of_freqs_across_read.txt".format(
                    k, conc_label.replace( " ", "_" ),
                    ordered_kmers_description_fnames.replace(" ", "_") )
            out_F = os.path.join( out_DIR, out_basename )
            with open( out_F, "w" ) as f_out:
                #### write a header line
                f_out.write("{0}: {1}\n".format(protein.replace("_", " "), conc_label ))
                f_out.write("\tKL Div(Uniform||Observed)\t\tlog2(Obs/Unif) at Pos 1\tPos. 2\n")

                #### Go through and write out all of the kmers
                for kmer in ordered_kmers_to_consider_L:
                    for other_kmer, KL, observed_over_unif_L in kmer_KLDiv_tuples_L:
                        if (other_kmer == kmer):
                            f_out.write("\n{0}\t{1:.4g}\t\t".format( kmer, KL ))
                            for ratio in observed_over_unif_L:
                                f_out.write( "{0:.3f}\t".format( ratio ) )

            ### </ Make a table of KL div by decreasing R of sig. R kmers > ###

            #### try to get the Z-score threshold for sig. enrichment
            kmers_type_annot = "Enriched $k$mers"
            try:
                num_std = int( float(ordered_kmers_description_fnames.split("std")[0]) )
                kmers_type_annot += (" (Z-score $\geq${})".format( num_std ))
            except ValueError:
                if (ordered_kmers_description_fnames == "least"):
                    kmers_type_annot = "Least Enriched $k$mers"
                else:
                    kmers_type_annot = ordered_kmers_description_fnames

            #### Make the plot
            returned_fig_D = RBNS_plots.make_rectangular_heatmap_plot_RBNS_freqs(
                    kmers_to_KL_Div_D,
                    kmers_to_log2_Obs_over_Exp_L_D,
                    control_kmers_to_KL_Div_D,
                    control_kmers_to_log2_Obs_over_Exp_L_D,
                    order_of_kmers_L = ordered_kmers_to_consider_L,
                    kmers_type_annot = kmers_type_annot,
                    title = "{0}: {1} {2}mer frequencies across reads".format(
                        protein.replace("_", " "), conc_label, k ),
                    colorbar_label = r"$log_2$(Observed / Uniform freq.)",
                    max_log2_val_colormap = max_log2_val_colormap )
            return_D[k]["fig"] = returned_fig_D["fig"]

    return return_D







def analyze_freqs_by_position_one_library(
        protein,
        main_DIR,
        conc_for_fastq,
        ks_L,
        make_output_Fs = True,
        num_controls = 20,
        max_log2_val_colormap = None):
    """
    -   Calculates the KL divergence of
        (Uniform across read || Observed freqs. across read) for each kmer,
        and outputs a .txt table with kmers in descending order of
        KL Divergence

    - INPUTs:
            - make_output_Fs:
                - If True, makes a .txt out file and a plot;
                - If Flase, doesn't make .txt/.pdf (this is used the first time
                    around to get the maximum absolute log2 value so that on
                    the second time around when plots are made, all the
                    colorbars can be coordinated togeter)
            - max_log2_val_colormap:
                - If passed in, the heatmap colorbar will go from
                    -max_log2_val_colormap to max_log2_val_colormap
    """
    return_D = {}

    if (conc_for_fastq == "input"):
        conc_label = "Input lib."
    else:
        conc_label = "{} nM lib.".format( conc_for_fastq )

    frequency_Ds_DIR = os.path.join( main_DIR, "frequency_Ds" )
    RBNS_utils.make_dir( frequency_Ds_DIR )

    #### go through each of the k's
    for k in ks_L:

        #### Load the previously pickled dictionary of kmer frequencies at each
        ####    position
        D_F = os.path.join( frequency_Ds_DIR,
                "{0}_{1}.{2}mer.frequencies.by_position.pkl".format(
                    protein,
                    conc_for_fastq,
                    k ) )
        with open( D_F ) as f:
            freqs_by_pos_kmer_D = pickle.load( f )
        num_kmers_per_rd = len( freqs_by_pos_kmer_D.keys() )

        #### A uniform distribution over all positions in the read
        uniform_L = [1./num_kmers_per_rd] * num_kmers_per_rd

        #### Now go through each of the kmers and get the KL divergence of
        ####    KLDiv( uniform || observed freqs. across read )
        kmer_KLDiv_tuples_L = []
        kmer_to_KLDiv_D = {}
        for kmer in RBNS_utils.yield_kmers( k ):
            obs_freqs_L = [freqs_by_pos_kmer_D[x][kmer] for x in range(
                num_kmers_per_rd )]
            sum_obs_freqs = sum( obs_freqs_L )
            #### Normalize the obs_freqs_L
            obs_freqs_L = [x/sum_obs_freqs for x in obs_freqs_L]

            #### Get the KL Divergence
            KL = RBNS_utils.KL_divergence(uniform_L, obs_freqs_L)
            kmer_to_KLDiv_D[kmer] = KL
            #### Also get the log2(OBSERVED/UNIFORM) at each position
            try:
                observed_over_unif_L = [math.log(obs_freqs_L[i] / uniform_L[i], 2)\
                        for i in range( num_kmers_per_rd )]
            except ValueError:
                observed_over_unif_L = []
                for i in range( num_kmers_per_rd ):
                    try:
                        observed_over_unif_L.append( math.log(obs_freqs_L[i] / uniform_L[i], 2) )
                    except ValueError:
                        observed_over_unif_L.append( 1. )

            kmer_KLDiv_tuples_L.append( (kmer, KL, observed_over_unif_L) )

        #### Sort the kmers by decreasing
        kmer_KLDiv_tuples_L.sort( key = lambda x: -1*x[1] )

        #### Get the mean KL divergence
        KL_divs = [x[1] for x in kmer_KLDiv_tuples_L]
        mean_KL, std_KL = RBNS_utils.mean_std( KL_divs )
        #### a 3 STD threshold for the "most unequal"
        three_STD_thresh = mean_KL + (3 * std_KL)
        neg1_STD_thresh = mean_KL - std_KL

        #### Go through and get the kmers & KL divergences for those that are
        ####    >= 3 STD
        three_STD_tuples_L = kmer_KLDiv_tuples_L[:30]
        #three_STD_tuples_L = [x for x in kmer_KLDiv_tuples_L if x[1] >=\
        #        three_STD_thresh ]
        #print "{0}".format( len(three_STD_tuples_L) )
        #### Go through each of the significant kmers
        sig_kmers_to_KL_Div_D = {}
        sig_kmers_to_log2_Obs_over_Exp_L_D = {}
        #### the maximum absolute log2 value plotted, so that all of the
        ####    libraries can have the same colormap scale
        max_abs_log2_plotted = 0.
        for kmer, KL, observed_over_unif_L in three_STD_tuples_L:
            sig_kmers_to_KL_Div_D[kmer] = KL
            sig_kmers_to_log2_Obs_over_Exp_L_D[kmer] = observed_over_unif_L
            #### update max_abs_log2_plotted
            max_abs_log2_plotted = max( max_abs_log2_plotted,
                    max( observed_over_unif_L ), -1*min( observed_over_unif_L ))

        #### Get "control" kmer distributions that have KL divergence below
        ####    the mean
        ctrl_STD_tuples_L = [x for x in kmer_KLDiv_tuples_L if x[1] <
                mean_KL ]
        control_delta = int(len( ctrl_STD_tuples_L ) / float( num_controls ))
        #### go through and get the 20 evently spaced controls to plot
        control_tuples_L = [ctrl_STD_tuples_L[control_delta*x] for x in
                range(num_controls)]
        #### Go through each of the control kmers
        low_kmers_to_KL_Div_D = {}
        low_kmers_to_log2_Obs_over_Exp_L_D = {}
        for kmer, KL, observed_over_unif_L in control_tuples_L:
            low_kmers_to_KL_Div_D[kmer] = KL
            low_kmers_to_log2_Obs_over_Exp_L_D[kmer] = observed_over_unif_L
            #### update max_abs_log2_plotted
            max_abs_log2_plotted = max( max_abs_log2_plotted,
                    max( observed_over_unif_L ), -1*min( observed_over_unif_L ))

        #### add the max_abs_log2_plotted to the return_D
        return_D[k] = {"max_abs_log2_plotted": max_abs_log2_plotted}

        if (make_output_Fs == True):
            #### Make the out_F
            out_DIR = os.path.join( main_DIR, "tables/by_position" )
            RBNS_utils.make_dir( out_DIR )

            out_basename = "{0}mers.{1}greatest_KL_div_of_freqs_across_read.txt".format(
                    k, conc_label.replace( " ", "_" ) )
            out_F = os.path.join( out_DIR, out_basename )
            with open( out_F, "w" ) as f_out:
                #### write a header line
                f_out.write("{0}: {1}\n".format(protein, conc_label ))
                f_out.write("\tKL Div(Uniform||Observed)\t\tlog2(Obs/Unif) at Pos 1\tPos. 2\n")

                #### Go through and write out all of the kmers
                for kmer, KL, observed_over_unif_L in kmer_KLDiv_tuples_L:
                    f_out.write("\n{0}\t{1:.4g}\t\t".format( kmer, KL ))
                    for ratio in observed_over_unif_L:
                        f_out.write( "{0:.3f}\t".format( ratio ) )

            return_D[k] = {"out_F": out_F,
                    "kmer_to_KLDiv_D": kmer_to_KLDiv_D}

            #### Make a plot using the helper function in
            ####    /helpers/python_helpers/plots.py
            returned_fig_D = RBNS_plots.make_rectangular_heatmap_plot_RBNS_freqs(
                    sig_kmers_to_KL_Div_D,
                    sig_kmers_to_log2_Obs_over_Exp_L_D,
                    low_kmers_to_KL_Div_D,
                    low_kmers_to_log2_Obs_over_Exp_L_D,
                    title = "{0}: {1} {2}mer frequencies across reads".format(
                        protein.replace( "_", " " ), conc_label, k ),
                    colorbar_label = r"$log_2$(Observed / Uniform freq.)",
                    max_log2_val_colormap = max_log2_val_colormap)
            return_D[k]["fig"] = returned_fig_D["fig"]

    return return_D







def plot_abs_ratio_of_kmers_at_each_position_relative_to_input_lib(
        ordered_kmers_to_consider_Ls_by_k_D,
        bottom_kmers_to_consider_Ls_by_k_D,
        protein,
        main_DIR,
        conc_for_fastqs_L,
        ks_L,
        ordered_kmers_description_fnames = "",
        make_output_Fs = True,
        max_log2_val_colormap = None ):
    """
    - Makes a plot of the enrichment of the specified kmers in
        ordered_kmers_to_consider_Ls_by_k_D (typically Z-score > 3 or 2)
        as well as the bottom kmers, where enrichment at each position
        is defined as the freq. at that position /
        the average input freq ( averaged over all position )
    """
    assert( conc_for_fastqs_L[0] == "input" )

    return_D = {}

    frequency_Ds_DIR = os.path.join( main_DIR, "frequency_Ds" )

    #### go through each of the k's
    for k in ks_L:

        max_abs_log2_plotted = 0.
        freqs_by_conc_pos_kmer_D = {}

        abs_freqs_figs_L = []

        ordered_kmers_to_consider_L = ordered_kmers_to_consider_Ls_by_k_D[k]
        bottom_kmers_to_consider_L = bottom_kmers_to_consider_Ls_by_k_D[k]
        if ( len( ordered_kmers_to_consider_L ) == 0 ):
            return_D[k] = {}
            continue

        #### Go through each of the concentrations to load the freq Ds
        for conc_for_fastq in conc_for_fastqs_L:

            #### Load the previously pickled dictionary of kmer frequencies at each
            ####    position
            D_F = os.path.join( main_DIR,
                    "frequency_Ds/{0}_{1}.{2}mer.frequencies.by_position.pkl".format(
                        protein, conc_for_fastq, k ) )
            with open( D_F ) as f:
                freqs_by_pos_kmer_D = pickle.load( f )
            freqs_by_conc_pos_kmer_D[conc_for_fastq] = freqs_by_pos_kmer_D
            num_kmers_per_rd = len( freqs_by_pos_kmer_D.keys() )

        #### Go through each of the concentrations
        for conc_for_fastq in conc_for_fastqs_L:

            if (conc_for_fastq == "input"):
                conc_label = "Input lib."
            else:
                conc_label = "{} nM lib.".format( conc_for_fastq )

            input_freqs_by_pos_kmer_D = freqs_by_conc_pos_kmer_D['input']
            thisconc_freqs_by_pos_kmer_D = freqs_by_conc_pos_kmer_D[conc_for_fastq]

            #### Now go through each of the ordered kmers
            sig_kmers_to_overall_log2R_D = {}
            sig_kmers_to_log2_R_L_D = {}
            bottom_kmers_to_overall_log2R_D = {}
            bottom_kmers_to_log2_R_L_D = {}
            for kmers_L, kmers_to_overall_log2R_D, kmers_to_log2_R_L_D in\
                    [(ordered_kmers_to_consider_L,
                        sig_kmers_to_overall_log2R_D,
                        sig_kmers_to_log2_R_L_D),
                     (bottom_kmers_to_consider_L,
                         bottom_kmers_to_overall_log2R_D,
                         bottom_kmers_to_log2_R_L_D)]:
                for kmer in kmers_L:

                    #### Get the average freq. of this kmer in the input (averaging
                    ####    over all start positions)
                    input_freqs_L = [input_freqs_by_pos_kmer_D[idx][kmer] for\
                            idx in range( num_kmers_per_rd )]
                    avg_input_freq = sum( input_freqs_L ) / num_kmers_per_rd

                    #### Now go through and
                    thisconc_freqs_L = [thisconc_freqs_by_pos_kmer_D[idx][kmer] for\
                            idx in range( num_kmers_per_rd )]
                    avg_thisconc_freq = sum( thisconc_freqs_L ) / num_kmers_per_rd

                    overall_R = avg_thisconc_freq / avg_input_freq
                    kmers_to_overall_log2R_D[kmer] = math.log( overall_R, 2. )

                    #### Go through each position
                    log2_R_L = []
                    for pos, this_conc_freq in enumerate( thisconc_freqs_L ):
                        log2_R_L.append( math.log( this_conc_freq / avg_input_freq, 2. ) )
                    max_abs_log2_plotted = max( max_abs_log2_plotted,
                            max( [abs(x) for x in log2_R_L] ) )
                    kmers_to_log2_R_L_D[kmer] = log2_R_L

            ##### Make a plot with this set of kmers
            kmers_type_annot = "Enriched $k$mers"
            try:
                num_std = int( float(ordered_kmers_description_fnames.split("std")[0]) )
                kmers_type_annot += (" (Z-score $\geq${})".format( num_std ))
            except ValueError:
                pass
            returned_D = RBNS_plots.make_rectangular_heatmap_plot_RBNS_freqs(
                    sig_kmers_to_overall_log2R_D,
                    sig_kmers_to_log2_R_L_D,
                    bottom_kmers_to_overall_log2R_D,
                    bottom_kmers_to_log2_R_L_D,
                    kmers_type_annot = kmers_type_annot,
                    title = "{0} {1}".format( protein, conc_label ),
                    colorbar_label = r"$\log_2$(freq / avg. input freq.)",
                    order_of_kmers_L = ordered_kmers_to_consider_L,
                    order_of_control_kmers_L = bottom_kmers_to_consider_L,
                    annot_above_KL_heatmap_bar = r"$R$",
                    right_bar_color_by = "R",
                    control_kmers_type_annot = "Least enriched {}mers".format( k ),
                    log2_colormap = "PuOr",
                    KLdiv_colormap = "PuOr",
                    max_log2_val_colormap = max_log2_val_colormap )
            abs_freqs_figs_L.append( returned_D['fig'] )


        #### add the max_abs_log2_plotted to the return_D
        return_D[k] = {"max_abs_log2_plotted": max_abs_log2_plotted,
                "abs_freqs_figs_L": abs_freqs_figs_L }


        if make_output_Fs:

            #### Make the out_F
            out_DIR = os.path.join( main_DIR, "tables/by_position" )
            RBNS_utils.make_dir( out_DIR )

            #### Make a table of KL div by decreasing R of sig. R kmers
            out_basename = "{0}.{1}mers.sig_R_{2}.abs_freqs_across_read.pdf".format(
                    protein, k, ordered_kmers_description_fnames.replace(" ", "_") )
            out_F = os.path.join( out_DIR, out_basename )

            RBNS_plots.plot_multiple_inOnePDF(
                    abs_freqs_figs_L,
                    out_F )

    return return_D







