#!/usr/bin/env python
import os, sys
import itertools
import subprocess
import pprint

import RBNS_utils
import RBNS_settings

nts_L = ["A", "C", "G", "U"]



###############################################################################
######################### < CALLED BY RBNS_logos.py > #########################


def make_logos_from_kmer_offset_R_T_Ls_by_classnum_D(
        kmer_offset_R_T_Ls_by_classnum_D,
        out_motifs_DIR,
        protein_name = "",
        total_kmers_in_aligned_F = 10000,
        freq_to_include_in_logo = 0.25 ):
        #freq_to_include_in_logo = 0.0 ):
    """
    - INPUT:
        kmer_PWMs_by_k_classnum_D: see, e.g., the
        add_kmer_to_existing_kmer_PWMs_by_classnum_D_or_start_new_class()
        function in ~/RBNS_pipeline/RBNS_logos.py:

        kmer_PWMs_by_classnum_D[class_num] = {
            "PWM_D": PWM_D,
            "founding_kmer": kmer,
            "excess_R": excess_R,
            "kmer_offset_R_tuples_L":
                kmer_offset_R_tuples_L,
            "score_B_by_pos_D": score_B_by_pos_D}
    """
    os.system( "mkdir -p {}".format( out_motifs_DIR ) )

    motif_classes_L = kmer_offset_R_T_Ls_by_classnum_D.keys()
    motif_classes_L.sort()

    prob_logos_by_classnum_D = {}
    seq_logos_by_classnum_D = {}
    pruned_pos_to_include_by_classnum_D = {}

    aligned_kmers_F = ""

    for motif_class_num in motif_classes_L:

        print "motif_class_num: {}".format( motif_class_num )

        kmer_offset_R_tuples_L = kmer_offset_R_T_Ls_by_classnum_D[motif_class_num]
        pprint.pprint( kmer_offset_R_tuples_L )
        total_excess_R = sum( x[2] for x in kmer_offset_R_tuples_L )

        k = len( kmer_offset_R_tuples_L[0][0] )

        #### The starter basename for the files
        starting_basename = "{0}_{1}mer_logo{2}".format( protein_name, k, motif_class_num )
        aligned_kmers_F = os.path.join( out_motifs_DIR,
                starting_basename + ".aligned_kmers.txt" )

        #### Get the excess offset, each of which is repeated before (e.g., -2, -1)
        ####    and after (e.g., 5, 6 for k = 5)
        offsets_L = list(set([x[1] for x in kmer_offset_R_tuples_L]) )
        offsets_L.sort()
        print "offsets_L:"
        pprint.pprint( offsets_L )
        #### The max_offset (NOT the max. of absolute values - e.g., if it's
        ####    [-2, -1, 0, 1], the max_offset would be 1
        max_offset = max( offsets_L )
        print "max_offset: {}".format( max_offset )

        if ( len( offsets_L ) == 1 ):
            total_length = abs( offsets_L[0] ) + k
        else:
            first_offset = offsets_L[0]
            last_offset = offsets_L[-1]
            #### If the first is 0, the last must be positive
            if ( first_offset == 0 ):
                total_length = offsets_L[-1] + k
            #### If the last is 0, the first must be negative
            elif ( last_offset == 0 ):
                total_length = abs(offsets_L[0]) + k
            else:
                #### If they're both the same sign
                if ( ( first_offset * last_offset ) > 0 ):
                    total_length = max( [abs(x) for x in offsets_L] ) + k
                #### If they're the opposite sign
                else:
                    total_length = abs( offsets_L[-1] ) + abs( offsets_L[0] ) + k
        print "total_length: {}".format( total_length )

        # now convert the counts into numbers that sum up to (approximately)
        # total_kmers_in_txt_file
        num_blanks_by_pos_D = {}
        for pos in range( total_length ):
            num_blanks_by_pos_D[pos] = 0.

        with open( aligned_kmers_F, "w" ) as f:
            for kmer, offset, excess_R in kmer_offset_R_tuples_L:

                total_this_kmer = int(excess_R * total_kmers_in_aligned_F /\
                        total_excess_R )
                # add any leading or lagging -
                #### e.g., if the max_offset is 1 but this offset is -2, there
                ####    should be 3 before; if the max_offset is 1, min_offset is -1
                ####    (for 5mers so that total length is 7), for:
                ####        1: there should be 2 before
                ####        -1: there should be 0 before
                num_before = total_length - max_offset + offset - k
                print "{0}: num_before {1}".format( kmer, num_before )
                #num_before = max_offset + offset
                num_after = total_length - len(kmer) - num_before
                print "\tnum_before {0}".format( num_after )
                string_for_kmer_w_dashes = "-"*num_before + kmer.replace("T", "U") + "-"*num_after

                for i in range( total_this_kmer ):
                    str_to_write_out = string_for_kmer_w_dashes.replace("-", nts_L[i%4] )
                    f.write( str_to_write_out + "\n" )
                for pos, char in enumerate( string_for_kmer_w_dashes ):
                    if (char == "-"):
                        num_blanks_by_pos_D[pos] += total_this_kmer

        pos_to_include_L = []
        for pos in range( total_length ):
            freq_nts_this_pos = 1 - (num_blanks_by_pos_D[pos] /\
                    total_kmers_in_aligned_F)
            if (freq_nts_this_pos >= freq_to_include_in_logo ):
                print "Pos. {0} accepted w/ freq_nts_this_pos: {1}".format(
                       pos, freq_nts_this_pos )
                pos_to_include_L.append( pos )
            else:
                print "Pos. {0} REJECTED w/ freq_nts_this_pos: {1}".format(
                       pos, freq_nts_this_pos )
        pos_to_include_L.sort()
        pruned_pos_to_include_by_classnum_D[motif_class_num] = pos_to_include_L
        pruned_idx_start = pos_to_include_L[0]

        #### Make a pruned_aligned_kmers_F
        pruned_aligned_kmers_F = os.path.join( os.path.dirname( aligned_kmers_F ),
            "pruned_" + os.path.basename( aligned_kmers_F ))
        pruned_f = open( pruned_aligned_kmers_F, "w" )
        for line in open( aligned_kmers_F ):
            pruned_f.write( line.strip()[pruned_idx_start:
                (pruned_idx_start+len(pos_to_include_L))] + "\n" )
        pruned_f.close()

        #### Make the sequence & probability logos for this motif_class_num
        returned_D = make_seqlogo_from_aligned_kmers_L(
            out_motifs_DIR,
            starting_basename,
            pruned_aligned_kmers_F )
        prob_logos_by_classnum_D[motif_class_num] =\
                returned_D["prob_logo_out_F"]
        seq_logos_by_classnum_D[motif_class_num] =\
                returned_D["seq_logo_out_F"]

    return_D = {"prob_logos_by_classnum_D": prob_logos_by_classnum_D,
            "seq_logos_by_classnum_D": seq_logos_by_classnum_D,
            "pruned_pos_to_include_by_classnum_D": pruned_pos_to_include_by_classnum_D,
            "aligned_kmers_F": aligned_kmers_F }
    pprint.pprint( pruned_pos_to_include_by_classnum_D )
    return return_D






def make_seqlogo_from_aligned_kmers_L(
        out_DIR,
        basename_start,
        aligned_F,
        title = "",
        starting_pos = None):
    """
    INPUT:
        - out_startname, the starter path to the files that will be produced
        ( should be like ".../PTB2_set0" ("_seqlogo.png" or similar will be
        appended to seq_logo_out_startname)
        - pwm_f, a position weight matrix file that was created by the
          make_PWM_from_dict() function above
    - See http://weblogo.threeplusone.com/manual.html#CLI for options
    """
    for scale_width in [True, False]:
        ###################### < SEQ LOGO, no axes > #######################
        if scale_width is False:
            seq_logo_out_F = os.path.join( out_DIR, basename_start + ".seq_logo.noaxes.png" )
        else:
            seq_logo_out_F = os.path.join( out_DIR, basename_start +\
                    ".seq_logo.cols_scaled_by_num_nt.noaxes.png" )
        cmd_seqlogo = ['weblogo']
        #### if starting_pos other than 1 was passed in, add it
        if (starting_pos != None):
            cmd_seqlogo += ['--first-index {0}'.format( starting_pos )]
        cmd_seqlogo += [
           "-U 'bits'",
           "--composition 'equiprobable'",
           "--weight 0",
           "--format png",
           "-n 200",
           "--stack-width 40",
           "--sequence-type rna",
           "--resolution 900",
           "--size large",
           "--errorbars NO",
           "--show-yaxis NO",
           "--show-xaxis NO",
           "--fineprint ''",
           "--scale-width YES",
           "--title '{0}'".format(title),
           "--color-scheme classic",
           "<",
           aligned_F,
           ">",
           seq_logo_out_F
            ]
        #### If the width SHOULDN'T be scaled
        if scale_width is False:
            cmd_seqlogo = cmd_seqlogo[:-4] + ["--scale-width", "NO"] + cmd_seqlogo[-4:]

        p = subprocess.Popen(' '.join( cmd_seqlogo ), shell=True)
        p.wait()
        ###################### </ SEQ LOGO, no axes > ######################


        ########################### < SEQ LOGO > ###########################
        if scale_width is False:
            seq_logo_out_F = os.path.join( out_DIR, basename_start + ".seq_logo.png" )
        else:
            seq_logo_out_F = os.path.join( out_DIR, basename_start +\
                    ".seq_logo.cols_scaled_by_num_nt.png" )
        cmd_seqlogo = ['weblogo']
        #### if starting_pos other than 1 was passed in, add it
        if (starting_pos != None):
            cmd_seqlogo += ['--first-index {0}'.format( starting_pos )]
        cmd_seqlogo += [
           "-U 'bits'",
           "--composition 'equiprobable'",
           "--weight 0",
           "--format png",
           "-n 200",
           "--stack-width 40",
           "--sequence-type rna",
           "--resolution 900",
           "--size large",
           "--errorbars NO",
           "--fineprint ''",
           "--scale-width YES",
           "--title '{0}'".format(title),
           "--color-scheme classic",
           "<",
           aligned_F,
           ">",
           seq_logo_out_F
            ]
        #### If the width SHOULDN'T be scaled
        if scale_width is False:
            cmd_seqlogo = cmd_seqlogo[:-4] + ["--scale-width", "NO"] + cmd_seqlogo[-4:]

        p = subprocess.Popen(' '.join( cmd_seqlogo ), shell=True)
        p.wait()
        ########################### </ SEQ LOGO > ##########################


        ########################### < PROB. LOGO > #########################
        #### Get the probability logo
        if scale_width is False:
            prob_logo_out_F = os.path.join( out_DIR, basename_start + ".prob_logo.png" )
        else:
            prob_logo_out_F = os.path.join( out_DIR, basename_start +\
                    ".prob_logo.cols_scaled_by_num_nt.png" )

        cmd_probslogo = ['weblogo']
        #### if starting_pos other than 1 was passed in, add it
        if (starting_pos != None):
            cmd_probslogo += ['--first-index {0}'.format( starting_pos )]
        cmd_probslogo += [
           '-U',
           "'probability'",
           "--composition",
           "'equiprobable'",
           "--format png",
           "--sequence-type rna",
           "--resolution 900",
           "--size large",
           "--errorbars NO",
           "--fineprint ''",
           "--title '{0}'".format(title),
           "--color-scheme classic",
           "<",
           aligned_F,
           ">",
           prob_logo_out_F
            ]
        #### If the width SHOULDN'T be scaled
        if scale_width is False:
            cmd_probslogo = cmd_probslogo[:-4] + ["--scale-width", "NO"] + cmd_probslogo[-4:]
        p = subprocess.Popen(' '.join( cmd_probslogo ), shell=True)
        p.wait()
        ########################## </ PROB. LOGO > #########################

        return_D = {"seq_logo_out_F": seq_logo_out_F,
                "prob_logo_out_F": prob_logo_out_F}
    return return_D




######################### </ CALLED BY RBNS_logos.py > ########################
###############################################################################







def get_kmer_excessR_regardlessofmotifnum_from_F(
        kmers_F ):
    """
    - Returns a list of tuples of: [(kmer, excess_R), ...],
        not including the offset of each kmer or which motif_class it was from
    """
    kmer_excessR_T_L = []
    with open( kmers_F ) as f:
        for line in f:
            try:
                kmer_w_flanking, excess_R = line.strip().split(" ")
                excess_R = float( excess_R )
                kmer = kmer_w_flanking.strip("_").rstrip("_")
                kmer_excessR_T_L.append( ( kmer, excess_R ) )
            except:
                pass
    kmer_excessR_T_L.sort( key = lambda x: -1 * x[1] )
    return kmer_excessR_T_L






def get_kmer_excessR_of_kmers_in_both(
        kmer_excessR_T_L_1,
        kmer_excessR_T_L_2 ):
    """
    - Give two tuple lists like:
        [("GCATG", 1.6), ...],
        returns a list similar to that including ONLY the kmers that are in
        both, averaging the R's
    """
    kmer_excessR_T_L_in_both = []

    for kmer_1, excess_R_1 in kmer_excessR_T_L_1:
        for kmer_2, excess_R_2 in kmer_excessR_T_L_2:
            if ( kmer_1 == kmer_2 ):
                avg_R = ( excess_R_1 + excess_R_2 ) / 2.
                kmer_excessR_T_L_in_both.append( ( kmer_1, avg_R ) )

    kmer_excessR_T_L_in_both.sort( key = lambda x: -1 * x[1] )
    return kmer_excessR_T_L_in_both




def get_in_both_kmers_F_from_2_indiv_notrequiring_same_logonum(
        kmers_for_logo_F_1,
        kmers_for_logo_F_2 ):
    """
    - Given 2 kmers_for_logo_F, makes a composite one that includes kmers
        in both, NOT requiring that they be in the same motif #'s in the
        two kmers_for_logo_F's
    - Returns: kmer_offset_R_T_Ls_by_classnum_D
    """
    #### Return lists of tuples like:
    kmer_R_T_L_1 = get_kmer_excessR_regardlessofmotifnum_from_F(
            kmers_for_logo_F_1 )
    kmer_R_T_L_2 = get_kmer_excessR_regardlessofmotifnum_from_F(
            kmers_for_logo_F_2 )
    ##### Get the kmers & avg. excess_R for kmers that occur in BOTH
    kmer_R_T_L_in_both = get_kmer_excessR_of_kmers_in_both(
            kmer_R_T_L_1,
            kmer_R_T_L_2 )

    #### Now that we have the kmers in both and their average excess_R,
    ####    align them (and split them up into motif classes)
    kmer_offset_R_T_Ls_by_classnum_D = {}
    #### Go through all kmers
    for kmer_to_align, excess_R in kmer_R_T_L_in_both:
        returned_D = add_kmer_to_existing_motif_class_or_start_new(
                kmer_to_align,
                excess_R,
                kmer_offset_R_T_Ls_by_classnum_D )
        #### Get the updated kmer_offset_R_T_Ls_by_classnum_D
        kmer_offset_R_T_Ls_by_classnum_D = returned_D['kmer_offset_R_T_Ls_by_classnum_D']

    return kmer_offset_R_T_Ls_by_classnum_D







def make_final_file_aligned_kmers_weights_per_motif_class(
        kmer_offset_R_T_Ls_by_classnum_D,
        out_aligned_kmers_F ):
    """
    - After each of the "half" logos are made, this makes a succinct .txt file
        of the aligned kmers and excess weights that went into making each
        motif class
    - The out_aligned_kmers_F can be manually adjusted if necessary (i.e.,
        if alignments were not optimal),
        and then passed into the make_logos_from_manually_aligned_file()
        function in ~/python_lib/RBNS_logo_helpers.py to make each of the individual
        motifs as well as the bar plot with the final percentage assigned to
        each motif


    INPUT:

    - kmer_offset_R_T_Ls_by_classnum_D is like:
        {{0: [ ('ATACA', 0, 0.7520439188674823),
                ('ACACA', 2, 0.7199040161693568),
                ('CAATA', -2, 0.5747841800013938),
                ('CTACA', 0, 0.5363634377737869),
                ('TTACA', 0, 0.5255695916563865)],
    - out_aligned_kmers_F:
        the .txt output file that will be made

    """
    out_f = open( out_aligned_kmers_F, 'w' )

    for motif_num in range( len( kmer_offset_R_T_Ls_by_classnum_D ) ):

        kmer_offset_R_tuples_L = kmer_offset_R_T_Ls_by_classnum_D[motif_num]

        #### Get the lowest (most negative) offset, which is where all
        ####    kmers will start)
        lowest_offset = min( [x[1] for x in kmer_offset_R_tuples_L] )
        max_positions_after_0 = max(
                [len(x[0]) + x[1] for x in kmer_offset_R_tuples_L] )
        total_num_aligned_chars = abs( lowest_offset ) + max_positions_after_0

        out_f.write( "#" * 33 + " < motif {} > ".format( motif_num + 1 )\
                + "#" * 33 + "\n" )

        for kmer, offset, excess_R in kmer_offset_R_tuples_L:
            num_underscore_before = offset - lowest_offset
            num_underscore_after = total_num_aligned_chars - len( kmer ) -\
                    num_underscore_before
            out_f.write( "{0}{1}{2} {3:.3f}\n".format(
                "_"*num_underscore_before, kmer, "_"*num_underscore_after,
                    excess_R ) )

        out_f.write( "\n" )

    #### At the end of the file, make an "end motif" line (the word "motif"
    ####    signals there was a previous motif)
    out_f.write( "#"*32 + " < end motifs > " + "#"*32 )

    out_f.close()





def make_final_file_aligned_kmers_weights_per_motif_class_only_pruned_positions(
        unpruned_kmers_for_logo_F,
        pruned_pos_to_include_by_classnum_D ):
    """
    - Give an original, unpruned unpruned_kmers_for_logo_F, makes another
        kmers_for_logo_F that ONLY includes the pruned positions for each
        motif class_num, which are included in pruned_pos_to_include_by_classnum_D

    - Called by the make_logos_from_kmer_offset_R_T_Ls_by_classnum_D() function
        above

    6/2/16
    """
    pruned_kmers_for_logo_F = unpruned_kmers_for_logo_F.split(".txt")[0] +\
            ".pruned.txt"
    out_f = open( pruned_kmers_for_logo_F, "w" )

    motif_class = -1
    pruned_pos_to_start = ""
    pruned_pos_to_end = ""
    with open( unpruned_kmers_for_logo_F ) as f:
        for line in f:
            #### If it's the starting indicator line (it will contain "motif ")
            if ( line.find( "motif " ) != -1 ):
                motif_class += 1
                pruned_pos_to_include_L = pruned_pos_to_include_by_classnum_D[motif_class]
                pruned_pos_to_start = pruned_pos_to_include_L[0]
                pruned_pos_to_end = pruned_pos_to_start + len( pruned_pos_to_include_L )
                out_f.write( line )
            else:
                try:
                    kmer_w_flank_underscore, excess_R = line.strip().split(" ")
                    pruned_kmer_w_flank_underscore = kmer_w_flank_underscore[
                            pruned_pos_to_start:pruned_pos_to_end]
                    out_f.write( "{0} {1}\n".format(
                        pruned_kmer_w_flank_underscore, excess_R ) )

                except:
                    out_f.write( line )

    out_f.close()
    return pruned_kmers_for_logo_F





def add_kmer_to_existing_motif_class_or_start_new(
        kmer,
        excess_R,
        kmer_offset_R_T_Ls_by_classnum_D ):
    """
    - For a kmer & its excess_R, will try to:
        1. Add it to an existing motif_class (in kmer_offset_R_T_Ls_by_classnum_D)
        2. Start a new motif_class if that's not possible

    - Returns:        return_D =
        {"kmer_offset_R_T_Ls_by_classnum_D": kmer_offset_R_T_Ls_by_classnum_D,
         "motif_class_added_to": motif_class_added_to,
         "offset_used": offset_used}

    5/30/16
    """
    possible_comps_by_foundingk_alignk_D = {
            4: {4: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0"]},

            5: {5: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1"],
                4: ["side1_mismatch0", "side0_mismatch1", ]},

            6: {6: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1", "side0_mismatch2", "side2_mismatch1"],
                5: ["side1_mismatch0", "side0_mismatch1", "side1_mismatch1",
                    "side2_mismatch0"],
                4: ["side1_mismatch0", "side0_mismatch1"]},

            7: {7: ["side1_mismatch0", "side2_mismatch0", "side0_mismatch1",
                "side1_mismatch1", "side3_mismatch0", "side0_mismatch2",
                "side2_mismatch1", "side1_mismatch2"],
                6: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1", "side2_mismatch1", "side3_mismatch0",
                    "side0_mismatch2"],
                5: ["side1_mismatch0", "side0_mismatch1", "side1_mismatch1",
                    "side2_mismatch0"],
                4: ["side1_mismatch0", "side0_mismatch1"]},
            8: {8: ["side0_mismatch1",
                    "side0_mismatch2",
                    "side1_mismatch0",
                    "side1_mismatch1",
                    "side1_mismatch2",
                    "side2_mismatch0",
                    "side2_mismatch1",
                    "side3_mismatch0",
                    "side3_mismatch1",
                    "side4_mismatch0"],
                7: ["side1_mismatch0",
                    "side2_mismatch0",
                    "side0_mismatch1",
                    "side1_mismatch1",
                    "side3_mismatch0",
                    "side0_mismatch2",
                    "side2_mismatch1",
                    "side1_mismatch2"],
                6: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1", "side2_mismatch1", "side3_mismatch0",
                    "side0_mismatch2"],
                5: ["side1_mismatch0", "side0_mismatch1", "side1_mismatch1",
                    "side2_mismatch0"],
                4: ["side1_mismatch0", "side0_mismatch1"]} }

    score_D = {"match": 1,
        "mismatch": -1,
        "side": -0.7}

    name_score_T_Ls_by_alignk_D = {
        8: [],
        7: [],
        6: [],
        5: [],
        4: []}

    #### Go through and score all of the possibilities in
    ####    possible_comps_by_foundingk_alignk_D
    for founding_k, founding_k_D in possible_comps_by_foundingk_alignk_D.iteritems():
        for align_k, L in founding_k_D.iteritems():
            for alignment in L:
                side = int( alignment.split("side")[-1][0] )
                mismatch = int( alignment[-1] )
                match = align_k - side - mismatch

                score = match*score_D["match"] +\
                    mismatch * score_D["mismatch"] +\
                    side * score_D["side"]
                name = "{0}_{1}".format( founding_k, alignment )
                name_score_T_Ls_by_alignk_D[align_k].append(
                        (name, score) )

    #### Go through all of the align_ks, sort by descending score
    pref_order_Ls_by_alignk_D = {}
    for align_k, name_score_T_L in name_score_T_Ls_by_alignk_D.iteritems():
        name_score_T_L.sort( key = lambda x: -1* x[1] )
        pref_order_Ls_by_alignk_D[align_k] = [x[0] for x in name_score_T_L]

    existing_motif_class_nums_L = range( len( kmer_offset_R_T_Ls_by_classnum_D ) )
    #### the length of the kmer to be aligned
    align_k = len( kmer )

    #### Go through each of the existing motif classes &
    matches_to_existing_L = []
    pos_of_best_match_in_orderedL = 100
    best_existing_motif_class = None
    best_offset = None
    for existing_motif_class in existing_motif_class_nums_L:

        founding_kmer = kmer_offset_R_T_Ls_by_classnum_D[existing_motif_class][0][0]

        #### If the founding_kmer is a homopolymer, deal with it differently
        if ( is_homopolymer( founding_kmer ) == True ):

            #### Determine if this is a homopolymer off-by-1; if so, include
            ####    it with the minor_nt at the center position
            homo_off_by_1_returned_D = is_homopolymer_w_1_intervening(
                    kmer )
            if (homo_off_by_1_returned_D["is_homopolymer_w_1_intervening"] == False):
                continue

            #### Make sure that the homopolymer nt is the same between the
            ####    founding_kmer and kmer
            if (founding_kmer[0] != homo_off_by_1_returned_D["major_nt"]):
                continue

            #### If it's made it this far, kmer is a homopolymer off-by-1; if
            ####    NO add'l kmers have been added to this motif class, add
            ####    this kmer to have offset 0
            if (len( kmer_offset_R_T_Ls_by_classnum_D[existing_motif_class] ) == 1):

                pos_of_best_match_in_orderedL = -1
                best_existing_motif_class = existing_motif_class
                best_offset = 0

            #### Otherwise, if an off-by-1 homopolymer has already been added
            ####    to this motif class, match up the non-homopolymer nt
            ####    positions
            else:
                non_homo_already_added_kmer = kmer_offset_R_T_Ls_by_classnum_D[existing_motif_class][1][0]
                num_nt_before_non_homo_of_already_added = is_homopolymer_w_1_intervening(
                        non_homo_already_added_kmer )["pos_of_minor_nt"]
                num_nt_before_non_homo_this_kmer = is_homopolymer_w_1_intervening(
                        kmer )["pos_of_minor_nt"]
                offset_to_align_nonhomo = num_nt_before_non_homo_of_already_added -\
                        num_nt_before_non_homo_this_kmer

                pos_of_best_match_in_orderedL = -1
                best_existing_motif_class = existing_motif_class
                best_offset = offset_to_align_nonhomo


        #### If the founding_kmer is NOT a homopolymer, do the normal alignment
        else:
            #### best_match_returned_D has keys:
            ####    best_match_offset (is None if there's no match)
            ####    best_match_side (the number of unaligned positions; None if
            ####        theres no match)
            ####    best_match: like "side0_mismatch1"
            best_match_returned_D = get_best_match_of_kmer_to_foundingkmer(
                    kmer,
                    founding_kmer,
                    possible_comps_by_foundingk_alignk_D )

            if type( best_match_returned_D["best_match_offset"] ) is int:
                foundingk_match_name = "{0}_{1}".format(
                        len(founding_kmer), best_match_returned_D["best_match"] )
                try:
                    pos_of_this_match_in_orderedL = pref_order_Ls_by_alignk_D[
                            align_k].index( foundingk_match_name )
                    if (pos_of_this_match_in_orderedL < pos_of_best_match_in_orderedL):
                        pos_of_best_match_in_orderedL = pos_of_this_match_in_orderedL
                        best_existing_motif_class = existing_motif_class
                        best_offset = best_match_returned_D["best_match_offset"]
                except ValueError:
                    pass

    #### If the kmer aligned to any existing motif class
    if (type( best_existing_motif_class ) is int):
        kmer_offset_R_T_Ls_by_classnum_D[best_existing_motif_class].append(
                ( kmer, best_offset, excess_R ) )
        offset_used = best_offset
        motif_class_added_to = best_existing_motif_class
    #### Otherwise, start a new class
    else:
        motif_class_added_to = len( kmer_offset_R_T_Ls_by_classnum_D )
        kmer_offset_R_T_Ls_by_classnum_D[ motif_class_added_to ] =\
                [(kmer, 0, excess_R)]
        offset_used = 0

    return_D = {"kmer_offset_R_T_Ls_by_classnum_D":
            kmer_offset_R_T_Ls_by_classnum_D,
        "motif_class_added_to": motif_class_added_to,
            "offset_used": offset_used}
    return return_D





###############################################################################
######################### < CALLED BY RBNS_main.py > ##########################


def make_logos_from_PWM_D(
        PWM_D,
        output_DIR,
        identifier_name,
        title = "",
        T = False,
        first_index = 0 ):
    """
    - INPUT:
        - output_DIR, where a temporary PWM_file will be made and the output
          logos will go
        - identifier_name, a unique file name start that will be used to create
          the PWM_file and SeqLogo
        - title: optionally, the title above the SeqLogo
    """
    PWM_F = os.path.join(output_DIR, identifier_name + ".PWM")
    #### Make the PWM file and logos using the helper function below
    make_PWM_file_from_D( PWM_F, PWM_D )
    # make the SeqLogo using the just-created PWM file
    out_fstart = os.path.join( output_DIR, identifier_name )
    seq_logo_out_F, prob_logo_out_F = make_seqlogo_from_PWM_F(
            out_fstart,
            PWM_F,
            title = title,
            first_index = first_index )
    return seq_logo_out_F, prob_logo_out_F








def make_PWM_file_from_D(
        PWM_out_F,
        PWM_D ):
    """
    - A helper function called by make_logos_from_PWM_D() above
    """
    pos_L = PWM_D.keys()
    pos_L.sort()
    letters_L = PWM_D[pos_L[0]].keys()
    letters_L.sort()
    #### RNA logos must have "U" in the 2nd row header line
    letters_L_header = ["A","C","G","U"]
    with open(PWM_out_F, 'w') as f:
        f.write('ID test\n')
        f.write("\t".join(["PO"] + letters_L_header) + "\n")
        for position in pos_L:
            f.write(str(position))
            for letter in letters_L:
                f.write("\t" + str(PWM_D[position][letter]))
            f.write("\n")






def make_seqlogo_from_PWM_F(
        out_fstart,
        PWM_F,
        title = "",
        first_index = 0 ):
    """
    INPUT:
        - out_fstart, the starter path to the files that will be produced
        ( should be like ".../PTB2_set0" ("_seqlogo.png" or similar will be
        appended to seq_logo_out_fstart)
        - PWM_f, a position weight matrix file that was created by the
          make_PWM_from_dict() function above
    - See http://weblogo.threeplusone.com/manual.html#CLI for options
    """
    for scale_width in [True, False]:
        ############################# < PROB LOGO > ###########################
        if scale_width is False:
            prob_logo_out_F = out_fstart + "_probs_logo.png"
        else:
            prob_logo_out_F = out_fstart + "_probs_logo.cols_scaled_by_num_nt.png"
        cmd_probslogo = [
           'weblogo',
           '-U',
           "'probability'",
           "--composition",
           "'equiprobable'",
           "--format png",
           "--sequence-type rna",
           "--resolution 900",
           "--size large",
           "--errorbars NO",
           "--fineprint ''",
           "--first-index {}".format( starting_nt_idx ),
           '--title "{0}"'.format(title),
           "--color-scheme classic",
           "<",
           PWM_F,
           ">",
           prob_logo_out_F
            ]
        #### If the width SHOULDN'T be scaled
        if scale_width is False:
            cmd_probslogo = cmd_probslogo[:-4] + ["--scale-width", "NO"] + cmd_probslogo[-4:]
        cmd = ' '.join(cmd_probslogo)
        FNULL = open(os.devnull, 'w')
        p = subprocess.Popen( cmd, shell = True,
                stdout=FNULL, stderr=subprocess.STDOUT )
        p.wait()

        ############################# < SEQ LOGO > ############################
        if scale_width is False:
            seq_logo_out_F = out_fstart + "_seq_logo.png"
        else:
            seq_logo_out_F = out_fstart + "_seq_logo.cols_scaled_by_num_nt.png"
        cmd_seqlogo = [
           'weblogo',
           "-U 'bits'",
           "--composition 'equiprobable'",
           "--weight 0",
           "--format png",
           "-n 200",
           "--stack-width 40",
           "--sequence-type rna",
           "--resolution 900",
           "--size large",
           "--errorbars NO",
           "--fineprint ''",
           "--first-index {}".format( starting_nt_idx ),
           '--title "{0}"'.format(title),
           "--color-scheme classic",
           "<",
           PWM_F,
           ">",
           seq_logo_out_F
            ]
        #### If the width SHOULDN'T be scaled
        if scale_width is False:
            cmd_seqlogo = cmd_seqlogo[:-4] + ["--scale-width", "NO"] + cmd_seqlogo[-4:]
        cmd = ' '.join( cmd_seqlogo )
        p = subprocess.Popen( cmd, shell = True,
                stdout=FNULL, stderr=subprocess.STDOUT )
        p.wait()
        FNULL.close()
    return seq_logo_out_F, prob_logo_out_F


######################### </ CALLED BY RBNS_main.py > #########################
###############################################################################





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



def get_best_match_of_kmer_to_foundingkmer(
        kmer_to_align,
        founding_kmer,
        possible_comps_by_foundingk_alignk_D ):
    """
    - Will try to align the kmer_to_align to the founding_kmer (trying all
        possible sliding combinations), with the number of mismatches
        allowed specified by possible_comps_by_foundingk_alignk_D

    - Returns:
        return_D = {"best_match_offset": best_match_offset,
            "best_match_side": best_match_side,
            "best_match": best_match}
    5/30/16
    """
    best_match_offset = None
    best_match_side = None
    best_match = None
    best_match_position_in_allowedmatchesL = 100
    ####    allowed_matches_L is like:
    ####        ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0"],
    ####        where side is the number of unaligned (overhang) positions,
    ####        and mismatch is the # of mismatched among the aligned positions
    ####    - This list is ordered from best -> worst, so if multiple offsets
    ####        are in the list, we'll use the one that has the lowest
    ####        best_match_position_in_allowedmatchesL
    allowed_matches_L = possible_comps_by_foundingk_alignk_D\
            [len(founding_kmer)][len(kmer_to_align)]
    sides_allowed_L = [int(x.split("side")[-1][0]) for x in allowed_matches_L]

    for offset in range( -3, len( founding_kmer ) ):

        #### Get the # of nt hanging of the "side" of the founding kmer
        if (offset < 0):
            side = abs( offset )
            pos_to_align = len( kmer_to_align ) - side
            kmer_to_align_for_mismatches = kmer_to_align[(-1*pos_to_align):]
            founding_kmer_for_mismatches = founding_kmer[:pos_to_align]

        else:
            side = max( offset + len( kmer_to_align ) - len( founding_kmer ), 0 )

            if (side == 0):
                kmer_to_align_for_mismatches = kmer_to_align
            else:
                kmer_to_align_for_mismatches = kmer_to_align[:-1*side]

            founding_kmer_for_mismatches = founding_kmer[
                    offset:offset+len( kmer_to_align_for_mismatches )]

        if side not in sides_allowed_L:
            continue

        #### Get the # of mismatches between the kmer_to_align_for_mismatches
        ####    and founding_kmer_for_mismatches
        mismatches = RBNS_utils.hamming_distance(
                kmer_to_align_for_mismatches,
                founding_kmer_for_mismatches )
        #### See if this side/mismatch combination is allowed
        this_match = "side{0}_mismatch{1}".format( side, mismatches )

        if (this_match in allowed_matches_L):
            pos_in_allowedmatchesL = allowed_matches_L.index( this_match )
            #### If this offset is better than (i.e., has a lover index in
            ####    allowed_matches_L) the previous best one, record it
            if (pos_in_allowedmatchesL < best_match_position_in_allowedmatchesL):
                best_match_position_in_allowedmatchesL = pos_in_allowedmatchesL
                best_match_offset = offset
                best_match_side = side
                best_match = this_match

    return_D = {"best_match_offset": best_match_offset,
            "best_match_side": best_match_side,
            "best_match": best_match}

    return return_D











