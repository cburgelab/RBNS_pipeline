#!/usr/bin/env python
import pprint
import os, sys
import cPickle as pickle
import datetime
import math
import scipy
import glob

import pylab as pl
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
from matplotlib.ticker import FuncFormatter
from scipy.stats import gaussian_kde
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
from matplotlib.backends.backend_pdf import PdfPages

import RBNS_utils
#import RBNS_logo_helpers


#### To turn on LaTeX text; if you do not have an installation of LaTeX on
####    your path (e.g., .../texlive/2013/bin/x86_64-linux), can comment this out
try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
except:
    pass

light_to_dark_D = {1: [5],
                   2: [7, 3],
                   3: [8, 5, 2],
                   4: [9, 6, 4, 1],
                   5: [9, 7, 5, 3, 1],
                   6: [9, 7, 6, 4, 2, 0],
                   7: [9, 7, 6, 4, 3, 1, 0],
                   8: [9, 7, 6, 5, 4, 3, 2, 0],
                   9: [9, 8, 7, 6, 5, 4, 3, 2, 0],
                   10: [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                   11: [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
                   12: [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]}
dark_to_light_D = { 1: [5],
                    2: [3, 7],
                    3: [2, 5, 8],
                    4: [1, 3, 6, 9],
                    5: [1, 2, 5, 7, 9],
                    6: [0, 2, 4, 6, 7, 9],
                    7: [0, 1, 3, 4, 6, 7, 9],
                    8: [0, 2, 3, 4, 5, 6, 7, 9],
                    9: [0, 2, 3, 4, 5, 6, 7, 8, 9],
                    10: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                    11: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}



shades_of_orange_D = {
        ### dark
        0: "#4e3300",
        1: "#764c00",
        2: "#9d6600",
        3: "#c47f00",
        4: "#eb9800",
        5: "#ffac14",
        6: "#ffba3b",
        7: "#ffc862",
        8: "#ffd589",
        9: "#ffe3b1",
        10: "#fff8eb"}
        ### light

shades_of_red_D = {
        ### dark
        0: "#4e0000",
        1: "#760000",
        2: "#9d0000",
        3: "#c40000",
        4: "#eb0000",
        5: "#ff1414",
        6: "#ff3b3b",
        7: "#ff8989",
        8: "#ffb1b1",
        9: "#ffd8d8",
        10: "#ffecec"}
        ### light


shades_of_blue_D = {
        0: "#00004e",
        1: "#000089",
        2: "#0000b1",
        3: "#0000d8",
        4: "#0000ff",
        5: "#2727ff",
        6: "#6262ff",
        7: "#8989ff",
        8: "#c4c4ff",
        9: "#ebebff",
        10: "#ebebff",
        11: "w"}


nt_colors_D = {
        "A": "#008000",
        "C": "#0000FF",
        "G": "#FFA500",
        "T": "#FF0000",
        "A_input": "#94ff94",
        "C_input": "#9d9dff",
        "G_input": "#ffd589",
        "T_input": "#ff9d9d" }


rainbow_colors_L = ["#E05B51", "#EF6C00", "#FDD835", "#8BC34A", "#2196F3",
        "#d0a2ff", "#ff99ab", "#b5eef3", "#c47d33", "#960e29", "#808080"]


default_scatter_colors_L = ["#0066F7", "#FE8C0C", "#F71800", "#00C611",
        "#7400c6", "#ff4a6a", "#95461d"]


legend_locations_D = {"upper_right": 1,
        "upper_left": 2,
        "lower_left": 3,
        "lower_right": 4,
        "right": 5,
        "center_left": 6,
        "center_right": 7,
        "lower_center": 8,
        "upper_center": 9,
        "center": 10}

nt_colors_by_nt_D = {
        "A": "#008000",
        "C": "#0000FF",
        "G": "#FFA500",
        "T": "#FF0000",
        "U": "#FF0000",
        "C+G": "#805380",
        "A+T": "#804000",
        "A+U": "#804000"}





def plot_multiple_inOnePDF(
        figs_L,
        output_PDF_F,
        also_png = True,
        fig_type = "pdf" ):
    """
    INPUT:
        - figs_L: a list of fig from previously created matplotlib plots
        - output_PDF_F: an absolute path of the combined PDF to be created
    """
    print "Saving to {}".format( output_PDF_F )
    output_DIR = os.path.dirname( output_PDF_F )
    output_basename = os.path.basename( output_PDF_F ).rsplit(".", 1)[0]


    if (fig_type == "pdf"):
        pdf = PdfPages( output_PDF_F )
    for fig_num, fig in enumerate( figs_L ):

        if (fig_type == "pdf"):
            pdf.savefig( fig )

        if also_png:
            #### Also save this as a .png
            png_F = os.path.join( output_DIR,
                    output_basename + "_{}.png".format( fig_num ) )
            print "Saving to: {}".format( png_F )
            fig.savefig( png_F, format = 'png', dpi = 200, rasterized = False )

    if (fig_type == "pdf"):
        pdf.close()






def make_rectangular_heatmap_plot(
        checkerboard_D,
        annots_for_checkerboard_D = None,
        label_key_1 = "",
        label_key_2 = "",
        label_keys_1 = True,
        label_keys_2 = True,
        title = "",
        colorbar_label = "",
        order_of_keys_1_L = "ordinal",
        order_of_keys_2_L = "ordinal",
        override_keys_1_labels_L = None,
        override_keys_2_labels_L = None,
        max_value_for_colormap = None,
        min_value_for_colormap = None,
        colormap = "RedBlue",
        colormap_center = None,
        zero_color = None,
        horiz_lines_L = [],
        horiz_line_kw_args_D = {},
        vert_lines_L = [],
        vert_line_kw_args_D = {},
        lower_text_color = "#ffffff",
        upper_text_color = "k",
        out_F = None,
        nan_color = "#808080",
        special_colors_D = {},
        gridlines_color = False,
        height_per_row = "default",
        title_fontsize = 18,
        annot_fontsize = 14,
        label_1_fontsize = 18,
        label_2_fontsize = 18,
        label_fontsize = 18,
        label_1_rotation = 0,
        bottom_of_main_heatmap = 0.11,
        top_of_main_heatmap = 0.92 ):
    """
    - Makes a rectangular (generally square) plot with each sub-square shaded
        according to its value in checkerboard_D
    - INPUTs:
        - checkerboard_D:
            - a dictionary with format checkerboard_D[key_1][key_2]
              - key_1 will be plotted along the x-axis with label_key_1
              - key_2 will be plotted along the y-axis with label_key_2
        - colormap:
            - One of the colormaps below

        - label_fontsize: The font of the main x- and y-axis captions
        - label_1_fontsize/label_2_fontsize: the size of the x- and y-axis
            (label_1 and label_2) ticks (if label_keys_1/2 are True)
        - annot_fontsize: The size of the labels INSIDE each square of the
            heatmap

        - horiz_lines_L and vert_lines_L: should be between 0 & 1 (i.e.,
            fractional of how far across the axes to go)

    7/15/15
    """
    #pprint.pprint( checkerboard_D )

    from matplotlib.colors import rgb2hex

    #### GET THE COLORMAP
    if (colormap == "RedBlue"):
        cmap = pl.cm.seismic
    elif (colormap == "YellowGreenBlue"):
        cmap = pl.cm.YlGnBu
    elif (colormap == "BrownGreen"):
        cmap = pl.cm.BrBG
    elif (colormap == "Reds"):
        cmap = pl.cm.Reds
    elif (colormap == "Greens"):
        cmap = pl.cm.Greens
    elif (colormap == "Blues"):
        cmap = pl.cm.Blues
    elif (colormap == "Purples"):
        cmap = pl.cm.Purples
    elif (colormap == "Greys"):
        cmap = pl.cm.Greys
    else:
        print "COLORMAP {} NOT RECOGNIZED".format( colormap )

    #### go through and get the keys
    if (order_of_keys_1_L == "ordinal"):
        keys_1_L = checkerboard_D.keys()
        keys_1_L.sort()
    else:
        keys_1_L = order_of_keys_1_L
    #### a list of strings, which will be used as the x-ticks
    str_keys_1_L = [str(x).replace( "_", "\_" ) for x in keys_1_L]

    if (order_of_keys_2_L == "ordinal"):
        keys_2_L = checkerboard_D[keys_1_L[0]].keys()
        keys_2_L.sort()
    else:
        keys_2_L = order_of_keys_2_L
    #### a list of strings, which will be used as the y-ticks
    str_keys_2_L = [str(x).replace( "_", "\_" ) for x in keys_2_L]

    #### Go through and get the MAX and MIN values in checkerboard_D
    min_val = 1000000
    max_val = -1000000
    for key_1 in keys_1_L:
        for key_2 in keys_2_L:
            val = checkerboard_D[key_1][key_2]
            min_val = min( min_val, val )
            max_val = max( max_val, val )
    if (max_value_for_colormap != None):
        #max_val = min( max_val, max_value_for_colormap )
        max_val = max_value_for_colormap
    if (min_value_for_colormap != None):
        #min_val = max( min_val, min_value_for_colormap )
        min_val = min_value_for_colormap

    #### if a colormap_center was provided, use that as the center
    if (colormap_center != None):
        #### get the farthest from the center
        max_dist_from_center = max(abs( max_val - colormap_center ),
                abs( colormap_center - min_val ))
        colormap_bottom = colormap_center - max_dist_from_center
        colormap_top = colormap_center + max_dist_from_center
    #### if a colormap_center was NOT provided, simply use the max_ and min_vals
    else:
        colormap_bottom = min_val
        colormap_top = max_val
    norm = mpl.colors.Normalize( colormap_bottom, colormap_top )

    ####    Make a figure and a main axis for the heat map
    fig = plt.figure()
    #gs = GridSpec( 1, 1, bottom = 0.18, left = 0.18, right = 0.82)
    #ax = fig.add_subplot( gs[0,0] )
    #if (height_per_row == None):
    ax = fig.add_subplot(111)
    #fig.add_axes( [0.15, 0.1, 0.66, 0.8] )
    #ax = plt.gca()
    ax_colorbar_limits_L = [
            0.83,
            bottom_of_main_heatmap,
            0.02,
            top_of_main_heatmap - bottom_of_main_heatmap]
    if (height_per_row != "default"):
        top_of_fig = 0.94
        total_height = height_per_row * len( keys_2_L )
        current_limits_L = ax.get_position()
        axis_limits_L = [0.2, top_of_fig - total_height, 0.65, total_height]
        axis_limits_Bbox = matplotlib.transforms.Bbox(
                array([[axis_limits_L[0], axis_limits_L[1]],
                    [axis_limits_L[2], axis_limits_L[3]]]))
        ax.set_position( axis_limits_Bbox )
        #ax = fig.add_axes( axis_limits_L )
        ax_colorbar_limits_L = [axis_limits_L[0] + axis_limits_L[2] + 0.01,
                axis_limits_L[1], 0.02, axis_limits_L[3]]

    ##### go through each of the squares
    for key_1_num, key_1 in enumerate( keys_1_L ):
        for key_2_num, key_2 in enumerate( keys_2_L ):
            val = checkerboard_D[key_1][key_2]
            ##### get the color for this square
            #print val
            #print "bottom: {}".format( colormap_bottom )
            #print "top: {}".format( colormap_top )
            #### if val is "nan", change it to 0
            if (math.isnan(val) == True):
                color = nan_color
                zero_to_one_val = 0.5
            elif val in special_colors_D:
                color = special_colors_D[val]
                zero_to_one_val = 0.5
            else:
                zero_to_one_val = get_zero_to_one_val(
                        val,
                        colormap_bottom,
                        colormap_top )
                color = rgb2hex(cmap( zero_to_one_val ))

            if (zero_color != None) and ( float(val) == 0. ):
                color = zero_color

            #### if the zero_to_one_val is LESS than ~0.4, that means the
            ####    background square is dark blue, so the text on top should
            ####    be white
            if (zero_to_one_val <= 0.4):
                text_color = lower_text_color
            else:
                text_color = upper_text_color
            #### the the coordinates of the bottom left corner, on a scale of
            ####    (0, 1)
            left_edge = key_1_num / float(len( keys_1_L ))
            bottom_edge = key_2_num / float(len( keys_2_L ))
            center_x = (key_1_num + 0.5) / float(len( keys_1_L ))
            center_y = (key_2_num + 0.5) / float(len( keys_2_L ))
            plt.gca().add_patch(Rectangle(
                    (left_edge, bottom_edge),
                    1./len( keys_1_L ),
                    1./len( keys_2_L ) , color=color ))
            #### If there is an annotation to plot on top of this square
            if (annots_for_checkerboard_D != None):
                annot = annots_for_checkerboard_D[key_1][key_2]
                plt.annotate(
                    r'{0}'.format( annot.replace("_", "\_") ),
                    (center_x, center_y),
                    ha = 'center',
                    va = 'center',
                    fontsize = annot_fontsize,
                    color = text_color,
                    xycoords = 'axes fraction')

    for vert_line in vert_lines_L:
        ax.axvline( x = vert_line,
                **vert_line_kw_args_D )
    for horiz_line in horiz_lines_L:
        ax.axhline( y = horiz_line,
                **horiz_line_kw_args_D )

    axes().set_aspect('equal')

    #ax.set_title(r"{}".format( title.replace("_", "\_") ),
    ax.set_title(r"{}".format( title.replace("_", "\_") ),
            fontsize=title_fontsize)
    #ax.set_xlabel(r"{}".format( label_key_1.replace("_", "\_") ),
    ax.set_xlabel(r"{}".format( label_key_1 ),
            fontsize = label_fontsize)
    ax.set_ylabel(r"{}".format( label_key_2 ),
            fontsize = label_fontsize)

    #### SET THE x-ticks
    xticks_L = [(i+0.5)/len( keys_1_L ) for i in range(len(keys_1_L))]
    if (label_keys_1 == True):
        if (label_1_rotation == 0):
            ha_xticklabel = 'center'
        else:
            ha_xticklabel = 'right'
        censored_xticks_L = []
        censored_tick_lables_L = []
        for idx, xtick_label in enumerate( str_keys_1_L ):
            if ( len( xtick_label.strip() ) > 0 ):
                censored_xticks_L.append( xticks_L[idx] )
                censored_tick_lables_L.append( xtick_label )
        ax.set_xticks( censored_xticks_L )
        if override_keys_1_labels_L is not None:
            censored_tick_lables_L = override_keys_1_labels_L
        ax.set_xticklabels(
                censored_tick_lables_L,
                rotation = label_1_rotation,
                ha = ha_xticklabel )
        ax.tick_params(axis='x', labelsize = label_1_fontsize)
    else:
        ax.set_xticklabels( [""] * len(str_keys_1_L))
        ax.tick_params(axis='x', labelsize = 4 )
        ax.set_xticks( xticks_L )

    #### SET THE y-ticks
    yticks_L = [(i+0.5)/len( keys_2_L ) for i in range(len(keys_2_L))]
    ax.set_yticks( yticks_L )
    if (label_keys_2 == True):
        if override_keys_2_labels_L is not None:
            str_keys_2_L = override_keys_2_labels_L
        ax.set_yticklabels( str_keys_2_L ) #, size=16, rotation=45, ha='right')
        ax.tick_params(axis='y', labelsize = label_2_fontsize)
    else:
        ax.set_yticklabels( [""] * len(str_keys_2_L))
        ax.tick_params(axis='y', labelsize = 4 )

    ##### If the gridlines should be included
    if (gridlines_color != False):
        for i in range( 1, len(keys_1_L) ):
            plt.axvline( x = i / float(len( keys_1_L )),
                    color = gridlines_color )
        for i in range( 1, len(keys_2_L) ):
            plt.axhline( y = i / float(len( keys_2_L )),
                    color = gridlines_color )

    ## Plot the colorbar to the left of the main figure
    axcolor = fig.add_axes( ax_colorbar_limits_L )
    cb = mpl.colorbar.ColorbarBase(
        axcolor,
        cmap = cmap,
        norm = norm,
        orientation='vertical')
    #cb = plt.colorbar()
    #cb.set_label(r'{}'.format( colorbar_label.replace("_", "\_") ), size = 18)
    cb.set_label(r'{}'.format( colorbar_label ), size = 18)
    plt.gca().tick_params(axis='y', labelsize=18)

    #plt.gca().tight_layout()
    plt.gcf().subplots_adjust( bottom = bottom_of_main_heatmap )
    plt.gcf().subplots_adjust( top = top_of_main_heatmap )

    return_D = {"fig": fig}
    #### save the plot if an out_F was provided
    if (out_F != None):

        fig.savefig( out_F )
        print "Saved to: {}".format( out_F )
        png_out_F = out_F.replace( "pdf", "png" )
        fig.savefig( png_out_F )
        print "and: {}".format( png_out_F )
        return_D["out_F"] = out_F
        return_D["png_out_F"] = png_out_F
    return return_D






def make_rectangular_heatmap_plot_RBNS_freqs(
        sig_kmers_to_KL_Div_D,
        sig_kmers_to_log2_Obs_over_Exp_L_D,
        #### These are the same, but for randomly chosen "controls" from
        ####    the bottom half of the KL Div. distribution
        ctrl_kmers_to_KL_Div_D,
        ctrl_kmers_to_log2_Obs_over_Exp_L_D,
        kmers_type_annot = "$k$mers with KL Div. Z-score $\geq3$",
        control_kmers_type_annot = r"Example control $k$mers with bottom 50\% KL Div.",
        annot_above_KL_heatmap_bar = "KL Div.",
        order_of_kmers_L = "decreasing KL div",
        order_of_control_kmers_L = "decreasing KL div",
        title = "",
        colorbar_label = "",
        log2_colormap = "RedBlue",
        KLdiv_colormap = "BuGn",
        right_bar_color_by = "KL",
        out_F = None,
        title_fontsize = 14,
        max_log2_val_colormap = None):
    """
    - Makes a plot of kmers log2(enrichment) of frequency relative to uniform
        background across the positions of the reads

    - INPUTs:

        - sig_kmers_to_KL_Div_D and ctrl_kmers_to_KL_Div_D:
            dictionaries of the kmer to its KL divergence from uniformly
            distributed over the read

        - sig_kmers_to_log2_Obs_over_Exp_L_D and ctrl_kmers_to_log2_Obs_over_Exp_L_D:
            dictionaries of the kmer to a list of its log2(R) of Observed over
            Expected Uniform freq. (each element in the list represents a
            read start position)

        - kmers_type_annot:
            - the annotation that will go about the main heatmap plot,
                explaining e.g. whether the kmers shown are all those with
                Z-score >= 3, or whether they're the sig. enriched kmers

        - order_of_kmers_L:
            - if it == "decreasing KL div", will plot them by decreasing KL Div.;
                otherwise, can pass in a list of kmers and it will plot them in
                the order (e.g., to plot kmers in the same order in which
                they're significantly enriched)

            - right_bar_color_by: either "KL" or "R"

    - Originally made to plot enrichments of kmer frequencies at each position
        over those expected from a uniform distribution - called by
        analyze_freqs_by_position_one_barcodes_ordered_kmers_to_consider() in
        /net/utr/data/atf/pfreese/other_RBNS/motif_alignments/scripts/rbns_motifs.py

    """
    from matplotlib.colors import rgb2hex

    assert( right_bar_color_by in ["KL", "R"] )

    #### GET THE COLORMAPS
    if (log2_colormap == "RedBlue"):
        log2_cmap = pl.cm.seismic
    elif (log2_colormap == "YellowGreenBlue"):
        log2_cmap = pl.cm.YlGnBu
    elif (log2_colormap == "BrownGreen"):
        log2_cmap = pl.cm.BrBG
    elif (log2_colormap == "RdPu"):
        log2_cmap = pl.cm.RdPu
    elif (log2_colormap == "PuOr"):
        log2_cmap = pl.cm.PuOr
    else:
        print "COLORMAP {} NOT RECOGNIZED".format( log2_colormap )

    if (KLdiv_colormap == "RedBlue"):
        KLdiv_cmap = pl.cm.seismic
    elif (KLdiv_colormap == "YellowGreenBlue"):
        KLdiv_cmap = pl.cm.YlGnBu
    elif (KLdiv_colormap == "BrownGreen"):
        KLdiv_cmap = pl.cm.BrBG
    elif (KLdiv_colormap == "RdPu"):
        KLdiv_cmap = pl.cm.RdPu
    elif (KLdiv_colormap == "BuGn"):
        KLdiv_cmap = pl.cm.BuGn
    elif (KLdiv_colormap == "PuOr"):
        KLdiv_cmap = pl.cm.PuOr
    else:
        print "COLORMAP {} NOT RECOGNIZED".format( KLdiv_colormap )

    #### Go through and get the MAX and MIN values in checkerboard_D
    greatest_abs_log2 = 0
    for log2_L in sig_kmers_to_log2_Obs_over_Exp_L_D.values():
        num_kmers_each_read = len( log2_L )
        for log2 in log2_L:
            greatest_abs_log2 = max( greatest_abs_log2, abs(log2) )
    #### if ax_log2_val_colormap was not passed in, use the greatest_abs_log2
    ####    values
    if (max_log2_val_colormap == None):
        colormap_top = greatest_abs_log2
        colormap_bottom = -1*greatest_abs_log2
    else:
        colormap_top = max_log2_val_colormap
        colormap_bottom = -1 * max_log2_val_colormap

    norm_log2 = mpl.colors.Normalize(colormap_bottom, colormap_top)

    #### Get the order in which the sig. kmers should be plotted (from top ->
    ####    bottom, in decreasing KL divergence)
    sig_kmers_tuples_L = [(kmer, sig_kmers_to_KL_Div_D[kmer]) for kmer in\
            sig_kmers_to_KL_Div_D]
    sig_kmers_tuples_L.sort( key = lambda x: -1*x[1] )
    if (order_of_kmers_L == "decreasing KL div"):
        sig_kmers_L = [x[0] for x in sig_kmers_tuples_L]
    else:
        sig_kmers_L = order_of_kmers_L
    num_sig_kmers = len( sig_kmers_L )
    #### The same thing, but for the "control" kmers
    if ( order_of_control_kmers_L == "decreasing KL div" ):
        ctrl_kmers_tuples_L = [(kmer, ctrl_kmers_to_KL_Div_D[kmer]) for kmer in\
                ctrl_kmers_to_KL_Div_D]
        ctrl_kmers_tuples_L.sort( key = lambda x: -1*x[1] )
        ctrl_kmers_L = [x[0] for x in ctrl_kmers_tuples_L]
    else:
        assert( type( order_of_control_kmers_L ) is list )
        ctrl_kmers_L = order_of_control_kmers_L
    num_ctrl_kmers = len( ctrl_kmers_L )

    #### the amount of space between the "significant" kmers axis and the
    ####    "control" kmers axis
    y_space_between_ctrl_and_sig = 0.08

    ####    Make a figure for the
    fig = plt.figure()
    #### the boundaries of the "control" kmers log2 axis
    ctrl_log2_ax_bottom = 0.08
    #### the boundaries of the log2 axis for the sig. kmers
    log2_ax_top = 0.9


    #### The amount of space for the sig. + control kmers (0.85, which is the
    ####    entire 0.9 range minus 0.05 that will separate the sig. & control
    ####    kmers
    yspace_per_kmer = (log2_ax_top - ctrl_log2_ax_bottom - y_space_between_ctrl_and_sig) /\
            (num_sig_kmers + num_ctrl_kmers)
    log2_ax_bottom = log2_ax_top - (yspace_per_kmer * num_sig_kmers)
    log2_ax_left = 0.24
    log2_ax_right = 0.76

    #### The amount of x-axis space for each of the read positions
    width_per_pos = (log2_ax_right - log2_ax_left) / num_kmers_each_read

    #### The amount of space between the main heatmap & the KL Div bar, and the
    ####    KL Div bar and the heatmap scale
    width_between_scale_bars = 0.03

    #### Annotate the desired kmers
    num_kmers_to_writeout = 30
    top_kmer_y = 0.93
    bottom_kmer_y = 0.03
    delta_y_kmer = (top_kmer_y - bottom_kmer_y) / num_kmers_to_writeout
    axlog2 = fig.add_axes( [log2_ax_left, log2_ax_bottom,
        log2_ax_right - log2_ax_left, log2_ax_top - log2_ax_bottom] )

    #### Label the title
    plt.annotate(r"{}".format( title.replace("_", "\_") ),
        (0.5, 0.96),
        ha = 'center',
        va = 'center',
        fontsize = title_fontsize,
        xycoords = 'figure fraction')
    #### Label the "significant" and "control" kmers
    plt.annotate( control_kmers_type_annot,
        ((log2_ax_right + log2_ax_left)/2., ctrl_log2_ax_bottom+(yspace_per_kmer*num_ctrl_kmers)+0.005),
        ha = 'center',
        va = 'bottom',
        fontsize = 12,
        xycoords = 'figure fraction')
    plt.annotate(r"{}".format( kmers_type_annot ),
        ((log2_ax_right + log2_ax_left)/2., log2_ax_top+0.005),
        ha = 'center',
        va = 'bottom',
        fontsize = 12,
        xycoords = 'figure fraction')
    ###########################################################################
    ######### Max a heatmap of the log2(Obs / Expected) at each position for
    ########    the "significant" kmers
    for kmer_num, kmer in enumerate( sig_kmers_L ):
        if (kmer_num < num_kmers_to_writeout):
            y_pos_annot = top_kmer_y - (kmer_num*delta_y_kmer)
            #### This version goes from 0 -> 1
            y_pos_heatmap_box = 1. - (kmer_num+0.5)/num_sig_kmers
            #### This converts the 0 to 1 version to, e.g.,
            ####    [log2_ax_bottom, log2_ax_top]
            y_pos_fig = (log2_ax_top-log2_ax_bottom)*y_pos_heatmap_box +\
                    log2_ax_bottom
            if (kmer_num % 2) == 0:
                line_width = 0.2
            else:
                line_width = 0.0
            plt.annotate(
                kmer,
                #xy = (log2_ax_left - 2*width_per_pos, y_pos_fig),
                xy = (log2_ax_left , y_pos_fig),
                xytext = ((log2_ax_left - (2*width_per_pos))/2., y_pos_annot),
                ha = 'center',
                va = 'center',
                fontsize = 12,
                textcoords='figure fraction',
                xycoords = 'figure fraction',
                arrowprops = dict(
                    facecolor='black',
                    shrink=0.,
                    frac = 0.,
                    headwidth = 0.,
                    width = 0.0,
                    lw = line_width))
                    #relpos = (1, 0.5)))

        bottom_edge = 1. - (kmer_num+1.)/num_sig_kmers
        #y_bottom_this_kmer = log2_ax_top - (yspace_per_kmer*(kmer_num + 1))
        #### Go through each of the positions
        for read_pos in range( num_kmers_each_read ):
            left_this_kmer = log2_ax_left + (read_pos * width_per_pos)
            left_edge = read_pos / float( num_kmers_each_read )

            val = sig_kmers_to_log2_Obs_over_Exp_L_D[kmer][read_pos]

            zero_to_one_val = get_zero_to_one_val(
                    val,
                    colormap_bottom,
                    colormap_top )
            color = rgb2hex(log2_cmap( zero_to_one_val ))
            #### the the coordinates of the bottom left corner, on a scale of
            ####    (0, 1)
            gca().add_patch(Rectangle(
                    (left_edge, bottom_edge),
                    1./num_kmers_each_read,
                    1./num_sig_kmers,
                    color=color ))

    #### SET THE x-ticks
    if (num_kmers_each_read >= 20):
        xticks_L = [(i+0.5)/num_kmers_each_read for i in range(num_kmers_each_read)\
                if ((i % 2 ) == 1)]
        str_xaxis_keys_L = [str(i) for i in range( 1, num_kmers_each_read+1 )\
                if ((i % 2 ) == 0)]
    else:
        xticks_L = [(i+0.5)/num_kmers_each_read for i in range(num_kmers_each_read)]
        str_xaxis_keys_L = [str(x) for x in range( 1, num_kmers_each_read+1 )]
    axlog2.set_xticks( xticks_L )
    axlog2.set_xticklabels( str_xaxis_keys_L )

    ##### ERASE THE y-ticks
    axlog2.set_yticks( [] )
    axlog2.tick_params(labelleft = True)

    axlog2.tick_params(axis='x', labelsize=12)
    ###########################################################################


    ###########################################################################
    #### Make a heatmap showing the KL Divergence for the sig. kmers
    #### Annotate that this is the KL Divergence (or other)
    plt.annotate( annot_above_KL_heatmap_bar,
        (log2_ax_right + 1.5*width_between_scale_bars, log2_ax_top + 0.01),
        ha = 'center',
        va = 'bottom',
        fontsize = 12,
        xycoords = 'figure fraction')

    KL_colormap_top = 0.
    for kmer in sig_kmers_L:
        KL_colormap_top = max( KL_colormap_top, sig_kmers_to_KL_Div_D[kmer] )
    for value in ctrl_kmers_to_KL_Div_D.values():
        KL_colormap_top = max( KL_colormap_top, value )
    KL_ax_left = log2_ax_right + width_between_scale_bars
    axKL = fig.add_axes( [KL_ax_left, log2_ax_bottom,
        width_between_scale_bars, log2_ax_top - log2_ax_bottom] )


    ##### go through each of the squares
    for kmer_num, kmer in enumerate( sig_kmers_L ):
        bottom_edge = 1. - (kmer_num+1.)/num_sig_kmers

        val = sig_kmers_to_KL_Div_D[kmer]

        if ( right_bar_color_by == "KL" ):
            zero_to_one_val = get_zero_to_one_val(
                    val,
                    0,
                    KL_colormap_top )
        if ( right_bar_color_by == "R" ):
            zero_to_one_val = get_zero_to_one_val(
                    val,
                    colormap_bottom,
                    colormap_top )
        color = rgb2hex(KLdiv_cmap( zero_to_one_val ))
        #### the the coordinates of the bottom left corner, on a scale of
        ####    (0, 1)
        gca().add_patch(Rectangle(
                (0., bottom_edge),
                1.,
                1./num_sig_kmers,
                color=color ))

    #### ERASE THE x-ticks and y-ticks
    axKL.set_xticks( [] )
    axKL.set_xticklabels([])
    axKL.set_yticks( [] )
    axKL.set_yticklabels([])
    ###########################################################################



    ###########################################################################
    ######### Make a heatmap of the log2(Obs / Expected) at each position for
    ########    the "control" kmers
    axlog2_ctrl = fig.add_axes( [log2_ax_left, ctrl_log2_ax_bottom,
        log2_ax_right - log2_ax_left, yspace_per_kmer*num_ctrl_kmers] )
    for kmer_num, kmer in enumerate( ctrl_kmers_L ):
        bottom_edge = 1. - (kmer_num+1.)/num_ctrl_kmers
        #y_bottom_this_kmer = log2_ax_top - (yspace_per_kmer*(kmer_num + 1))
        #### Go through each of the positions
        for read_pos in range( num_kmers_each_read ):
            left_this_kmer = log2_ax_left + (read_pos * width_per_pos)
            left_edge = read_pos / float( num_kmers_each_read )

            val = ctrl_kmers_to_log2_Obs_over_Exp_L_D[kmer][read_pos]

            zero_to_one_val = get_zero_to_one_val(
                    val,
                    colormap_bottom,
                    colormap_top )
            color = rgb2hex(log2_cmap( zero_to_one_val ))
            #### the the coordinates of the bottom left corner, on a scale of
            ####    (0, 1)
            gca().add_patch(Rectangle(
                    (left_edge, bottom_edge),
                    1./num_kmers_each_read,
                    1./num_ctrl_kmers,
                    color=color ))

    #### SET THE x-ticks
    axlog2_ctrl.set_xticks( xticks_L )
    axlog2_ctrl.set_xticklabels( str_xaxis_keys_L ) #, size=16, rotation=45, ha='right')

    ##### ERASE THE y-ticks
    axlog2_ctrl.set_yticks( [] )
    axlog2_ctrl.tick_params(labelleft = True)

    axlog2_ctrl.tick_params(axis='x', labelsize=12)
    ###########################################################################


    ###########################################################################
    #### Make a heatmap showing the KL Divergence for the "control" kmers
    axKL_ctrl = fig.add_axes( [KL_ax_left, ctrl_log2_ax_bottom,
        width_between_scale_bars, yspace_per_kmer*num_ctrl_kmers] )

    ##### go through each of the squares
    for kmer_num, kmer in enumerate( ctrl_kmers_L ):
        bottom_edge = 1. - (kmer_num+1.)/num_ctrl_kmers

        val = ctrl_kmers_to_KL_Div_D[kmer]

        if ( right_bar_color_by == "KL" ):
            zero_to_one_val = get_zero_to_one_val(
                    val,
                    0,
                    KL_colormap_top )
        if ( right_bar_color_by == "R" ):
            zero_to_one_val = get_zero_to_one_val(
                    val,
                    colormap_bottom,
                    colormap_top )
        color = rgb2hex(KLdiv_cmap( zero_to_one_val ))
        #### the the coordinates of the bottom left corner, on a scale of
        ####    (0, 1)
        gca().add_patch(Rectangle(
                (0., bottom_edge),
                1.,
                1./num_ctrl_kmers,
                color=color ))

    #### ERASE THE x-ticks and y-ticks
    axKL_ctrl.set_xticks( [] )
    axKL_ctrl.set_xticklabels([])
    axKL_ctrl.set_yticks( [] )
    axKL_ctrl.set_yticklabels([])
    ###########################################################################

    ###########################################################################
    ## Plot the colorbar to the right of the main figure
    axcolor = fig.add_axes( [log2_ax_right + 3*width_between_scale_bars,
        ctrl_log2_ax_bottom, 0.02, log2_ax_top - ctrl_log2_ax_bottom] )
    cb = mpl.colorbar.ColorbarBase(
        axcolor,
        cmap = log2_cmap,
        norm=norm_log2,
        orientation='vertical')
    cb.set_label(r'{}'.format( colorbar_label ), size = 16)
    gca().tick_params(axis='y', labelsize=16)
    ###########################################################################

    return_D = {"fig": fig}
    #### save the plot if an out_F was provided
    if (out_F != None):
        fig.savefig( out_F )
        return_D["out_F"] = out_F
    return return_D






def scatter_nt_freqs_by_read_position(
        freqs_by_pos_D,
        conc_annot,
        protein_name_for_plotting = "",
        y_axis_min = "from_data",
        y_axis_max = "from_data",
        RNA = True ):
    """
    - Makes a scatter of the nt frequencies by position within the read, and
        returns the figure
    """

    num_nt = len( freqs_by_pos_D.keys() )

    fig = plt.figure()
    ax = fig.gca()
    #### plot a dotted horizontal line at 25%
    ax.axhline( y = 25.,
            color='k', linestyle='--', linewidth = 2)

    xs_L = [x + 1 for x in range( num_nt )]
    A_percentages_L = [freqs_by_pos_D[pos]["A"] * 100. for pos in range(num_nt)]
    C_percentages_L = [freqs_by_pos_D[pos]["C"] * 100. for pos in range(num_nt)]
    G_percentages_L = [freqs_by_pos_D[pos]["G"] * 100. for pos in range(num_nt)]
    T_percentages_L = [freqs_by_pos_D[pos]["T"] * 100. for pos in range(num_nt)]

    min_val = min(
        A_percentages_L + C_percentages_L + G_percentages_L + T_percentages_L )
    max_val = max(
        A_percentages_L + C_percentages_L + G_percentages_L + T_percentages_L )
    y_range = max_val - min_val

    h1 = plt.scatter( xs_L, A_percentages_L, c = nt_colors_D["A"],
            edgecolor = nt_colors_D["A"] )
    h2 = plt.scatter( xs_L, C_percentages_L, c = nt_colors_D["C"],
            edgecolor = nt_colors_D["C"] )
    h3 = plt.scatter( xs_L, G_percentages_L, c = nt_colors_D["G"],
            edgecolor = nt_colors_D["G"] )
    h4 = plt.scatter( xs_L, T_percentages_L, c = nt_colors_D["T"],
            edgecolor = nt_colors_D["T"] )

    x_max = 1.25 * num_nt
    ax.set_xlim(0, x_max)
    if ( y_axis_min == "from_data" ) and ( y_axis_max == "from_data" ):
        ax.set_ylim( min_val - (0.25 * y_range), max_val + (0.25 * y_range))
    else:
        ax.set_ylim( y_axis_min, y_axis_max )

    if (RNA == True):
        ax.legend( [h1,h2,h3,h4], ["A", "C", "G", "U"], loc = 5,\
            fancybox = True, fontsize = 18, borderaxespad = 0)
    else:
        ax.legend( [h1,h2,h3,h4], ["A", "C", "G", "T"], loc = 5,\
            fancybox = True, fontsize = 18, borderaxespad = 0)

    ax.set_xlabel(r"Read Position", fontsize = 18)
    #ax.set_ylabel( r"nt %", fontsize = 18 )
    ax.set_title( "{0} Nucleotide frequencies, {1}".format(
        protein_name_for_plotting.replace("_", "\_"), conc_annot ), fontsize = 18 )

    ax.tick_params( axis = 'x', labelsize = 18 )
    ax.tick_params( axis = 'y', labelsize = 18 )

    #### Add % signs to the y-axis labels
    formatter = FuncFormatter( axis_labels_to_percent )
    plt.gca().yaxis.set_major_formatter( formatter )
    plt.show()

    return fig



def make_rectangular_heatmap_NT_freq_across_read_all_libs(
        freqs_by_pos_D_by_conc_annot_D,
        conc_annots_L,
        title = "",
        colorbar_label = r"$log_2$(Observed / Uniform freq.)",
        log2_colormap = "RedBlue",
        KLdiv_colormap = "BuGn",
        out_F = None,
        title_fontsize = 16,
        conc_annot_fontsize = 8,
        max_log2_val_colormap = None ):
    """
    - Makes a plot of kmers log2(enrichment) of frequency relative to uniform
        background across the positions of the reads for ALL submitted RBNS
        experiments

        - called by make_meta_plots_of_kmer_freqs_by_position_in_submitted_exps()
            in ~/python_lib/RBNS_helpers/RBNS_kmers_by_position.py

    - INPUTs:

        - RBP_to_PD_KL_Div_D and ctrl_kmers_to_KL_Div_D:
            dictionaries of the RBNS exp. to its KL divergence from uniformly
            distributed over the read

        - RBP_to_PD_log2_Obs_over_Exp_L_D and ctrl_kmers_to_log2_Obs_over_Exp_L_D:
            dictionaries of the kmer to a list of its log2(R) of Observed over
            Expected Uniform freq. (each element in the list represents a
            read start position)


        - RBP_to_yaxis_annot_D: what should be plotted for each experiment
            along the y-axis (e.g., "RBFOX2: UGCAUG")

    - Originally made to plot enrichments of kmer frequencies at each position
        over those expected from a uniform distribution - called by
        analyze_freqs_by_position_one_barcodes_ordered_kmers_to_consider() in
        /net/utr/data/atf/pfreese/other_RBNS/motif_alignments/scripts/rbns_motifs.py

    """
    from matplotlib.colors import rgb2hex

    #### GET THE COLORMAPS
    if (log2_colormap == "RedBlue"):
        log2_cmap = pl.cm.seismic
    elif (log2_colormap == "YellowGreenBlue"):
        log2_cmap = pl.cm.YlGnBu
    elif (log2_colormap == "BrownGreen"):
        log2_cmap = pl.cm.BrBG
    elif (log2_colormap == "RdPu"):
        log2_cmap = pl.cm.RdPu
    else:
        print "COLORMAP {} NOT RECOGNIZED".format( log2_colormap )

    if (KLdiv_colormap == "RedBlue"):
        KLdiv_cmap = pl.cm.seismic
    elif (KLdiv_colormap == "YellowGreenBlue"):
        KLdiv_cmap = pl.cm.YlGnBu
    elif (KLdiv_colormap == "BrownGreen"):
        KLdiv_cmap = pl.cm.BrBG
    elif (KLdiv_colormap == "RdPu"):
        KLdiv_cmap = pl.cm.RdPu
    elif (KLdiv_colormap == "BuGn"):
        KLdiv_cmap = pl.cm.BuGn
    else:
        print "COLORMAP {} NOT RECOGNIZED".format( KLdiv_colormap )

    random_len = len( freqs_by_pos_D_by_conc_annot_D['Input'] )

    max_abs_val = 0.
    unnorm_vals_L_by_nt_concannot_D = {}
    for nt in ["T", "G", "C", "A"]:
        unnorm_vals_L_by_nt_concannot_D[nt] = {}
        for conc_annot in conc_annots_L:
            freqs_L = [freqs_by_pos_D_by_conc_annot_D[conc_annot][pos][nt] for\
                    pos in range( random_len )]
            mean_freq_at_each_pos = sum( freqs_L ) / random_len
            log_2_ratio_L = [math.log( freqs_L[x] / mean_freq_at_each_pos, 2. )\
                    for x in range( random_len )]
            max_abs_val = max( max_abs_val, max( [abs(x) for x in log_2_ratio_L] ) )
            unnorm_vals_L_by_nt_concannot_D[nt][conc_annot] = log_2_ratio_L
    #colormap_bottom = -1 * max_abs_val
    #colormap_top = max_abs_val
    colormap_bottom = -0.5
    colormap_top = 0.5

    norm_log2 = mpl.colors.Normalize(colormap_bottom, colormap_top)


    fig = plt.figure()

    plt.annotate(r"{}".format( title ),
        (0.5, 0.95),
        ha = 'center',
        va = 'bottom',
        fontsize = title_fontsize,
        xycoords = 'figure fraction')

    gca().set_xticks( [] )
    gca().set_xticklabels([])
    gca().set_yticks( [] )
    gca().set_yticklabels([])
    gca().spines['top'].set_color('none')
    gca().spines['bottom'].set_color('none')
    gca().spines['left'].set_color('none')
    gca().spines['right'].set_color('none')
    gca().tick_params(labelcolor="w", top='off',left="off",right="off",\
                                            bottom="off")

    ax_left = 0.1
    ax_right = 0.85
    max_ax_top = 0.
    orig_bottom = 0.04
    #### The amount of x-axis space for each of the read positions
    width_per_pos = ( ax_right - ax_left) / random_len
    for nt_idx, nt in enumerate( ["T", "G", "C", "A"] ):

        ax_bottom = ( 0.23 * nt_idx ) + orig_bottom
        ax_top = ax_bottom + 0.2

        ### The middle of this axis
        plt.annotate( nt.replace( "T", "U" ),
            (ax_left - 0.08, (ax_bottom + ax_top) / 2.),
            ha = 'center',
            va = 'center',
            fontsize = 12,
            xycoords = 'figure fraction')

        ax_this_nt = fig.add_axes( [ax_left, ax_bottom,
            ax_right - ax_left, ax_top - ax_bottom] )

        ###########################################################################
        ######### Max a heatmap of the log2(Obs / Expected) at each position for
        ########    the "significant" kmers
        for lib_idx, conc_annot in enumerate( conc_annots_L ):

            bottom_edge = 1. - ( lib_idx + 1.) / len( conc_annots_L )
            #### Go through each of the positions
            for read_pos in range( random_len ):

                left_this_pos = ax_left + (read_pos * width_per_pos)
                left_edge = read_pos / float( random_len )

                val = unnorm_vals_L_by_nt_concannot_D[nt][conc_annot][read_pos]

                zero_to_one_val = get_zero_to_one_val(
                        val,
                        colormap_bottom,
                        colormap_top )
                color = rgb2hex(log2_cmap( zero_to_one_val ))
                #### the the coordinates of the bottom left corner, on a scale of
                ####    (0, 1)
                gca().add_patch(Rectangle(
                        (left_edge, bottom_edge),
                        1. / random_len,
                        1. / len( conc_annots_L ),
                        color = color ))

        #### SET THE x-ticks
        if ( random_len  > 20 ):
            xticks_L = [(i+0.5)/random_len for i in range(random_len)\
                    if ((i % 2 ) == 1)]
            str_xaxis_keys_L = [str(i) for i in range( 1, random_len + 1 )\
                    if ((i % 2 ) == 0)]
        else:
            xticks_L = [(i+0.5)/random_len for i in range( random_len )]
            str_xaxis_keys_L = [str(x) for x in range( 1, random_len + 1 )]
        ax_this_nt.set_xticks( xticks_L )
        ax_this_nt.set_xticklabels( str_xaxis_keys_L )

        ##### Make the x-ticks from top to bottom
        yticks_L = [(i+0.5)/len( conc_annots_L ) for i in range(len(conc_annots_L))][::-1]
        ax_this_nt.set_yticks( yticks_L )
        ax_this_nt.tick_params( labelleft = True )
        ax_this_nt.set_yticklabels( conc_annots_L )
        if ( nt_idx == 0 ):
            ax_this_nt.set_xlabel( "Read position" )

        ax_this_nt.tick_params(axis='x', labelsize = 10 )
        ax_this_nt.tick_params(axis='y', labelsize = conc_annot_fontsize )

    ###########################################################################
    ## Plot the colorbar to the right of the main figure
    axcolor = fig.add_axes( [ ax_right + 0.02,
        orig_bottom, 0.02, ax_top - orig_bottom] )
    cb = mpl.colorbar.ColorbarBase(
        axcolor,
        cmap = log2_cmap,
        norm = norm_log2,
        orientation='vertical')
    cb.set_label(r'{}'.format( colorbar_label ), size = 16)
    gca().tick_params(axis='y', labelsize=16)
    ###########################################################################

    return_D = {'fig': fig}

    return return_D




def plot_stacked_bargraph_with_logos(
        nums_L,
        SeqLogos_L,
        out_F,
        title = "",
        xlab = "",
        ylab = "",
        colors_for_bargraph = "default",
        collapse_less_than_x_perc_into_Other = 5. ):
    """
    - Plots a basic stacked histogram with the SeqLogos_L, each having
        nums_L weight, to out_F

    - collapse_less_than_x_perc_into_Other:
        - if None: won't collapse any categories
        - Otherwise, a float between 0 and 100 (e.g., 5.). Individual
            categories less than this will be collapsed into an "other"
            category with no SeqLogo
    """
    import PIL
    from PIL import Image
    from matplotlib.backends.backend_pdf import PdfPages

    nums_L = [float(x) for x in nums_L]
    #### Normalize nums_L so they all sum to 100.
    nums_L = [x * 100./sum( nums_L ) for x in nums_L]

    #### Combine any nums that are individually less than 5% into 1 composite
    ####    "other" category
    if (type(collapse_less_than_x_perc_into_Other) is float) or\
            (type(collapse_less_than_x_perc_into_Other) is int):
        total_other = 0.
        pared_nums_L = []
        pared_SeqLogos_L = []
        for num_idx, num in enumerate( nums_L ):
            if (num < collapse_less_than_x_perc_into_Other):
                total_other += num
            else:
                pared_nums_L.append( num )
                try:
                    pared_SeqLogos_L.append( SeqLogos_L[num_idx] )
                except IndexError:
                    pass
        #### If there was any "total_other" added, include it
        if (total_other > 0.):
            pared_nums_L.append( total_other )
        nums_L = pared_nums_L
        SeqLogos_L = pared_SeqLogos_L


    #### Get the colors to use based on the number of logos
    if (colors_for_bargraph == "default"):
        colors_to_use_L = ["#a6cee3", "#1f78b4", "#b2df8a","#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00"]

    #### If a list of colors to use were explicitly passed in
    elif (type( colors_for_bargraph ) is list):
        colors_to_use_L = colors_for_bargraph

    else:
        print "ERROR: colors_for_bargraph {} NOT RECOGNIZED".format(
                colors_for_bargraph )

    pdf = PdfPages(out_F)

    fig = plt.figure(figsize=(4,3), dpi=300)

    strt_ax = plt.gca()
    strt_ax.spines['top'].set_color('none')
    strt_ax.spines['bottom'].set_color('none')
    strt_ax.spines['left'].set_color('none')
    strt_ax.spines['right'].set_color('none')
    strt_ax.tick_params(top='off',left="off",right="off",\
            bottom="off")
    for which_ax in ["x", "y"]:
        plt.tick_params(
                axis = which_ax,
                which = 'both',
                left = 'off',
                right = 'off',
                bottom = 'off',
                top = 'off',
                labelbottom = 'off',
                labeltop = 'off',
                labelleft = 'off',
                labelright = 'off' )

    bg_ax = fig.add_axes([0, 0, 1, 1])
    bg_ax.set_xlim(0, 1.)
    bg_ax.set_ylim(0, 1.)
    bg_ax.spines['top'].set_color('none')
    bg_ax.spines['bottom'].set_color('none')
    bg_ax.spines['left'].set_color('none')
    bg_ax.spines['right'].set_color('none')
    bg_ax.tick_params(top='off',left="off",right="off",\
            bottom="off")
    for which_ax in ["x", "y"]:
        plt.tick_params(
                axis = which_ax,
                which = 'both',
                left = 'off',
                right = 'off',
                bottom = 'off',
                top = 'off',
                labelbottom = 'off',
                labeltop = 'off',
                labelleft = 'off',
                labelright = 'off' )
    if (title != ""):
        plt.annotate( title.replace("_", "\_")[:20], (0.5, 0.95), fontsize = 14,
                xycoords='figure fraction',
                verticalalignment = 'center',
                horizontalalignment = 'center')


    ##################### < the axis for the bar plot > #######################
    lft = 0.45
    btm = 0.05
    wdth = 0.1
    hgh = 0.85
    ax = fig.add_axes([lft, btm, wdth, hgh])
    ax.set_xlim(0, 1.)
    ax.set_ylim(0, 100)
    ax.set_xlabel(xlab, fontsize=13)
    ax.set_ylabel(ylab, fontsize=13)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')

    ##################### < / the axis for the bar plot > #####################


    num_seqLogos = len( SeqLogos_L )
    if (num_seqLogos <= 4):
        height_of_seqlogo = 0.43
        width_of_seqlogo = .39
        num_logos_right_of_bar = 2
        seqLogos_right_L = SeqLogos_L[:2]
        seqLogos_left_L = SeqLogos_L[2:]
    elif (num_seqLogos <= 6):
        height_of_seqlogo = 0.33333
        width_of_seqlogo = .2777
        num_logos_right_of_bar = 3
        seqLogos_right_L = SeqLogos_L[:3]
        seqLogos_left_L = SeqLogos_L[3:]
    else:
        height_of_seqlogo = 0.25
        width_of_seqlogo = .2083
        num_logos_right_of_bar = 4
        seqLogos_right_L = SeqLogos_L[:4]
        seqLogos_left_L = SeqLogos_L[4:8]

    # if the length of SeqLogos is less than the length of nums_L, the last
    # number is the "other signficant enriched" category
    if (len(SeqLogos_L) < len(nums_L)):
        other_significant = True
    else:
        other_significant = False

    handles_L = []
    cur_bottom = 100.
    bottoms_L = []
    for num, val in enumerate(nums_L):
        cur_bottom -= val
        bottoms_L.append( cur_bottom )
        color_idx = num % len( colors_to_use_L )
        hand = ax.bar(0.,
                val,
                width = 1.,
                color = colors_to_use_L[color_idx],
                bottom = cur_bottom)
        handles_L.append(hand)
    for num, cur_bottom in enumerate(bottoms_L):
        try:
            y_to_connect_line_to = (cur_bottom/100.) + (nums_L[num]/200.)
            y_to_connect_line_to = (y_to_connect_line_to*hgh) + btm
            if (num < num_logos_right_of_bar):
                x_left = .6
                y_bot = 1.-(height_of_seqlogo*(num+1))
                y_bot = (hgh*y_bot) + btm
                xs = [x_left, 0.55]
            elif ( num < 2*num_logos_right_of_bar ):
                x_left = 0.38 - width_of_seqlogo
                y_bot = 1.-(height_of_seqlogo*(num - num_logos_right_of_bar + 1))
                y_bot = (hgh*y_bot) + btm
                xs = [x_left + width_of_seqlogo, 0.45]
            # the ys for the for the line is the same for either left
            # or right SeqLogos
            ys = [y_bot+(height_of_seqlogo*hgh/2), y_to_connect_line_to]
            SeqLogo_F = SeqLogos_L[num]
            ax2 = fig.add_axes(
                    [x_left, y_bot, width_of_seqlogo, height_of_seqlogo*hgh],
                    zorder = 1 )
            ax2.spines['top'].set_color('none')
            ax2.spines['bottom'].set_color('none')
            ax2.spines['left'].set_color('none')
            ax2.spines['right'].set_color('none')
            ax2.tick_params(which = 'both', top='off',left="off",right="off",\
                    bottom="off", labelbottom = "off", labeltop = "off",
                    labelleft = "off", labelright= 'off')
            basewidth = 1000
            img = Image.open(SeqLogo_F)
            wpercent = (basewidth / float(img.size[0]))
            hsize = int((float(img.size[1]) * float(wpercent)))
            img2 = img.resize((basewidth, hsize), PIL.Image.ANTIALIAS)
            plt.imshow(img2, aspect = "equal", interpolation='none',
                    zorder = 1)
            bg_ax.plot( xs, ys, '-', linewidth = .6, color = "#808080",
                    zorder = 3)
        except IndexError: pass

    if (other_significant == True):
        x_left = .6
        y_to_draw_line_to = (nums_L[-1]/200.)
        y_to_draw_line_to = (hgh*y_to_draw_line_to) + btm
        plt.annotate("Other sig. enriched", (x_left+0.01, btm+0.02),
                fontsize=10, xycoords='figure fraction',
                verticalalignment = 'center', horizontalalignment = 'left')
        xs = [x_left, 0.55]
        ys = [btm+0.02, y_to_draw_line_to]
        bg_ax.plot( xs, ys, '-', linewidth = .6, color = "#808080", zorder = 3)

    for which_ax in ["x", "y"]:
        if (which_ax == "x"):
            zorder = 1
        else:
            zorder = 3
        ax.tick_params(
                axis = which_ax,
                which = 'both',
                bottom = 'off',
                top = 'off',
                labelbottom = 'off',
                zorder = zorder )

    pdf.savefig(fig, dpi = 900)
    pdf.close()

    # also save as an SVG
    out_svg_F = out_F.rsplit(".", 1)[0] + ".svg"
    fig.savefig( out_svg_F, dpi = 900, format = "svg" )

    # also save as a PNG
    out_png_F = out_F.rsplit(".", 1)[0] + ".png"
    fig.savefig( out_png_F, dpi = 900, format = "png" )

    return_D = {"fig": fig,
            "out_F": out_F,
            "out_svg_F": out_svg_F,
            "out_png_F": out_png_F}
    return return_D






###############################################################################
#################################### < UTILS > ################################


def get_zero_to_one_val(
            val,
            min_val,
            max_val ):
    try:
        zero_to_one_val = float( val - min_val ) / (max_val - min_val)
    except ZeroDivisionError:
        print "min_val = max_val = {}, so instead of ZeroDivisionError, returning 0".format(
                max_val )
        zero_to_one_val = 0.
    try:
        assert( zero_to_one_val <= 1. )
    except AssertionError:
        print "val: {}".format( val )
        print "min_val: {}".format( min_val )
        print "max_val: {}".format( max_val )


    try:
        assert( zero_to_one_val >= 0. )
    except AssertionError:
        print "val: {}".format( val )
        print "min_val: {}".format( min_val )
        print "max_val: {}".format( max_val )
    return zero_to_one_val


def axis_labels_to_percent(
        y,
        position,
        num_decimals = None ):
    """
    - Transforms tick axis lables to include a percentage sign. Taken from:
        http://matplotlib.org/examples/pylab_examples/histogram_percent_demo.html

    - To utilize this in a plot for Y-AXIS labels, for example, include:

        formatter = FuncFormatter( axis_labels_to_percent )
        plt.gca().yaxis.set_major_formatter( formatter )
        plt.show()

    """
    # Ignore the passed in position. This has the effect of scaling the
    # default tick locations.
    #s = str(100 * y)

    # The percent symbol needs escaping in latex
    try:
        if matplotlib.rcParams['text.usetex'] == True:
            pass
    except NameError:
        from pylab import matplotlib
    if matplotlib.rcParams['text.usetex'] == True:
        if (num_decimals == 0):
            return str(int(y)) + r'$\%$'
        elif (num_decimals == 1):
            return str(int(y)) + r'{0:.1f}$\%$'.format( y )
        else:
            return str(y) + r'$\%$'
    else:
        return str(y) + '%'

################################### </ UTILS > ################################
###############################################################################





