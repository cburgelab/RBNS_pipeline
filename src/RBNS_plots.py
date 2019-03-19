#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
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


#### If there's an installation of LaTeX on you path (e.g., .../texlive/2013/bin/x86_64-linux)
#try:
#    plt.rc( 'text', usetex = True )
#    plt.rc( 'font', family = 'serif' )
#except:
#    pass

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
    """

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
    ax = fig.add_subplot(111)
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
        ax_colorbar_limits_L = [axis_limits_L[0] + axis_limits_L[2] + 0.01,
                axis_limits_L[1], 0.02, axis_limits_L[3]]

    ##### go through each of the squares
    for key_1_num, key_1 in enumerate( keys_1_L ):
        for key_2_num, key_2 in enumerate( keys_2_L ):

            val = checkerboard_D[key_1][key_2]

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

    ax.set_title(r"{}".format( title.replace("_", "\_") ),
            fontsize=title_fontsize)
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
        ax.set_yticklabels( str_keys_2_L )
        ax.tick_params(axis='y', labelsize = label_2_fontsize)
    else:
        ax.set_yticklabels( [""] * len(str_keys_2_L))
        ax.tick_params(axis='y', labelsize = 4 )

    ##### If the gridlines should be included
    if ( gridlines_color != False ):
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
    cb.set_label(r'{}'.format( colorbar_label ), size = 18)
    plt.gca().tick_params(axis='y', labelsize=18)

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
    - Makes a plot of the log2(enrichment) of the 4 base frequencies relative to
        uniform background across the positions of the reads for ALL submitted
        RBNS experiments

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
                    bottom="off", labelbottom = "off", labeltop = "off",
                    labelleft = "off", labelright= 'off')

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
        seqlogo_proportions_L,
        seqLogo_Fs_L,
        out_F,
        title = "",
        xlab = "",
        ylab = "",
        colors_for_bargraph = "default",
        collapse_less_than_x_perc_into_Other = 5. ):
    """
    - Plots a basic stacked histogram with the seqLogo_Fs_L, each having
        seqlogo_proportions_L weight (the summer stepwise R - 1 of kmer in
        that logo), to out_F

    - collapse_less_than_x_perc_into_Other:
        - if None: won't collapse any categories
        - Otherwise, a float between 0 and 100 (e.g., 5.). Individual
            categories less than this will be collapsed into an "other"
            category with no SeqLogo
    """
    import PIL
    from PIL import Image
    from matplotlib.backends.backend_pdf import PdfPages

    seqlogo_proportions_L = [float(x) for x in seqlogo_proportions_L]
    #### Normalize seqlogo_proportions_L so they all sum to 100.
    seqlogo_proportions_L = [x * 100./sum( seqlogo_proportions_L )\
            for x in seqlogo_proportions_L]

    #### Combine any nums that are individually less than 5% into 1 composite
    ####    "other" category
    if (type(collapse_less_than_x_perc_into_Other) is float) or\
            (type(collapse_less_than_x_perc_into_Other) is int):
        total_other = 0.
        pared_nums_L = []
        pared_seqLogo_Fs_L = []
        for num_idx, num in enumerate( seqlogo_proportions_L ):
            if (num < collapse_less_than_x_perc_into_Other):
                total_other += num
            else:
                pared_nums_L.append( num )
                try:
                    pared_seqLogo_Fs_L.append( seqLogo_Fs_L[num_idx] )
                except IndexError:
                    pass
        #### If there was any "total_other" added, include it
        if (total_other > 0.):
            pared_nums_L.append( total_other )
        seqlogo_proportions_L = pared_nums_L
        seqLogo_Fs_L = pared_seqLogo_Fs_L


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

    fig = plt.figure( figsize=(4,3), dpi = 300 )

    strt_ax = plt.gca()
    strt_ax.spines['top'].set_color('none')
    strt_ax.spines['bottom'].set_color('none')
    strt_ax.spines['left'].set_color('none')
    strt_ax.spines['right'].set_color('none')
    strt_ax.tick_params( top='off', left="off", right="off", bottom = "off" )

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
    ##################### </ the axis for the bar plot > ######################

    num_seqLogos = len( seqLogo_Fs_L )
    if (num_seqLogos <= 4):
        height_of_seqlogo = 0.43
        width_of_seqlogo = .39
        num_logos_right_of_bar = 2
        seqLogos_right_L = seqLogo_Fs_L[:2]
        seqLogos_left_L = seqLogo_Fs_L[2:]
    elif (num_seqLogos <= 6):
        height_of_seqlogo = 0.33333
        width_of_seqlogo = .2777
        num_logos_right_of_bar = 3
        seqLogos_right_L = seqLogo_Fs_L[:3]
        seqLogos_left_L = seqLogo_Fs_L[3:]
    else:
        height_of_seqlogo = 0.25
        width_of_seqlogo = .2083
        num_logos_right_of_bar = 4
        seqLogos_right_L = seqLogo_Fs_L[:4]
        seqLogos_left_L = seqLogo_Fs_L[4:8]

    # if the length of SeqLogos is less than the length of seqlogo_proportions_L,
    # the last number is the "other signficant enriched" category
    if (len(seqLogo_Fs_L) < len(seqlogo_proportions_L)):
        other_significant = True
    else:
        other_significant = False

    handles_L = []
    cur_bottom = 100.
    bottoms_L = []
    for num, val in enumerate( seqlogo_proportions_L ):
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
            y_to_connect_line_to = (cur_bottom/100.) + (seqlogo_proportions_L[num]/200.)
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
            SeqLogo_F = seqLogo_Fs_L[num]
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

    if other_significant:
        x_left = .6
        y_to_draw_line_to = (seqlogo_proportions_L[-1]/200.)
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








def plot_Ppaired_ratio_of_motif(
        avg_fold_probs_by_motif_and_pos_D,
        annots_L,
        motifs_of_interest_L,
        pulldown_annot,
        out_F_start,
        read_len,
        title = "",
        plot_signif = True,
        skip_if_already_exists = True,
        dicts_DIR = "/net/nevermind/data/nm/RBNS_results/pulldown_libs/analyses/structure",
        include_legend = True ):
    """
    - Makes a plot of P_paired(pulldown)/P_paired(input) for each motif (kmer)
        in motifs_of_interest_L

    - avg_fold_probs_by_motif_and_pos_D has format:
        avg_fold_probs_by_motif_and_pos_D[lib_annot][motif_num][motif_pos]=\
                motif_avg_bp_prob

    - lib_annot is like: "input" or "80_nM"

    - out_F_start is a file start, to which "Ppairedanalysis_{kmer}.pdf" will be added for
        each of the kmers in motifs_of_interest_L

    - If dicts_DIR exists, it should have the mean & stdev. of each kmer in 0 nM
        libraries by position within the motif, to get significance - if these
        do not exist, it will be skipped
    """
    k = len( motifs_of_interest_L[0] )

    ##### The dictionary of folded input reads of a particular read_len to use
    ####    for the mean & st.dev. for determining significance
    mean_STD_Ppaired_L_by_pos_for_each_kmer_D = {}
    if plot_signif:
        Ds_DIR = os.path.join( dicts_DIR, "all_RBPs_0nM_fld_all{0}mers".format( k ) )
        for kmer in motifs_of_interest_L:
            #### Each D has format:
            ## {0: {'mean': 0.99960544423771491, 'std_ratio': 0.015726334560542742},
            ## 1: {'mean': 0.99460381791438968, 'std_ratio': 0.005669277527848831},
            ## 2: {'mean': 0.99608373137948758, 'std_ratio': 0.005302256957967216},
            ## 3: {'mean': 0.99783444433615498, 'std_ratio': 0.0051178830567742625},
            ## 4: {'mean': 0.99797767093026613, 'std_ratio': 0.0081676955702984725}}
            D_basename = "0nM.{0}.len_{1}.mean_and_std_by_pos_D.pkl".format(
                    kmer,
                    read_len )
            D_F = os.path.join( Ds_DIR, D_basename )
            if not os.path.exists( D_F ):
                plot_signif = False
            else:
                mean_STD_Ppaired_L_by_pos_for_each_kmer_D[kmer] =\
                    pickle.load( open( D_F ) )


    #### get the concentrations to plot (all except the input lib.)
    conc_keys_to_plot_L = []
    for annot in annots_L:
        if ( annot.find( "input" ) == -1 ):
            conc_keys_to_plot_L.append( annot )

    # sort the concentrations to plot in ascending order
    try:
        conc_keys_to_plot_L.sort( key = lambda x: int( x.split("_")[0] ) )
    except:
        conc_keys_to_plot_L.sort()
    num_concs_to_plot = len( conc_keys_to_plot_L )
    labels_to_plot_L = [x.replace( "_", " " ) for x in conc_keys_to_plot_L]

    # go through each of the kmers of interest, making a new plot for each one
    PDF_Fs_by_motifidx_D = {}
    PDF_Fs_by_motif_D = {}
    lower_upper_sig_by_kmer_pos_D = {}
    ratios_by_motif_pos_D = {}
    sigB_by_motif_pos_D = {}
    for motif_num, motif in enumerate( motifs_of_interest_L ):

        ratios_by_motif_pos_D[motif] = {}
        sigB_by_motif_pos_D[motif] = {}

        by_motif_D = {}
        lower_upper_sig_by_kmer_pos_D[motif] = {}

        out_F = out_F_start + "_Ppairedanalysis_ratio.{0}.pdf".format( motif )
        out_abs_F = out_F_start + "_Ppairedanalysis_abs.{0}.abs.pdf".format( motif )
        out_merged_F = out_F_start + "_Ppaired_analysis.{0}.pdf".format( motif )
        two_up_merged_F = out_F_start + "_Ppaired_analysis.{0}-2up.pdf".format( motif )
        if os.path.exists( out_merged_F ) and\
                skip_if_already_exists:
            continue

        out_Ds_DIR = os.path.join( os.path.dirname( out_F_start ), 'Ds' )
        os.system( "mkdir -p {}".format( out_Ds_DIR ) )
        out_D_F = os.path.join( out_Ds_DIR,
                "Ppaired_Ls_by_libannot_D.{0}.{1}.pkl".format( motif_num, motif ) )

        # a list of the 0 concentration ratios, which will be used to get the
        # st.dev. and mean for determining if the pulldown library is signif.
        zero_conc_ratios_L = []
        pulldown_conc_ratios_L = []

        abs_Ppaireds_by_motifpos_annot_D = {}

        # get the input avg. bp of prob.
        input_avg_prob_by_position_D = {}
        conc_avg_prob_by_position_D = {}
        L_of_Ls_to_plot = []
        for motif_pos in range(len( motif )):

            abs_Ppaireds_by_motifpos_annot_D[motif_pos] = {}
            ratios_to_plot_L = []
            # go through each of the concentrations to plot
            for lib_annot in conc_keys_to_plot_L:

                conc_avg_prob =\
                    avg_fold_probs_by_motif_and_pos_D[lib_annot][motif_num][motif_pos]
                try:
                    input_avg_prob =\
                        avg_fold_probs_by_motif_and_pos_D["input"][motif_num][motif_pos]
                except KeyError:
                    input_avg_prob = avg_fold_probs_by_motif_and_pos_D["{0}_input".format(lib_annot)][motif_num][motif_pos]
                abs_Ppaireds_by_motifpos_annot_D[motif_pos][lib_annot] = conc_avg_prob
                abs_Ppaireds_by_motifpos_annot_D[motif_pos]['input'] = input_avg_prob

                ratio = conc_avg_prob / input_avg_prob
                #print "Ratio for motif {0}, position {1}, library {2} is: {3}".format(
                #        motif, motif_pos, lib_annot, ratio)
                ratios_to_plot_L.append( ratio )
                if ( lib_annot == "0_nM" ):
                    zero_conc_ratios_L.append( ratio )
                elif ( lib_annot == pulldown_annot ):
                    pulldown_conc_ratios_L.append( ratio )
                    ratios_by_motif_pos_D[motif][motif_pos] = ratio
                try:
                    by_motif_D[lib_annot].append( ratio )
                except KeyError:
                    by_motif_D[lib_annot] = [ratio]

            L_of_Ls_to_plot.append( ratios_to_plot_L )

            #### If the +/- 2 std. ratios for significance should be calculated
            if plot_signif:
                mean = mean_STD_Ppaired_L_by_pos_for_each_kmer_D[motif][motif_pos]["mean"]
                std = mean_STD_Ppaired_L_by_pos_for_each_kmer_D[motif][motif_pos]["std_ratio"]
                ratio_below_signif = mean - (2*std)
                ratio_above_signif = mean + (2*std)
                lower_upper_sig_by_kmer_pos_D[motif][motif_pos] =\
                        {"lower_sig_ratio": ratio_below_signif,
                                "upper_sig_ratio": ratio_above_signif}

        if ( len( avg_fold_probs_by_motif_and_pos_D ) > 12 ):
            many_legend_annots = True
            figsize=(6, 4)
        else:
            many_legend_annots = False
            figsize = (7, 3)
        fig2 = plt.figure( dpi = 300, figsize = figsize )
        plt.tick_params(\
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are
                bottom='off',      # ticks along the bottom edge are
                top='off',         # ticks along the top edge are
                labelbottom='off')
        ax2 = fig2.add_axes( [0.05, 0.05, 0.9, 0.9] )
        ax2.spines['top'].set_color('none')
        ax2.spines['bottom'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.spines['right'].set_color('none')
        ax2.tick_params(which = 'both', top='off',left="off",right="off",\
                bottom="off", labelbottom = "off", labeltop = "off",
                labelleft = "off", labelright= 'off')
        # THE WIDTH OF ONE LETTER'S SPACE
        one_letter_width = 0.11
        # the LEFTMOST VERTICAL LINE (also the start of the main horiz. line)
        main_H_line_min = 0.1
        main_H_line_max = main_H_line_min + (one_letter_width * len( motif ))

        if ( title != "" ):
            plt.annotate( title, ((main_H_line_min+main_H_line_max)/2., 0.8),
                    fontsize=18,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'center')

        LFC_height = 0.2
        LFC_0 = 0.62

        # offset to set the colored bars
        offset_for_bars = one_letter_width*0.5

        in_bar_width = one_letter_width / (num_concs_to_plot + 1)
        max_ratio = 0.
        min_ratio = 20.
        for pos, nt in enumerate( motif ):
            enrich_handles_L = 0
            handles_L = []
            this_letter_center = main_H_line_min + ((pos+0.5)*one_letter_width)
            # get the enrichment of each of the 4 types
            for struc_num, ratio in enumerate(L_of_Ls_to_plot[pos]):

                # convert the ratio to a LFC
                try:
                    lfc = math.log( ratio, 2 )
                except ValueError:
                    print "VALUE ERROR FOR motif {0}, position {1}".format(
                            motif, pos)
                    lfc = 0.

                max_ratio = max(max_ratio, ratio)
                min_ratio = min(min_ratio, ratio)

                lft = this_letter_center + (((-0.5*num_concs_to_plot)+struc_num)*in_bar_width)
                btm = min( LFC_0, LFC_0+(lfc*LFC_height) )
                wdth = in_bar_width
                hgh = abs(lfc) * LFC_height
                axb = fig2.add_axes([lft, btm, wdth, hgh])
                axb.set_xlim(0, 1.)
                axb.set_ylim(0, 1.)
                axb.spines['top'].set_color('none')
                axb.spines['bottom'].set_color('none')
                axb.spines['left'].set_color('none')
                axb.spines['right'].set_color('none')
                #ax.tick_params(labelcolor="w", top='off',left="off",right="off",\
                        #        bottom="off")
                axb.tick_params(which = 'both', top='off',left="off",right="off",\
                    bottom="off", labelbottom = "off", labeltop = "off",
                    labelleft = "off", labelright= 'off')
                color_num_to_use = light_to_dark_D[num_concs_to_plot][struc_num]
                hand = axb.bar(0., 1., width = 1.,
                        #color = shades_of_red_D[color_num_to_use], edgecolor = 'k')
                        color = shades_of_blue_D[color_num_to_use], edgecolor = 'k')
                handles_L.append( hand )

        # add lines at 0 and +/-1 LFC
        axbg = fig2.add_axes( [0., 0., 1., 1.] )
        axbg.patch.set_visible(False)
        axbg.spines['top'].set_color('none')
        axbg.spines['bottom'].set_color('none')
        axbg.spines['left'].set_color('none')
        axbg.spines['right'].set_color('none')
        axbg.tick_params(which = 'both', top='off',left="off",right="off",\
                bottom="off", labelbottom = "off", labeltop = "off",
                labelleft = "off", labelright= 'off')


        potential_ratios_to_plot = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2,
                1.5, 1.7, 2, 2.5, 3, 5, 10]
        # get the lowest ratio to plot
        min_ratio_to_plot = 0.
        min_ratio_to_plot_idx = -1
        max_ratio_to_plot = 20.
        max_ratio_to_plot_idx = -1
        # get the
        for ratio_idx, ratio in enumerate( potential_ratios_to_plot ):
            if (ratio < min_ratio):
                min_ratio_to_plot = ratio
                min_ratio_to_plot_idx = ratio_idx
        for ratio_idx, ratio in enumerate( potential_ratios_to_plot ):
            if (ratio < max_ratio):
                pass
            else:
                max_ratio_to_plot = ratio
                max_ratio_to_plot_idx = ratio_idx
                break

        axbg.axhline( y=LFC_0, xmin = main_H_line_min, xmax = main_H_line_max,
                color='k', linestyle='-', linewidth = 2)
        plt.annotate("P$_{paired}$(RBNS)/", (0.03, LFC_0), fontsize=12,
                xycoords='figure fraction',
                verticalalignment = 'center', horizontalalignment = 'right',
                rotation = 90)
        plt.annotate("P$_{paired}$(input)", (0.053, LFC_0), fontsize=12,
                xycoords='figure fraction',
                verticalalignment = 'center', horizontalalignment = 'right',
                rotation = 90)
        # go through and annotate all of the desired ratios
        for ratio_idx in range( min_ratio_to_plot_idx, max_ratio_to_plot_idx + 1 ):
            ratio = potential_ratios_to_plot[ratio_idx]
            height_on_axis = LFC_0 + (LFC_height * math.log( ratio, 2 ))
            plt.annotate(str( ratio ), (main_H_line_min - 0.01, height_on_axis),
                    fontsize=12,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'right')
            plt.annotate(str( ratio ), (main_H_line_max + 0.01, height_on_axis),
                    fontsize=12,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'left')
            # also add mini-horizontal lines inside the plot
            axbg.axhline( y=height_on_axis,
                    xmin = main_H_line_min,
                    xmax = main_H_line_min + 0.01,
                    color='k', linestyle='-', linewidth = 1)
            axbg.axhline( y=height_on_axis,
                    xmin = main_H_line_max - 0.01,
                    xmax = main_H_line_max,
                    color='k', linestyle='-', linewidth = 1)
        # make a vertical line at the left and right ends
        min_v_line = LFC_0 + (LFC_height * math.log( min_ratio_to_plot, 2 )) + 0.002
        max_v_line = LFC_0 + (LFC_height * math.log( max_ratio_to_plot, 2 )) - 0.002
        axbg.axvline( x = main_H_line_min,
                ymin = min_v_line,
                ymax = max_v_line,
                color='k', linestyle='-', linewidth = 2)
        axbg.axvline( x = main_H_line_max,
                ymin = min_v_line,
                ymax = max_v_line,
                color='k', linestyle='-', linewidth = 2)

        if (include_legend == True):
            # if motif is below length 7, plot it to the right
            if many_legend_annots:
                legend_fontsize = 8
                box_height = 0.035 * num_concs_to_plot
                if (len(motif) < 7):
                    leg = axbg.legend( handles_L, labels_to_plot_L,
                            loc=(main_H_line_max + 0.05, LFC_0 - (box_height/2.)),
                            fancybox=True, fontsize = legend_fontsize, borderaxespad=0)
                else:
                    leg = axbg.legend( handles_L, labels_to_plot_L,
                            loc=(main_H_line_min + ((len(motif)/2.)-0.5)*one_letter_width,
                                max_v_line),
                            fancybox=True, fontsize=11, borderaxespad=0)
            else:
                legend_fontsize = 11
                box_height = 0.08 * num_concs_to_plot
                if (len(motif) < 7):
                    leg = axbg.legend( handles_L, labels_to_plot_L,
                            loc=(main_H_line_max + 0.05, LFC_0 - (box_height/2.)),
                            fancybox=True, fontsize = legend_fontsize, borderaxespad=0)
                else:
                    leg = axbg.legend( handles_L, labels_to_plot_L,
                            loc=(main_H_line_min + ((len(motif)/2.)-0.5)*one_letter_width,
                                max_v_line),
                            fancybox=True, fontsize=11, borderaxespad=0)


        # annotate the letters below the bars
        for pos, nt in enumerate( motif.replace("T", "U") ):
            this_letter_center = main_H_line_min + ((pos+0.5)*one_letter_width)
            plt.annotate(nt, (this_letter_center, min_v_line - 0.02), fontsize= 24,
                    xycoords='figure fraction',
                    verticalalignment = 'top', horizontalalignment = 'center')
            # annotate with a star if the pulldown ratio is significant
            if plot_signif:
                pulldown_ratio = pulldown_conc_ratios_L[pos]
                ratio_below_signif = lower_upper_sig_by_kmer_pos_D[motif][pos]["lower_sig_ratio"]
                ratio_above_signif = lower_upper_sig_by_kmer_pos_D[motif][pos]["upper_sig_ratio"]
                if (pulldown_ratio >= ratio_above_signif) or\
                        (pulldown_ratio <= ratio_below_signif):
                    sigB_by_motif_pos_D[motif][pos] = True
                    plt.annotate("*", (this_letter_center, min_v_line - 0.16), fontsize=12,
                            xycoords='figure fraction',
                            verticalalignment = 'center', horizontalalignment = 'center')
                else:
                    sigB_by_motif_pos_D[motif][pos] = False

        fig2.savefig( out_F )

        with open( out_D_F, "wb" ) as f:
            pickle.dump( by_motif_D, f )

        PDF_Fs_by_motifidx_D[motif_num] = out_F
        PDF_Fs_by_motif_D[motif] = out_F

        #######################################################################
        #######################################################################
        ####### < Now make a stacked barplot of the *absolute* Ppaired > ######
        num_concs_to_plot_inc_input = num_concs_to_plot + 1
        fig3 = plt.figure( dpi = 300, figsize = figsize )
        plt.tick_params(\
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are
                bottom='off',      # ticks along the bottom edge are
                top='off',         # ticks along the top edge are
                labelbottom='off')
        ax2 = fig3.add_axes( [0.05, 0.05, 0.9, 0.9] )
        ax2.spines['top'].set_color('none')
        ax2.spines['bottom'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.spines['right'].set_color('none')
        ax2.tick_params(which = 'both', top='off',left="off",right="off",\
                bottom="off", labelbottom = "off", labeltop = "off",
                labelleft = "off", labelright= 'off')
        if ( title != "" ):
            plt.annotate( title, ((main_H_line_min+main_H_line_max)/2., 0.85),
                    fontsize=18,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'center')
        # THE WIDTH OF ONE LETTER'S SPACE
        one_letter_width = 0.11
        # the LEFTMOST VERTICAL LINE (also the start of the main horiz. line)
        main_H_line_min = 0.1
        main_H_line_max = main_H_line_min + (one_letter_width * len( motif ))

        height_of_0Ppaired = .2
        height_of_1Ppaired = .8

        # offset to set the colored bars
        offset_for_bars = one_letter_width*0.5

        in_bar_width = one_letter_width / (num_concs_to_plot_inc_input + 1)
        max_ratio = 0.
        min_ratio = 20.
        for pos, nt in enumerate( motif ):

            enrich_handles_L = 0
            handles_L = []
            this_letter_center = main_H_line_min + ((pos+0.5)*one_letter_width)

            #### Go through all of the lib annots
            for struc_num, lib_annot in enumerate( ['input'] + conc_keys_to_plot_L ):

                Ppaired = abs_Ppaireds_by_motifpos_annot_D[pos][lib_annot]


                lft = this_letter_center + (((-0.5*num_concs_to_plot_inc_input)+struc_num)*in_bar_width)
                btm = height_of_0Ppaired
                wdth = in_bar_width
                hgh = (height_of_1Ppaired - height_of_0Ppaired) * Ppaired
                axb = fig3.add_axes([lft, btm, wdth, hgh])
                axb.set_xlim(0, 1.)
                axb.set_ylim(0, 1.)
                axb.spines['top'].set_color('none')
                axb.spines['bottom'].set_color('none')
                axb.spines['left'].set_color('none')
                axb.spines['right'].set_color('none')
                #ax.tick_params(labelcolor="w", top='off',left="off",right="off",\
                        #        bottom="off")
                axb.tick_params(which = 'both', top='off',left="off",right="off",\
                    bottom="off", labelbottom = "off", labeltop = "off",
                    labelleft = "off", labelright= 'off')

                if ( lib_annot == 'input' ):
                    color = '#808080'
                else:
                    color_num_to_use = light_to_dark_D[num_concs_to_plot][struc_num-1]
                    color = shades_of_orange_D[color_num_to_use]
                hand = axb.bar(0., 1., width = 1.,
                        #color = shades_of_red_D[color_num_to_use], edgecolor = 'k')
                        color = color, edgecolor = 'k')
                handles_L.append( hand )

        axbg = fig3.add_axes( [0., 0., 1., 1.] )
        axbg.patch.set_visible(False)
        axbg.spines['top'].set_color('none')
        axbg.spines['bottom'].set_color('none')
        axbg.spines['left'].set_color('none')
        axbg.spines['right'].set_color('none')
        axbg.tick_params(which = 'both', top='off',left="off",right="off",\
                bottom="off", labelbottom = "off", labeltop = "off",
                labelleft = "off", labelright= 'off')

        ##### Plot the bottom & right & left lines
        axbg.axhline( y = height_of_0Ppaired, xmin = main_H_line_min, xmax = main_H_line_max,
                color='k', linestyle='-', linewidth = 2)
        plt.annotate("P$_{paired}$", (0.053, (height_of_0Ppaired+height_of_1Ppaired)/2),
                fontsize=14,
                xycoords='figure fraction',
                verticalalignment = 'center', horizontalalignment = 'right',
                rotation = 90)
        axbg.axvline( x = main_H_line_min,
                ymin = height_of_0Ppaired,
                ymax = height_of_1Ppaired,
                color='k', linestyle='-', linewidth = 2)
        axbg.axvline( x = main_H_line_max,
                ymin = height_of_0Ppaired,
                ymax = height_of_1Ppaired,
                color='k', linestyle='-', linewidth = 2)

        for Ppaired in [0., 0.2, 0.4, 0.6, 0.8, 1.]:

            height_on_axis = height_of_0Ppaired + ( height_of_1Ppaired - height_of_0Ppaired ) * Ppaired

            plt.annotate( "{:.1f}".format( Ppaired ),
                    (main_H_line_min - 0.01, height_on_axis),
                    fontsize=12,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'right')
            plt.annotate( "{:.1f}".format( Ppaired ),
                    (main_H_line_max + 0.01, height_on_axis),
                    fontsize=12,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'left')
            # also add mini-horizontal lines inside the plot
            if ( Ppaired > 0. ):
                axbg.axhline( y=height_on_axis,
                        xmin = main_H_line_min,
                        xmax = main_H_line_min + 0.01,
                        color='k', linestyle='-', linewidth = 1)
                axbg.axhline( y=height_on_axis,
                        xmin = main_H_line_max - 0.01,
                        xmax = main_H_line_max,
                        color='k', linestyle='-', linewidth = 1)

        if (include_legend == True):
            legend_fontsize = 8
            box_height = 0.035 * num_concs_to_plot_inc_input
            if (len(motif) < 7):
                leg = axbg.legend( handles_L, ['Input'] + labels_to_plot_L,
                        loc = ( main_H_line_max + 0.05, height_of_0Ppaired + 0.12 ),
                        fancybox = True, fontsize = legend_fontsize + 2, borderaxespad=0 )

        # annotate the letters below the bars
        for pos, nt in enumerate( motif.replace("T", "U") ):
            this_letter_center = main_H_line_min + ((pos+0.5)*one_letter_width)
            plt.annotate(nt, (this_letter_center, height_of_0Ppaired - 0.02), fontsize= 24,
                    xycoords='figure fraction',
                    verticalalignment = 'top', horizontalalignment = 'center')
        fig3.savefig( out_abs_F )

        #### Merge the out_F & out_abs_F into the out_merged_F
        RBNS_utils.merge_PDFs_mult_on_page(
                [out_abs_F, out_F],
                out_merged_F,
                direc = 'vertical' )

        #### Move the -2up.pdf to the original merged F
        #os.system( "mv {0} {1}".format( two_up_merged_F, out_merged_F ) )

    return_D = {
            "ratios_by_motif_pos_D": ratios_by_motif_pos_D,
            "sigB_by_motif_pos_D": sigB_by_motif_pos_D,
            "PDF_Fs_by_motifidx_D": PDF_Fs_by_motifidx_D,
            "PDF_Fs_by_motif_D": PDF_Fs_by_motif_D }
    return return_D






def plot_enrichment_by_5_Ppaired_bins(
        enrichments_by_kmer_conc_bin_D,
        annots_L,
        motifs_of_interest_L,
        Ppaired_bin_upper_limits_L,
        out_F_start,
        read_len,
        conc_keys_to_plot_override_L = [],
        plot_signif = True,
        title = "",
        dicts_DIR = "/net/nevermind/data/nm/RBNS_results/pulldown_libs/analyses/structure",
        include_legend = True ):
    """
    - out_F_start is a file start, to which "enrichment_by_Ppairedbin_{kmer}.pdf" will be added for
        each of the kmers in motifs_of_interest_L

    - If dicts_DIR exists, it should have the mean & stdev. of each kmer in 0 nM
        libraries by position within the motif, to get significance - if these
        do not exist, it will be skipped
    """
    k = len( motifs_of_interest_L[0] )

    sig_Ds_by_kmer_D = {}
    if plot_signif:
        for kmer in motifs_of_interest_L:
            D_F = os.path.join( dicts_DIR, "all_RBPs_RbyPpairedBin_0nM_fld_all{}mers".format( k ),
                    "0nM.{0}.len_{1}.mean_and_std_of_adjbins_by_lowerbin_D.pkl".format(
                        kmer, read_len ) )
            if os.path.exists( D_F ):
                sig_Ds_by_kmer_D[kmer] = pickle.load( open( D_F ) )
            else:
                plot_signif = False

    motifs_of_interest_L = [x.replace("U", "T") for x in motifs_of_interest_L]
    Ppaired_bin_upper_limits_L.sort()

    #### If the concentrations were providen directly (e.g., if weighted_by_CG)
    if ( len( conc_keys_to_plot_override_L ) > 0 ):
        conc_keys_to_plot_L = conc_keys_to_plot_override_L
    # get the concentrations to plot (all except the input lib.)
    else:
        conc_keys_to_plot_L = []
        for annot in annots_L:
            if ( annot.find( "input") == -1 ):
                conc_keys_to_plot_L.append( annot )

    # sort the concentrations to plot in ascending order
    try:
        conc_keys_to_plot_L.sort( key = lambda x: float(x.split("_")[0]) )
    except ValueError:
        conc_keys_to_plot_L.sort()

    num_concs_to_plot = len( conc_keys_to_plot_L )
    num_Ppaired_bins = len( Ppaired_bin_upper_limits_L )

    sigB_by_kmer_conc_bin_D = {}
    for kmer in motifs_of_interest_L:
        sigB_by_kmer_conc_bin_D[kmer] = {}
        for conc in conc_keys_to_plot_L:
            sigB_by_kmer_conc_bin_D[kmer][conc] = {}

    # make the labels for the different Ppaired bins
    labels_to_plot_L = []
    lower_limit = 0
    lower_limit_string = ""
    for bin_num, upper_limit in enumerate(Ppaired_bin_upper_limits_L):
        upper_limit_string = "{0:.2f}".format(upper_limit).rstrip("0")
        if (bin_num == 0):
            label = "0-{0}".format( upper_limit_string )
        elif (bin_num == len(Ppaired_bin_upper_limits_L) - 1):
            label = "{0}-1".format( lower_limit_string )
        else:
            label = "{0}-{1}".format( lower_limit_string, upper_limit_string )
        lower_limit_string = upper_limit_string
        labels_to_plot_L.append( label )

    # go through each of the kmers of interest, making a new plot for each one
    for motif_num, kmer in enumerate( motifs_of_interest_L ):

        out_F = out_F_start + "_enrichment_by_Ppairedbin_{0}.pdf".format(
                kmer.replace("T","U") )
        if os.path.exists( out_F ):
            continue

        pulldown_conc_ratios_L = []

        # get the input avg. bp of prob.
        input_avg_prob_by_position_D = {}
        conc_avg_prob_by_position_D = {}
        L_of_Ls_to_plot = []

        max_enrich= 0.
        for conc in conc_keys_to_plot_L:
            ratios_to_plot_L = []
            for Ppaired_bin in range(len(Ppaired_bin_upper_limits_L)):
                enrich = enrichments_by_kmer_conc_bin_D[kmer][conc][Ppaired_bin]
                ratios_to_plot_L.append( enrich )
                max_enrich = max( max_enrich, enrich )
            L_of_Ls_to_plot.append( ratios_to_plot_L )

        # get the vertical enrichments to plot, based on the max. enrichment
        potential_upper_enriches_L = [2, 3, 4, 5, 8, 10, 15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 150, 200, 250, 300]
        upper_lim_found_B = False
        for upper_lim in potential_upper_enriches_L:
            if (upper_lim > max_enrich) and (upper_lim_found_B == False):
                upper_lim_found_B = True
                if (upper_lim == 2):
                    enriches_to_mark_L = [1,1.2,1.4,1.6,1.8,2.]
                elif (upper_lim == 3):
                    enriches_to_mark_L = [1,1.4,1.8,2.2,2.6,3]
                elif (upper_lim <= 8):
                    enriches_to_mark_L = range( 1, upper_lim + 1 )
                elif (upper_lim == 10):
                    enriches_to_mark_L = [1, 2, 4, 6, 8, 10]
                elif (upper_lim <= 50):
                    enriches_to_mark_L = [1] + range( 5, upper_lim+1, 5)
                elif (upper_lim <= 120):
                    enriches_to_mark_L = [1] + range( 10, upper_lim+1, 10)
                else:
                    enriches_to_mark_L = [1] + range( 50, upper_lim+1, 50)
                break

        fig2 = plt.figure( dpi=300 )
        plt.tick_params(\
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are
                bottom='off',      # ticks along the bottom edge are
                top='off',         # ticks along the top edge are
                labelbottom='off')
        #ax = fig.add_axes([.2,0,.6,1])
        #ax.spines['top'].set_color('none')
        #ax.spines['bottom'].set_color('none')
        #ax.tick_params(labelcolor="w", top='off',bottom="off")
        ax_for_lines = fig2.add_axes( [0., 0., 1., 1.], zorder = 2, alpha = 0. )
        #ax_for_lines.patch.set_visible(False)
        ax_for_lines.spines['top'].set_color('none')
        ax_for_lines.spines['bottom'].set_color('none')
        ax_for_lines.spines['left'].set_color('none')
        ax_for_lines.spines['right'].set_color('none')
        ax_for_lines.tick_params(which = 'both', top='off',left="off",right="off",\
                bottom="off", labelbottom = "off", labeltop = "off",
                labelleft = "off", labelright= 'off')

        axis_x_left = 0.05
        axis_width = 0.9
        axis_y_bottom = 0.05
        axis_height = 0.9

        #ax2 = fig2.add_axes(
        #        [axis_x_left, axis_y_bottom, axis_width, axis_height] )
        #ax2.spines['top'].set_color('none')
        #ax2.spines['bottom'].set_color('none')
        #ax2.spines['left'].set_color('none')
        #ax2.spines['right'].set_color('none')
        #ax2.tick_params(which = 'both', top='off',left="off",right="off",\
        #ax2.tick_params(which = 'both', top='off',\
                #bottom="off", labelbottom = "off", labeltop = "off",
                #labelleft = "off", labelright= 'off')
        # THE WIDTH OF ONE LETTER'S SPACE
        one_letter_width = axis_width / (len(conc_keys_to_plot_L) + 1.)
        # the LEFTMOST VERTICAL LINE (also the start of the main horiz. line)
        num_letters_for_buffer_width = 0.35
        main_H_line_min = 0.1
        main_H_line_max = main_H_line_min + (one_letter_width*(num_concs_to_plot + num_letters_for_buffer_width))

        # get the middle of the horizontal line, which will be used for the
        # center of the title
        middle_of_plot = (main_H_line_max + main_H_line_min)/2.

        enrich_0 = 0.2
        enrich_max = 0.8
        one_enrich_unit_height = (enrich_max - enrich_0)/enriches_to_mark_L[-1]
        # the middle of the plot vertically, for the center of the legend
        middle_plot_vertical = (enrich_0 + enrich_max)/2.

        # offset to set the colored bars
        offset_for_bars = one_letter_width*0.5

        in_bar_width = one_letter_width / (num_Ppaired_bins + 1)

            ###### Plot the significance over all 4 bins
            #if plot_signif:

            #    abs_delta_last_minus_first = abs( L_of_Ls_to_plot[pos][-1] - L_of_Ls_to_plot[pos][0] )
            #    #### Divide this by the overall enrichment so the scales are the same
            #    normed_abs_delta_last_minus_first = abs_delta_last_minus_first / avg_enrich_all_bins

            #    std_this_motif_all_bins = sig_Ds_by_kmer_D[kmer]['0_to_4']['std']
            #    if ( normed_abs_delta_last_minus_first >= 2 * std_this_motif_all_bins ):

            #        lft_of_conc = this_letter_center - ( 2.3 * in_bar_width )
            #        right_of_conc = this_letter_center + ( 2.3 * in_bar_width )

            #        ax_for_lines.axhline( y = enrich_0 - 0.03,
            #            xmin = lft_of_conc, xmax = right_of_conc,
            #            color = 'k', linestyle = '-', linewidth = 1 )
            #        ax_for_lines.axvline( x = lft_of_conc, ymin = enrich_0 - 0.028, ymax = enrich_0,
            #            color = 'k', linestyle = '-', linewidth = 1 )
            #        ax_for_lines.axvline( x = right_of_conc, ymin = enrich_0 - 0.028, ymax = enrich_0,
            #            color = 'k', linestyle = '-', linewidth = 1 )
            #        plt.annotate( "**", ( this_letter_center, enrich_0 - 0.045 ), fontsize=12,
            #                xycoords='figure fraction',
            #                verticalalignment = 'center', horizontalalignment = 'center' )
            #        sigB_by_kmer_conc_bin_D[kmer][conc]["0_to_4"] = True
            #    else:
            #        sigB_by_kmer_conc_bin_D[kmer][conc]["0_to_4"] = False

        enrich_0 = 0.2
        enrich_max = 0.8
        axbg = fig2.add_axes( [0., 0., 1., 1.], zorder = 2, alpha = 0. )
        #axbg.patch.set_visible(False)
        axbg.spines['top'].set_color('none')
        axbg.spines['bottom'].set_color('none')
        axbg.spines['left'].set_color('none')
        axbg.spines['right'].set_color('none')
        axbg.tick_params(which = 'both', top='off',left="off",right="off",\
                bottom="off", labelbottom = "off", labeltop = "off",
                labelleft = "off", labelright= 'off')
        # make a dotted horizontal line for the average enrichment for each
        # concentration
        avg_Rs_L = []
        for pos, conc in enumerate( conc_keys_to_plot_L ):
            this_letter_center = main_H_line_min + ((pos+0.5 +\
                (num_letters_for_buffer_width)/2.)*one_letter_width)
            avg_enrich = enrichments_by_kmer_conc_bin_D[kmer][conc]["overall"]
            avg_Rs_L.append( avg_enrich )
            #height_on_axis = enrich_0 + (avg_enrich*one_enrich_unit_height)
            #axbg.axhline( y = height_on_axis,
            #        xmin = this_letter_center - (0.5*one_letter_width),
            #        xmax = this_letter_center + (0.5*one_letter_width),
            #        color='k', linestyle='--', linewidth = 1, zorder = 3 )

        # the horizontal line for an enrichment of 0
        #axbg.axhline( y=LFC_0-LFC_height, xmin=0.065, xmax=0.065+0.11*(len(motif)),
        #        color='k', linestyle='--', linewidth = 1)
        #axbg.axhline( y=LFC_0+LFC_height, xmin=0.074, xmax=0.065+0.11*(len(motif)),
        #        color='k', linestyle='--', linewidth = 1)
        enrich_0 = 0.2
        enrich_max = 0.8
        plt.annotate("RBNS $R$", (0.053, (enrich_max+enrich_0)/2. ), fontsize=12,
                xycoords='figure fraction',
                verticalalignment = 'center', horizontalalignment = 'right',
                rotation = 90)
        # go through and annotate all of the desired ratios
        for enrich_to_mark in enriches_to_mark_L:
            enrich_0 = 0.2
            enrich_max = 0.8
            height_on_axis = enrich_0 + (enrich_to_mark*one_enrich_unit_height)
            plt.annotate(str( float(enrich_to_mark) ).rstrip("0").rstrip("."),
                    (main_H_line_min - 0.01, height_on_axis),
                    fontsize=12,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'right')
            plt.annotate(str( float(enrich_to_mark) ).rstrip("0").rstrip("."),
                    (main_H_line_max + 0.01, height_on_axis),
                    fontsize=12,
                    xycoords='figure fraction',
                    verticalalignment = 'center', horizontalalignment = 'left')
            # also add mini-horizontal lines inside the plot
            axbg.axhline( y=height_on_axis,
                    xmin = main_H_line_min,
                    xmax = main_H_line_min + 0.01,
                    color='k', linestyle='-', linewidth = 1)
            axbg.axhline( y=height_on_axis,
                    xmin = main_H_line_max - 0.01,
                    xmax = main_H_line_max,
                    color='k', linestyle='-', linewidth = 1)
        # make a vertical line at the left and right ends
        min_v_line = enrich_0
        max_v_line = enrich_max - 0.002
        axbg.axvline( x = main_H_line_min,
                ymin = min_v_line,
                ymax = max_v_line,
                color='k', linestyle='-', linewidth = 2)
        axbg.axvline( x = main_H_line_max,
                ymin = min_v_line,
                ymax = max_v_line,
                color='k', linestyle='-', linewidth = 2)
        axbg.axhline( y=min_v_line, xmin = main_H_line_min, xmax = main_H_line_max,
                color='k', linestyle='-', linewidth = 2)

        # annotate the concentrations
        for pos, conc in enumerate( conc_keys_to_plot_L ):
            this_letter_center = main_H_line_min +\
                ((pos+0.5+(num_letters_for_buffer_width)/2.)*one_letter_width)
            if ( len( conc_keys_to_plot_L ) in [7, 8] ):
                fontsize = 12
            elif ( len( conc_keys_to_plot_L ) > 8 ):
                fontsize = 11
            else:
                fontsize= 14
            plt.annotate( conc.replace( "_", " " ),
                    (this_letter_center, min_v_line - 0.05), fontsize = fontsize,
                    xycoords='figure fraction',
                    verticalalignment = 'top', horizontalalignment = 'center')

        # make a title
        if ( title != "" ):
            plt.annotate( title,
                (middle_of_plot, enrich_max + 0.07),
                xycoords='figure fraction',
                fontsize= 16,
                verticalalignment = 'top', horizontalalignment = 'center')


        ##### Plot the significance over all 4 bins
        if plot_signif:
            #### pos is the index for the protien concentration
            for pos, conc in enumerate( conc_keys_to_plot_L ):
                this_letter_center = main_H_line_min + ((pos+0.5 +\
                    (num_letters_for_buffer_width)/2.)*one_letter_width)
                ### The avg. enrichment over all bins
                avg_enrich_all_bins = enrichments_by_kmer_conc_bin_D[kmer][conc]["overall"]

                abs_delta_last_minus_first = abs( L_of_Ls_to_plot[pos][-1] - L_of_Ls_to_plot[pos][0] )
                #### Divide this by the overall enrichment so the scales are the same
                normed_abs_delta_last_minus_first = abs_delta_last_minus_first / avg_enrich_all_bins

                std_this_motif_all_bins = sig_Ds_by_kmer_D[kmer]['0_to_4']['std']
                if ( normed_abs_delta_last_minus_first >= 2 * std_this_motif_all_bins ):

                    lft_of_conc = this_letter_center - ( 2.3 * in_bar_width )
                    right_of_conc = this_letter_center + ( 2.3 * in_bar_width )

                    ax_for_lines.axhline( y = min_v_line - 0.03,
                        xmin = lft_of_conc, xmax = right_of_conc,
                        color = 'k', linestyle = '-', linewidth = 1 )
                    ax_for_lines.axvline( x = lft_of_conc, ymin = min_v_line - 0.028, ymax = enrich_0,
                        color = 'k', linestyle = '-', linewidth = 1 )
                    ax_for_lines.axvline( x = right_of_conc, ymin = min_v_line - 0.028, ymax = enrich_0,
                        color = 'k', linestyle = '-', linewidth = 1 )
                    plt.annotate( "**", ( this_letter_center, min_v_line - 0.045 ), fontsize=12,
                            xycoords='figure fraction',
                            verticalalignment = 'center', horizontalalignment = 'center' )
                    sigB_by_kmer_conc_bin_D[kmer][conc]["0_to_4"] = True
                else:
                    sigB_by_kmer_conc_bin_D[kmer][conc]["0_to_4"] = False

        ######### < moved down to here >
        #### pos is the index for the protien concentration
        for pos, conc in enumerate( conc_keys_to_plot_L ):
            enrich_handles_L = 0
            handles_L = []
            this_letter_center = main_H_line_min + ((pos+0.5 +\
                (num_letters_for_buffer_width)/2.)*one_letter_width)
            ### The avg. enrichment over all bins
            avg_enrich_all_bins = enrichments_by_kmer_conc_bin_D[kmer][conc]["overall"]
            #### Go through the 5 bins for this protein concentration
            for Ppaired_bin_idx, enrich in enumerate( L_of_Ls_to_plot[pos] ):

                lft = this_letter_center + (((-0.5*num_Ppaired_bins)+Ppaired_bin_idx)*in_bar_width)
                btm = enrich_0
                wdth = in_bar_width
                hgh = enrich*one_enrich_unit_height
                axb = fig2.add_axes([lft, btm, wdth, hgh], zorder = 2 )
                axb.set_xlim(0, 1.)
                axb.set_ylim(0, 1.)
                axb.spines['top'].set_color('none')
                axb.spines['bottom'].set_color('none')
                axb.spines['left'].set_color('none')
                axb.spines['right'].set_color('none')
                #ax.tick_params(labelcolor="w", top='off',left="off",right="off",\
                        #        bottom="off")
                axb.tick_params(which = 'both', top='off',left="off",right="off",\
                    bottom="off", labelbottom = "off", labeltop = "off",
                    labelleft = "off", labelright= 'off')
                color_num_to_use = light_to_dark_D[len(Ppaired_bin_upper_limits_L)][Ppaired_bin_idx]
                hand = axb.bar(0., 1., width = 1.,
                        #color = shades_of_red_D[color_num_to_use], edgecolor = 'k')
                        color = shades_of_red_D[color_num_to_use], edgecolor = 'k')
                handles_L.append( hand )
                #### If it's NOT the first bin, see if the SIGNFICANCE should be plotted
                if ( Ppaired_bin_idx > 0 ) and plot_signif:
                    prev_bin_enrich = L_of_Ls_to_plot[pos][Ppaired_bin_idx-1]
                    #### Get the mean between these two bins so the difference can
                    ####    be normalized (since the 0 nM std. are around 1 )
                    mean_these_bins = math.sqrt( enrich * prev_bin_enrich )
                    norm_abs_delta_this_bin_vs_prev = abs( enrich - prev_bin_enrich ) / mean_these_bins
                    #{0: {'mean': -0.00471339254492531, 'std': 0.035055744570428218},
                    #       1: {'mean': 0.0025747329424797849, 'std': 0.027230896635803935},
                    #       2: {'mean': 0.023499924740518468, 'std': 0.037633649955424724},
                    #       3: {'mean': 0.044438761269242379, 'std': 0.12016792840035895},
                    #       '0_to_4': {'mean': 0.065800026407315318, 'std': 0.1106762538396405}}
                    std_this_motif_and_bin = sig_Ds_by_kmer_D[kmer][Ppaired_bin_idx-1]['std']
                    if ( norm_abs_delta_this_bin_vs_prev >= (2 * std_this_motif_and_bin) ):
                        #### Make a line
                        lower_x = lft - (wdth * 0.35)
                        upper_x = lft + (wdth * 0.35)
                        ax_for_lines.axhline( y = enrich_0 - 0.005,
                            xmin = lower_x, xmax = upper_x,
                            color = 'k', linestyle = '-', linewidth = 1 )
                        ax_for_lines.axvline( x = lower_x, ymin = enrich_0 - 0.005, ymax = enrich_0,
                            color = 'k', linestyle = '-', linewidth = 1 )
                        ax_for_lines.axvline( x = upper_x, ymin = enrich_0 - 0.005, ymax = enrich_0,
                            color = 'k', linestyle = '-', linewidth = 1 )
                        plt.annotate( "*", ( lft, enrich_0 - 0.02 ), fontsize=12,
                                xycoords='figure fraction',
                                verticalalignment = 'center', horizontalalignment = 'center' )
                        #### The significance bin should be indexed by the
                        ####    lower bin idx (so subtract 1)
                        sigB_by_kmer_conc_bin_D[kmer][conc][Ppaired_bin_idx-1] = True
                    else:
                        sigB_by_kmer_conc_bin_D[kmer][conc][Ppaired_bin_idx-1] = False


        if (include_legend == True):
            #### Put the legend on the right if the first concentration has
            ####    a higher average R than the last concentration
            if ( avg_Rs_L[0] > avg_Rs_L[-1] ):
                loc = (axis_x_left + axis_width - 0.21, max_v_line - 0.2)
            else:
                loc = (axis_x_left + 0.08, max_v_line - 0.2)
            leg = axbg.legend( handles_L, labels_to_plot_L,
                loc = loc,
                fancybox=True, fontsize=12, borderaxespad=0 )
            leg.get_frame().set_alpha(0.5)

        fig2.savefig( out_F )

    return_D = { 'sigB_by_kmer_conc_bin_D': sigB_by_kmer_conc_bin_D }

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




def merge_4_PDFs_on_1_page(
        PDFs_L,
        out_F ):
    """
    - Merges 4 PDFs (absolute paths to PDFs in PDFs_L into
        out_F

        PDF_1   PDF_2
        PDF_3   PDF_4

    """
    from PyPDF2 import PdfFileReader, PdfFileMerger, PdfFileWriter
    from PyPDF2.generic import RectangleObject
    #from pdfnup import generateNup # commented out by MA/KK

    output = PdfFileWriter()

    for F in PDFs_L[:4]:
        print F
        input = PdfFileReader( file( F , "rb") )
        img = input.getPage(0)
        output.addPage( img )

    outputStream = file( out_F , "wb")
    output.write(outputStream)
    outputStream.close()

    #generateNup( out_F, 4 ) # commented out by MA/KK




################################### </ UTILS > ################################
###############################################################################





