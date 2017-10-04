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







def make_density_plots_of_multiple_histvaluesLs(
        xs_Ls_L,
        captions_L,
        #density_covar = 0.01,
        out_F = None,
        include_num_points_in_caption = True,
        lower_x_limit = 0.,
        upper_x_limit = "max",
        upper_y_limit = "max",
        xs_vert_lines_L = [0.],
        title = "",
        x_label = "",
        y_label = "",
        legend_location = 1, ### upper right
        title_fontsize = 18,
        label_fontsize = 16,
        legend_fontsize = 14,
        density_covar = 0.3,
        colors = "reds_light_to_dark"):
    """
    - For each list in xs_Ls_L, a Gaussian kernel density plot will be made

    - lower_x_limit, upper_x_limit; and upper_y_limit:
        - "max"
            1.1 * the max of any value in xs_Ls_L
        - "max_except_first_L"
            The maximum value in any list EXCEPT the first list (e.g., if the
                first list is 0 nM enrichments, but these should be ignored)
        - Otherwise, pass in a float directly
    """
    num_groups = len( xs_Ls_L )
    assert( len( captions_L ) == num_groups )

    if (colors == "reds_dark_to_light"):
        #### Get the list of hexadecimal colors
        colors_to_use_L = return_shades_of_color_L(
                num_groups,
                color = "red" )
    elif (colors == "reds_light_to_dark"):
        colors_to_use_L = return_shades_of_color_L(
                num_groups,
                color = "red",
                dark_to_light = "light_to_dark")
    elif (colors == "blues_dark_to_light"):
        colors_to_use_L = return_blues_colors_L(
                num_groups,
                color = "blue" )
    elif (colors == "blues_light_to_dark"):
        colors_to_use_L = return_reds_colors_L(
                num_groups,
                color = "blue",
                dark_to_light = "light_to_dark")
    elif (type( colors ) is list):
        assert( len( colors ) == num_groups )
        colors_to_use_L = colors

    ##### First, get the max and min x values
    min_x = lower_x_limit
    all_xs_L = []
    for xs_L in xs_Ls_L:
        all_xs_L += xs_L
    if (upper_x_limit == "max"):
        max_x = 1.2 * max( all_xs_L )
    else:
        max_x = upper_x_limit
    xs = np.arange( min_x , max_x, (max_x - min_x)/10001. )

    ##### plot the data
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #for vert_line in xs_vert_lines_L:
    #    plt.axvline( x = vert_line, linewidth=2, linestyle="--", color='#808080')

    for idx, vert_line in enumerate( xs_vert_lines_L ):
        #plt.axvline( x = vert_line, linewidth=2, linestyle="--", color=colors[idx])
        plt.axvline( x = vert_line, linewidth=2, linestyle="--", color = '#808080' )

    ##### turn each of the xs_L into a density
    if include_num_points_in_caption:
        mod_captions_L = []
    legend_handles_L = []
    max_y = 0.
    for line_num, xs_L in enumerate( xs_Ls_L ):

        density = gaussian_kde( xs_L )

        #### COVARIANCE
        density.covariance_factor = lambda : density_covar
        density._compute_covariance()

        #### fill in the y values using density class
        ys = density( xs )
        if (upper_y_limit == "max_except_first_L") and (line_num > 0):
            max_y = max( max_y, max( ys ) )
        elif (upper_y_limit == "max"):
            max_y = max( max_y, max( ys ) )
        hand_L = ax.plot( xs, ys, color = colors_to_use_L[line_num],
                linewidth = 3)
        legend_handles_L.append( hand_L[0] )

        #### If desired, add the number of points to the captions
        if include_num_points_in_caption:
            mod_captions_L.append( captions_L[line_num] + " ($n={0:,}$)".format( len(xs_L) ) )

    if include_num_points_in_caption:
        captions_to_use_L = mod_captions_L
    else:
        captions_to_use_L = captions_L

    #### Make a legend
    num_nonzero = len( [x for x in captions_L if x != ""] )
    if ( num_nonzero > 0 ):
        leg = ax.legend(
                legend_handles_L,
                captions_to_use_L,
                loc = legend_location,
                fancybox = True,
                fontsize = legend_fontsize,
                borderaxespad = 0)
        leg.get_frame().set_alpha(0.5)

    #### Title and Axis labels
    ax.set_title( title, fontsize = title_fontsize )
    ax.set_xlabel( x_label, fontsize = label_fontsize )
    ax.set_ylabel( y_label, fontsize = label_fontsize )
    ax.tick_params( axis='x', labelsize = label_fontsize )
    ax.tick_params( axis='y', labelsize = label_fontsize )
    ax.set_xlim( xs[0], xs[-1] )
    ax.set_ylim( 0, 1.1 * max_y )

    return_D = {"fig": fig}
    if (out_F != None):
        fig.savefig( out_F )
        return_D["out_F"] = out_F
    return return_D




def make_plot_of_Zscores_from_enrichments_txt_F(
        txt_F ):
    """
    - For a particular enrichments txt_F from the RBNS pipeline (e.g., the
        Z-score = 2 or 3 table),
        makes a heatmap of, for each kmer, its Z-scores over the concentrations

    - Meant to visually help discriminated "contaminated" or other odd kmers
        within an enriched list, to see that they have a different profile over
        the different concentrations vs. other kmers

    6/12/16
    """
    import numpy as np
    out_DIR = os.path.join( txt_F.split("/tables")[0], "plots/Z_score_over_concs" )
    os.system( "mkdir -p {}".format( out_DIR ) )

    out_basename = os.path.basename( txt_F ).split( ".txt" )[0] +\
            ".Zscore_of_R_over_concs.pdf"

    out_same_order_F = os.path.join( out_DIR, out_basename.split(".pdf")[0] +\
            ".same_order_as_top.txt" )
    out_diff_order_F = os.path.join( out_DIR, out_basename.split(".pdf")[0] +\
            ".diff_order_than_top.txt" )
    same_order_f = open( out_same_order_F, "w" )
    diff_order_f = open( out_diff_order_F, "w" )

    #### First get the concentrations from the header line
    ####    [SFPQ]  5 nM    20 nM   80 nM   320 nM  1300 nM
    with open( txt_F ) as f:
        for line in f:
            header_line_L = line.strip().split('\t')
            break

    protein_w_brackets = header_line_L[0]
    protein = protein_w_brackets[1:-1]
    concs_L = header_line_L[1:]

    checkerboard_D = {}
    for conc in concs_L:
        checkerboard_D[conc] = {}

    with open( txt_F ) as f:
        next( f )
        for line in f:
            line_L = line.strip().split('\t')

            Rs_L = [float(x) for x in line_L[1:]]
            Rs_concs_T_L = zip( Rs_L, concs_L )
            Rs_concs_T_L.sort( key = lambda x: -1 * x[0] )
            ordered_concs_L = [x[1] for x in Rs_concs_T_L]
            break

    #### Now go through all the (non-header) lines
    kmers_order_L = []
    with open( txt_F ) as f:
        next( f )
        for line in f:
            line_L = line.strip().split('\t')

            kmer = line_L[0].replace( "T", "U" )
            kmers_order_L.append( kmer )
            Z_scores_L = [float(x) for x in line_L[1:]]
            mean_Zscore = np.mean( Z_scores_L )
            std_Zscore = np.std( Z_scores_L )

            this_Rs_concs_T_L = zip( Z_scores_L, concs_L )
            this_Rs_concs_T_L.sort( key = lambda x: -1 * x[0] )
            this_ordered_concs_L = [x[1] for x in this_Rs_concs_T_L]
            if ( ordered_concs_L == this_ordered_concs_L ):
                same_order_f.write( line )
            else:
                diff_order_f.write( line )

            for conc_idx, conc in enumerate( concs_L ):
                checkerboard_D[conc][kmer] = (Z_scores_L[conc_idx] - mean_Zscore) / std_Zscore

    if ( len( kmers_order_L ) < 20 ):
        label_2_fontsize = 12
    elif ( len( kmers_order_L ) < 40 ):
        label_2_fontsize = 9
    else:
        label_2_fontsize = 7

    #### Make the plot
    out_F = os.path.join( out_DIR, out_basename )
    make_rectangular_heatmap_plot(
            checkerboard_D,
            label_key_1 = protein_w_brackets,
            label_key_2 = "Top {}mers".format( len( kmer ) ),
            order_of_keys_1_L = concs_L,
            order_of_keys_2_L = kmers_order_L[::-1],
            colorbar_label = r"Z-score of $R$ over {0} concentrations".format( len(concs_L) ),
            colormap_center = 0.,
            label_2_fontsize = label_2_fontsize,
            title = "{} Z-scores over concentrations".format( protein ),
            out_F = out_F )

    same_order_f.close()
    diff_order_f.close()



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








