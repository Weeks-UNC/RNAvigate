
import sys
import sys
import os
sys.path.append('{}/arcPlot/'.format(os.environ['HOME']))
sys.path.append('{}/RNATools/'.format(os.environ['HOME']))

import arcPlot as ap
import RNAtools2 as RNAtools
from pmanalysis import PairMap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

##########################COPYPASTA FROM SHAPEMAPPER2
from numpy import isnan, nan, sqrt
from numpy import nanpercentile as percentile
from math import ceil

import matplotlib as mp
mp.use('Agg')
mp.rcParams["font.sans-serif"].insert(0,"Arial")
mp.rcParams["font.family"] = "sans-serif"
mp.rcParams["pdf.fonttype"] = 42 # use TrueType fonts when exporting PDFs
                                 # (embeds most fonts - this is especially
                                 #  useful when opening in Adobe Illustrator)
mp.rcParams['xtick.direction'] = 'out'
mp.rcParams['ytick.direction'] = 'out'
mp.rcParams['legend.fontsize'] = 14
mp.rcParams['grid.color'] = ".8"
mp.rcParams['grid.linestyle'] = '-'
mp.rcParams['grid.linewidth'] = 1
mp.use('Agg')

from matplotlib.patches import Rectangle

rx_color = "red"
bg_color = "blue"
dc_color = "darkgoldenrod"
#####################################################

sns.set_style("ticks")
sns.set_context("talk")


#
# 0 Nucleotide
# 1 Sequence
# 2 Modified_mutations
# 3 Modified_read_depth
# 4 Modified_effective_depth
# 5 Modified_rate
# 6 Modified_off_target_mapped_depth
# 7 Modified_low_mapq_mapped_depth
# 8 Modified_primer_pair_1_mapped_depth
# 9 Untreated_mutations
# 10 Untreated_read_depth
# 11 Untreated_effective_depth
# 12 Untreated_rate
# 13 Untreated_off_target_mapped_depth
# 14 Untreated_low_mapq_mapped_depth
# 15 Untreated_primer_pair_1_mapped_depth
# 16 Denatured_mutations
# 17 Denatured_read_depth
# 18 Denatured_effective_depth
# 19 Denatured_rate
# 20 Denatured_off_target_mapped_depth
# 21 Denatured_low_mapq_mapped_depth
# 22 Denatured_primer_pair_1_mapped_depth
# 23 Reactivity_profile
# 24 Std_err
# 25 HQ_profile
# 26 HQ_stderr
# 27 Norm_profile
# 28 Norm_stderr

def addSeqBar(axis, profile, yvalue=-0.01):
    color_dict = {"A": "#f20000","U": "#f28f00","G": "#00509d","C": "#00c200"}
    for n in ['A','U','C','G']:
        mask = profile['Sequence']==n
        x = profile['Nucleotide'][mask_18S]
        axis.scatter(x, [-0.01 for x in x1], marker='.', c=color_dict[n])

def plotReactivities(axis, sample, label, column='Reactivity_profile'):
    axis.plot(sample['Nucleotide'])


def metric_abbreviate(num):
    suffixes = {3:'k',
                6:'M',
                9:"G"}
    s = str(num)
    # replace trailing zeros with metric abbreviation
    zero_count = len(s)-len(s.rstrip('0'))
    suffix = ''
    new_string = str(s)
    for num_zeros in sorted(suffixes.keys()):
        if num_zeros <= zero_count:
            suffix = suffixes[num_zeros]
            new_string = s[:-num_zeros]
    new_string = new_string+suffix
    return new_string


def render_profiles(sample, name, qc_pass):
    num = sample['Nucleotide']
    seq = sample['Sequence']
    reactivity = sample['Norm_profile']
    stderr = sample['Norm_stderr']
    rx_rates = sample['Modified_rate']
    bg_rates = sample['Untreated_rate']
    dc_rates = sample['Denatured_rate']
    rx_depth = sample['Modified_effective_depth']
    bg_depth = sample['Untreated_effective_depth']
    dc_depth = sample['Denatured_effective_depth']
    rx_simple_depth = sample['Modified_read_depth']
    bg_simple_depth = sample['Untreated_read_depth']
    dc_simple_depth = sample['Denatured_read_depth']
    # FIXME: break each panel into a separate func
    if rx_rates is None and bg_rates is None and dc_rates is None:
        no_mapped = True
    else:
        no_mapped = False

    legend_labels = []
    if rx_rates is not None:
        legend_labels.append("Modified")
        rx_err = np.sqrt(rx_rates) / np.sqrt(rx_depth)
    else:
        rx_rates = np.zeros((len(rx_depth),))
        rx_err = np.zeros((len(rx_depth),))
    if bg_rates is not None:
        legend_labels.append("Untreated")
        bg_err = np.sqrt(bg_rates) / np.sqrt(bg_depth)
    else:
        bg_rates = np.zeros((len(rx_depth),))
        bg_err = np.zeros((len(rx_depth),))
    if dc_rates is not None:
        legend_labels.append("Denatured")
        dc_err = np.sqrt(dc_rates) / np.sqrt(dc_depth)
    else:
        dc_rates = np.zeros((len(rx_depth),))
        dc_err = np.zeros((len(rx_depth),))

    # Add a zeroeth nuc so axis numbering works correctly
    # There's probably a better way to do this
    num = np.append(0, num)
    if reactivity is not None:
        reactivity = np.append(0, reactivity)
        stderr = np.append(0, stderr)
    rx_depth = np.append(0, rx_depth)
    bg_depth = np.append(0, bg_depth)
    dc_depth = np.append(0, dc_depth)
    rx_simple_depth = np.append(0, rx_simple_depth)
    bg_simple_depth = np.append(0, bg_simple_depth)
    dc_simple_depth = np.append(0, dc_simple_depth)
    rx_rates = np.append(0, rx_rates)
    bg_rates = np.append(0, bg_rates)
    dc_rates = np.append(0, dc_rates)
    rx_err = np.append(0, rx_err)
    bg_err = np.append(0, bg_err)
    dc_err = np.append(0, dc_err)

    if reactivity is not None:
        orange_thresh = 0.4
        red_thresh = 0.85

        gray_vals = []
        gray_nums = []
        gray_errs = []
        black_vals = []
        black_nums = []
        black_errs = []
        orange_vals = []
        orange_nums = []
        orange_errs = []
        red_vals = []
        red_nums = []
        red_errs = []
        for i in range(len(reactivity)):
            if np.isnan(reactivity[i]):
                gray_vals.append(-1)
                gray_nums.append(num[i])
                gray_errs.append(0)
            elif reactivity[i] < orange_thresh:
                black_vals.append(reactivity[i])
                black_nums.append(num[i])
                black_errs.append(stderr[i])
            elif reactivity[i] < red_thresh:
                orange_vals.append(reactivity[i])
                orange_nums.append(num[i])
                orange_errs.append(stderr[i])
            else:
                red_vals.append(reactivity[i])
                red_nums.append(num[i])
                red_errs.append(stderr[i])

    yMin, ymax = (-0.5, 4)
    left_inches = 0.9
    right_inches = 0.4
    sp_width = len(num)*0.032
    fig_width = max(7,sp_width+left_inches+right_inches)
    fig = plt.figure(figsize=(fig_width,8))

    left_percent = left_inches/fig_width
    right_percent = 1-right_inches/fig_width
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(313)
    ax3 = plt.subplot(312)
    plt.subplots_adjust(hspace=0.5, left=left_percent,right=right_percent,top=0.94)

    near_black = (0,0,1/255.0)

    if reactivity is not None:
        if len(gray_nums)>0:
            ax1.bar(gray_nums,gray_vals,
                     align="center",
                     width=1.05, color="0.80", edgecolor="0.80",linewidth=0.0)
        if len(black_nums)>0:
            ax1.bar(black_nums,black_vals,
                     align="center",
                     width=1.05, color="black", edgecolor="black",linewidth=0.0,
                     yerr=black_errs,ecolor=near_black,capsize=1)
        if len(orange_nums)>0:
            ax1.bar(orange_nums,orange_vals,
                     align="center",
                     width=1.05, color="orange",edgecolor="orange",linewidth=0.0,
                     yerr=orange_errs,ecolor=near_black,capsize=1)
        if len(red_nums)>0:
            ax1.bar(red_nums,red_vals,
                     align="center",
                     width=1.05,color="red",edgecolor="red",linewidth=0.0,
                     yerr=red_errs,ecolor=near_black,capsize=1)

    #print("title: "+name)
    ax1title = ax1.set_title(name, horizontalalignment="left", fontsize=16)
    x,y = ax1title.get_position()
    ax1title.set_position((0,y))
    ax1.set_ylim(yMin,ymax)
    ax1.set_xlim(1,len(num))
    #ax1.set_yticks(fontsize=9)

    if not qc_pass:
        msg = "Note: possible data quality issue - see log file"
        return_flag = False
        if no_mapped:
            msg = "ERROR: no reads mapped to this RNA"
            return_flag = True
        txt = plt.text(30,573,
                 msg,
                 ha='left', va='top',
                 fontsize=16, color='red',
                 transform=mp.transforms.IdentityTransform())
        if return_flag:
            plt.show()
            return

    #tickNums = range(num[0]+10,num[-1]+1,10)
    #tickPos = range(num[0]+9,num[-1],10)
    #ax1.set_xticks(tickPos,tickNums,fontsize=9,rotation=30)
    #ax1.set_xticks(fontsize=9)

    ax1.yaxis.grid(True)
    ax1.set_axisbelow(True)

    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    for loc, spine in ax1.spines.items():
        if loc == 'bottom':
            spine.set_position(('outward', 6))  # move outward (down) 6 pts
            spine.set_smart_bounds(True)
    for loc, spine in ax1.spines.items():
        if loc == 'left':
            spine.set_position(('outward', 6))  # move outward (left) 6 pts
            spine.set_smart_bounds(True)

    # need to add labels after moving spines, otherwise they will disappear
    ax1xlabel = ax1.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax1xlabel.get_position()
    ax1xlabel.set_position((0,y))
    ax1ylabel = ax1.set_ylabel("Shape Reactivity", horizontalalignment="left", fontsize=14)
    x,y = ax1ylabel.get_position()
    ax1ylabel.set_position((x,0))

    if reactivity is not None:
        # add a SHAPE colorbar to the vertical axis
        # uses a little transformation magic to place correctly
        inv = ax1.transData.inverted()
        for loc, spine in ax1.spines.items():
            if loc == 'left':
                trans = spine.get_transform()
        pt = trans.transform_point([0,0])
        pt2 = inv.transform_point(pt)
        rectX = pt2[0]
        ptA = (0,0)
        ptB = (6,0)
        ptA2 = inv.transform_point(ptA)
        ptB2 = inv.transform_point(ptB)
        rectW = ptB2[0]-ptA2[0]
        rect = Rectangle((rectX,-0.5), rectW, orange_thresh+0.5, facecolor="black", edgecolor="none")
        ax1.add_patch(rect)
        rect.set_clip_on(False)
        rect = Rectangle((rectX,orange_thresh), rectW, red_thresh-orange_thresh, facecolor="orange", edgecolor="none")
        ax1.add_patch(rect)
        rect.set_clip_on(False)
        rect = Rectangle((rectX,red_thresh), rectW, 4-red_thresh, facecolor="red", edgecolor="none")
        ax1.add_patch(rect)
        rect.set_clip_on(False)

    ax1.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax1.get_yaxis().tick_left()

    ax1.tick_params(axis='y',which='minor',left=False)
    #ax1.tick_params(axis='x',which='minor')

    ax1.minorticks_on()

    yticks = ax1.get_yticks()
    stripped_ticks = [str(val).rstrip('0').rstrip('.') for val in yticks]
    ax1.set_yticklabels(stripped_ticks)

    for line in ax1.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax1.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax1.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)


    # put nuc sequence below axis
    font_prop = mp.font_manager.FontProperties(family = "monospace", style="normal",weight="bold",size="4.5")
    for i in range(seq.shape[0]):
        nuc = seq[i]
        if nuc == "T":
            nuc = "U"
        color_dict = {"A": "#f20000",
                      "U": "#f28f00",
                      "G": "#00509d",
                      "C": "#00c200"}
        if nuc in color_dict:
            col = color_dict[nuc]
        elif nuc.upper() in color_dict:
            col = color_dict[nuc.upper()]
        else:
            col = "black"
        ax1.annotate(nuc, xy=(i+1, -0.67),fontproperties=font_prop,color=col,annotation_clip=False, horizontalalignment="center")

    handles = []
    h, = ax2.plot(num, rx_simple_depth, linewidth = 1.5, color=rx_color, alpha=1.0)
    ax2.plot(num, rx_depth, linewidth = 1.0, color=rx_color, alpha=0.3)
    handles.append(h)
    h, = ax2.plot(num, bg_simple_depth, linewidth = 1.5, color=bg_color, alpha=1.0)
    ax2.plot(num, bg_depth, linewidth = 1.0, color=bg_color, alpha=0.3)
    handles.append(h)
    h, = ax2.plot(num, dc_simple_depth, linewidth = 1.5, color=dc_color, alpha=1.0)
    ax2.plot(num, dc_depth, linewidth = 1.0, color=dc_color, alpha=0.3)
    handles.append(h)
    ax2.set_xlim(1,len(num))
    #ax2.legend(["+Reagent","Background","Denatured"], bbox_to_anchor=(1.1,1.1))
    leg = ax2.legend(handles, legend_labels, loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    xmin, xmax, ymin, ymax = ax2.axis()
    ax2.set_ylim(0,ymax)
    #ax2.set_yscale('log')
    #ax2.set_yscale('symlog')# useful, but disabled because of a matplotlib/pyparsing bug
    ax2xlabel = ax2.set_xlabel("Nucleotide\n(note: effective read depths shown in lighter colors)", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax2xlabel.get_position()
    ax2xlabel.set_position((0,y))

    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax2.get_yaxis().tick_left()

    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='minor',left=False)
    #ax2.tick_params(axis='x',which='minor')

    #xlabels = ["%.2f"%v for v in xticks]
    #ax3.set_xticks(xticks)
    #ax3.set_xticklabels(xlabels,rotation = -45, horizontalalignment='left')

    yticks = [int(y) for y in ax2.get_yticks()]
    formatted_ticks = []
    for val in yticks:
        formatted_ticks.append(metric_abbreviate(val))
    ax2.set_yticklabels(formatted_ticks)

    for line in ax2.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax2.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax2.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax2.yaxis.grid(True)
    ax2.set_axisbelow(True)

    ax2ylabel = ax2.set_ylabel("Read depth", horizontalalignment="left", fontsize=14)
    x, y = ax2ylabel.get_position()
    ax2ylabel.set_position((x, 0))
    # tried to make an offset, smaller font note about effective depths,
    # but couldn't get positioning/transforms to work properly.
    # For now just putting in xaxis label

    # choose a decent range for axis, excluding high-background positions
    good_indices = []
    for i in range(len(bg_rates)):
        if bg_rates[i]<=0.05 or np.isnan(bg_rates[i]):
            good_indices.append(i)
    temp_rates = [rx_rates[i] for i in good_indices]
    near_top_rate = percentile(temp_rates,98.0)
    maxes = [0.32,0.16,0.08,0.04,0.02,0.01]
    ymax = maxes[0]
    for i in range(len(maxes)):
        if near_top_rate<maxes[i]:
            ymax = maxes[i]

    rx_upper = rx_rates + rx_err
    rx_lower = rx_rates - rx_err
    bg_upper = bg_rates + bg_err
    bg_lower = bg_rates - bg_err
    dc_upper = dc_rates + dc_err
    dc_lower = dc_rates - dc_err

    ax3xlabel = ax3.set_xlabel("Nucleotide", horizontalalignment="left", fontsize=14, labelpad=0)
    x,y = ax3xlabel.get_position()
    ax3xlabel.set_position((0,y))
    ax3ylabel = ax3.set_ylabel("Mutation rate (%)", horizontalalignment="left", fontsize=14)
    x,y = ax3ylabel.get_position()
    ax3ylabel.set_position((x,0))

    ax3.plot(num, rx_rates, zorder=3, color=rx_color, linewidth=1.5)
    ax3.plot(num, bg_rates, zorder=2, color=bg_color, linewidth=1.5)
    ax3.plot(num, dc_rates, zorder=2, color=dc_color, linewidth=1.5)
    ax3.fill_between(num, rx_lower, rx_upper, edgecolor="none", alpha=0.5, facecolor=rx_color)
    ax3.fill_between(num, bg_lower, bg_upper, edgecolor="none", alpha=0.5, facecolor=bg_color)
    ax3.fill_between(num, dc_lower, dc_upper, edgecolor="none", alpha=0.5, facecolor=dc_color)
    ax3.legend(legend_labels, loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    ax3.set_xlim((1,len(rx_rates)))
    ax3.set_ylim((0,ymax))

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax3.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax3.get_yaxis().tick_left()

    ax3.minorticks_on()
    ax3.tick_params(axis='y',which='minor',left=False)

    ticks = [x*100 for x in ax3.get_yticks()]
    ax3.set_yticklabels([str(val).rstrip('0').rstrip('.') for val in ticks])

    for line in ax3.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in ax3.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in ax3.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    ax3.yaxis.grid(True)
    ax3.set_axisbelow(True)

    # TODO: add a tick for the first nuc - can't seem to add one without screwing
    # up all the other ticks

    plt.show()


def arcPlot_pairmap(title, fa, ct, shapefile, pm):
    ctStructNum = 0
    filterNC = True
    pairmapAll = False
    aplot = ap.ArcPlot(title=title, fasta=fa)
    aplot.addCT(RNAtools.CT(ct, structNum=ctStructNum, filterNC=filterNC),
                panel=1)
    aplot.addPairMap(PairMap(pm), panel=-1, plotall=pairmapAll)
    aplot.readProfile(shapefile)
    aplot.writePlot('show')


def plotIgnoredCorrs(title, file):
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    x, y = [], []
    with open(file) as f:
        for line in f.readlines():
            if line.startswith('Pair '):
                pair = line.split()[1].strip('()').split(',')
                x.extend(pair)
                y.extend(pair[::-1])
    ax.scatter(x, y, marker='s', s=1)
    ax.set(title=title+": Background Correlations")
    return x, y
