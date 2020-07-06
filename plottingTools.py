
from matplotlib.patches import Rectangle
import matplotlib as mp
from numpy import nanpercentile as percentile
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# COPYPASTA FROM SHAPEMAPPER2

mp.use('Agg')
mp.rcParams["font.sans-serif"].insert(0, "Arial")
mp.rcParams["font.family"] = "sans-serif"
mp.rcParams["pdf.fonttype"] = 42  # use TrueType fonts when exporting PDFs
# (embeds most fonts - this is especially
#  useful when opening in Adobe Illustrator)
mp.rcParams['xtick.direction'] = 'out'
mp.rcParams['ytick.direction'] = 'out'
mp.rcParams['legend.fontsize'] = 14
mp.rcParams['grid.color'] = ".8"
mp.rcParams['grid.linestyle'] = '-'
mp.rcParams['grid.linewidth'] = 1
mp.use('Agg')


rx_color = "red"
bg_color = "blue"
dc_color = "darkgoldenrod"
#####################################################

sns.set_style("ticks")
sns.set_context("talk")


# profile.txt header names and indices
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


def plotSkyline(axis, profile, label=None, column='Reactivity_profile'):
    x = []
    y = []
    for n, r in zip(profile['Nucleotide'], profile[column]):
        x.extend([n - 0.5, n + 0.5])
        y.extend([r, r])
    axis.plot(x, y, label=label)


def addSeqBar(axis, profile, yvalue=-0.017):
    font_prop = mp.font_manager.FontProperties(
        family="monospace", style="normal", weight="bold", size="12")
    color_dict = {"A": "#f20000", "U": "#f28f00",
                  "G": "#00509d", "C": "#00c200"}
    for i, row in profile.iterrows():
        col = color_dict[row['Sequence'].upper()]
        axis.annotate(row['Sequence'], xy=(i + 1, yvalue),
                      fontproperties=font_prop,
                      color=col, annotation_clip=False,
                      horizontalalignment="center")


def metric_abbreviate(num):
    suffixes = {3: 'k',
                6: 'M',
                9: "G"}
    s = str(num)
    # replace trailing zeros with metric abbreviation
    zero_count = len(s) - len(s.rstrip('0'))
    suffix = ''
    new_string = str(s)
    for num_zeros in sorted(suffixes.keys()):
        if num_zeros <= zero_count:
            suffix = suffixes[num_zeros]
            new_string = s[:-num_zeros]
    new_string = new_string + suffix
    return new_string


def getWidth(sample):
    left_inches = 0.9
    right_inches = 0.4
    sp_width = len(sample['Nucleotide']) * 0.1
    fig_width = max(7, sp_width + left_inches + right_inches)
    return fig_width


def plotProfile(axis, sample, name):
    yMin, ymax = (-0.5, 4)
    near_black = (0, 0, 1 / 255.0)
    orange_thresh = 0.4
    red_thresh = 0.85
    cindex = np.zeros(len(sample['Norm_profile']), dtype=int)
    cindex[np.logical_not(np.isnan(sample['Norm_profile']))] += 1
    cindex[sample["Norm_profile"] > orange_thresh] += 1
    cindex[sample['Norm_profile'] > red_thresh] += 1
    colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
    profile = sample['Norm_profile'].copy()
    profile[np.isnan(profile)] = -1
    axis.bar(sample['Nucleotide'], profile, align="center",
             width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
             yerr=sample['Norm_stderr'], ecolor=near_black, capsize=1)

    axis.set_title(name, fontsize=16)
    axis.set_ylim(yMin, ymax)
    axis.set_xlim(1, len(sample['Nucleotide']))

    axis.yaxis.grid(True)
    axis.set_axisbelow(True)
    axis.set_xlabel("Nucleotide", fontsize=14, labelpad=0)
    axis.set_ylabel("Shape Reactivity", fontsize=14)

    # add a SHAPE colorbar to the vertical axis
    # uses a little transformation magic to place correctly
    inv = axis.transData.inverted()
    for loc, spine in axis.spines.items():
        if loc == 'left':
            trans = spine.get_transform()
    pt = trans.transform_point([0, 0])
    pt2 = inv.transform_point(pt)
    rectX = pt2[0]
    ptA = (0, 0)
    ptB = (6, 0)
    ptA2 = inv.transform_point(ptA)
    ptB2 = inv.transform_point(ptB)
    rectW = ptB2[0] - ptA2[0]
    rect = Rectangle((rectX, -0.5), rectW, orange_thresh +
                     0.5, facecolor="black", edgecolor="none")
    axis.add_patch(rect)
    rect.set_clip_on(False)
    rect = Rectangle((rectX, orange_thresh), rectW, red_thresh -
                     orange_thresh, facecolor="orange", edgecolor="none")
    axis.add_patch(rect)
    rect.set_clip_on(False)
    rect = Rectangle((rectX, red_thresh), rectW, 4 -
                     red_thresh, facecolor="red", edgecolor="none")
    axis.add_patch(rect)
    rect.set_clip_on(False)

    axis.get_xaxis().tick_bottom()   # remove unneeded ticks
    axis.get_yaxis().tick_left()

    axis.tick_params(axis='y', which='minor', left=False)

    axis.minorticks_on()

    yticks = axis.get_yticks()
    stripped_ticks = [str(int(val)) for val in yticks]
    axis.set_yticklabels(stripped_ticks)

    for line in axis.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in axis.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in axis.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    addSeqBar(axis, sample, yvalue=-0.67)  # put nuc sequence below axis


def plotDepth(axis, sample):
    axis.plot(sample['Nucleotide'], sample['Modified_read_depth'],
              linewidth=1.5, color=rx_color, alpha=1.0, label="Modified")
    axis.plot(sample['Nucleotide'], sample['Untreated_read_depth'],
              linewidth=1.5, color=bg_color, alpha=1.0, label="Untreated")
    axis.plot(sample['Nucleotide'], sample['Denatured_read_depth'],
              linewidth=1.5, color=dc_color, alpha=1.0, label="Denatured")
    axis.set_xlim(1, len(sample['Nucleotide']))
    axis.legend(loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    axis.plot(sample['Nucleotide'], sample['Modified_effective_depth'],
              linewidth=1.0, color=rx_color, alpha=0.3)
    axis.plot(sample['Nucleotide'], sample['Untreated_effective_depth'],
              linewidth=1.0, color=bg_color, alpha=0.3)
    axis.plot(sample['Nucleotide'], sample['Denatured_effective_depth'],
              linewidth=1.0, color=dc_color, alpha=0.3)

    xmin, xmax, ymin, ymax = axis.axis()
    axis.set_ylim(0, ymax)
    axis.set_xlabel("Nucleotide\n" +
                    "(note: effective read depths shown in lighter colors)",
                    fontsize=14, labelpad=0)

    axis.minorticks_on()
    axis.tick_params(axis='y', which='minor', left=False)

    yticks = [int(y) for y in axis.get_yticks()]
    formatted_ticks = []
    for val in yticks:
        formatted_ticks.append(metric_abbreviate(val))
    axis.set_yticklabels(formatted_ticks)

    for line in axis.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in axis.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in axis.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    axis.yaxis.grid(True)
    axis.set_axisbelow(True)

    axis.set_ylabel("Read depth", fontsize=14)
    # tried to make an offset, smaller font note about effective depths,
    # but couldn't get positioning/transforms to work properly.
    # For now just putting in xaxis label


def plotMutationRates(axis, sample):
    # choose a decent range for axis, excluding high-background positions
    temp_rates = sample['Modified_rate'][sample['Untreated_rate'] <= 0.05]
    near_top_rate = percentile(temp_rates, 98.0)
    maxes = np.array([0.32, 0.16, 0.08, 0.04, 0.02, 0.01])
    ymax = np.amin(maxes[maxes > near_top_rate])

    rx_err = sample['Modified_rate'] / sample['Modified_effective_depth']
    rx_upper = sample['Modified_rate'] + rx_err
    rx_lower = sample['Modified_rate'] - rx_err
    bg_err = sample['Untreated_rate'] / sample['Untreated_effective_depth']
    bg_upper = sample['Untreated_rate'] + bg_err
    bg_lower = sample['Untreated_rate'] - bg_err
    dc_err = sample['Denatured_rate'] / sample['Denatured_effective_depth']
    dc_upper = sample['Denatured_rate'] + dc_err
    dc_lower = sample['Denatured_rate'] - dc_err

    axis.set_xlabel("Nucleotide", fontsize=14)
    axis.set_ylabel("Mutation rate (%)", fontsize=14)

    axis.plot(sample['Nucleotide'], sample['Modified_rate'], zorder=3,
              color=rx_color, linewidth=1.5, label='Modified')
    axis.plot(sample['Nucleotide'], sample['Untreated_rate'], zorder=2,
              color=bg_color, linewidth=1.5, label='Untreated')
    axis.plot(sample['Nucleotide'], sample['Denatured_rate'], zorder=2,
              color=dc_color, linewidth=1.5)
    axis.fill_between(sample['Nucleotide'], rx_lower, rx_upper,
                      edgecolor="none", alpha=0.5, facecolor=rx_color)
    axis.fill_between(sample['Nucleotide'], bg_lower, bg_upper,
                      edgecolor="none", alpha=0.5, facecolor=bg_color)
    axis.fill_between(sample['Nucleotide'], dc_lower, dc_upper,
                      edgecolor="none", alpha=0.5, facecolor=dc_color)
    axis.legend(loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
    axis.set_xlim((1, len(sample['Modified_rate'])))
    axis.set_ylim((0, ymax))

    axis.minorticks_on()
    axis.tick_params(axis='y', which='minor', left=False)

    ticks = [x * 100 for x in axis.get_yticks()]
    axis.set_yticklabels([str(int(val)) for val in ticks])

    for line in axis.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(1)

    for line in axis.get_xticklines():
        line.set_markersize(7)
        line.set_markeredgewidth(2)

    for line in axis.xaxis.get_ticklines(minor=True):
        line.set_markersize(5)
        line.set_markeredgewidth(1)

    axis.yaxis.grid(True)
    axis.set_axisbelow(True)


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
    ax.set(title=title + ": Background Correlations")
    return x, y
