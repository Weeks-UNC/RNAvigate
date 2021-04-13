#!/usr/bin/env python
import matplotlib as mp
import seaborn as sns
import pandas as pd
from scipy import stats

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


def plotRegression(ax, p1, p2, ctfile=None):
    """Plots scatter plot of reactivity profile vs. reactivity profile and
    computes regression metrics R^2 and slope, which are annotated on the axis.
    If a ctfile is provided, paired and unpaired nucleotides will have distinct
    colors.

    Args:
        ax (pyplot axis): axis on which scatter plot appears
        p1 (numpy array): array containing reactivity profile values
        p2 (numpy array): second array containing reactivity profile values
        ctfile (str, optional): path to ct file. Defaults to 'None'.
    """
    ax.plot([0, 1], [0, 1], color='black')
    gradient, _, r_value, _, _ = stats.linregress(p1, p2)
    ax.text(0.1, 0.8,
            'R^2 = {:.2f}\nslope = {:.2f}'.format(r_value**2, gradient),
            transform=ax.transAxes)
    if ctfile is not None:
        ct = pd.read_csv(ctfile, sep='\s+', usecols=[4], names=['j'], header=1)
        paired = ct.j != 0
        unpaired = ct.j == 0
        ax.scatter(p1[paired], p2[paired], label="Paired")
        ax.scatter(p1[unpaired], p2[unpaired], label="Unpaired")
    else:
        ax.scatter(p1, p2)


def readHistograms(logfile):
    """Parses the shapemapper_log.txt file for read length and mutations per
    molecule histograms. Requires '--per-read-histograms' flag when running
    shapemapper

    Args:
        logfile (string): path to shapemapper_log.txt file

    Returns:
        dataFrame: a dataFrame containing [Read_length, Mutation_count,
            Modified_read_length, Modified_mutations_per_molecule,
            Untreated_read_length, Untreated_mutations_per_molecule]
    """
    with open(logfile, 'r') as f:
        flist = list(f)
        for i, line in enumerate(flist):
            if line.startswith("  |MutationCounter_Modified"):
                modlength = []
                for x in flist[i+6:i+27]:
                    modlength.append(float(x.strip().split('\t')[1]))
                modmuts = []
                for x in flist[i+32:i+53]:
                    modmuts.append(float(x.strip().split('\t')[1]))
            if line.startswith("  |MutationCounter_Untreated"):
                untlength = []
                for x in flist[i+6:i+27]:
                    untlength.append(float(x.strip().split('\t')[1]))
                untmuts = []
                for x in flist[i+32:i+53]:
                    untmuts.append(float(x.strip().split('\t')[1]))
    data = {'Read_length': ['(0,49)', '(50,99)', '(100,149)', '(150,199)',
                            '(200,249)', '(250,299)', '(300,349)', '(350,399)',
                            '(400,449)', '(450,499)', '(500,549)', '(550,599)',
                            '(600,649)', '(650,699)', '(700,749)', '(750,799)',
                            '(800,849)', '(850,899)', '(900,949)', '(950,999)',
                            '>1000'],
            'Mutation_count': list(range(21)),
            'Modified_read_length': modlength,
            'Modified_mutations_per_molecule': modmuts,
            'Untreated_read_length': untlength,
            'Untreated_mutations_per_molecule': untmuts}
    return pd.DataFrame(data)


def plotSkyline(axis, profile, label=None, column='Reactivity_profile'):
    """Creates a skyline plot on the given axis from the profile and column name
    passed. Label is the sample name that should appear on the legend.

    Args:
        axis (pyplot axis): axis on which to plot skyline
        profile (dataFrame): dataFrame containing profile information
        label (string, optional): Label given to plot on legend. Defaults to None.
        column (str, optional): Column to plot. Defaults to 'Reactivity_profile'.
    """
    x = []
    y = []
    # converts standard plot to skyline plot.
    for n, r in zip(profile['Nucleotide'], profile[column]):
        x.extend([n - 0.5, n + 0.5])
        y.extend([r, r])
    axis.plot(x, y, label=label)


def addSeqBar(axis, profile, yvalue=0.005):
    """Takes a profile and adds a colored sequence bar along the bottom of the
    given axis. ylim may need to be adjusted to accomodate. ylim must be set
    before calling this function.

    Args:
        axis (pyplot axis): axis to which sequence bar will be added
        profile (Pandas dataFrame): dataFrame containing profile data
        yvalue (float, optional): y data value at which sequence bar is added.
            Defaults to 0.005. (barely above x-axis)
    """
    # set font style and colors for each nucleotide
    font_prop = mp.font_manager.FontProperties(
        family="monospace", style="normal", weight="bold", size="12")
    color_dict = {"A": "#f20000", "U": "#f28f00",
                  "G": "#00509d", "C": "#00c200"}
    ymin, ymax = axis.get_ylim()
    yvalue = (ymax-ymin)*yvalue + ymin
    for i, seq in profile["Sequence"]:
        col = color_dict[seq.upper()]
        axis.annotate(seq, xy=(i + 1, yvalue), xycoords='data'
                      fontproperties=font_prop,
                      color=col, horizontalalignment="center")


def getWidth(sample):
    """Takes a dataFrame object and returns the appropriate plot width to pass
    to figsize based on the length of your RNA. For example:
    fig, ax = plt.subplots(1, figsize=(pt.getWidth(profile), 7))

    Args:
        sample (dataFrame object): Pandas dataFrame containing profile info.

    Returns:
        float: appropriate figure width for the size of RNA.
    """
    left_inches = 0.9
    right_inches = 0.4
    sp_width = len(sample['Nucleotide']) * 0.1
    fig_width = max(7, sp_width + left_inches + right_inches)
    return fig_width


def plotBMprofiles(ax, reactivities):
    """Reads in a BM file and plots a skyline of each population on a
    single axis. Population percentages are included in the legend.
    Also adds sequence bar along the bottom and sets the width of x axis.

    Args:
        ax (pyplot axis): ax to which skylines will be plotted
        reactivities (string): path to BM file with population reactivities
    """

    # read in 2 line header
    with open(reactivities) as inf:
        header1 = inf.readline().strip().split()
        header2 = inf.readline().strip().split()
    # number of components
    self.components = int(header1[0])
    # population percentage of each component
    self.p = header2[1:]
    # build column names for reading in BM file
    colnames = ["Nucleotide", "Sequence"]
    for i in range(self.components):
        colnames.append("nReact"+str(i))
        colnames.append("Raw"+str(i))
        colnames.append("blank"+str(i))
    colnames.append("Background")
    # read in BM file
    self.reactivities = pd.read_csv(reactivities, sep='\t', header=2,
                                    names=colnames)
    # Add skylines, seqbar, and legend to axis. Set axis width.
    for i in range(self.components):
        pt.plotSkyline(axis, self.reactivities,
                       label="{}: {}".format(i, self.p[i]),
                       column="nReact{}".format(i))
    pt.addSeqBar(axis, self.reactivities, yvalue=-0.1)
    axis.legend(title="Component: Population", loc=1)
    axis.set_xlim(0, len(self.reactivities))
