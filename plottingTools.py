#!/usr/bin/env python
import matplotlib as mp
import seaborn as sns
import pandas as pd
from scipy import stats
import numpy as np
import RNAtools2 as rna

sns.set_style("ticks")
sns.set_context("talk")


class shapeMapperLog():

    def __init__(self, logfile=None):
        if logfile is not None:
            self.readHistograms(logfile)

    def readHistograms(self, logfile):
        """Parses the shapemapper_log.txt file for read length and mutations per
        molecule histograms. Requires '--per-read-histograms' flag when running
        shapemapper. Stored as dataFrame at self.histogram, columns are:
        [Read_length, Mutation_count, Modified_read_length,
        Modified_mutations_per_molecule, Untreated_read_length,
        Untreated_mutations_per_molecule]

        Args:
            logfile (string): path to shapemapper_log.txt file
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
        data = {'Read_length': ['0-49', '50-99', '100-149', '150-199',
                                '200-249', '250-299', '300-349', '350-399',
                                '400-449', '450-499', '500-549', '550-599',
                                '600-649', '650-699', '700-749', '750-799',
                                '800-849', '850-899', '900-949', '950-999',
                                '>1000'],
                'Mutation_count': list(range(21)),
                'Modified_read_length': modlength,
                'Modified_mutations_per_molecule': modmuts,
                'Untreated_read_length': untlength,
                'Untreated_mutations_per_molecule': untmuts}
        self.histogram = pd.DataFrame(data)

    def plotMutsPerMol(self, ax, sample="Modified", upper_limit=10):
        x = self.histogram['Mutation_count'][:upper_limit]
        y = self.histogram[sample+'_mutations_per_molecule'][:upper_limit]
        ax.plot(x, y)

    def setMutsPerMol(self, ax, labels=[]):
        ax.legend(labels=labels,
                  title="Samples")
        ax.set(xlabel="Mutations per molecule",
               ylabel="Percentage of Reads",
               title='Mutations per molecule distribution')

    def plotReadLength(self, ax, sample="Modified", upper_limit=10, n=1, of=1):
        width = 0.8/of
        x = np.arange(upper_limit) - 0.4 - (width/2) + (width*n)
        y = self.histogram[sample+'_read_length'][:upper_limit]
        ax.bar(x, y, width)

    def setReadLength(self, ax, labels=[], upper_limit=10):
        ax.legend(labels=labels,
                  title="Samples")
        ax.set(xticks=range(upper_limit),
               xlabel='Read Length',
               ylabel='Percentage of Reads',
               title='Read length distribution')
        ax.set_xticklabels(self.histogram["Read_length"][:upper_limit],
                           rotation=45, ha="right")


class ReactivityProfile():

    def __init__(self, profile=None, ctfile=None, sample=None):
        self.sample = sample
        if profile is not None:
            self.profile = pd.read_csv(profile, sep='\t')
            self.length = len(self.profile)
        if ctfile is not None:
            self.ct = rna.CT(ctfile)

    def plotSkyline(self, axis, column='Reactivity_profile', label=None):
        """Creates a skyline plot on the given axis from the profile and column name
        passed. Label is the sample name that should appear on the legend.

        Args:
            axis (pyplot axis): axis on which to plot skyline
            profile (dataFrame): dataFrame containing profile information
            label (string, optional): Label given to plot on legend. Defaults to None.
            column (str, optional): Column to plot. Defaults to 'Reactivity_profile'.
        """
        x = [-0.5]
        y = [0]
        # converts standard plot to skyline plot.
        for n, r in zip(self.profile['Nucleotide'], self.profile[column]):
            x.extend([n - 0.5, n + 0.5])
            y.extend([r, r])
        x.append(x[-1])
        y.append(0)
        if label is None:
            label = self.sample
        axis.plot(x, y, label=label)

    def addSeqBar(self, axis, yvalue=0.005):
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
        for i, seq in enumerate(self.profile["Sequence"]):
            col = color_dict[seq.upper()]
            axis.annotate(seq, xy=(i + 1, yvalue), xycoords='data',
                          fontproperties=font_prop,
                          color=col, horizontalalignment="center")

    def figsize(self, rows, cols):
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
        sp_width = self.length * 0.1
        fig_width = max(7, sp_width + left_inches + right_inches)
        return (fig_width*cols, 6*rows)

    def plotBMprofiles(self, axis, reactivityfile):
        """Reads in a BM file and plots a skyline of each population on a
        single axis. Population percentages are included in the legend.
        Also adds sequence bar along the bottom and sets the width of x axis.

        Args:
            ax (pyplot axis): ax to which skylines will be plotted
            reactivities (string): path to BM file with population reactivities
        """
        # read in 2 line header
        with open(reactivityfile) as inf:
            header1 = inf.readline().strip().split()
            header2 = inf.readline().strip().split()
        # number of components
        components = int(header1[0])
        # population percentage of each component
        p = header2[1:]
        # build column names for reading in BM file
        colnames = ["Nucleotide", "Sequence"]
        for i in range(components):
            colnames.append("nReact"+str(i))
            colnames.append("Raw"+str(i))
            colnames.append("blank"+str(i))
        colnames.append("Background")
        # read in BM file
        self.profile = pd.read_csv(reactivityfile, sep='\t', header=2,
                                   names=colnames)
        # Add skylines, seqbar, and legend to axis. Set axis width.
        for i in range(components):
            self.plotSkyline(axis, label="{}: {}".format(i, p[i]),
                             column="nReact{}".format(i))
        self.addSeqBar(axis)
        axis.legend(title="Component: Population", loc=1)
        axis.set_xlim(0, len(self.profile))

    def plotRegression(self, ax, comp_profile, ctfile=None, column="Reactivity_profile"):
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
        p1 = self.profile[column].copy()
        p2 = comp_profile.profile[column].copy()

        ax.plot([0, 1], [0, 1], color='black')
        notNans = ~np.isnan(p1) & ~np.isnan(p2)
        p1 = p1[notNans]
        p2 = p2[notNans]
        gradient, _, r_value, _, _ = stats.linregress(p1, p2)
        ax.text(0.1, 0.8,
                'R^2 = {:.2f}\nslope = {:.2f}'.format(r_value**2, gradient),
                transform=ax.transAxes)
        if ctfile is not None:
            ct = pd.read_csv(ctfile, sep='\s+',
                             usecols=[4], names=['j'], header=0)
            paired = ct.j != 0
            unpaired = ct.j == 0
            paired = paired[notNans]
            unpaired = unpaired[notNans]
            ax.scatter(p1[paired], p2[paired], label="Paired")
            ax.scatter(p1[unpaired], p2[unpaired], label="Unpaired")
        else:
            ax.scatter(p1, p2)

    def boxplot(self, ax, samples=None):
        cols = ["Modified_rate", "Untreated_rate"]
        if samples is not None:
            profs = [self.profile[cols].assign(Sample=1)]
            for i, d in enumerate(samples):
                profs.append(d.profile[cols].assign(Sample=i+2))
            data = pd.concat(profs)
        else:
            data = self.profile[["Modified_rate", "Untreated_rate"]]
            data = data.assign(Location=1)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        data.head()
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v')
        ax.set(yscale='log', ylim=(0.0005, 0.5), ylabel="Mutation Rate")
