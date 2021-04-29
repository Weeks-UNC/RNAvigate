#!/usr/bin/env python

# general python packages
import matplotlib.colors as mc
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import pandas as pd
from scipy import stats
import numpy as np
from numpy import nanpercentile as percentile

# packages from github.com/Weeks-UNC
import RNAtools2 as rna
import pmanalysis as pm


def create_code_button():
    from IPython.display import display, HTML
    display(HTML('''<script>
                 code_show=true;
                 function code_toggle() {
                 if (code_show) {$('div.input').hide();}
                 else {$('div.input').show();}
                 code_show = !code_show
                 }
                 $( document ).ready(code_toggle);
                 </script>
                 <form action="javascript:code_toggle()">
                 <input type="submit" value="Hide/show raw code.">
                 </form>'''))


# STYLE SHEET
###############################################################################
sns.set_context("talk")
sns.set_style("ticks")
colors = [
    '#a100ffff',  # Purple
    '#edc600ff',  # Yellow
    '#0092edff',  # Blue
    '#ff8300ff',  # Orange
    '#ff48e9ff',  # Pink
    '#3fd125ff'  # Green
]
sns.set_palette(colors)

# COPIED FROM SHAPEMAPPER2 some of this might be inappropriately applied
# to all plots
# TODO: look into passing dict to mp.rc()
###############################################################################
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
###############################################################################


class Sample():

    def __init__(self,
                 sample=None,
                 profile=None,
                 ctfile=None,
                 structure=None,
                 logfile=None,
                 ringfile=None,
                 pairfile=None,
                 dance_reactivities=None,
                 structure_cassettes=False,
                 clip=[0, 0]):
        if structure_cassettes:
            self.clip = [14, -43]
        else:
            self.clip = clip
        self.sample = sample
        if profile is not None:
            self.profile = pd.read_csv(profile, sep='\t')
            self.length = len(self.profile)
        if ctfile is not None:
            self.ct = rna.CT(ctfile)
        if ringfile is not None:
            self.read_ring(ringfile)
        if pairfile is not None:
            self.read_pair(pairfile)
        if structure is not None:
            self.structure_type = structure.split('.')[-1]
            if self.structure_type == 'xrna':
                self.read_xrna(structure)
            elif self.structure_type == 'varna':
                self.read_varna(structure)
            elif self.structure_type == 'nsd':
                self.read_nsd(structure)
            elif self.structure_type == 'cte':
                self.read_cte(structure)
            else:
                print(f'{self.structure_type} not a valid structure file')
        if logfile is not None:
            self.read_log(logfile)
        if dance_reactivities is not None:
            self.read_dance_reactivities(dance_reactivities)

###############################################################################
# Parsing data files
#   structure files
#       read_xrna
#       read_cte
#       read_nsd
#   MaP data files
#       read_rings
#       read_pairs
#       read_log
#       read_dance_reactivities
#   Future:
#       read_varna
#       read_jumpdeletions
#       etc.
###############################################################################

    def read_xrna(self, xrnafile):
        basepairs, sequence, xcoordinates, ycoordinates = [], '', [], []
        with open(xrnafile, 'r') as file:
            for line in file.readlines():
                if any(line.startswith(nt) for nt in 'GUAC'):
                    nucdata = line.split(' ')
                    sequence += nucdata[0]
                    xcoordinates.append(float(nucdata[1]))
                    ycoordinates.append(float(nucdata[2]))
                elif line.startswith("<BasePairs"):
                    helix = line.split(' ')
                    i = int(helix[1].split('=')[1].strip("'"))
                    j = int(helix[3].split('=')[1].strip("'"))
                    length = int(helix[2].split('=')[1].strip("'"))
                    for x in range(length):
                        basepairs.append([i+x, j-x])
        self.structure_length = len(sequence)
        self.nucleotides = np.arange(self.length)+1
        self.basepairs = np.array(basepairs)
        self.sequence = sequence
        self.xcoordinates = np.array(xcoordinates)
        self.ycoordinates = np.array(ycoordinates)

    def read_cte(self, ctefile):
        names = ['nuc', 'seq', 'pair', 'xcoords', 'ycoords']
        ct = pd.read_csv(ctefile, sep=r'\s+', usecols=[0, 1, 4, 8, 10],
                         names=names, header=0)
        self.sequence = ''.join(list(ct['seq']))
        self.structure_length = len(self.sequence)
        self.nucleotides = [int(nuc) for nuc in ct.nuc]
        self.xcoordinates = np.array([float(x) for x in ct.xcoords])
        self.ycoordinates = np.array([float(y) for y in ct.ycoords])
        ct = ct[ct.nuc < ct.pair]
        self.pairs = [[int(i), int(j)] for i, j in zip(ct.nuc, ct.pair)]

    def read_nsd(self, nsdfile):
        basepairs, nucleotides, sequence, xcoordinates, ycoordinates = [], [], '', [], []
        with open(nsdfile, 'r') as file:
            item = ""
            for line in file.readlines():
                line = line.strip().split(' ')
                if "Strand:[" in line:
                    item = "strand"
                    continue
                elif "]" in line:
                    item = ""
                elif "Pairs:[" in line:
                    item = "pairs"
                    continue
                if item == "strand":
                    nucs = [item.split(':') for item in line]
                    nuc = int(nucs[1][1])
                    sequence += nucs[2][1]
                    xcoordinates.append(float(nucs[4][1]))
                    ycoordinates.append(float(nucs[5][1]))
                elif item == "pairs":
                    pairs = line[1].split(":")[1:]
                    basepair = [int(nuc.strip('"')) for nuc in pairs]
                    basepairs.append(basepair)
        self.sequence = sequence
        self.structure_length = len(self.sequence)
        self.nucleotides = np.array(nucleotides)
        self.basepairs = np.array(basepairs)
        self.xcoordinates = np.array(xcoordinates)
        self.ycoordinates = np.array(ycoordinates)

    def read_ring(self, ringfile):
        with open(ringfile, 'r') as file:
            header = file.readline().split('\t')
        self.rings_window = int(header[1].split('=')[1])
        self.rings = pd.read_csv(ringfile, sep='\t', header=1)

    def read_pair(self, pairfile):
        with open(pairfile, 'r') as file:
            header = file.readline().split('\t')
        self.pairs_window = int(header[1].split('=')[1])
        self.pairs = pd.read_csv(pairfile, sep='\t', header=1)

    def read_log(self, logfile):
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
        self.log = pd.DataFrame(data)

    def read_dance_reactivities(self, reactivityfile):
        # read in 2 line header
        with open(reactivityfile) as inf:
            header1 = inf.readline().strip().split()
            header2 = inf.readline().strip().split()
        # number of components
        self.dance_components = int(header1[0])
        # population percentage of each component
        self.dance_percents = header2[1:]
        # build column names for reading in BM file
        colnames = ["Nucleotide", "Sequence"]
        for i in range(self.dance_components):
            colnames.append("nReact"+str(i))
            colnames.append("Raw"+str(i))
            colnames.append("blank"+str(i))
        colnames.append("Background")
        # read in BM file
        dance_reactivities = pd.read_csv(reactivityfile, sep='\t',
                                         header=2, names=colnames)
        self.profile = self.profile.join(dance_reactivities, rsuffix='_dance')

###############################################################################
# Skyline plotting functions
#     get_skyline_figsize
#     plot_skyline
#     plot_sequence
#     make_skyline
#     make_dance_skyline
#     Future:
#         set_skyline
###############################################################################

    def get_skyline_figsize(self, rows, cols):
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

    def plot_skyline(self, axis, column='Reactivity_profile', label=None):
        """Creates a skyline plot on the given axis from the profile and column
        passed. Label is the sample name that should appear on the legend.

        Args:
            axis (pyplot axis): axis on which to plot skyline
            label (string, optional): Label for legend.
                Default: None.
            column (str, optional): Column to plot.
                Default: Reactivity_profile.
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

    def plot_sequence(self, axis, yvalue=0.005):
        """Adds a colored sequence bar along the bottom of the given axis.
        ylim may need to be adjusted to accomodate. ylim must be set
        before calling this function.

        Args:
            axis (pyplot axis): axis to which sequence bar will be added
            yvalue (float, optional): y data value at which sequence bar is added.
                Default: 0.005. (barely above x-axis)
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

    def make_skyline(self, column="Reactivity_profile"):
        fig, ax = plt.subplots(1, figsize=self.get_skyline_figsize(1, 1))
        self.plot_skyline(ax)
        self.plot_sequence(ax)
        ax.set(title="Raw Reactivity Profile",
               xlim=[0, self.length],
               xticks=range(0, self.length, 20))
        ax.legend(title="Samples")

    def make_dance_skyline(self):
        fig, ax = plt.subplots(1, figsize=self.get_skyline_figsize(1, 1))
        for i in range(self.dance_components):
            self.plot_skyline(ax, label=f"{i}: {self.dance_percents[i]}",
                              column=f"nReact{i}")
        self.plot_sequence(ax)
        ax.legend(title="Component: Population", loc=1)
        ax.set(xlim=(0, self.length),
               title="DANCE-MaP Reactivities",
               xlabel="Nucleotide",
               ylabel="Profile")

###############################################################################
# ShapeMapper Log plotting functions
#     plot_log_MutsPerMol
#     set_log_MutsPerMol
#     make_log_MutsPerMol
#     plot_log_ReadLength
#     set_log_ReadLength
#     make_log_ReadLength
#     make_log_qc
###############################################################################

    def plot_log_MutsPerMol(self, ax, sample="Modified", upper_limit=10):
        x = self.log['Mutation_count'][:upper_limit]
        y = self.log[sample+'_mutations_per_molecule'][:upper_limit]
        ax.plot(x, y, label=self.sample+": "+sample)

    def set_log_MutsPerMol(self, ax):
        ax.legend(title="Samples")
        ax.set(xlabel="Mutations per molecule",
               ylabel="Percentage of Reads",
               title='Mutations per molecule distribution')

    def make_log_MutsPerMol(self, ax):
        self.plot_log_MutsPerMol(ax, sample="Modified")
        self.plot_log_MutsPerMol(ax, sample="Untreated")
        self.set_log_MutsPerMol(ax)

    def plot_log_ReadLength(self, ax, sample="Modified", upper_limit=10, n=1, of=1):
        width = 0.8/of
        x = np.arange(upper_limit) - 0.4 - (width/2) + (width*n)
        y = self.log[sample+'_read_length'][:upper_limit]
        ax.bar(x, y, width, label=self.sample+": "+sample)

    def set_log_ReadLength(self, ax, upper_limit=10):
        ax.legend(title="Samples")
        ax.set(xticks=range(upper_limit),
               xlabel='Read Length',
               ylabel='Percentage of Reads',
               title='Read length distribution')
        ax.set_xticklabels(self.log["Read_length"][:upper_limit],
                           rotation=45, ha='right', label=self.sample)

    def make_log_ReadLength(self, ax):
        self.plot_log_ReadLength(ax, sample="Modified", n=1, of=2)
        self.plot_log_ReadLength(ax, sample="Untreated", n=2, of=2)
        self.set_log_ReadLength(ax)

    def make_log_qc(self):
        fig, ax = plt.subplots(1, 3, figsize=(21, 7))
        self.make_log_MutsPerMol(ax[0])
        self.make_log_ReadLength(ax[1])
        self.boxplot(ax[2])
        plt.tight_layout()

###############################################################################
# Arc Plot plotting functions
#     _add_arc
#     get_ap_figsize
#     set_ap
#     plot_ap_ct
#     plot_ap_ctcompare
#     plot_ap_profile
#     plot_ap_pairs
#     make_ap
#     Future:
#         plot_ap_probabilities
#         plot_ap_rings
###############################################################################

    def _add_arc(self, ax, i, j, window, color, alpha, panel):
        """internal function used to add arcs based on data

        Args:
            ax (pyplot axis): axis to which arc will be added
            i (int): leftmost position of left side
            j (int): leftmost position of right side
            window (int): number of nucleotides included in arc
            color (color): color of arc
            alpha (float): transparency value [0-1]
            panel (str): "top" or "bottom" panel of arc plot
        """
        if panel == "top":
            center = ((i+j)/2., 0)
            theta1 = 0
            theta2 = 180
        if panel == "bottom":
            center = ((i+j+window-1)/2., -2)
            theta1 = 180
            theta2 = 360
        radius = 0.5+(j+window-1-i)/2.
        arc = mp.patches.Wedge(center, radius, theta1, theta2, color=color,
                               alpha=alpha, width=window, ec='none')
        ax.add_patch(arc)

    def get_ap_figsize(self, rows, cols):
        dim = self.length * 0.1 + 1
        return (dim*cols, dim*rows)

    def set_ap(self, ax):
        """Sets the aspect ratio to equal, removes y-axis, and sets x-axis at zero

        Args:
            ax (pyplot axis): axis to be set up for arc plot
        """
        ax.set_aspect('equal')
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        pairs = self.ct.pairList()
        height = max([abs(i-j) for i, j in pairs])/2
        ax.set(xlim=(0, self.length),
               ylim=(-height-5, height+1))

    def plot_ap_ct(self, ax):
        """Adds arc plot representation of the ct structure to given axis

        Args:
            ax (pyplot axis): axis to which structure will be added
            ctpath (str): path to ct file
        """
        ct_pairs = self.ct.pairList()
        for i, j in ct_pairs:
            self._add_arc(ax, i, j, 1, 'black', 0.7, 'top')

    def plot_ap_ctcompare(self, ax, compctpath):
        """Adds structure comparison arc plot to the given axis

        Args:
            ax (pyplot axis): axis to which structure comparison is added
            ctpath1 (str): path to first ct file
            ctpath2 (str): path to second ct file
        """
        ct1 = self.ct.pairList()
        ct2 = rna.CT(compctpath).pairList()
        shared = ct1.union(ct2)
        ref = ct1.difference(ct2)
        comp = ct2.difference(ct1)
        sharedcolor = (150/255., 150/255., 150/255.)
        refcolor = (38/255., 202/255., 145/255.)
        compcolor = (153/255., 0.0, 1.0)
        for i, j in comp:
            self._add_arc(ax, i, j, 1, compcolor, 0.7, 'top')
        for i, j in ref:
            self._add_arc(ax, i, j, 1, refcolor, 0.7, 'top')
        for i, j in shared:
            self._add_arc(ax, i, j, 1, sharedcolor, 0.7, 'top')

    def plot_ap_profile(self, ax):
        """Adds bar graph of normalized reactivity profile to given axis

        Args:
            ax (pyplot axis): axis to which reactivity profile is added
        """
        near_black = (0, 0, 1 / 255.0)
        orange_thresh = 0.4
        red_thresh = 0.85
        values = self.profile['Norm_profile']
        cindex = np.zeros(len(values), dtype=int)
        # where values are not NaNs, add 1 to color index array
        cindex[np.array(np.logical_not(np.isnan(values)), dtype=bool)] += 1
        # where values excede orange threshold (0.4), add 1 to color index array
        cindex[np.array(values > orange_thresh, dtype=bool)] += 1
        # where values excede red threshold (0.85), add 1 to color index array
        cindex[np.array(values > red_thresh, dtype=bool)] += 1
        # create color map array based on cindex
        colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
        ax.bar(self.profile['Nucleotide'], values*5, align="center",
               width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
               yerr=self.profile['Norm_stderr'], ecolor=near_black, capsize=1)
        self.plot_sequence(ax, yvalue=0.5)

    def plot_ap_pairs(self, ax, all=False):
        """add PAIR-MaP data to the given axis

        Args:
            ax (pyplot axis): axis to which PAIR data is added
            pairmappath (str): path to PairMapper data file
            window (int, optional): Window size used in pairmapper. Defaults to 3.
        """
        window = self.pairs_window
        pm = self.pairs
        primary = pm[pm['Class'] == 1]
        primary = list(zip(primary['i'], primary['j']))
        secondary = pm[pm['Class'] == 2]
        secondary = list(zip(secondary['i'], secondary['j']))
        other = pm[pm['Class'] == 0]
        other = list(zip(other['i'], other['j']))
        for i, j in secondary:
            self._add_arc(ax, i, j, window,
                          (30/255., 194/255., 1.), 0.2, 'bottom')
        for i, j in primary:
            self._add_arc(ax, i, j, window, (0., 0., 243/255.), 0.7, 'bottom')

    def make_ap(self, type=None):
        fig, ax = plt.subplots(1, figsize=self.get_ap_figsize(1, 1))
        self.set_ap(ax)
        self.plot_ap_ct(ax)
        self.plot_ap_profile(ax)
        if type == "pairs":
            self.plot_ap_pairs(ax)
        self.plot_sequence(ax, yvalue=0.5)
        ax.annotate(self.sample, xy=(0.1, 0.9),
                    xycoords="axes fraction", fontsize=60, ha='center')

###############################################################################
# Secondary Structure graph plotting functions
#     get_ss_figsize
#     set_ss
#     plot_ss_structure
#     set_ss_seqcolors
#     setColorsByProfile
#     setColorsByNucleotide
#     setColorsByPosition
#     plot_ss_sequence
#     filterRings
#     plot_ss_rings
#     plot_ss_positions
#     make_ss
###############################################################################

    def get_ss_figsize(self, rows, cols):
        xmin = min(self.xcoordinates)
        xmax = max(self.xcoordinates)
        ymin = min(self.ycoordinates)
        ymax = max(self.ycoordinates)
        scale = 0.025
        return ((xmax-xmin)*scale*cols, (ymax-ymin)*scale*rows)

    def set_ss(self, ax):
        ax.set_aspect('equal')
        ax.axis('off')

    def plot_ss_structure(self, ax):
        for pair in self.basepairs:
            xcoords = [self.xcoordinates[pair[0]-1],
                       self.xcoordinates[pair[1]-1]]
            ycoords = [self.ycoordinates[pair[0]-1],
                       self.ycoordinates[pair[1]-1]]
            ax.plot(xcoords, ycoords, color='grey',
                    linestyle=(0, (1, 1)), zorder=0)
        ax.plot(self.xcoordinates, self.ycoordinates, color='grey', zorder=0)

    def set_ss_seqcolors(self):
        def whiteOrBlack(color):
            """excepts any valid mpl color and returns white or black,
            which ever has higher contrast"""
            try:
                c = mc.cnames[color]
            except:
                c = color
            if mc.rgb_to_hsv(mc.to_rgb(c))[2] < 0.179:
                return 'white'
            else:
                return 'black'
        self.letter_colors = np.apply_along_axis(whiteOrBlack, 0, self.colors)

    def setColorsByProfile(self, clip=False):
        if clip:
            start, end = self.clip
        else:
            start, end = 0, -1
        cmap = ['grey', 'black', 'orange', 'red']
        bins = [0, 0.4, 0.85]
        profcolors = []
        for x in self.profile["Norm_profile"][start:end]:
            profcolors.append(sum([b < x for b in bins]))
        self.colors = np.array([cmap[val] for val in profcolors])

    def setColorsByNucleotide(self):
        original_color_dict = {"A": "#f20000", "U": "#f28f00",
                               "G": "#00509d", "C": "#00c200"}
        color_dict = {"A": "#366ef0", "U": "#9bb9ff",
                      "G": "#f04c4c", "C": "#ffa77c"}
        self.colors = np.array([color_dict[seq] for seq in self.sequence])

    def setColorsByPosition(self):
        cmap = plt.get_cmap('viridis')
        self.colors = np.array([cmap(n/self.length)
                                for n in range(self.length)])

    def plot_ss_sequence(self, ax, colorby=None, markers='o'):
        if colorby == "profile":
            self.setColorsByProfile(clip=True)
        elif colorby == "sequence":
            self.setColorsByNucleotide()
        elif colorby == "position":
            self.setColorsByPosition()
        elif isinstance(colorby, list) and len(colorby) == self.length:
            self.colors = np.array(colorby)
        else:
            print("Invalid colorby: choices are profile, sequence, or a list" +
                  " of length equal to RNA, defaulting to 'sequence'.")
            self.setColorsByNucleotide()
        for nuc in "GUAC":
            mask = [nt == nuc for nt in self.sequence]
            xcoords = self.xcoordinates[mask]
            ycoords = self.ycoordinates[mask]
            if markers == 'sequence':
                marker = "$"+nuc+"}$"
            else:
                marker = markers
            ax.scatter(xcoords, ycoords, marker=marker,
                       c=self.colors[mask])

    def filterRings(self, cdAbove=None, cdBelow=None, statistic=None, zscore=None,
                    ctfile=None, clip=False):
        self.rings_mask = np.ones(len(self.rings), dtype=bool)
        if clip:
            start = self.clip[0]
            end = self.structure_length + self.clip[0] + 1
            i = self.rings["i"].values
            j = self.rings["j"].values
            iinrange = (start < i) & (i < end)
            jinrange = (start < j) & (j < end)
            self.rings_mask = self.rings_mask & iinrange & jinrange
            self.rings_offset = 14
        if cdAbove is not None or cdBelow is not None:
            pairs = [tuple(pair) for pair in self.basepairs]
            ct = RNA.CT()
            ct.pair2CT(pairs=basepairs, seq=self.sequence)
            mask = []
            for i, j in zip(self.rings["i"], self.rings["j"]):
                if isinstance(cdAbove, int):
                    mask.append(ct.contactDistance(i, j) > cdAbove)
                if isinstance(cdBelow, int):
                    mask.append(ct.contactDistance(i, j) < cdBelow)
            mask = np.array(mask, dtype=bool)
            self.rings_mask = self.rings_mask & mask
        if statistic is not None:
            mask = self.rings["Statistic"] > statistic
            self.rings_mask = self.rings_mask & mask
        if zscore is not None:
            mask = self.rings["Zij"] > zscore
            self.rings_mask = self.rings_mask & mask

    def plot_ss_rings(self, ax, statistic="Statistic", bins=None):
        # Blue to light gray to red.
        cmap = plt.get_cmap("coolwarm")

        if statistic == 'z':
            if bins is None:
                bins = [-50, -5, -1, 0, 1, 5, 50]
            else:
                bins = [-1e10, -bins[1], -bins[0], 0, bins[0], bins[1], 1e10]
        else:
            if bins is None:
                bins = [-1e10, -100, -20, 0, 20, 100, 1e10]
            else:
                bins = [-1e-10, -bins[1], -bins[0], 0, bins[0], bins[1], 1e10]
        # rings are plotted from lowest to highest
        filtered = self.rings[self.rings_mask].sort_values(by=['Statistic'])
        filtered["i"] -= self.rings_offset
        filtered["j"] -= self.rings_offset
        for i, j, stat, sign in zip(filtered['i'], filtered['j'], filtered[statistic], filtered["+/-"]):
            xcoords = [self.xcoordinates[i+1],
                       self.xcoordinates[j+1+self.rings_window]]
            ycoords = [self.ycoordinates[i+1],
                       self.ycoordinates[j+1+self.rings_window]]
            stat = stat*sign
            if stat < bins[1]:
                stat = 0.0
            elif stat > bins[5]:
                stat = 1.0
            else:
                stat = (stat - bins[1]) / (bins[5] - bins[1])
            ax.plot(xcoords, ycoords,
                    color=cmap(stat), alpha=0.5)

    def plot_ss_positions(self, ax, spacing=20):
        for i in range(0, self.structure_length, spacing):
            ax.annotate(i+1, xy=(self.xcoordinates[i], self.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords='data', fontsize=16,
                        bbox=dict(fc="white", ec='none', pad=0.3))

    def make_ss(self, setPlot=True, ss=True, positions=True,
                rings=True, sequence=True, colorby='profile',
                markers='o'):
        fig, ax = plt.subplots(1, figsize=self.get_ss_figsize(1, 1))
        self.filterRings(clip=True)
        if setPlot is True:
            self.set_ss(ax)
        if ss is True:
            self.plot_ss_structure(ax)
        if sequence is True:
            self.plot_ss_sequence(ax, colorby=colorby, markers=markers)
        if positions is True:
            self.plot_ss_positions(ax)
        if hasattr(self, "rings") and rings is True:
            self.plot_ss_rings(ax)

###############################################################################
# Standard ShapeMapper plotting functions
#     plot_sm_profile
#     plot_sm_depth
#     plot_sm_mutationRate
#     metric_abbreviate
#     make_shapemapper
###############################################################################

    def plot_sm_profile(self, axis):
        yMin, ymax = (-0.5, 4)
        near_black = (0, 0, 1 / 255.0)
        orange_thresh = 0.4
        red_thresh = 0.85
        self.setColorsByProfile()
        sample = self.profile["Norm_profile"].copy()
        sample[np.isnan(sample)] = -1
        axis.bar(self.profile['Nucleotide'], sample, align="center",
                 width=1.05, color=self.colors, edgecolor=self.colors,
                 linewidth=0.0, yerr=self.profile['Norm_stderr'],
                 ecolor=near_black, capsize=1)
        axis.set_title(self.sample, fontsize=16)
        axis.set_ylim(yMin, ymax)
        axis.set_xlim(1, self.length)
        axis.yaxis.grid(True)
        axis.set_axisbelow(True)
        axis.set_xlabel("Nucleotide", fontsize=14, labelpad=0)
        axis.set_ylabel("Shape Reactivity", fontsize=14)
        # add a SHAPE colorbar to the vertical axis
        # uses a little transformation magic to place correctly
        inv = axis.transData.inverted()
        for loc, spine in list(axis.spines.items()):
            if loc == 'left':
                trans = spine.get_transform()
        tp = trans.transform_point([0, 0])
        tp2 = inv.transform_point(tp)
        rectX = tp2[0]
        tpA = (0, 0)
        tpB = (6, 0)
        tpA2 = inv.transform_point(tpA)
        tpB2 = inv.transform_point(tpB)
        rectW = tpB2[0] - tpA2[0]
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

        # put nuc sequence below axis
        self.plot_sequence(axis)

    def plot_sm_depth(self, axis):
        sample = self.profile
        axis.plot(sample['Nucleotide'], sample['Modified_read_depth'],
                  linewidth=1.5, color=rx_color, alpha=1.0, label="Modified")
        axis.plot(sample['Nucleotide'], sample['Untreated_read_depth'],
                  linewidth=1.5, color=bg_color, alpha=1.0, label="Untreated")
        axis.plot(sample['Nucleotide'], sample['Denatured_read_depth'],
                  linewidth=1.5, color=dc_color, alpha=1.0, label="Denatured")
        axis.set_xlim(1, self.length)
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
            formatted_ticks.append(self.metric_abbreviate(val))
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

    def plot_sm_rates(self, axis):
        sample = self.profile
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

    def metric_abbreviate(self, num):
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

    def make_shapemapper(self):
        fig, ax = plt.subplots(3, 1, figsize=self.get_skyline_figsize(3, 1))
        self.plot_sm_profile(ax[0])
        self.plot_sm_rates(ax[1])
        self.plot_sm_depth(ax[2])

##############################################################################
# CONSTRUCTION ZONE
#     plotRegression (doesn't make sense without container class)
#     boxplot (broken for multiple samples)
##############################################################################

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
            data = data.assign(Sample=1)
            xticklabels = [self.sample]
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        data.head()
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v')
        # import MaPplotlib as MaP
        ax.set(yscale='log',
               ylim=(0.0005, 0.5),
               ylabel="Mutation Rate",
               xticklabels=xticklabels)
