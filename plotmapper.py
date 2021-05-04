#!/usr/bin/env python

# general python packages
import matplotlib as mp
import matplotlib.colors as mc
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

# packages from github.com/Weeks-UNC
import RNAtools2 as rna


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

# COPIED FROM SHAPEMAPPER2
# some of this might be inappropriately applied to all plots
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


rx_color = "red"
bg_color = "blue"
dc_color = "darkgoldenrod"
###############################################################################


class Sample():

    def __init__(self,
                 sample=None,
                 profile=None,
                 ct=None,
                 compct=None,
                 ss=None,
                 log=None,
                 rings=None,
                 pairs=None,
                 dance_reactivities=None):
        self.paths = {"profile": profile,
                      "ct": ct,
                      "comptct": compct,
                      "ss": ss,
                      "log": log,
                      "rings": rings,
                      "pairs": pairs,
                      "dance_reactivities": dance_reactivities}
        self.sample = sample
        if profile is not None:
            self.profile = pd.read_csv(profile, sep='\t')
            self.profile_length = len(self.profile)
            self.profile_sequence = ''.join(self.profile["Sequence"].values)
        if ct is not None:
            self.ct = rna.CT(ct)
            self.ct_length = len(self.ct.seq)
            self.ct_sequence = ''.join(self.ct.seq)
        if compct is not None:
            self.compct = rna.CT(compct)
            self.compct_length = len(self.compct.seq)
            self.compct_sequence = ''.join(self.compct.seq)
        if rings is not None:
            self.read_ring(rings)
        if pairs is not None:
            self.read_pair(pairs)
        if ss is not None:
            # get and check file extension, then read
            self.ss_type = ss.split('.')[-1].lower()
            valid_type = self.ss_type in ['xrna', 'varna', 'nsd', 'cte']
            message = f"stucture file type {self.ss_type} not supported"
            assert valid_type, message
            read_file = getattr(self, f"read_{self.ss_type}")
            read_file(ss)
        if log is not None:
            self.read_log(log)
        if dance_reactivities is not None:
            self.read_dance_reactivities(dance_reactivities)

###############################################################################
# Parsing data files
# TODO: make read_nsd function less sucky. I think it's YAML.
#   structure files
#       read_xrna
#       read_varna
#       read_cte
#       read_nsd Kludgey AF
#   MaP data files
#       read_ring
#       read_pair
#       read_log
#       read_dance_reactivities
#   Future:
#       read_jumpdeletions
#       etc.
###############################################################################

    def read_xrna(self, xrnafile):
        import xml.etree.ElementTree as xmlet
        tree = xmlet.parse(xrnafile)
        root = tree.getroot()
        # extract sequence, x and y coordinates
        nucList = root.findall('./Complex/RNAMolecule/')
        nucLists = []
        for i in nucList:
            if i.tag == 'NucListData':
                nucLists.append(i)
        sequence, xcoords, ycoords = '', [], []
        for nt in nucLists[0].text.split('\n'):
            if nt == '':
                continue
            line = nt.split()
            sequence += line[0]
            xcoords.append(float(line[1]))
            ycoords.append(float(line[2]))
        # extract pairing information
        basepairs = []
        for helix in root.findall('./Complex/RNAMolecule/BasePairs'):
            i_outter = int(helix.get('nucID'))
            j_outter = int(helix.get('bpNucID'))
            length = int(helix.get('length'))
            helix_list = [(i_outter+nt, j_outter-nt) for nt in range(length)]
            basepairs.extend(helix_list)
        # make expected arrays
        self.ss_length = len(sequence)
        self.nucleotides = np.arange(self.ss_length)+1
        self.basepairs = np.array(basepairs)
        self.ss_sequence = sequence
        self.xcoordinates = np.array(xcoords)/20
        self.ycoordinates = np.array(ycoords)/20

    def read_varna(self, varnafile):
        import xml.etree.ElementTree as xmlet
        tree = xmlet.parse(varnafile)
        root = tree.getroot()
        # extract sequence, y and x coordinates
        sequence, xcoords, ycoords = "", [], []
        for nt in root.findall('./RNA/bases/nt'):
            base = nt.find('base').text
            sequence += base
            for i in nt:
                if i.get('r') == 'pos':
                    xcoords.append(float(i.get('x')))
                    ycoords.append(float(i.get('y')))
        # extract pairing information
        basepairs = []
        for pair in root.findall('./RNA/BPs/bp'):
            i = int(pair.get('part5'))+1
            j = int(pair.get('part3'))+1
            basepairs.append((i, j))
        self.ss_length = len(sequence)
        self.nucleotides = np.arange(self.ss_length)
        self.basepairs = basepairs
        self.ss_sequence = sequence
        self.xcoordinates = np.array(xcoords)/65
        self.ycoordinates = np.array(ycoords)/65

    def read_cte(self, ctefile):
        names = ['nuc', 'seq', 'pair', 'xcoords', 'ycoords']
        ct = pd.read_csv(ctefile, sep=r'\s+', usecols=[0, 1, 4, 8, 10],
                         names=names, header=0)
        self.ss_sequence = ''.join(list(ct['seq']))
        self.ss_length = len(self.ss_sequence)
        self.nucleotides = [int(nuc) for nuc in ct.nuc]
        self.xcoordinates = np.array([float(x) for x in ct.xcoords])/30.5
        self.ycoordinates = np.array([float(y) for y in ct.ycoords])/30.5
        ct = ct[ct.nuc < ct.pair]
        self.basepairs = [[int(i), int(j)] for i, j in zip(ct.nuc, ct.pair)]

    def read_nsd(self, nsdfile):
        basepairs, sequence, xcoordinates, ycoordinates = [], '', [], []
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
                    nt = [item.split(':') for item in line]
                    sequence += nt[2][1]
                    xcoordinates.append(float(nt[4][1]))
                    ycoordinates.append(float(nt[5][1]))
                elif item == "pairs":
                    pairs = line[1].split(":")[1:]
                    basepair = [int(nuc.strip('"')) for nuc in pairs]
                    basepairs.append(basepair)
        self.ss_sequence = sequence
        self.ss_length = len(self.ss_sequence)
        self.basepairs = np.array(basepairs)
        self.xcoordinates = np.array(xcoordinates)/30.5
        self.ycoordinates = np.array(ycoordinates)/30.5

    def read_ring(self, rings):
        with open(rings, 'r') as file:
            header = file.readline().split('\t')
        self.rings_window = int(header[1].split('=')[1])
        self.rings = pd.read_csv(rings, sep='\t', header=1)

    def read_pair(self, pairs):
        with open(pairs, 'r') as file:
            header = file.readline().split('\t')
        self.pairs_window = int(header[1].split('=')[1])
        self.pairs = pd.read_csv(pairs, sep='\t', header=1)

    def read_log(self, log):
        with open(log, 'r') as f:
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

    def get_skyline_figsize(self, rows, columns):
        """Pass this function call to figsize within pyplot.subplots() to set
        an appropriate figure width and height for skyline plots.

        Args:
            rows (int): number of rows in pyplot figure
            cols (int): number of columns in pyplot figure

        Returns:
            tuple: (width, height) appropriate figsize for pyplot figure
        """
        left_inches = 0.9
        right_inches = 0.4
        ax_width = self.profile_length * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        return (fig_width*columns, fig_height*rows)

    def plot_skyline(self, axis, column='Reactivity_profile', label=None):
        """Plots a skyline on the given axis from the profile file and column
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

    def plot_sequence(self, axis, type, yvalue=0.005):
        """Adds a colored sequence bar along the bottom of the given axis.
        ylim must be set before calling this function.

        Args:
            axis (pyplot axis): axis to which sequence bar will be added
            yvalue (float, optional): y value as a fraction of the ylimits set
                at which sequence bar is added.
                Default: 0.005. (barely above x-axis)
        """
        # set font style and colors for each nucleotide
        font_prop = mp.font_manager.FontProperties(
            family="monospace", style="normal", weight="bold", size="12")
        color_dict = {"A": "#f20000", "U": "#f28f00",
                      "G": "#00509d", "C": "#00c200"}
        # transform yvalue to a y-axis data value
        ymin, ymax = axis.get_ylim()
        yvalue = (ymax-ymin)*yvalue + ymin
        sequence = getattr(self, f"{type}_sequence")
        for i, seq in enumerate(sequence):
            col = color_dict[seq.upper()]
            axis.annotate(seq, xy=(i + 1, yvalue), xycoords='data',
                          fontproperties=font_prop,
                          color=col, horizontalalignment="center")

    def set_skyline(self, ax, title="Raw Reactivity Profile"):
        ax.set(title=title,
               xlim=[0, self.profile_length],
               xticks=range(0, self.profile_length, 20))
        ax.set_xticks(range(0, self.profile_length, 5), minor=True)
        ax.legend(title="Samples")

    def make_skyline(self, ax=None, column="Reactivity_profile"):
        """Creates a skyline figure, including sequence, title, legend, and
        axis labels.

        Args:
            column (str, optional): Name of column from profile.txt to plot.
                Defaults to "Reactivity_profile".
        """
        if ax is None:
            fig, ax = plt.subplots(1, figsize=self.get_skyline_figsize(1, 1))
        self.plot_skyline(ax, column=column)
        self.plot_sequence(ax, type="profile")
        self.set_skyline(ax, title=column.replace("_", " "))

    def make_dance_skyline(self, ax=None):
        """Creates a skyline figure representing the component reactivities of
        the ensemble.
        """
        if ax is None:
            fig, ax = plt.subplots(1, figsize=self.get_skyline_figsize(1, 1))
        for i in range(self.dance_components):
            self.plot_skyline(ax, label=f"{i}: {self.dance_percents[i]}",
                              column=f"nReact{i}")
        self.plot_sequence(ax, type="profile")
        ax.legend(title="Component: Population", loc=1)
        ax.set(xlim=(0, self.profile_length),
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
        """Adds a line plot of the distribution of mutations per molecule to
        the given axis.

        Args:
            ax (pyploy axis): axis on which to add plot
            sample (str, optional): Options are "Modified" or "Untreated".
                    Defaults to "Modified".
            upper_limit (int, optional): Upper limit to mutations per molecule.
                    Defaults to 10.
        """
        x = self.log['Mutation_count'][:upper_limit]
        y = self.log[sample+'_mutations_per_molecule'][:upper_limit]
        ax.plot(x, y, label=self.sample+": "+sample)

    def set_log_MutsPerMol(self, ax):
        """Adds a title, axis labels, and legend appropriate for mutations per
        molecule to the given axis

        Args:
            ax (pyplot axis): axis to set
        """
        ax.legend(title="Samples")
        ax.set(xlabel="Mutations per molecule",
               ylabel="Percentage of Reads",
               title='Mutations per molecule distribution')

    def make_log_MutsPerMol(self, ax):
        """Adds modified and unmodified samples, and sets title, labels, and
        legend for the given axis.

        Args:
            ax (pyplot axis): axis on which to make plot
        """
        self.plot_log_MutsPerMol(ax, sample="Modified")
        self.plot_log_MutsPerMol(ax, sample="Untreated")
        self.set_log_MutsPerMol(ax)

    def plot_log_ReadLength(self, ax, sample="Modified", upper_limit=10,
                            n=1, of=1):
        """adds the read length distribution of the sample to the given axis.

        Args:
            ax (pyplot axis): axis on which to add data
            sample (str, optional): Options are "Modified" or "Untreated".
                    Defaults to "Modified".
            upper_limit (int, optional): upper limit of read length to plot.
                    Defaults to 10. (bin 450-499)
            n (int, optional): Sample number out of...
                    Defaults to 1.
            of (int, optional): Total number of samples to be added to plot.
                    Defaults to 1.
        """
        width = 0.8/of
        x = np.arange(upper_limit) - 0.4 - (width/2) + (width*n)
        y = self.log[sample+'_read_length'][:upper_limit]
        ax.bar(x, y, width, label=self.sample+": "+sample)

    def set_log_ReadLength(self, ax, upper_limit=10):
        """Adds title, axis labels, bin labels, and legend appropriate for read
        depth to the axis.

        Args:
            ax (pyplot axis): axis to be set
            upper_limit (int, optional): Number of bin labels to add.
                    Defaults to 10.
        """
        ax.legend(title="Samples")
        ax.set(xticks=range(upper_limit),
               xlabel='Read Length',
               ylabel='Percentage of Reads',
               title='Read length distribution')
        ax.set_xticklabels(self.log["Read_length"][:upper_limit],
                           rotation=45, ha='right', label=self.sample)

    def make_log_ReadLength(self, ax):
        """Adds modified and untreated readlengths to the given axis, and sets
        labels, titles, and legend.

        Args:
            ax (pyplot axis): axis on which to add data
        """
        self.plot_log_ReadLength(ax, sample="Modified", n=1, of=2)
        self.plot_log_ReadLength(ax, sample="Untreated", n=2, of=2)
        self.set_log_ReadLength(ax)

    def make_log_qc(self):
        """Creates a three panel plot with read depth, mutations per molecule,
        and boxplot of reactivities.
        """
        fig, ax = plt.subplots(1, 3, figsize=(21, 7))
        self.make_log_MutsPerMol(ax[0])
        self.make_log_ReadLength(ax[1])
        self.plot_boxplot(ax[2])
        plt.tight_layout()

###############################################################################
# Arc Plot plotting functions
#     _add_arc
#     get_ap_figsize
#     set_ap
#     plot_ap_ct
#     plot_ap_ctcompare
#     plot_ap_profile
#     plot_ap_rings
#     plot_ap_pairs
#     make_ap
#     Future:
#         plot_ap_probabilities
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
        """Pass this function call to figsize within pyplot.subplots() to set
        an appropriate figure width and height for arc plots.

        Args:
            rows (int): number of rows in pyplot figure
            cols (int): number of columns in pyplot figure

        Returns:
            tuple: (width, height) appropriate figsize for pyplot figure
        """
        dim = self.profile_length * 0.1 + 1
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
        ax.set(xlim=(0, self.profile_length),
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

    def plot_ap_ctcompare(self, ax, panel='top'):
        """Adds structure comparison arc plot to the given axis

        Args:
            ax (pyplot axis): axis to which structure comparison is added
            ctpath1 (str): path to first ct file
            ctpath2 (str): path to second ct file
        """
        ct1 = set(self.ct.pairList())
        ct2 = set(self.compct.pairList())
        shared = ct1.union(ct2)
        ref = ct1.difference(ct2)
        comp = ct2.difference(ct1)
        sharedcolor = (150/255., 150/255., 150/255.)
        refcolor = (38/255., 202/255., 145/255.)
        compcolor = (153/255., 0.0, 1.0)
        for i, j in comp:
            self._add_arc(ax, i, j, 1, compcolor, 0.7, panel)
        for i, j in ref:
            self._add_arc(ax, i, j, 1, refcolor, 0.7, panel)
        for i, j in shared:
            self._add_arc(ax, i, j, 1, sharedcolor, 0.7, panel)
        handles = [mp.patches.Patch(color=sharedcolor, alpha=0.7),
                   mp.patches.Patch(color=refcolor, alpha=0.7),
                   mp.patches.Patch(color=compcolor, alpha=0.7)]
        labels = ["Shared pairs", "Reference only", "Compared only"]
        ax.legend(title="top panel", handles=handles, labels=labels,
                  loc="upper right")

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
        # where values excede orange threshold, add 1 to color index array
        cindex[np.array(values > orange_thresh, dtype=bool)] += 1
        # where values excede red threshold (0.85), add 1 to color index array
        cindex[np.array(values > red_thresh, dtype=bool)] += 1
        # create color map array based on cindex
        colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
        ax.bar(self.profile['Nucleotide'], values*5, align="center",
               width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
               yerr=self.profile['Norm_stderr'], ecolor=near_black, capsize=1)
        self.plot_sequence(ax, 'ct', yvalue=0.5)

    def plot_ap_rings(self, ax, statistic="Statistic", bins=None, **kwargs):
        """Adds ring correlation information to the arc plot

        Args:
            ax (pyplot axis): axis on which to add data
            statistic (str, optional): Options are 'z' or 'Statistic'.
                    Defaults to "Statistic".
            bins (list of 2 ints, optional): statistic bins used for coloring.
                    Defaults to None (5, 50 for 'z', 20, 100 for 'Statistic').
        """
        # Blue to light gray to red.
        self.filter_rings("ct", **kwargs)
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
        filtered = self.rings[self.rings['mask']].sort_values(by=['Statistic'])
        for i, j, stat, sign in zip(filtered['i_offset'], filtered['j_offset'],
                                    filtered[statistic], filtered["+/-"]):
            stat = stat*sign
            if stat < bins[1]:
                stat = 0.0
            elif stat > bins[5]:
                stat = 1.0
            else:
                stat = (stat - bins[1]) / (bins[5] - bins[1])
            self._add_arc(ax, i, j, self.rings_window, cmap(stat), 0.5,
                          'bottom')

    def plot_ap_pairs(self, ax, all=False):
        """add PAIR-MaP data to the given axis

        Args:
            ax (pyplot axis): axis to which PAIR data is added
            pairmappath (str): path to PairMapper data file
            window (int, optional): Window size used in pairmapper.
                    Defaults to 3.
        """
        window = self.pairs_window
        pairs = self.pairs
        primary = pairs[pairs['Class'] == 1]
        primary = list(zip(primary['i'], primary['j']))
        secondary = pairs[pairs['Class'] == 2]
        secondary = list(zip(secondary['i'], secondary['j']))
        for i, j in secondary:
            self._add_arc(ax, i, j, window,
                          (30/255., 194/255., 1.), 0.2, 'bottom')
        for i, j in primary:
            self._add_arc(ax, i, j, window, (0., 0., 243/255.), 0.7, 'bottom')
        if all:
            other = pairs[pairs['Class'] == 0]
            other = list(zip(other['i'], other['j']))
            for i, j in other:
                self._add_arc(ax, i, j, window, 'gray', 0.5, 'bottom')

    def make_ap(self, ax=None, type=None, ctcompare=True):
        """Creates figure with arc plot. Includes ct, profile, sequence, and
        data type passed to type. [rings, pairs, probabilities]

        Args:
            type (str, optional): datatype for bottom panel.
                    options: 'rings', 'pairs', 'probabilities'.
                    Defaults to None.
        """
        if ax is None:
            fig, ax = plt.subplots(1, figsize=self.get_ap_figsize(1, 1))
        self.set_ap(ax)
        if hasattr(self, "compct") and ctcompare is True:
            self.plot_ap_ctcompare(ax)
        else:
            self.plot_ap_ct(ax)
        self.plot_ap_profile(ax)
        if type is not None:
            try:
                getattr(self, f"plot_ap_{type}")(ax)
            except AttributeError:
                print(f"make_ap doesn't support type: {type}")
        self.plot_sequence(ax, "ct", yvalue=0.5)
        ax.annotate(self.sample, xy=(0.1, 0.9),
                    xycoords="axes fraction", fontsize=60, ha='center')

###############################################################################
# Secondary Structure graph plotting functions
#     get_ss_figsize
#     set_ss
#     plot_ss_structure
#     set_ss_seqcolors
#     plot_ss_sequence
#     plot_ss_rings
#     plot_ss_positions
#     make_ss
###############################################################################

    def get_ss_figsize(self, rows, columns):
        """Returns a tuple of the width and height for a secondary structure
        graph, based on the number of rows and columns. Pass this function call
        to the figsize argument of plt.subplots.

        Args:
            rows (int): number of rows in figure
            columns (int): number of columns in plot

        Returns:
            [type]: [description]
        """
        xmin = min(self.xcoordinates)
        xmax = max(self.xcoordinates)
        ymin = min(self.ycoordinates)
        ymax = max(self.ycoordinates)
        scale = 0.5
        width = (xmax-xmin)*scale
        height = (ymax-ymin)*scale
        return (width*columns, height*rows)

    def set_ss(self, ax):
        """Turn off axes and axis border and set aspect ratio to equal

        Args:
            ax (pyplot axis): axis to set
        """
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(self.sample, fontsize=30)

    def plot_ss_structure(self, ax):
        """Plots the skeleton of the secondary structure, with base pairs as
        dotted lines

        Args:
            ax (pyplot axis): axis on which to plot
        """
        for pair in self.basepairs:
            xcoords = [self.xcoordinates[pair[0]-1],
                       self.xcoordinates[pair[1]-1]]
            ycoords = [self.ycoordinates[pair[0]-1],
                       self.ycoordinates[pair[1]-1]]
            ax.plot(xcoords, ycoords, color='grey',
                    linestyle=(0, (1, 1)), zorder=0)
        ax.plot(self.xcoordinates, self.ycoordinates, color='grey', zorder=0)

    def get_ss_seqcolors(self, colors):
        """creates an array of white and black that contrasts well with the
        given color array
        """
        def whiteOrBlack(color):
            """excepts any valid mpl color and returns white or black,
            which ever has higher contrast"""
            try:
                c = mc.cnames[color]
            except KeyError:
                c = color
            if mc.rgb_to_hsv(mc.to_rgb(c))[2] < 0.179:
                return 'white'
            else:
                return 'black'
        letter_colors = np.apply_along_axis(whiteOrBlack, 0, colors)
        return letter_colors

    def plot_ss_sequence(self, ax, colorby='sequence', markers='o'):
        """plots the location of nucleotides on top of the secondary structure
        skeleton.

        Args:
            ax (pyplot axis): axis on which to plot
            colorby (str or list, optional): Options are 'sequence', 'profile',
                    'position' or a lis of valid matplotlib colors.
                    Defaults to None.
            markers (str, optional): Options are a matplotlib marker type or
                    'sequence', which uses the sequence letters as markers.
                    Defaults to 'o' (a filled circle).
        """
        if isinstance(colorby, list) and len(colorby) == self.ss_length:
            self.colors = np.array(colorby)
        try:
            colors = getattr(self, f"get_colorby_{colorby}")("ss")
        except AttributeError:
            print("Invalid colorby: choices are profile, sequence, position" +
                  " or a list of length equal to structure sequence, using" +
                  " sequence.")
            self.get_colorby_sequence("ss")
        if markers == "sequence":
            for nuc in "GUAC":
                mask = [nt == nuc for nt in self.ss_sequence]
                xcoords = self.xcoordinates[mask]
                ycoords = self.ycoordinates[mask]
                marker = "$"+nuc+"}$"
                marker = markers
                ax.scatter(xcoords, ycoords, marker=marker,
                           c=colors[mask])
        else:
            ax.scatter(self.xcoordinates, self.ycoordinates, marker=markers,
                       c=colors)

    def plot_ss_rings(self, ax, statistic="Statistic", bins=None, **kwargs):
        """Adds ring correlations as colored lines to the secondary structure
        graph.

        Args:
            ax (pyplot axis): axis on which to add data
            statistic (str, optional): Options are 'z' or 'Statistic'.
                    Defaults to "Statistic".
            bins (list of 2 ints, optional): statistic bins used for coloring.
                    Defaults to None (5, 50 for 'z', 20, 100 for 'Statistic').
        """
        # Blue to light gray to red.
        self.filter_rings("ss", **kwargs)
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
        filtered = self.rings[self.rings['mask']].sort_values(by=['Statistic'])
        for i, j, stat, sign in zip(filtered['i_offset'], filtered['j_offset'],
                                    filtered[statistic], filtered["+/-"]):
            xcoords = [self.xcoordinates[i-1],
                       self.xcoordinates[j-1+self.rings_window]]
            ycoords = [self.ycoordinates[i-1],
                       self.ycoordinates[j-1+self.rings_window]]
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
        """adds position labels to the secondary structure graph

        Args:
            ax (pyplot axis): axis on which to add labels
            spacing (int, optional): labels every nth nucleotide.
                    Defaults to 20.
        """
        for i in range(0, self.ss_length, spacing):
            ax.annotate(i+1, xy=(self.xcoordinates[i], self.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords='data', fontsize=16,
                        bbox=dict(fc="white", ec='none', pad=0.3))

    def make_ss(self, ax=None, ss=True, positions=True, rings=True, sequence=True,
                colorby='profile', markers='o'):
        """Creates a full secondary structure graph with data plotted

        Args:
            ss (bool, optional): whether to plot skeleton.
                    Defaults to True.
            positions (bool, optional): whether to label positions.
                    Defaults to True.
            rings (bool, optional): whether to plot rings.
                    Defaults to True.
            sequence (bool, optional): whether to add nucleotide markers.
                    Defaults to True.
            colorby (str, optional): method used to color nucleotides.
                    Defaults to 'profile'.
            markers (str, optional): marker type to use for nucleotides.
                    Defaults to 'o'.
        """
        if ax is None:
            fig, ax = plt.subplots(1, figsize=self.get_ss_figsize(1, 1))
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
#     plot_sm_rates
#     metric_abbreviate
#     make_shapemapper
###############################################################################

    def plot_sm_profile(self, axis):
        """Plots classic ShapeMapper normalized reactivity on the given axis

        Args:
            axis (pyplot axis): axis on which to add plot
        """
        yMin, ymax = (-0.5, 4)
        near_black = (0, 0, 1 / 255.0)
        orange_thresh = 0.4
        red_thresh = 0.85
        colors = self.get_colorby_profile("profile")
        sample = self.profile["Norm_profile"].copy()
        sample[np.isnan(sample)] = -1
        axis.bar(self.profile['Nucleotide'], sample, align="center",
                 width=1.05, color=colors, edgecolor=colors,
                 linewidth=0.0, yerr=self.profile['Norm_stderr'],
                 ecolor=near_black, capsize=1)
        axis.set_title(self.sample, fontsize=16)
        axis.set_ylim(yMin, ymax)
        axis.set_xlim(1, self.profile_length)
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
        self.plot_sequence(axis, type="profile")

    def plot_sm_depth(self, axis):
        """Plots classic ShapeMapper read depth on the given axis

        Args:
            axis (pyplot axis): axis on which to add plot
        """
        sample = self.profile
        axis.plot(sample['Nucleotide'], sample['Modified_read_depth'],
                  linewidth=1.5, color=rx_color, alpha=1.0, label="Modified")
        axis.plot(sample['Nucleotide'], sample['Untreated_read_depth'],
                  linewidth=1.5, color=bg_color, alpha=1.0, label="Untreated")
        axis.plot(sample['Nucleotide'], sample['Denatured_read_depth'],
                  linewidth=1.5, color=dc_color, alpha=1.0, label="Denatured")
        axis.set_xlim(1, self.profile_length)
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
                        "NOTE: effective read depths shown in lighter colors",
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
        """Plots classic ShapeMapper mutation rates on the given axis

        Args:
            axis (pyplot axis): axis on which to add plot
        """
        sample = self.profile
        # choose a decent range for axis, excluding high-background positions
        temp_rates = sample['Modified_rate'][sample['Untreated_rate'] <= 0.05]
        near_top_rate = np.nanpercentile(temp_rates, 98.0)
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
        """takes a large number and applies an appropriate abbreviation

        Args:
            num (int): number to be abbreviated

        Returns:
            str: abbreviated number
        """
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
        """Creates a figure with the three classic Shapemapper plots.
        """
        fig, ax = plt.subplots(3, 1, figsize=self.get_skyline_figsize(3, 1))
        self.plot_sm_profile(ax[0])
        self.plot_sm_rates(ax[1])
        self.plot_sm_depth(ax[2])

##############################################################################
# CONSTRUCTION ZONE
#     writeCTE
#     plot_regression
#     get_boxplot_data
#     plot_boxplot
##############################################################################

    def get_pairs_sens_PPV(self):
        "Returns sensitivity and PPV for pair data to the ct structure"
        import pmanalysis as pma
        pm = pma.PairMap(self.paths["pairs"])
        ct = self.ct.copy()
        ct.filterNC()
        ct.filterSingleton()
        p, s = pm.ppvsens_duplex(ct, ptype=1, exact=False)
        return p, s

    def writeCTE(self, outputPath):
        """writes the current structure out to CTE format for Structure Editor.

        Args:
            outputPath (string): path to output cte file to be created
        """
        pairs = [tuple(pair) for pair in self.pairs]
        ct = RNA.CT()
        ct.pair2CT(pairs=pairs, seq=self.sequence)
        # set scaling factors based on data source.
        xscale = {'xrna': 1.525, 'varna': 0.469,
                  'nsd': 1, 'cte': 1}[self.ss_type]
        yscale = {'xrna': -1.525, 'varna': 0.469,
                  'nsd': 1, 'cte': 1}[self.ss_type]
        ctlen = len(ct.num)
        w = open(outputPath, 'w')
        w.write('{0:6d} {1}\n'.format(ctlen, ct.name))
        line = ('{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}' +
                ' ;! X: {5:5d} Y: {6:5d}\n')
        for i in range(ctlen):
            xcoord = round(xscale*self.xcoordinates[i])
            ycoord = round(yscale*self.ycoordinates[i])
            num, seq, cti = ct.num[i], ct.seq[i], ct.ct[i]
            # no nums after ctlen, resets to zero
            nextnum = (num+1) % (ctlen+1)
            cols = [num, seq, num-1, nextnum, cti, xcoord, ycoord]
            w.write(line.format(*cols))
        w.close()

    def plot_regression(self, comp_sample, ax=None, colorby="ct",
                        column="Reactivity_profile"):
        """Plots scatter plot of reactivity profile vs. reactivity profile and
        computes regression metrics R^2 and slope, which are annotated on the
        axis. Can color scatter plot by nucleotide or by paired status.

        Args:
            ax (pyplot axis): axis on which scatter plot appears
            comp_sample (plotmapper sample): sample to be compared
            colorby (str, optional): How to color the scatter plot. Options are
                "nucleotide" or "ct".
                Defaults to None
            column (str, optional): column from profile to be compared
                Defaults to "Reactivity_profile".
        """
        if ax is None:
            ax = plt.gca()
        p1 = self.profile[column].copy()
        p2 = comp_sample.profile[column].copy()

        ax.plot([0, 1], [0, 1], color='black')
        notNans = ~np.isnan(p1) & ~np.isnan(p2)
        p1 = p1[notNans]
        p2 = p2[notNans]
        gradient, _, r_value, _, _ = stats.linregress(p1, p2)
        ax.text(0.1, 0.8, f'R^2 = {r_value**2:.2f}\nslope = {gradient:.2f}',
                transform=ax.transAxes)
        if colorby == "ct":
            paired_list = self.ct.pairedResidueList()
            paired_mask = np.zeros(self.profile_length, dtype=bool)
            for i in paired_list:
                paired_mask[i] = True
            paired_mask = paired_mask[notNans]
            unpaired_mask = ~paired_mask
            ax.scatter(p1[paired_mask], p2[paired_mask], label="Paired")
            ax.scatter(p1[unpaired_mask], p2[unpaired_mask], label="Unpaired")
        if colorby == "nucleotide":
            for nuc in "GUAC":
                sequence = self.profile["Sequence"][notNans]
                nuc_mask = [nt == nuc for nt in sequence]
                ax.scatter(p1[nuc_mask], p2[nuc_mask], label=nuc)
        else:
            ax.scatter(p1, p2)
        s1, s2 = self.sample, comp_sample.sample
        ax.set(xscale='log',
               xlim=[0.0001, 0.3],
               xlabel=s1,
               yscale='log',
               ylim=[0.0001, 0.3],
               ylabel=s2,
               title=f'{s1} vs. {s2}: {column}')
        ax.legend(title=colorby)

    def get_boxplot_data(self, sample=1,
                         cols=["Modified_rate", "Untreated_rate"]):
        """Returns slice of profile with additional ID column "Sample" for easy
        use with sns.boxplot

        Args:
            sample (int, optional): x-value for this sample's boxplot.
                Defaults to 1.
            cols (list, optional): list of columns from profile to use.
                Defaults to ["Modified_rate", "Untreated_rate"].

        Returns:
            DataFrame: slice of profile formatted for sns.boxplot
        """
        data = self.profile[cols].copy()
        data = data.assign(Sample=sample)
        return data

    def plot_boxplot(self, ax, other_samples=[]):
        """Plots a boxplot for self and other samples on the given axis.

        Args:
            ax (pyplot axis): axis on which to add plot
            other_samples (list, optional): list of other MaP Samples to plot.
                Defaults to empty list.
        """
        data = [self.get_boxplot_data(sample=1)]
        xticklabels = [self.sample]
        for i, sample in enumerate(other_samples):
            data.append(sample.get_boxplot_data(sample=i+2))
            xticklabels.append(sample.sample)
        data = pd.concat(data)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        data.head()
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v')
        ax.set(yscale='log',
               ylim=(0.0005, 0.5),
               ylabel="Mutation Rate",
               xticklabels=xticklabels)

###############################################################################
# CONSTRUCTION ZONE
# Heatmap plotting functions
# TODO: This code needs to be reworked.
#     plotCorrs2D
#     getBitmapLegendHandles
#     getCmaps
#     get_bitmap
#
# ! Example usage, because you will never figure this out on your own.
# # Create the plot
# fig, ax = plt.subplots(1, figsize=(10, 10))
# # Make and plot the bitmaps
# ct = example.get_bitmap('ct')
# mask = example.get_bitmap('fasta')
# pairmap = example.get_bitmap('rings')
# example.plotCorrs2D(axis=ax, pairmap=pairmap, ct=ct, mask=mask)
# # Get legend handles and create legend and title
# handles = example.getBitmapLegendHandles()
# handles = [handles['Primary'], handles['Secondary'],
#            handles['Helices'], handles['No Data']]
# ax.legend(handles=handles, loc='lower left')
# ax.set_title("Example PAIR-MaP", y=1.07, fontsize=40)
###############################################################################

    def plotCorrs2D(self, axis, pairmap=0, allcorrs=0, bg_corrs=0, ct=0,
                    mask=0):
        ct_cmap, mask_cmap, pm_cmap = self.getCmaps()
        if type(allcorrs) is np.ndarray:
            axis.imshow(allcorrs, cmap='bwr')
        if type(bg_corrs) is np.ndarray:
            axis.imshow(bg_corrs, cmap='gray')
        if type(pairmap) is np.ndarray:
            axis.imshow(pairmap, cmap=pm_cmap)
        if type(ct) is np.ndarray:
            axis.imshow(ct, cmap=ct_cmap, interpolation='bicubic')
        if type(mask) is np.ndarray:
            axis.imshow(mask, cmap=mask_cmap)

    def getBitmapLegendHandles(self):
        handles = {}
        handles['Primary'] = mp.patches.Patch(color='red', label='Primary')
        handles['Secondary'] = mp.patches.Patch(
            color='blue', label='Secondary')
        handles['No Data'] = mp.patches.Patch(color='green', label='No Data',
                                              alpha=0.2)
        handles['Helices'] = mp.patches.Patch(edgecolor='green', fill=False,
                                              label='Known Helices')
        handles['BG Correlations'] = mp.patches.Patch(color='black',
                                                      label='BG correlations')
        return handles

    def getCmaps(self):
        cmap = plt.get_cmap('Greens')
        ct_cmap = cmap(np.arange(cmap.N))
        ct_cmap[:, -1] = np.concatenate((np.linspace(0, 1, cmap.N/2),
                                         np.linspace(1, 0, cmap.N/2)),
                                        axis=None)
        ct_cmap = mp.colors.ListedColormap(ct_cmap)

        mask_cmap = cmap(np.arange(cmap.N))
        mask_cmap[:, -1] = np.linspace(0, 0.2, cmap.N)
        mask_cmap = mp.colors.ListedColormap(mask_cmap)

        N = 256
        vals = np.ones((N, 4))
        vals[:, 0] = np.concatenate((np.linspace(1, 0, N/2),
                                     np.linspace(0, 1, N/2)), axis=None)
        vals[:, 1] = np.concatenate((np.linspace(1, 0, N/2),
                                     np.linspace(0, 0, N/2)), axis=None)
        vals[:, 2] = np.concatenate((np.linspace(1, 1, N/2),
                                     np.linspace(1, 0, N/2)), axis=None)
        pm_cmap = mp.colors.ListedColormap(vals)
        return ct_cmap, mask_cmap, pm_cmap

    def get_bitmap(self, type, window=1, limit=200):
        valid_types = ['fasta', 'ct', 'rings', 'log', 'pairs']
        assert type in valid_types, f'Please provide valid type: {valid_types}'

        size = self.profile_length
        if type == 'fasta':
            seq = self.profile['Sequence'].copy()
            bitmap = np.zeros((size, size))
            for i in range(size):
                for j in range(size):
                    if (i > j) or seq[i].islower() or seq[j].islower():
                        bitmap[i, j] = 1
        elif type == 'ct':
            bitmap = np.ones((size, size))
            basepairs = self.ct.pairList()
            for x in [-1, 0, 1]:
                for y in [-1, 0, 1]:
                    for i, j in basepairs:
                        bitmap[i+x-1, j+y-1] = 0
        elif type == 'rings':
            bitmap = np.zeros((size, size))
            corrs = self.rings
            for i in range(window*window):
                x = i % window
                y = int(i/window)
                bitmap[corrs.i+x, corrs.j+y] += corrs.Statistic*corrs['+/-']
            bitmap[bitmap > limit] = limit
            bitmap[bitmap < -limit] = -limit
            bitmap = bitmap/limit
            bitmap[0, 0] = 1
            bitmap[0, 1] = -1
            bitmap += 1
            bitmap /= 2
        elif type == 'log':
            bitmap = np.ones((size, size))
            with open(type) as f:
                for line in f.readlines():
                    if line.startswith('Pair '):
                        pair = line.split()[1].strip('()').split(',')
                        bitmap[pair[0], pair[1]] = 0
        elif type == 'pairs':
            bitmap = np.zeros((size, size))
            pairs = self.pairs
            primary = pairs[pairs['Class'] == 1]
            secondary = pairs[pairs['Class'] == 2]
            for i in range(9):
                x = i % 3
                y = int(i/3)
                bitmap[secondary['i']+x, secondary['j']+y] = 0.5
            for i in range(9):
                x = i % 3
                y = int(i/3)
                bitmap[primary['i']+x, primary['j']+y] = 1
        return bitmap

###############################################################################
# filters, clipping, paddings, and coloring functions
# TODO:      filter_pairs
#   get_clip_pad
#   filter_rings
#   get_colorby_profile
#   get_colorby_sequence
#   get_colorby_position
#   Future:
#       filter_ct?
#       filter_profile?
###############################################################################

    def get_alignment(self, fit_this, to_that):
        from Bio.pairwise2 import align, format_alignment
        this_sequence = getattr(self, f"{fit_this}_sequence").upper()
        that_sequence = getattr(self, f"{to_that}_sequence").upper()
        alignment = align.globalxs(this_sequence, that_sequence, -1, -0.1,
                                   penalize_end_gaps=False)
        print(format_alignment(*alignment[0]))

    def get_clip_pad(self, fit_this, to_that):
        attribute = f"{fit_this}_{to_that}_clip_pad"
        if hasattr(self, attribute):
            return getattr(self, attribute)
        this_sequence = getattr(self, f"{fit_this}_sequence")
        that_sequence = getattr(self, f"{to_that}_sequence")
        this_length = getattr(self, f"{fit_this}_length")
        that_length = getattr(self, f"{to_that}_length")
        # find the best match
        str1, str2 = this_sequence, that_sequence
        if len(str1) < len(str2):
            str1, str2 = str2, str1
        len1 = len(str1)
        len2 = len(str2)
        scores = []
        for i in range(len1 - len2 + 1):
            substr1 = str1[i:]
            count = 0
            for nt in range(len2):
                if substr1[nt] == str2[nt]:
                    count += 1
            score = count/len2
            scores.append(score)
        best_match = max(scores)
        good_match = 0.9
        bad_match_message = f"{fit_this}:{to_that} sequences don't match well."
        assert best_match > good_match, bad_match_message
        # assign clip and pad and return
        index = scores.index(best_match)
        if this_length == that_length:
            clip = [0, that_length]
            pad = [0, 0]
        elif this_length > that_length:
            clip = [index, index + that_length]
            pad = [0, 0]
        else:
            clip = [0, that_length]
            pad = [index-1, that_length - this_length - index + 1]
        setattr(self, attribute, [clip, pad])
        return [clip, pad]

    def filter_rings(self, type, cdAbove=None, cdBelow=None, statistic=None,
                     zscore=None, ct=None):
        clip_pad = self.get_clip_pad("profile", type)
        self.rings['mask'] = np.ones(len(self.rings), dtype=bool)
        self.rings['i_offset'] = self.rings['i'].copy
        self.rings['j_offset'] = self.rings['j'].copy
        if clip_pad is not None:
            start = clip_pad[0][0]
            end = clip_pad[0][1]
            i = self.rings["i"].values
            j = self.rings["j"].values
            iinrange = (start < i) & (i < end)
            jinrange = (start < j) & (j < end)
            self.rings['mask'] = self.rings['mask'] & iinrange & jinrange
            offset = clip_pad[1][0]-start
            self.rings['i_offset'] = i+offset
            self.rings['j_offset'] = j+offset
        if cdAbove is not None or cdBelow is not None:
            if ct == 'ss':
                basepairs = [tuple(pair) for pair in self.basepairs]
                ct = rna.CT()
                ct.pair2CT(pairs=basepairs, seq=self.ss_sequence)
            elif ct == 'ct':
                ct = self.ct
            cd_mask = []
            for mask, i, j in zip(self.rings["i_offset"],
                                  self.rings["j_offset"],
                                  self.rings['mask']):
                if isinstance(cdAbove, int):
                    mask.append(ct.contactDistance(i, j) > cdAbove)
                if isinstance(cdBelow, int):
                    mask.append(ct.contactDistance(i, j) < cdBelow)
            cd_mask = np.array(cd_mask, dtype=bool)
            self.rings['mask'] = self.rings['mask'] & cd_mask
        if statistic is not None:
            mask = self.rings["Statistic"] > statistic
            self.rings['mask'] = self.rings['mask'] & mask
        if zscore is not None:
            mask = self.rings["Zij"] > zscore
            self.rings['mask'] = self.rings['mask'] & mask

    def get_colorby_profile(self, type):
        clip_pad = self.get_clip_pad("profile", type)
        if clip_pad is not None:
            start, end = clip_pad[0]
        else:
            start, end = 0, -1
        cmap = ['grey', 'black', 'orange', 'red']
        bins = [0, 0.4, 0.85]
        profcolors = []
        for x in self.profile["Norm_profile"][start:end]:
            profcolors.append(sum([b < x for b in bins]))
        if clip_pad is not None:
            profcolors = ['grey']*clip_pad[1][0] + profcolors
            profcolors += ['grey']*clip_pad[1][1]
        colors = np.array([cmap[val] for val in profcolors])
        return colors

    def get_colorby_sequence(self, type, colors='new'):
        nuc_colors = {"old": {"A": "#f20000", "U": "#f28f00",
                              "G": "#00509d", "C": "#00c200"},
                      "new": {"A": "#366ef0", "U": "#9bb9ff",
                              "G": "#f04c4c", "C": "#ffa77c"}}[colors]
        seq = getattr(self, f"{type}_sequence")
        colors = np.array([nuc_colors[nuc.upper()] for nuc in seq])
        return colors

    def get_colorby_position(self, type, cmap='rainbow'):
        length = getattr(self, f"{type}_length")
        cmap = plt.get_cmap(cmap)
        colors = np.array([cmap(n/length) for n in range(length)])
        return colors

###############################################################################
# Plotting functions that accept a list of samples
#   array_ap
#   array_ss
#
###############################################################################


def array_ap(samples=[], rows_cols=None, **kwargs):
    if rows_cols is None:
        if len(samples) > 9:
            rows = 3
            columns = round(len(samples)/3)
        elif len(samples) > 2:
            rows = 2
            columns = round(len(samples)/2)
        else:
            rows = 1
            columns = len(samples)
    fig = plt.figure(figsize=samples[0].get_ap_figsize(rows, columns))
    rcn = rows*100 + columns*10
    axes = [fig.add_subplot(rcn+i+1) for i in range(len(samples))]
    for i, sample in enumerate(samples):
        sample.make_ap(ax=axes[i], **kwargs)


def array_qc(samples=[]):
    fig = plt.figure(figsize=(14, 14))
    ax1 = fig.add_subplot(421)
    ax2 = fig.add_subplot(422)
    ax3 = fig.add_subplot(423)
    ax4 = fig.add_subplot(424)
    ax5 = fig.add_subplot(212)
    num_samples = len(samples)
    for i, sample in enumerate(samples):
        sample.plot_log_MutsPerMol(ax1, sample="Untreated")
        sample.plot_log_MutsPerMol(ax2, sample="Modified")
        sample.plot_log_ReadLength(
            ax3, sample="Untreated", n=i, of=num_samples)
        sample.plot_log_ReadLength(ax4, sample="Modified", n=i, of=num_samples)
    samples[0].set_log_MutsPerMol(ax1)
    samples[0].set_log_MutsPerMol(ax2)
    samples[0].set_log_ReadLength(ax3)
    samples[0].set_log_ReadLength(ax4)
    samples[0].plot_boxplot(ax5, samples[1:])
    ax5.set_title("boxplot")
    plt.tight_layout()


def array_skyline(samples, **kwargs):
    fig, ax = plt.subplots(1, figsize=samples[0].get_skyline_figsize(1, 1))
    for sample in samples[:-1]:
        sample.plot_skyline(ax)
    samples[-1].make_skyline(ax, **kwargs)


def array_ss(samples, **kwargs):
    fig, ax = plt.subplots(1, 4, figsize=samples[0].get_ss_figsize(1, 4))
    for i, sample in enumerate(samples):
        sample.make_ss(ax[i], **kwargs)

###############################################################################
# Parsing arguments
#
#
###############################################################################
