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
import math
import os.path

# Special python packages
try:
    import Bio.PDB
    from Bio import SeqIO
except ModuleNotFoundError:
    print("Biopython package is missing. Parsing PDB or deletions won't work.")
try:
    import py3Dmol
except ModuleNotFoundError:
    print("py3Dmol package is missing. Plotting 3D structures will not work.")

# scripts in JNBTools
import RNAtools3 as rna


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
    '#0092edff',  # Blue
    '#ff8300ff',  # Orange
    '#a100ffff',  # Purple
    '#edc600ff',  # Yellow
    '#ff48e9ff',  # Pink
    '#3fd125ff'  # Green
]
sns.set_palette(colors)

# Hard coded defaults
###############################################################################


def get_nt_color(nt, colors="new"):
    nt_color = {"old": {"A": "#f20000",  # red
                        "U": "#f28f00",  # yellow
                        "G": "#00509d",  # blue
                        "C": "#00c200"},  # green
                "new": {"A": "#366ef0",  # blue
                        "U": "#9bb9ff",  # light blue
                        "G": "#f04c4c",  # red
                        "C": "#ffa77c"}  # light red
                }[colors][nt]
    return nt_color


def get_default_fill(metric):
    fill = {'Class': -1,
            'Statistic': 0.0,
            'Zij': 0.0,
            'Metric': 0.0,
            'Distance': 1000,
            'Percentile': 0.0}[metric]
    return fill


def get_default_cmap(metric):
    cmap = {'Class': mp.colors.ListedColormap([[0.3, 0.3, 0.3, 0.2],
                                               [0.0, 0.0, 0.95, 0.6],
                                               [0.12, 0.76, 1.0, 0.6]]),
            'Statistic': 'bwr',
            'Zij': 'bwr',
            'Metric': 'YlGnBu',
            'Distance': 'jet',
            'Percentile': 'YlGnBu',
            'Probability': 'rainbow_r'
            }[metric]
    cmap = plt.get_cmap(cmap)
    cmap = cmap(np.arange(cmap.N))
    cmap[:, -1] = np.full((len(cmap)), 0.6)  # set default alpha to 0.6
    if metric == 'Class':
        cmap[0, -1] = 0.2  # alpha of non primary and secondary pairs to 0.2
    if metric == 'Distance':
        # set color of max distance and no data distances to gray
        cmap[-1, :] = np.array([80/255., 80/255., 80/255., 0.2])
    cmap = mp.colors.ListedColormap(cmap)
    return cmap


def get_default_metric(ij_data):
    metric = {'rings': 'Statistic',
              'pairs': 'Class',
              'deletions': 'Percentile',
              'probs': 'Probability'
              }[ij_data]
    return metric


def get_default_min_max(metric):
    min_max = {'Percentile': [0.98, 1.0],
               "Statistic": [-100, 100],
               'Zij': [-50, 50],
               "Class": [0, 2],
               "Metric": [0, 0.001],
               'Distance': [10, 80],
               'Probability': [0, 1]
               }[metric]
    return min_max


def view_colormap(ij_data=None, metric=None, ticks=None, values=None,
                  title=None, cmap=None):
    """Given an ij_data (ij data) will display a colorbar for the default
    values (metric, cmap, min_max).

    Args:
        ij_data (str, optional): string matching an ij data type.
            Options are "rings", "pairs" or "deletions".
            Defaults to None.
        metric (str, optional): string matching column name of ij data.
            Default determined by get_default_metric.
        ticks (list, optional): locations to add ticks. scale is 0-10.
            Defaults to [0.5, 0.95], or [10/6, 30/6/, 50/6] for "Class" metric.
        title (str, optional): string for title of colorbar.
            Defaults to "{ij_data}: {metric}"
        cmap (str, optional): string matching a valid matplotlib colormap.
            Default determined by get_default_cmap.
    """
    if metric is None:
        metric = get_default_metric(ij_data)
    if ticks is None:
        if metric == "Class":
            ticks = [10/6, 30/6, 50/6]
        else:
            ticks = [0, 2, 4, 6, 8, 10]
    if values is None:
        if metric == "Class":
            values = ['Complementary', 'Primary', 'Secondary']
        else:
            mn, mx = get_default_min_max(metric)
            values = [f"{mn + ((mx-mn)/5)*i:.2}" for i in range(6)]
    if title is None:
        title = f"{ij_data.capitalize()}: {metric.lower()}"
    if cmap is None:
        cmap = get_default_cmap(metric)
    else:
        cmap = plt.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))

    _, ax = plt.subplots(1, figsize=(6, 2))
    ax.imshow([colors], extent=[0, 10, 0, 1])
    ax.set_title(title)
    ax.set_xticks(ticks)
    ax.set_xticklabels(values)
    ax.set_yticks([])


# COPIED FROM SHAPEMAPPER2
# some of this might be inappropriately applied to all plots
# TODO: look into passing dict to mp.rc()
###############################################################################
mp.rcParams["font.sans-serif"].insert(0, "Arial")
shapemapper_style = {"font.family": "sans-serif",
                     "pdf.fonttype": 42,
                     # use TrueType fonts when exporting PDFs
                     # (embeds most fonts - this is especially
                     #  useful when opening in Adobe Illustrator)
                     'xtick.direction': 'out',
                     'ytick.direction': 'out',
                     'legend.fontsize': 14,
                     'grid.color': ".8",
                     'grid.linestyle': '-',
                     'grid.linewidth': 1}

rx_color = "red"
bg_color = "blue"
dc_color = "darkgoldenrod"
###############################################################################


class Sample():

    def __init__(self,
                 sample=None,
                 fasta=None,
                 profile=None,
                 ct=None,
                 compct=None,
                 ss=None,
                 log=None,
                 rings=None,
                 deletions=None,
                 pairs=None,
                 pdb=None,
                 probs=None,
                 dance_prefix=None):
        self.paths = {"fasta": fasta,
                      "profile": profile,
                      "ct": ct,
                      "comptct": compct,
                      "ss": ss,
                      "log": log,
                      "rings": rings,
                      "deletions": deletions,
                      "pairs": pairs,
                      "pdb": pdb,
                      "probs": probs}
        self.sample = sample

        self.sequence = {}  # stores sequences:
        # profile == rings == pairs == dance
        # fasta == deletions
        # ct, ss, 3d can be different
        self.length = {}  # stores length of all sequences
        self.ij_data = {}  # stores rings, pairs, deletions as dataFrames
        self.window = {}  # each ij_data value needs a window
        self.header = {}  # each ij_data value needs a header
        self.profile = {}  # TODO: add shape, frag, rnp as key:value pairs
        self.structure = {"ct": {},
                          "ss": {},
                          "pdb": {}}  # add ct's, ss's, pdb's
        if profile is not None:
            self.read_profile(profile)
        if ct is not None:
            self.read_ct("ct", ct)
        if compct is not None:
            self.read_ct("compct", compct)
        if rings is not None:
            assert hasattr(self, "profile"), "Rings plotting requires profile."
            self.read_rings(rings)
        if pairs is not None:
            assert hasattr(self, "profile"), "Pairs plotting requires profile."
            self.read_pairs(pairs)
        if probs is not None:
            assert hasattr(self, "profile"), "Probabilities require profile."
            self.read_probs(probs)
        if deletions is not None:
            assert fasta is not None, "Deletions plotting requires fasta"
            self.read_deletions(deletions, fasta)
        if ss is not None:
            self.read_ss("ss", ss)
        if log is not None:
            self.read_log(log)
        if pdb is not None:
            self.read_pdb(pdb)
        if dance_prefix is not None:
            self.init_dance(dance_prefix)

###############################################################################
# Parsing data files
# TODO: make read_nsd function less sucky. I think it's YAML.
#   structure files
#       read_pdb
#       read_ss          Note: nsd portion is kludgey AF should use yaml
#   MaP data files
#       read_rings
#       read_deletions
#       read_pairs
#       read_log
#       read_dance_reactivities
###############################################################################

    def read_ct(self, ct_name, ct):
        setattr(self, ct_name, rna.CT(ct))
        self.sequence[ct_name] = ''.join(self.ct.seq)
        self.length[ct_name] = len(self.sequence["ct"])

    def read_pdb(self, pdb):
        parser = Bio.PDB.PDBParser()
        self.pdb = parser.get_structure('RNA', pdb)
        self.sequence["pdb"] = ''
        with open(pdb) as file:
            for line in file.readlines():
                line = [field.strip() for field in line.split()]
                if line[0] == "SEQRES":
                    self.sequence["pdb"] += ''.join(line[4:])
        self.pdb_validres = []
        for res in self.pdb[0]["A"].get_residues():
            res_id = res.get_id()
            if res_id[0] == " ":
                self.pdb_validres.append(res_id[1])
        self.length["pdb"] = len(self.sequence["pdb"])

    def read_ss(self, ss_name, ss):
        # get and check file extension, then read
        self.ss_type = ss.split('.')[-1].lower()
        valid_type = self.ss_type in ['xrna', 'varna', 'nsd', 'cte']
        message = f"stucture file type {self.ss_type} not supported"
        assert valid_type, message
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        if self.ss_type == "xrna":
            tree = xmlet.parse(ss)
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
                helix_list = [(i_outter+nt, j_outter-nt)
                              for nt in range(length)]
                basepairs.extend(helix_list)
            # make expected arrays
            xcoords = np.array(xcoords)
            ycoords = np.array(ycoords)
        elif self.ss_type == "varna":
            tree = xmlet.parse(ss)
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
            xcoords = np.array(xcoords)
            ycoords = np.array(ycoords)
        elif self.ss_type == "cte":
            names = ['nuc', 'seq', 'pair', 'xcoords', 'ycoords']
            ct = pd.read_csv(ss, sep=r'\s+', usecols=[0, 1, 4, 8, 10],
                             names=names, header=0)
            sequence = ''.join(list(ct['seq']))
            xcoords = np.array([float(x) for x in ct.xcoords])
            ycoords = np.array([float(y) for y in ct.ycoords])
            ct = ct[ct.nuc < ct.pair]
            basepairs = [[int(i), int(j)] for i, j in zip(ct.nuc, ct.pair)]
        elif self.ss_type == "nsd":
            basepairs, sequence, xcoords, ycoords = [], '', [], []
            with open(ss, 'r') as file:
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
                        xcoords.append(float(nt[4][1]))
                        ycoords.append(float(nt[5][1]))
                    elif item == "pairs":
                        pairs = line[1].split(":")[1:]
                        basepair = [int(nuc.strip('"')) for nuc in pairs]
                        basepairs.append(basepair)
            xcoords = np.array(xcoords)
            ycoords = np.array(ycoords)
        # store attributes
        coord_scale_factor = {"xrna": 1/20,
                              "varna": 1/65,
                              "cte": 1/30.5,
                              "nsd": 1/30.5}[self.ss_type]
        self.sequence[ss_name] = sequence.upper().replace("T", "U")
        self.length[ss_name] = len(sequence)
        self.basepairs = basepairs
        self.xcoordinates = xcoords*coord_scale_factor
        self.ycoordinates = ycoords*coord_scale_factor

    def read_profile(self, profile):
        self.profile = pd.read_csv(profile, sep='\t')
        sequence = ''.join(self.profile["Sequence"].values)
        self.sequence["profile"] = sequence.upper().replace("T", "U")
        self.length["profile"] = len(self.sequence["profile"])

    def read_probs(self, probs):
        with open(probs, 'r') as file:
            self.header["probs"] = file.readline()
        self.window["probs"] = 1
        data = pd.read_csv(probs, sep='\t', header=1)
        data["Probability"] = 10 ** (-data["-log10(Probability)"])
        self.ij_data["probs"] = data
        self.sequence["probs"] = self.sequence["profile"]
        self.length["probs"] = self.length["profile"]

    def read_rings(self, rings):
        with open(rings, 'r') as file:
            self.header["rings"] = file.readline()
        split_header = self.header["rings"].split('\t')
        window = split_header[1].split('=')[1]
        self.window["rings"] = int(window)
        self.ij_data["rings"] = pd.read_csv(rings, sep='\t', header=1)
        self.sequence["rings"] = self.sequence["profile"]
        self.length["rings"] = self.length["profile"]

    def read_deletions(self, deletions, fasta):
        fasta = list(SeqIO.parse(open(fasta), 'fasta'))
        self.deletions_gene = fasta[0].id
        self.sequence["deletions"] = str(fasta[0].seq).upper().replace("T",
                                                                       "U")
        self.length["deletions"] = len(self.sequence["deletions"])
        column_names = ['Gene', 'i', 'j', 'Metric']
        data = pd.read_csv(deletions, sep='\t', names=column_names, header=0)
        data["Percentile"] = data['Metric'].rank(method='max', pct=True)
        self.ij_data["deletions"] = data
        self.window["deletions"] = 1

    def read_pairs(self, pairs):
        with open(pairs, 'r') as file:
            self.header["pairs"] = file.readline()
        split_header = self.header["pairs"].split('\t')
        window = split_header[1].split('=')[1]
        self.window["pairs"] = int(window)
        self.ij_data["pairs"] = pd.read_csv(pairs, sep='\t', header=1)
        self.sequence["pairs"] = self.sequence["profile"]
        self.length["pairs"] = self.length["profile"]

    def read_log(self, log):
        with open(log, 'r') as f:
            flist = list(f)
            log_format_test = 0
            for i, line in enumerate(flist):
                if line.startswith("  |MutationCounter_Modified"):
                    log_format_test += 1
                    modlength = []
                    for x in flist[i+6:i+27]:
                        modlength.append(float(x.strip().split('\t')[1]))
                    modmuts = []
                    for x in flist[i+32:i+53]:
                        modmuts.append(float(x.strip().split('\t')[1]))
                if line.startswith("  |MutationCounter_Untreated"):
                    log_format_test += 1
                    untlength = []
                    for x in flist[i+6:i+27]:
                        untlength.append(float(x.strip().split('\t')[1]))
                    untmuts = []
                    for x in flist[i+32:i+53]:
                        untmuts.append(float(x.strip().split('\t')[1]))
        message = ("Histogram data missing from log file. Requires" +
                   " --per-read-histogram flag when running ShapeMapper.")
        assert log_format_test >= 2, message
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

    def init_dance(self, prefix):
        reactivityfile = f"{prefix}-reactivities.txt"
        # read in 2 line header
        with open(reactivityfile) as inf:
            header1 = inf.readline().strip().split()
            header2 = inf.readline().strip().split()
        # number of components
        self.dance_components = int(header1[0])
        # population percentage of each component
        self.dance_percents = header2[1:]
        # dance is a list containing one sample for each component
        self.dance = [Sample() for _ in range(self.dance_components)]
        # build column names for reading in BM file
        for i, sample in enumerate(self.dance):
            sample.sample = f"{i} - {self.dance_percents[i]}"
            sample.paths = {"profile": reactivityfile,
                            "rings": f"{prefix}-{i}-rings.txt",
                            "pairs": f"{prefix}-{i}-pairmap.txt",
                            "ct": [f"{prefix}-{i}.f.ct",  # if using --pk
                                   f"{prefix}-{i}.ct"]}  # if regular fold used
            # read in "profile" from reactivities
            colnames = ["Nucleotide", "Sequence", "Reactivity_profile",
                        "Reactivity", "Background"]
            col_offset = 3 * i
            last_col = 3 * self.dance_components + 2
            columns = [0, 1, 2 + col_offset, 3 + col_offset, last_col]
            sample.profile = pd.read_csv(reactivityfile, sep='\t', header=2,
                                         names=colnames, usecols=columns)
            sample.sequence["profile"] = self.sequence['profile']
            sample.length["profile"] = len(sample.sequence["profile"])
            # read in other attributes
            if os.path.isfile(sample.paths["rings"]):
                sample.read_rings(sample.paths["rings"])
            if os.path.isfile(sample.paths["pairs"]):
                sample.read_pairs(sample.paths["pairs"])
            # possible that these both exist
            for ct_file in sample.paths["ct"]:
                if os.path.isfile(ct_file):
                    sample.read_ct("ct", ct_file)
                    sample.paths["ct"] = ct_file

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
        ax_width = self.length["profile"] * 0.1
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

    def plot_sequence(self, axis, sequence, yvalue=0.005):
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
        sequence = self.sequence[sequence]
        for i, seq in enumerate(sequence):
            col = color_dict[seq.upper()]
            axis.annotate(seq, xy=(i + 1, yvalue), xycoords='data',
                          fontproperties=font_prop,
                          color=col, horizontalalignment="center")

    def set_skyline(self, ax, title="Raw Reactivity Profile"):
        ax.set(title=title,
               xlim=[0, self.length["profile"]],
               xticks=range(0, self.length["profile"], 20))
        ax.set_xticks(range(0, self.length["profile"], 5), minor=True)
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
        self.plot_sequence(ax, "profile")
        self.set_skyline(ax, title=column.replace("_", " "))

    def make_dance_skyline(self, ax=None):
        """Creates a skyline figure representing the component reactivities of
        the ensemble.
        """
        if ax is None:
            fig, ax = plt.subplots(1, figsize=self.get_skyline_figsize(1, 1))
        for i, sample in enumerate(self.dance):
            sample.plot_skyline(ax, label=f"{i} - {self.dance_percents[i]}")
        self.plot_sequence(ax, sequence="profile")
        ax.legend(title="Component: Population", loc=1)
        ax.set(xlim=(0, self.length["profile"]),
               title=f"{self.sample}: DANCE-MaP Reactivities",
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
#     get_boxplot_data
#     plot_boxplot
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
# Arc Plot plotting functions
#     add_arc
#     get_ap_figsize
#     set_ap
#     plot_ap_ct
#     plot_ap_ctcompare
#     plot_ap_profile
#     plot_ap_data
#     make_ap
#     Future:
#         plot_ap_probabilities
###############################################################################

    def add_arc(self, ax, i, j, window, color, panel):
        """internal function used to add arcs based on data

        Args:
            ax (pyplot axis): axis to which arc will be added
            i (int): leftmost position of left side
            j (int): leftmost position of right side
            window (int): number of nucleotides included in arc
            color (color): color of arc
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
                               width=window, ec='none')
        ax.add_patch(arc)

    def get_ap_figsize(self, rows, cols, sequence="ct"):
        """Pass this function call to figsize within pyplot.subplots() to set
        an appropriate figure width and height for arc plots.

        Args:
            rows (int): number of rows in pyplot figure
            cols (int): number of columns in pyplot figure

        Returns:
            tuple: (width, height) appropriate figsize for pyplot figure
        """
        dim = len(self.sequence[sequence]) * 0.1 + 1
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
        ax.set(xlim=(0, self.length["ct"]),
               ylim=(-height-5, height+1))

    def plot_ap_ct(self, ax, panel='top'):
        """Adds arc plot representation of the ct structure to given axis

        Args:
            ax (pyplot axis): axis to which structure will be added
            ctpath (str): path to ct file
        """
        ct_pairs = self.ct.pairList()
        for i, j in ct_pairs:
            self.add_arc(ax, i, j, 1, (0.1, 0.1, 0.1, 0.7), panel)

    def plot_ap_ctcompare(self, ax, panel='top'):
        """Adds structure comparison arc plot to the given axis

        Args:
            ax (pyplot axis): axis to which structure comparison is added
        """
        ct1 = set(self.ct.pairList())
        ct2 = set(self.compct.pairList())
        shared = ct1.union(ct2)
        ref = ct1.difference(ct2)
        comp = ct2.difference(ct1)
        sharedcolor = (150/255., 150/255., 150/255., 0.7)
        refcolor = (38/255., 202/255., 145/255., 0.7)
        compcolor = (153/255., 0.0, 1.0, 0.7)
        for i, j in comp:
            self.add_arc(ax, i, j, 1, compcolor, panel)
        for i, j in ref:
            self.add_arc(ax, i, j, 1, refcolor, panel)
        for i, j in shared:
            self.add_arc(ax, i, j, 1, sharedcolor, panel)
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

    def plot_ap_data(self, ax, ij_data, metric=None, **kwargs):
        self.filter_ij_data(ij_data, "ct", **kwargs)
        ij_colors = self.get_ij_colors(ij_data, metric)
        window = self.window[ij_data]
        for i, j, color in zip(*ij_colors):
            self.add_arc(ax, i, j, window, color, "bottom")

        # TODO: create invisible axis and create colorbar
        # from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        # axins1 = inset_axes(ax, width="5%", height="5%", loc='lower left')
        # im1 = axins1.imshow([min_max, min_max], cmap=cmap)
        # axins2 = inset_axes(ax, width="5%", height="50%", loc='lower right')
        # fig.colorbar(im1, cax=axins2, ticks=min_max)
        # axins2.set_title(metric)
        # axins1.set_visible(False)

    def make_ap(self, ax=None, ij_data=None, metric=None, ctcompare=True,
                **kwargs):
        """Creates figure with arc plot. Includes ct, profile, sequence, and
        data type passed to data. [rings, pairs, probabilities]

        Args:
            data (str, optional): datatype for bottom panel.
                    options: 'rings', 'pairs', 'probabilities'.
                    Defaults to None.
        """
        if ax is None:
            figsize = self.get_ap_figsize(1, 1)
            fig, ax = plt.subplots(1, figsize=figsize)
        self.set_ap(ax)
        if hasattr(self, "compct") and ctcompare is True:
            self.plot_ap_ctcompare(ax)
        else:
            self.plot_ap_ct(ax)
        if hasattr(self, "profile"):
            self.plot_ap_profile(ax)
        if metric == "Distance":
            self.set_3d_distances(ij_data)
        if ij_data is not None:
            self.plot_ap_data(ax, ij_data, metric, **kwargs)
        self.plot_sequence(ax, "ct", yvalue=0.5)
        ax.annotate(self.sample, xy=(0.1, 0.9),
                    xycoords="axes fraction", fontsize=60, ha='center')
        return ax

    def make_dance_ap(self, ij_data=None, metric=None, **kwargs):
        rows, cols = get_rows_columns(self.dance_components)
        figsize = self.dance[0].get_ap_figsize(rows, cols)
        _, axes = plt.subplots(rows, cols, figsize=figsize, squeeze=False)
        for i, dance in enumerate(self.dance):
            row = i // cols
            col = i % rows
            ax = axes[col, row]
            self.set_ap(ax)
            dance.plot_ap_ct(ax)
            if hasattr(dance, "profile"):
                dance.plot_ap_profile(ax)
            if ij_data is not None:
                dance.plot_ap_data(ax, ij_data, metric, **kwargs)
            dance.plot_sequence(ax, "ct", yvalue=0.5)
            ax.annotate(dance.sample, xy=(0.1, 0.9),
                        xycoords="axes fraction", fontsize=60, ha='center')

###############################################################################
# Secondary Structure graph plotting functions
#     get_ss_figsize
#     set_ss
#     plot_ss_structure
#     get_ss_seqcolors
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
        if isinstance(colorby, list) and len(colorby) == self.length["ss"]:
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
                mask = [nt == nuc for nt in self.sequence["ss"]]
                xcoords = self.xcoordinates[mask]
                ycoords = self.ycoordinates[mask]
                marker = "$"+nuc+"}$"
                marker = markers
                ax.scatter(xcoords, ycoords, marker=marker,
                           c=colors[mask])
        else:
            ax.scatter(self.xcoordinates, self.ycoordinates, marker=markers,
                       c=colors)

    def add_ss_lines(self, ax, i, j, color, linewidth=1.5):
        xi = self.xcoordinates[i-1]
        yi = self.ycoordinates[i-1]
        xj = self.xcoordinates[j-1]
        yj = self.ycoordinates[j-1]
        ax.plot([xi, xj], [yi, yj], color=color, linewidth=linewidth)

    def plot_ss_data(self, ax, ij_data, metric=None, all_pairs=False,
                     **kwargs):
        self.filter_ij_data(ij_data, "ss", all_pairs=all_pairs, **kwargs)
        ij_colors = self.get_ij_colors(ij_data, metric)
        window = self.window[ij_data]
        for i, j, color in zip(*ij_colors):
            if window == 1:
                self.add_ss_lines(ax, i, j, color)
            else:
                for w in range(window):
                    self.add_ss_lines(ax, i+w, j+window-1 -
                                      w, color, linewidth=6)

    def plot_ss_positions(self, ax, spacing=20):
        """adds position labels to the secondary structure graph

        Args:
            ax (pyplot axis): axis on which to add labels
            spacing (int, optional): labels every nth nucleotide.
                    Defaults to 20.
        """
        for i in range(0, self.length["ss"], spacing):
            ax.annotate(i+1, xy=(self.xcoordinates[i], self.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords='data', fontsize=16,
                        bbox=dict(fc="white", ec='none', pad=0.3))

    def make_ss(self, ax=None, ij_data=None, metric=None, positions=True,
                colorby=None, markers='o', all_pairs=False, **kwargs):
        """Creates a full secondary structure graph with data plotted

        Args:
            ax (pyplot axis, option): axis to plot
            ij_data (str, optional): string matching the ij_data to plot.
                    Options: rings, pairs, deletions.
            metric (str, optional): string matching a column name of the
                    ij_data dataframe. Defaults defined at top of file.
            positions (bool, optional): whether to label positions.
                    Defaults to True.
            colorby (str, optional): method used to color nucleotides.
                    Options: sequence, profile, positions
                    Defaults to 'profile',
                    If profile is missing, sequence is used.
            markers (str, optional): marker type to use for nucleotides.
                    Defaults to 'o'. (dots)
        """
        if ax is None:
            _, ax = plt.subplots(1, figsize=self.get_ss_figsize(1, 1))
        self.set_ss(ax)
        self.plot_ss_structure(ax)
        if colorby is None:
            if hasattr(self, "profile"):
                colorby = "profile"
            else:
                colorby = "sequence"
        self.plot_ss_sequence(ax, colorby=colorby, markers=markers)
        if positions is True:
            self.plot_ss_positions(ax)
        if metric == "Distance":
            self.set_3d_distances(ij_data)
        if ij_data is not None:
            self.plot_ss_data(ax, ij_data, metric, all_pairs, **kwargs)
        return ax

###############################################################################
# Standard ShapeMapper plotting functions
#     plot_sm_profile
#     plot_sm_depth
#     plot_sm_rates
#     metric_abbreviate
#     make_shapemapper
###############################################################################

    def plot_sm_profile(self, axis=None):
        """Plots classic ShapeMapper normalized reactivity on the given axis

        Args:
            axis (pyplot axis): axis on which to add plot
        """
        if axis is None:
            _, axis = plt.subplots(figsize=self.get_skyline_figsize(1, 1))
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
        axis.set_xlim(1, self.length["profile"])
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
        axis.set(yticks=yticks,
                 yticklabels=stripped_ticks)

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
        self.plot_sequence(axis, "profile")

    def plot_sm_depth(self, axis=None):
        """Plots classic ShapeMapper read depth on the given axis

        Args:
            axis (pyplot axis): axis on which to add plot
        """
        if axis is None:
            _, axis = plt.subplots(figsize=self.get_skyline_figsize(1, 1))
        sample = self.profile
        axis.plot(sample['Nucleotide'], sample['Modified_read_depth'],
                  linewidth=1.5, color=rx_color, alpha=1.0, label="Modified")
        axis.plot(sample['Nucleotide'], sample['Untreated_read_depth'],
                  linewidth=1.5, color=bg_color, alpha=1.0, label="Untreated")
        axis.plot(sample['Nucleotide'], sample['Denatured_read_depth'],
                  linewidth=1.5, color=dc_color, alpha=1.0, label="Denatured")
        axis.set_xlim(1, self.length["profile"])
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
        formatted_ticks = [self.metric_abbreviate(val) for val in yticks]
        axis.set(yticks=yticks, yticklabels=formatted_ticks)
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

    def plot_sm_rates(self, axis=None):
        """Plots classic ShapeMapper mutation rates on the given axis

        Args:
            axis (pyplot axis): axis on which to add plot
        """
        if axis is None:
            _, axis = plt.subplots(figsize=self.get_skyline_figsize(1, 1))
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
        yticks = axis.get_yticks()
        yticklabels = [str(int(x*100)) for x in yticks]
        axis.set(yticks=yticks, yticklabels=yticklabels)
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
# CONSTRUCTION ZONE: pdb figures
##############################################################################

    def get_xyz_coord(self, nt, atom="O2'"):
        xyz = [float(c) for c in self.pdb[0]["A"][int(nt)][atom].get_coord()]
        return xyz

    def get_3d_distance(self, i, j):
        valid = [nt.get_id()[1] for nt in self.pdb[0]["A"].get_residues()]
        if i in valid and j in valid:
            xi, yi, zi = self.get_xyz_coord(i)
            xj, yj, zj = self.get_xyz_coord(j)
            distance = ((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)**0.5
        else:
            distance = 1000
        return distance

    def set_3d_distances(self, ij_data):
        self.filter_ij_data(ij_data, "pdb")
        data = self.ij_data[ij_data].copy()
        if "Distance" not in data:
            distances = []
            for _, i, j in data[["i_offset", "j_offset"]].itertuples():
                distances.append(self.get_3d_distance(i, j))
            data["Distance"] = distances
        self.ij_data[ij_data] = data

    def add_3d_lines(self, view, i, j, color, viewer=None):
        xi, yi, zi = self.get_xyz_coord(i)
        xj, yj, zj = self.get_xyz_coord(j)
        cylinder_specs = {"start": {"x": xi, "y": yi, "z": zi},
                          "end":  {"x": xj, "y": yj, "z": zj},
                          "radius": 0.5,
                          "fromCap": 2,
                          "toCap": 2,
                          "color": color}
        if viewer is not None:
            view.addCylinder(cylinder_specs, viewer=viewer)
        else:
            view.addCylinder(cylinder_specs)

    def plot_3d_data(self, view, ij_data, metric=None, viewer=None,
                     **kwargs):
        self.filter_ij_data(ij_data, "pdb", **kwargs)
        ij_colors = self.get_ij_colors(ij_data, metric)
        for i, j, color in zip(*ij_colors):
            color = mp.colors.rgb2hex(color)
            color = "0x"+color[1:]
            window = self.window[ij_data]
            for w in range(window):
                io = i+w
                jo = j+window-1-w
                if io in self.pdb_validres and jo in self.pdb_validres:
                    self.add_3d_lines(view, io, jo, color, viewer)

    def set_3d_colors(self, view, colorby, viewer=None):
        colorby_function = getattr(self, "get_colorby_"+colorby)
        colors = colorby_function("pdb")
        color_selector = {}
        for res in self.pdb.get_residues():
            res = res.get_id()
            color = colors[res[1]-1]
            if color in color_selector.keys():
                color_selector[color].append(res[1])
            else:
                color_selector[color] = [res[1]]
        for color in color_selector.keys():
            selector = {'resi': color_selector[color]}
            style = {"cartoon": {"color": color, "opacity": 0.8}}
            if viewer is None:
                view.setStyle(selector, style)
            else:
                view.setStyle(selector, style, viewer=viewer)
        selector = {'resi': self.pdb_validres, 'invert': 'true'}
        style = {"cross": {"hidden": "true"}}
        if viewer is None:
            view.setStyle(selector, style)
        else:
            view.setStyle(selector, style, viewer=viewer)

    def set_3d_view(self, view):
        with open(self.paths["pdb"], 'r') as pdb_file:
            pdb_str = pdb_file.read()
        view.addModel(pdb_str, 'pdb')
        view.setStyle({'chain': 'A'}, {"cartoon": {'color': 'grey'}})
        view.zoomTo()
        return view

    def make_3d(self, view=None, viewer=None, ij_data=None, metric=None,
                colorby="sequence", **kwargs):
        if view is None:
            view = py3Dmol.view()
            view = self.set_3d_view(view)
        if hasattr(self, "profile"):
            colorby = "profile"
        self.set_3d_colors(view, colorby, viewer)
        if metric == "Distance":
            self.set_3d_distances(ij_data)
        if ij_data is not None:
            self.plot_3d_data(view, ij_data, metric, viewer, **kwargs)
        return view

##############################################################################
# CONSTRUCTION ZONE
#     get_pairs_sens_PPV
#     writeCTE
#     plot_regression
##############################################################################

    def get_pairs_sens_PPV(self, ct="ct"):
        "Returns sensitivity and PPV for pair data to the ct structure"
        import pmanalysis as pma
        pm = pma.PairMap(self.paths["pairs"])
        ct = getattr(self, ct).copy()
        ct.filterNC()
        ct.filterSingleton()
        p, s = pm.ppvsens_duplex(ct, ptype=1, exact=False)
        return p, s

    def writeCTE(self, outputPath):
        """writes the current structure out to CTE format for Structure Editor.

        Args:
            outputPath (string): path to output cte file to be created
        """
        pairs = [tuple(pair) for pair in self.basepairs]
        ct = rna.CT()
        ct.pair2CT(pairs=pairs, seq=self.sequence["ss"])
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

    def make_regression(self, comp_sample, ax=None, colorby="ct",
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
            _, ax = plt.subplots(1, 1, figsize=(7, 7))
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
            paired_mask = np.zeros(self.length["profile"], dtype=bool)
            for i in paired_list:
                paired_mask[i] = True
            paired_mask = paired_mask[notNans]
            unpaired_mask = ~paired_mask
            ax.scatter(p1[paired_mask], p2[paired_mask], label="Paired")
            ax.scatter(p1[unpaired_mask], p2[unpaired_mask], label="Unpaired")
        elif colorby == "nucleotide":
            for nuc in "GUAC":
                sequence = self.profile["Sequence"][notNans]
                nuc_mask = [nt == nuc for nt in sequence]
                color = get_nt_color(nuc)
                ax.scatter(p1[nuc_mask], p2[nuc_mask], label=nuc, color=color)
        else:
            ax.scatter(p1, p2)
        s1, s2 = self.sample, comp_sample.sample
        ax.set(xscale='log',
               xlim=[0.00001, 0.3],
               xlabel=s1,
               yscale='log',
               ylim=[0.00001, 0.3],
               ylabel=s2,
               title=f'{s1} vs. {s2}: {column}')
        ax.legend(title=colorby, loc="lower right")

###############################################################################
# Heatmap plotting functions
#     get_distance_matrix
#     plot_contour_distances
#     plot_heatmap_data
#     make_heatmap
#     Future:
#       set_heatmap
###############################################################################

    def get_distance_matrix(self, structure):
        if not hasattr(self, "distances"):
            self.distances = {}
        if structure in self.distances.keys():
            return self.distances[structure]
        if structure == "ct":
            matrix = self.ct.get_distance_matrix()
            self.distances["ct"] = matrix
        elif structure == "pdb":
            length = self.length["pdb"]
            fill = get_default_fill('Distance')
            matrix = np.full([length, length], fill)
            for i in range(length):
                for j in range(length):
                    if matrix[i, j] == fill:
                        matrix[i, j] = self.get_3d_distance(i+1, j+1)
                        matrix[j, i] = matrix[i, j]
            self.distances["pdb"] = matrix
        return matrix

    def plot_contour_distances(self, ax, structure, ij_data, levels):
        clip_pad = self.get_clip_pad(structure, ij_data)
        length = self.length[structure]+sum(clip_pad[1])
        fill = get_default_fill('Distance')
        distances = np.full([length, length], fill)
        start = clip_pad[1][0]
        end = length - clip_pad[1][1]
        distances[start:end, start:end] = self.get_distance_matrix(structure)
        if levels is None:
            levels = {"ct": [5], "pdb": [20, 500]}[structure]
        cmap = mp.colors.LinearSegmentedColormap.from_list('contours',
                                                           ['black', 'gray'])
        x_y = list(range(1, length+1))
        ax.contour(x_y, x_y, distances, levels=levels, cmap=cmap,
                   linewidths=0.5)

    def plot_heatmap_data(self, ax, ij_data, structure, metric=None):
        clip_pad = self.get_clip_pad(ij_data, structure)
        data = self.ij_data[ij_data].copy()
        if metric is None:
            metric = get_default_metric(ij_data)
        columns = ["i", "j", metric]
        if ij_data == "rings":
            data = data[columns+["+/-"]]
            data[metric] = data[metric]*data["+/-"]
            data = data[columns]
        else:
            data = data[columns]
        data[["i", "j"]] += clip_pad[1][0]
        length = self.length[ij_data] + sum(clip_pad[1])
        fill = get_default_fill(metric)
        data_im = np.full([length, length], fill)
        window = self.window[ij_data]
        for _, i, j, value in data.itertuples():
            data_im[i-1:i-1+window, j-1:j-1+window] = value
            data_im[j-1:j-1+window, i-1:i-1+window] = value
        min_max = get_default_min_max(metric)
        if metric == "Percentile":
            min_max = [0, 1]
            cmap = plt.get_cmap("jet")
            cmap = cmap(np.linspace(0, 1, 256))
            cmap[:, -1] = 0.6
        elif metric == "Class":
            data_im = data_im + 1
            min_max = [0, 3]
            cmap = get_default_cmap(metric)
            cmap = cmap(np.arange(cmap.N))
            no_data = np.array([1, 1, 1, 0])
            cmap = np.vstack((no_data, cmap))
        else:
            cmap = get_default_cmap(metric)
            cmap = cmap(np.arange(cmap.N))
        if ij_data != "rings":
            cmap[0, :] = [1, 1, 1, 0]
        cmap = mp.colors.ListedColormap(cmap)
        ax.imshow(data_im, cmap=cmap, vmin=min_max[0], vmax=min_max[1],
                  interpolation='none')

    def make_heatmap(self, ij_data, structure, metric=None, ax=None,
                     levels=None):
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(10, 10))
        # self.set_heatmap(ax)
        self.plot_heatmap_data(ax, ij_data, structure, metric)
        self.plot_contour_distances(ax, structure, ij_data, levels)

###############################################################################
# i,j data filters and coloring functions
#   get_cd_mask
#   set_mask_offset
#   filter_ij_data
#   get_ij_colors
###############################################################################
    def filter_dance_rings(self, filterneg=True, cdfilter=15, sigfilter=20,
                           ssfilter=True):
        ctlist = [dance.ct for dance in self.dance]
        ringlist = [dance.ij_data["rings"].copy() for dance in self.dance]
        for index, rings in enumerate(ringlist):
            rings[['i_offset', 'j_offset']] = rings[['i', 'j']]
            mask = []
            for _, row in rings.iterrows():
                i, j, G, sign = row[['i', 'j', 'Statistic', '+/-']]
                true_so_far = True
                if ssfilter and true_so_far:
                    true_so_far = (ctlist[index].ct[i-1]
                                   == 0 and ctlist[index].ct[j-1] == 0)
                if filterneg and true_so_far:
                    true_so_far = sign == 1
                if sigfilter is not None and true_so_far:
                    true_so_far = G > sigfilter
                if cdfilter is not None and true_so_far:
                    for ct in ctlist:
                        if not true_so_far:
                            break
                        true_so_far = ct.contactDistance(i, j) >= cdfilter
                mask.append(true_so_far)
            rings['mask'] = mask
            self.dance[index].ij_data["rings"] = rings

    def get_cd_mask(self, ij_data, fit_to_sequence, cdAbove=None, cdBelow=None,
                    ss_only=False, ds_only=False, paired_only=False):
        """creates a mask for filtering i,j data based on contact distance.

        Args:
            data (str): string matching ij data attribute to be masked
            fit_to_sequence (string): the structure to derive contact distance
            cdAbove (int, optional): lower limit of contact distance
                    Defaults to None.
            cdBelow (int, optional): upper limit of contact distance
                    Defaults to None.

        Returns:
            boolean list: mask to be applied to data
        """
        if fit_to_sequence == 'ss':
            basepairs = [tuple(pair) for pair in self.basepairs]
            ct = rna.CT()
            ct.pair2CT(pairs=basepairs, seq=self.sequence["ss"])
        elif fit_to_sequence == 'ct':
            ct = self.ct
        ij_data = self.ij_data[ij_data]
        mask = []
        i_j_keep = ["i_offset", "j_offset", "mask"]
        for _, i, j, keep in ij_data[i_j_keep].itertuples():
            true_so_far = keep
            if paired_only and true_so_far:
                true_so_far = ct.ct[i-1] == j
            if ss_only and true_so_far:
                true_so_far = (ct.ct[i-1] == 0) and (ct.ct[j-1] == 0)
            if ds_only and true_so_far:
                true_so_far = (ct.ct[i-1] != 0) and (ct.ct[j-1] != 0)
            if cdAbove is not None and true_so_far:
                true_so_far = ct.contactDistance(i, j) > cdAbove
            if cdBelow is not None and true_so_far:
                true_so_far = ct.contactDistance(i, j) < cdBelow
            mask.append(true_so_far)
        return mask

    def get_profile_mask(self, ij_data, profAbove=None, profBelow=None):
        clip, pad = self.get_clip_pad(ij_data, "profile")
        data = self.ij_data[ij_data].copy()
        norm_prof = self.profile["Norm_profile"]
        mask = []
        for _, i, j in data[["i", "j"]].itertuples():
            if not (clip[1] > i > clip[0]) or not (clip[1] > j > clip[0]):
                mask.append(False)
                continue
            prof_i = norm_prof[i-1+pad[0]]
            prof_j = norm_prof[j-1+pad[0]]
            if profAbove is not None:
                above = (prof_i >= profAbove) and (prof_j >= profAbove)
            else:
                above = True
            if profBelow is not None:
                below = (prof_i <= profBelow) and (prof_j <= profBelow)
            else:
                below = True
            mask.append(above and below)
        return mask

    def set_mask_offset(self, ij_data, fit_to_sequence):
        """Creates 3 new columns of the attribute dataFrame: mask, i_offset,
        and j_offset. These are the clipped and padded values to fit to the
        given sequence.

        Args:
            attribute (str): string matching i, j data attribute.
            fit_to_sequence (str): string matching sequence to fit data to.
        """
        clip_pad = self.get_clip_pad(ij_data, fit_to_sequence)
        start = clip_pad[0][0]
        end = clip_pad[0][1]
        data = self.ij_data[ij_data]
        i = data["i"].values
        j = data["j"].values
        iinrange = (start < i) & (i < end)
        jinrange = (start < j) & (j < end)
        data['mask'] = np.ones(len(data), dtype=bool)
        data['mask'] = data['mask'] & iinrange & jinrange
        offset = clip_pad[1][0]-start
        data['i_offset'] = i+offset
        data['j_offset'] = j+offset

    def filter_ij_data(self, ij_data, fit_to_sequence, cdAbove=None,
                       cdBelow=None, ss_only=False, ds_only=False,
                       profAbove=None, profBelow=None, all_pairs=False,
                       **kwargs):
        """Apply given filters to i,j data attribute.

        Args:
            attribute (str): string matching i,j data attribute
            fit_to_sequence (str): string matching attribute to fit to
            cdAbove (int, optional): lower bound of contact distance
                    Defaults to None.
            cdBelow (int, optional): upper bound of contact distance
                    Defaults to None.
            **kwargs:
                Keywords should be column names of the given attribute.
                Arguments should be value of the lower bound for that column.
        """
        data = self.ij_data[ij_data]
        self.set_mask_offset(ij_data, fit_to_sequence)
        if fit_to_sequence == 'pdb':
            mask = []
            for i, j in zip(data["i_offset"], data["j_offset"]):
                mask.append(i in self.pdb_validres and j in self.pdb_validres)
            mask = np.array(mask, dtype=bool)
            data["mask"] = data["mask"] & mask
        if profAbove is not None or profBelow is not None:
            prof_mask = self.get_profile_mask(ij_data, profAbove, profBelow)
            data["mask"] = data["mask"] & prof_mask
        if cdAbove is not None or cdBelow is not None or ss_only or ds_only:
            cd_mask = self.get_cd_mask(ij_data, fit_to_sequence,
                                       cdAbove, cdBelow, ss_only, ds_only)
            data['mask'] = data['mask'] & cd_mask
        if not all_pairs and ij_data == 'pairs':
            data["mask"] = data["mask"] & (data["Class"] != 0)
        if ij_data == 'probs':
            data['mask'] = data['mask'] & (data["Probability"] >= 0.03)
        for key in kwargs.keys():
            try:
                mask = data[key] > kwargs[key]
                data["mask"] = data["mask"] & mask
            except KeyError:
                print(f"{key} is not a valid column of {ij_data} dataFrame")

    def get_ij_colors(self, ij_data, metric=None, min_max=None):
        """get lists for i_offset, j_offset, and mpl colors for the given
        i,j data attribute and metric.

        Args:
            attribute (str): string matching i,j data attribute
            metric (str, optional): column of i,j data attribute or "distance
                    Defaults to None.
            min_max (list of 2 ints, optional): min and max value for coloring
                    Defaults to None.

        Returns:
            3 lists: offset i and j values, and mpl colors
        """
        if metric is None:
            metric = get_default_metric(ij_data)
        if min_max is None:
            min_max = get_default_min_max(metric)
        columns = ["i_offset", "j_offset"]
        data = self.ij_data[ij_data]
        if ij_data == 'rings':
            columns.extend([metric, "+/-"])
            data = data[data["mask"]][columns].copy()
            data[metric] = data[metric]*data["+/-"]
        else:
            columns.append(metric)
            data = data[data["mask"]][columns].copy()
        if metric != 'Class':  # Normalize to between 0 and 1
            data.loc[data[metric] < min_max[0], metric] = min_max[0]
            data.loc[data[metric] > min_max[1], metric] = min_max[1]
            ascending = metric not in ['Distance']  # high distance is bad
            data = data.sort_values(by=metric, ascending=ascending)
            data[metric] = (data[metric]-min_max[0]) / (min_max[1]-min_max[0])
        else:  # Class data have a weird order
            data = pd.concat([data[data[metric] == 0],
                              data[data[metric] == 2],
                              data[data[metric] == 1]])
        i = data["i_offset"].values
        j = data["j_offset"].values
        cmap = get_default_cmap(metric)
        colors = cmap(data[metric].values)
        return i, j, colors

    def print_new_ij_file(self, ij_data, outfile=None, **kwargs):
        self.filter_ij_data(ij_data, "ct", **kwargs)
        if ij_data in self.header.keys():
            header = self.header[ij_data]
        else:
            header = ''
        data = self.ij_data[ij_data][self.ij_data[ij_data]["mask"]].copy()
        columns = list(data.columns)
        exclude_columns = ["i_offset", "j_offset",
                           "mask", "Distance", "Percentile"]
        for col in exclude_columns:
            if col in columns:
                columns.remove(col)
        csv = data.to_csv(columns=columns, sep='\t', index=False,
                          line_terminator='\n')
        if outfile is not None:
            with open(outfile, 'w') as out:
                out.write(header)
                out.write(csv)
        else:
            print(header, csv)

###############################################################################
# filters, clipping, paddings, and coloring functions
#   get_alignment
#   get_clip_pad
#   get_colorby_profile
#   get_colorby_sequence
#   get_colorby_position
###############################################################################

    def get_alignment(self, fit_this, to_that):
        """Might be used instead of clip/pad in the future. Would allow more
        flexibility for fitting data from different, but similar sequences.

        Args:
            fit_this (str): string matching attribute to be filtered
            to_that (str): string matching attribute to be fit to
        """
        from Bio.pairwise2 import align, format_alignment
        this_sequence = self.sequence[fit_this].upper()
        that_sequence = self.sequence[to_that].upper()
        alignment = align.globalxs(this_sequence, that_sequence, -1, -0.1,
                                   penalize_end_gaps=False)
        print(format_alignment(*alignment[0]))

    def get_clip_pad(self, fit_this, to_that):
        """Takes two attributes as strings, and returns a 5' and 3' clip and
        pad for the first to fit the data to the second.

        Args:
            fit_this (str): String matching attribute to be clipped/padded
            to_that (str): String matching attribute to be fit to

        Returns:
            [2x2 list]: 5' and 3' clip and 5' and 3' pad
        """
        if not hasattr(self, "clip_pad"):
            self.clip_pad = {}
        if fit_this not in self.clip_pad.keys():
            self.clip_pad[fit_this] = {}
        if to_that in self.clip_pad[fit_this].keys():
            return self.clip_pad[fit_this][to_that]
        this_sequence = self.sequence[fit_this]
        that_sequence = self.sequence[to_that]
        this_length = self.length[fit_this]
        that_length = self.length[to_that]
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
        # Assertion error if sequences don't match well
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
        self.clip_pad[fit_this][to_that] = [clip, pad]
        return [clip, pad]

    def get_colorby_profile(self, sequence):
        """Returns list of colors representing reactivity profile, matched to
        the given sequence.

        Args:
            sequence (str): string matching attribute with sequence to fit to.

        Returns:
            list of mpl colors: colors representing reactivity profile data
        """
        clip_pad = self.get_clip_pad("profile", sequence)
        if clip_pad is not None:
            start, end = clip_pad[0]
        else:
            start, end = 0, -1
        cmap = ['gray', 'black', 'orange', 'red']
        bins = [0, 0.4, 0.85]
        profcolors = []
        for x in self.profile["Norm_profile"][start:end]:
            profcolors.append(sum([b < x for b in bins]))
        if clip_pad is not None:
            profcolors = ['gray']*clip_pad[1][0] + profcolors
            profcolors += ['gray']*clip_pad[1][1]
        colors = np.array([cmap[val] for val in profcolors])
        return colors

    def get_colorby_sequence(self, sequence, colors='new'):
        """Returns list of mpl colors representing sequence.

        Args:
            sequence (str): string matching attribute with sequence to color
            colors (str, optional): Options are "new" or "old".
                    "new": A=blue, U=light blue, G=red, C=light red
                    "old": A=red, U=yellow, G=blue, C=green
                    Defaults to 'new'.

        Returns:
            list of mpl colors: color list representing sequence
        """
        seq = self.sequence[sequence]
        colors = np.array([get_nt_color(nt.upper(), colors) for nt in seq])
        return colors

    def get_colorby_position(self, sequence, cmap='rainbow'):
        """Returns list of mpl colors that spans the rainbow. Fits length of
        given sequence.

        Args:
            sequence (str): string matching sequence to fit colorlist to
            cmap (str, optional): mpl colormap to use. Defaults to 'rainbow'.

        Returns:
            list of mpl colors: spectrum of colors with same length as sequence
        """
        length = len(self.sequence[sequence])
        cmap = plt.get_cmap(cmap)
        colors = np.array([cmap(n/length) for n in range(length)])
        return colors

###############################################################################
# Plotting functions that accept a list of samples
#   array_qc
#   array_skyline
#   array_ap
#   array_ss
#   array_3d
###############################################################################


def get_rows_columns(number_of_samples, rows=None, cols=None):
    if isinstance(rows, int) and cols is None:
        cols = math.ceil(number_of_samples / rows)
    elif isinstance(cols, int) and rows is None:
        rows = math.ceil(number_of_samples / cols)
    elif number_of_samples < 10:
        rows, cols = [(0, 0), (1, 1), (1, 2), (1, 3), (2, 2),  # 0-4 samples
                      (2, 3), (2, 3), (3, 3), (3, 3), (3, 3)  # 5-9 samples
                      ][number_of_samples]
    else:
        cols = 4
        rows = math.ceil(number_of_samples / cols)
    return rows, cols


def array_ap(samples=[], rows=None, cols=None, **kwargs):
    rows, cols = get_rows_columns(len(samples), rows, cols)
    figsize = samples[0].get_ap_figsize(rows, cols)
    gs_kw = {'rows': rows, 'cols': cols, 'hspace': 0.01, 'wspace': 0.1}
    _, axes = plt.subplots(rows, cols, figsize=figsize, gridspec_kw=gs_kw)
    for i, sample in enumerate(samples):
        row = i // cols
        col = i % rows
        ax = axes[row, col]
        sample.make_ap(ax=ax, **kwargs)


def array_qc(samples=[]):
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(241)
    ax2 = fig.add_subplot(242)
    ax3 = fig.add_subplot(243)
    ax4 = fig.add_subplot(244)
    ax5 = fig.add_subplot(212)
    total = len(samples)
    for i, sample in enumerate(samples):
        sample.plot_log_MutsPerMol(ax1, sample="Untreated")
        sample.plot_log_MutsPerMol(ax2, sample="Modified")
        sample.plot_log_ReadLength(ax3, sample="Untreated", n=i, of=total)
        sample.plot_log_ReadLength(ax4, sample="Modified", n=i, of=total)
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


def array_3d(samples, rows=None, cols=None, **kwargs):
    """Given a list of samples, creates a grid of py3Dmol views. Kwargs passed
    to samples[0].make_3d(). Samples must all use the same structure.

    Args:
        samples (List of Sample objects): Sample objects from which to display.
        rows (int, optional): number of rows for figure.
            Default determined by get_rows_columns.
        cols (int, optional): number of columns for figure.
            Default determined by get_rows_columns.

    Returns:
        [type]: [description]
    """
    rows, cols = get_rows_columns(len(samples), rows, cols)
    view = py3Dmol.view(viewergrid=(rows, cols),
                        width=400*rows,
                        height=400*cols)
    view = samples[0].set_3d_view(view)
    for i, sample in enumerate(samples):
        row = i // cols
        col = i % rows
        view = sample.make_3d(view, (row, col), **kwargs)
    return view
