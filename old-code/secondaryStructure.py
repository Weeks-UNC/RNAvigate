#!/usr/bin/env python

# TODO: read in varna/xrna/nsd more elegantly
# TODO: Add no coloring option for sequence
# TODO: Add option for circles with sequence inside
# TODO: RING data ordering could be more sophisticated.

import matplotlib.colors as mc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import RNAtools2 as RNA
import numpy as np

sns.set_style("ticks")
sns.set_context("talk")


class SecondaryStructure():
    """A class for holding secondary structure coordinates and pairing data
    from nsd, cte, or xrna files, and plotting various data on top of it.
    """

    def __init__(self, structure_file=None, ringfile=None, profile=None,
                 structureCassettes=False):
        self.sequence = np.array([])
        self.length = 0
        self.nucleotide = np.array([])
        self.xcoordinates = np.array([])
        self.ycoordinates = np.array([])
        self.pairs = np.array([])
        self.colors = np.array([])
        self.structureCassettes = structureCassettes
        self.ss_filetype = structure_file.split('.')[-1]
        if self.ss_filetype == 'nsd':
            self.readNsdFile(structure_file)
        elif self.ss_filetype == 'cte':
            self.readCteFile(ctstructure_fileefile)
        elif self.ss_filetype == 'xrna':
            self.readXrnaFile(structure_file)
        else:
            "{} file type not recognized".format(self.ss_filetype)
        if profile is not None:
            self.readProfile(profile)
        if ringfile is not None:
            self.readRingFile(ringfile)

###############################################################################
# Manipulating Structure Coordinates
###############################################################################

    def readXrnaFile(self, xrnafile):
        pairs = []
        with open(xrnafile, 'r') as file:
            for line in file.readlines():
                if line.startswith("G") or line.startswith("U") or line.startswith("A") or line.startswith("C"):
                    nucdata = line.split(' ')
                    self.sequence = np.append(self.sequence, nucdata[0])
                    self.xcoordinates = np.append(
                        self.xcoordinates, float(nucdata[1]))
                    self.ycoordinates = np.append(
                        self.ycoordinates, float(nucdata[2]))
                elif line.startswith("<BasePairs"):
                    helix = line.split(' ')
                    i = int(helix[1].split('=')[1].strip("'"))
                    j = int(helix[3].split('=')[1].strip("'"))
                    length = int(helix[2].split('=')[1].strip("'"))
                    for x in range(length):
                        pairs.append([i+x, j-x])
        self.pairs = np.array(pairs)
        self.length = len(self.sequence)

    def readCteFile(self, ctefile):
        names = ['nuc', 'seq', 'pair', 'xcoords', 'ycoords']
        ct = pd.read_csv(ctefile, sep=r'\s+', usecols=[0, 1, 4, 8, 10],
                         names=names, header=0)
        self.sequence = list(ct['seq'])
        self.length = len(self.sequence)
        self.nucleotide = [int(nuc) for nuc in ct['nuc']]
        self.xcoordinates = [float(xcoord) for xcoord in ct['xcoords']]
        self.ycoordinates = [float(ycoord) for ycoord in ct['ycoords']]
        ct = ct[ct.nuc < ct.pair]
        self.pairs = [[int(x), int(y)] for x, y in zip(ct.nuc, ct.pair)]

    def readNsdFile(self, nsdfile):
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
                    base = nucs[2][1]
                    x = float(nucs[4][1])
                    y = float(nucs[5][1])
                    np.append(self.sequence, base)
                    np.append(self.nucleotide, nuc)
                    np.append(self.xcoordinates, x)
                    np.append(self.ycoordinates, y)
                    self.length += 1
                elif item == "pairs":
                    pairs = line[1].split(":")[1:]
                    pair = [int(nuc.strip('"')) for nuc in pairs]
                    np.append(self.pairs, [pair])

    def writeCTE(self, outputPath):
        """writes the current structure out to CTE format for Structure Editor.

        Args:
            outputPath (string): path to output cte file to be created
        """
        pairs = [tuple(pair) for pair in self.pairs]
        ct = RNA.CT()
        ct.pair2CT(pairs=pairs, seq=self.sequence)
        # set scaling factors based on data source.
        # TODO: figure out varna scaling factor
        xscale = {'xrna': 1.5, 'varna': 1, 'nsd': 1, 'cte': 1}[self.type]
        yscale = {'xrna': -1.5, 'varna': 1, 'nsd': 1, 'cte': 1}[self.type]
        ctlen = len(ct.num)
        w = open(outputPath, 'w')
        w.write('{0:6d} {1}\n'.format(ctlen, ct.name))
        line = '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d} ;! X: {5:5d} Y: {6:5d}\n'
        for i in range(ctlen):
            xcoord = int(xscale*self.xcoordinates[i])
            ycoord = int(yscale*self.ycoordinates[i])
            num, seq, cti = ct.num[i], ct.seq[i], ct.ct[i]
            # no nums after ctlen, resets to zero
            nextnum = (num+1) % (ctlen+1)
            cols = [num, seq, num-1, nextnum, cti, xcoord, ycoord]
            w.write(line.format(*cols))
        w.close()

    def plotSS(self, ax):
        for pair in self.pairs:
            xcoords = [self.xcoordinates[pair[0]-1],
                       self.xcoordinates[pair[1]-1]]
            ycoords = [self.ycoordinates[pair[0]-1],
                       self.ycoordinates[pair[1]-1]]
            ax.plot(xcoords, ycoords, color='grey',
                    linestyle=(0, (1, 1)), zorder=0)
        ax.plot(self.xcoordinates, self.ycoordinates, color='grey', zorder=0)

    def plotPositions(self, ax, spacing=20):
        for i in range(0, self.length, spacing):
            ax.annotate(i+1, xy=(self.xcoordinates[i], self.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords='data', fontsize=16,
                        bbox=dict(facecolor="white", edgecolor='none', pad=0.3))

###############################################################################
# Manipulating RInG data
###############################################################################

    def readRingFile(self, ringfile):
        with open(ringfile, 'r') as file:
            header = file.readline().split('\t')
        self.rings_window = int(header[1].split('=')[1])
        self.rings = pd.read_csv(ringfile, sep='\t', header=1)
        if self.structureCassettes:
            start = 14
            end = self.length + 15
            mask = []
            for i, j in zip(self.rings["i"], self.rings["j"]):
                iinrange = start < i < end
                jinrange = start < j < end
                mask.append(iinrange and jinrange)
            self.rings = self.rings[mask]
            self.rings["i"] -= start
            self.rings["j"] -= start
        self.resetRingsFilter()

    def resetRingsFilter(self):
        self.rings_filtered = self.rings.copy()

    def filterRings(self, cdAbove=None, cdBelow=None, statistic=None, zscore=None,
                    ctfile=None):
        self.resetRingsFilter()
        if cdAbove is not None or cdBelow is not None:
            pairs = [tuple(pair) for pair in self.pairs]
            ct = RNA.CT()
            ct.pair2CT(pairs=pairs, seq=self.sequence)
            mask = []
            for i, j in zip(self.rings_filtered["i"], self.rings_filtered["j"]):
                if isinstance(cdAbove, int):
                    mask.append(ct.contactDistance(i, j) > cdAbove)
                if isinstance(cdBelow, int):
                    mask.append(ct.contactDistance(i, j) < cdBelow)
            self.rings_filtered = self.rings_filtered[mask]
        if statistic is not None:
            mask = []
            for stat in self.rings_filtered["Statistic"]:
                mask.append(stat > statistic)
            self.rings_filtered = self.rings_filtered[mask]
        if zscore is not None:
            mask = []
            for zij in self.rings_filtered["Zij"]:
                mask.append(zij > zscore)
            self.rings_filtered = self.rings_filtered[mask]

    def plotRings(self, ax, statistic="Statistic", bins=None):
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
        self.rings_filtered.sort_values(by=['Statistic'], inplace=True)

        for i, j, stat, sign in zip(self.rings_filtered['i'], self.rings_filtered['j'], self.rings_filtered[statistic], self.rings_filtered["+/-"]):
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

###############################################################################
# Manipulating Profiles
###############################################################################

    def readProfile(self, profile):
        self.profile = pd.read_csv(profile, sep='\t')
        if self.structureCassettes:
            start = 14
            end = self.length + 15
            mask = [start < i < end for i in self.profile['Nucleotide']]
        self.profile = self.profile[mask]

    def setColorsByProfile(self):
        cmap = ['grey', 'black', 'orange', 'red']
        bins = [0, 0.4, 0.85]
        profcolors = []
        for x in self.profile["Norm_profile"]:
            profcolors.append(sum([b < x for b in bins]))
        self.colors = np.array([cmap[val] for val in profcolors])

###############################################################################
# General Plotting Functions
###############################################################################

    def setPlot(self, ax):
        ax.set_aspect('equal')
        ax.axis('off')

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

    def setLetterColors(self):
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

    def plotSequence(self, ax, colorby=None, markers='o'):
        if colorby == "profile":
            self.setColorsByProfile()
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

    def makePlot(self, ax=None, setPlot=True, ss=True, positions=True, rings=True,
                 sequence=True, colorby='sequence', markers='o'):
        if ax is None:
            ax = plt.gca()
        if setPlot == True:
            self.setPlot(ax)
        if ss == True:
            self.plotSS(ax)
        if sequence == True:
            self.plotSequence(ax, colorby=colorby, markers=markers)
        if positions == True:
            self.plotPositions(ax)
        if hasattr(self, "rings") and rings == True:
            self.plotRings(ax)

    def figsize(self, rows, cols):
        xmin = min(self.xcoordinates)
        xmax = max(self.xcoordinates)
        ymin = min(self.ycoordinates)
        ymax = max(self.ycoordinates)
        scale = 0.025
        return ((xmax-xmin)*scale*cols, (ymax-ymin)*scale*rows)