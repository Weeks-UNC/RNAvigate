#!/usr/bin/env python
import pmanalysis as pma
import RNAtools2 as RNA
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plottingTools as pt


def addArc(ax, i, j, window, color, alpha, panel):
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
    arc = mpl.patches.Wedge(center, radius, theta1, theta2, color=color,
                            alpha=alpha, width=window, ec='none')
    ax.add_patch(arc)


def getCTPairs(ctpath):
    """From a ct file, extract the i,j pairs and return as a list of tuples

    Args:
        ctpath (str): path to ct file

    Returns:
        list of tuples: list of i,j pairs, example: [(1, 20), (2, 19), (3, 18)]
    """
    names = ['i', 'j']
    ct = pd.read_csv(ctpath, sep=r'\s+', usecols=[0, 4], names=names, header=0)
    ct = ct[ct.j > ct.i]
    return {tuple(pair) for pair in zip(ct.i, ct.j)}


def addCT(ax, ctpath):
    """Adds arc plot representation of the ct structure to given axis

    Args:
        ax (pyplot axis): axis to which structure will be added
        ctpath (str): path to ct file
    """
    ct_pairs = getCTPairs(ctpath)
    for i, j in ct_pairs:
        addArc(ax, i, j, 1, 'black', 0.7, 'top')


def addCTCompare(ax, ctpath1, ctpath2):
    """Adds structure comparison arc plot to the given axis

    Args:
        ax (pyplot axis): axis to which structure comparison is added
        ctpath1 (str): path to first ct file
        ctpath2 (str): path to second ct file
    """
    ct1 = getCTPairs(ctpath1)
    ct2 = getCTPairs(ctpath2)
    shared = ct1.union(ct2)
    ref = ct1.difference(ct2)
    comp = ct2.difference(ct1)
    sharedcolor = (150/255., 150/255., 150/255.)
    refcolor = (38/255., 202/255., 145/255.)
    compcolor = (153/255., 0.0, 1.0)
    for i, j in comp:
        addArc(ax, i, j, 1, compcolor, 0.7, 'top')
    for i, j in ref:
        addArc(ax, i, j, 1, refcolor, 0.7, 'top')
    for i, j in shared:
        addArc(ax, i, j, 1, sharedcolor, 0.7, 'top')


def addProfile(ax, profilepath):
    """Adds bar graph of normalize reactivity profile to given axis

    Args:
        ax (pyplot axis): axis to which reactivity profile is added
        profilepath (str): path to profile.txt file
    """
    profile = pt.ReactivityProfile(profilepath)
    near_black = (0, 0, 1 / 255.0)
    orange_thresh = 0.4
    red_thresh = 0.85
    values = profile.profile['Norm_profile']
    cindex = np.zeros(len(values), dtype=int)
    # where values are not NaNs, add 1 to color index array
    cindex[np.array(np.logical_not(np.isnan(values)), dtype=bool)] += 1
    # where values excede orange threshold (0.4), add 1 to color index array
    cindex[np.array(values > orange_thresh, dtype=bool)] += 1
    # where values excede red threshold (0.85), add 1 to color index array
    cindex[np.array(values > red_thresh, dtype=bool)] += 1
    # create color map array based on cindex
    colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
    ax.bar(profile.profile['Nucleotide'], values*5, align="center",
           width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
           yerr=profile.profile['Norm_stderr'], ecolor=near_black, capsize=1)
    profile.addSeqBar(ax, yvalue=0.5)


def addPairmap(ax, pairmappath, window=3):
    """add PAIR-MaP data to the given axis

    Args:
        ax (pyplot axis): axis to which PAIR data is added
        pairmappath (str): path to PairMapper data file
        window (int, optional): Window size used in pairmapper. Defaults to 3.
    """
    pm = pd.read_csv(pairmappath, sep='\t', header=1)
    primary = pm[pm['Class'] == 1]
    primary = list(zip(primary['i'], primary['j']))
    secondary = pm[pm['Class'] == 2]
    secondary = list(zip(secondary['i'], secondary['j']))
    for i, j in secondary:
        addArc(ax, i, j, window, (30/255., 194/255., 1.), 0.2, 'bottom')
    for i, j in primary:
        addArc(ax, i, j, window, (0., 0., 243/255.), 0.7, 'bottom')


def setArcPlot(ax):
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


def plotCorrs2D(axis, pairmap=0, allcorrs=0, bg_corrs=0, ct=0, mask=0):
    ct_cmap, mask_cmap, pm_cmap = getCmaps()
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


def getBitmapLegendHandles():
    handles = {}
    handles['Primary'] = mpl.patches.Patch(color='red', label='Primary')
    handles['Secondary'] = mpl.patches.Patch(color='blue', label='Secondary')
    handles['No Data'] = mpl.patches.Patch(color='green', label='No Data',
                                           alpha=0.2)
    handles['Helices'] = mpl.patches.Patch(edgecolor='green', fill=False,
                                           label='Known Helices')
    handles['BG Correlations'] = mpl.patches.Patch(color='black',
                                                   label='BG correlations')
    return handles


def getCmaps():
    cmap = plt.get_cmap('Greens')
    ct_cmap = cmap(np.arange(cmap.N))
    ct_cmap[:, -1] = np.concatenate((np.linspace(0, 1, cmap.N/2),
                                     np.linspace(1, 0, cmap.N/2)), axis=None)
    ct_cmap = mpl.colors.ListedColormap(ct_cmap)

    mask_cmap = cmap(np.arange(cmap.N))
    mask_cmap[:, -1] = np.linspace(0, 0.2, cmap.N)
    mask_cmap = mpl.colors.ListedColormap(mask_cmap)

    N = 256
    vals = np.ones((N, 4))
    vals[:, 0] = np.concatenate((np.linspace(1, 0, N/2),
                                 np.linspace(0, 1, N/2)), axis=None)
    vals[:, 1] = np.concatenate((np.linspace(1, 0, N/2),
                                 np.linspace(0, 0, N/2)), axis=None)
    vals[:, 2] = np.concatenate((np.linspace(1, 1, N/2),
                                 np.linspace(1, 0, N/2)), axis=None)
    pm_cmap = mpl.colors.ListedColormap(vals)
    return ct_cmap, mask_cmap, pm_cmap


def toBitmap(filepath, ctpath='none', type='none', window=1, limit=200):
    valid_types = ['fasta', 'ct', 'corr', 'log', 'pairmap']
    if type not in valid_types:
        print(('Please provide a valid type: ' + ', '.join(valid_types)))
        return
    elif type == 'fasta':
        with open(filepath) as f:
            f.readline()
            seq = ""
            for line in f.readlines():
                seq += line.strip()
        size = len(seq)
        bitmap = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                if (i > j) or seq[i].islower() or seq[j].islower():
                    bitmap[i, j] = 1
    elif type == 'ct':
        size = ctLength(filepath)
        bitmap = np.ones((size, size))
        names = ['i', 'j']
        ct = pd.read_csv(filepath, sep=r'\s+', names=names,
                         header=0, usecols=[0, 4])
        ct = ct[ct.i < ct.j]
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                bitmap[ct['i']+x, ct['j']+y] = 0
    elif type == 'corr':
        size = ctLength(ctpath)
        bitmap = np.zeros((size, size))
        corrs = pd.read_csv(filepath, sep='\t', header=1)
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
        size = ctLength(ctpath)
        bitmap = np.ones((size, size))
        with open(filepath) as f:
            for line in f.readlines():
                if line.startswith('Pair '):
                    pair = line.split()[1].strip('()').split(',')
                    bitmap[pair[0], pair[1]] = 0
    elif type == 'pairmap':
        size = ctLength(ctpath)
        bitmap = np.zeros((size, size))
        pairs = pd.read_csv(filepath, sep='\t', header=1)
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


def ctLength(ctpath):
    with open(ctpath) as f:
        line1 = f.readline()
        length = int(line1.strip().split()[0])
    return length


def pairmapSensPPV(pairpath, ctpath):
    pm = pma.PairMap(pairpath)
    ct = RNA.CT(ctpath, filterNC=True, filterSingle=True)
    p, s = pm.ppvsens_duplex(ct, ptype=1, exact=False)
    return p, s