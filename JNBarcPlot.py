import pmanalysis as pma
import RNAtools2 as RNA
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plottingTools as pt


def addArc(ax, i, j, window, color, alpha, panel):
    if panel == "top":
        center = ((i+j)/2, 0)
        theta1 = 0
        theta2 = 180
    if panel == "bottom":
        center = ((i+j)/2, -3)
        theta1 = 180
        theta2 = 360
    radius = 1+(j-i)/2
    arc = mpl.patches.Wedge(center, radius, theta1, theta2, color=color,
                            alpha=alpha, width=window, ec='none')
    ax.addPath(arc)


def getCTPairs(ctpath):
    names = ['i', 'j']
    ct = pd.read_csv(ctpath, sep='\s+', usecols=[0, 4], names=names, header=1)
    ct = ct[ct.j > ct.i]
    return {tuple(pair) for pair in zip(ct.i, ct.j)}


def addCT(ax, ctpath):
    ct_pairs = getCTPairs(ctpath)
    for i, j in ct_pairs:
        addArc(ax, i, j, 1, 'black', 0.7, 'top')


def addCTCompare(ax, ctpath1, ctpath2):
    ct1 = getCTPairs(ctpath1)
    ct2 = getCTPairs(ctpath1)
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
    profile = pd.read_csv(profilepath, sep='\t')
    near_black = (0, 0, 1 / 255.0)
    orange_thresh = 0.4
    red_thresh = 0.85
    cindex = np.zeros(len(profile['Norm_profile']), dtype=int)
    cindex[np.array(np.logical_not(np.isnan(profile['Norm_profile'])),
           dtype=bool)] += 1
    cindex[np.array(profile["Norm_profile"] > orange_thresh, dtype=bool)] += 1
    cindex[np.array(profile['Norm_profile'] > red_thresh, dtype=bool)] += 1
    colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
    ax.bar(profile['Nucleotide'], profile['Norm_profile']*5, align="center",
           width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
           yerr=profile['Norm_stderr'], ecolor=near_black, capsize=1)
    pt.addSeqBar(ax, profile, yvalue=-1.5)


def addPairmap(ax, pairmappath, window=3):
    pm = pd.read_csv(pairmappath, sep='\t', header=1)
    primary = pm[pm['Class'] == 1]
    primary = zip(primary['i'], primary['j'])
    secondary = pm[pm['Class'] == 2]
    secondary = zip(secondary['i'], secondary['j'])
    for i, j in secondary:
        addArc(ax, i, j, window, (30/255, 194/255, 1), 0.4, 'bottom')
    for i, j in primary:
        addArc(ax, i, j, window, (0, 0, 243/255), 0.8, 'bottom')


def setArcPlot(ax):
    ax.set_aspect('equal')
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')


def plotCorrs2D(axis, pairmap=None, allcorrs=None, bg_corrs=None, ct=None,
                mask=None):
    ct_cmap, mask_cmap, pm_cmap = getCmaps()
    if type(allcorrs) is np.ndarray:
        axis.imshow(allcorrs, cmap='bwr', interpolation='bicubic')
    if type(bg_corrs) is np.ndarray:
        axis.imshow(bg_corrs, cmap='gray', interpolation='bicubic')
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


def fastaToMask(fastafile):
    with open(fastafile) as f:
        f.readline()
        seq = ""
        for line in f.readlines():
            seq += line.strip()
    size = len(seq)
    mask = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if (i > j) or seq[i].islower() or seq[j].islower():
                mask[i, j] = 1
    return mask


def ctLength(ctfile):
    with open(ctfile) as f:
        line1 = f.readline()
        length = int(line1.strip().split()[0])
        return length


def ignoredToBitmap(logfile, ctfile):
    size = ctLength(ctfile)
    bitmap = np.ones((size, size))
    with open(logfile) as f:
        for line in f.readlines():
            if line.startswith('Pair '):
                pair = line.split()[1].strip('()').split(',')
                bitmap[pair[0], pair[1]] = 0
    return bitmap


def pairmapSensPPV(pairfile, ctfile):
    pm = pma.PairMap(pairfile)
    ct = RNA.CT(ctfile, filterNC=True, filterSingle=True)
    p,s = pm.ppvsens_duplex(ct, ptype=1, exact=False)


def corrsToBitmap(corrfile, ctfile, window=1, limit=100):
    size = ctLength(ctfile)
    bitmap = np.zeros((size, size))
    corrs = pd.read_csv(corrfile, sep='\t', header=1)
    for i in range(window*window):
        x = i % window
        y = int(i/window)
        bitmap[corrs['i']+x, corrs['j']+y] += corrs['Statistic']*corrs['+/-']
    bitmap[bitmap > limit] = limit
    bitmap[bitmap < -limit] = -limit
    bitmap = bitmap/limit
    bitmap[0, 0] = 1
    bitmap[0, 1] = -1
    bitmap += 1
    bitmap /= 2
    return bitmap


def ctToBitmap(ctfile):
    size = ctLength(ctfile)
    bitmap = np.ones((size, size))
    names = ['i', 'j']
    ct = pd.read_csv(ctfile, sep='\s+', names=names,
                     header=0, usecols=[0, 4])
    ct = ct[(ct['j'] != 0) & (ct['i'] < ct['j'])]
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            bitmap[ct['i']+x, ct['j']+y] = 0
    return bitmap


def pairmapToBitmap(pairfile, ctfile):
    size = ctLength(ctfile)
    bitmap = np.zeros((size, size))
    pairs = pd.read_csv(pairfile, sep='\t', header=1)
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
