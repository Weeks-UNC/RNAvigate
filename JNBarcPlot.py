import arcPlot as ap
import RNAtools2 as RNA
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plottingTools as pt


def pairmapPlot(ax, ctfile, pairmapfile, profilefile):
    ct = RNA.CT(ctfile)
    ct_pairs = ct.pairList()
    pm = pd.read_csv(pairmapfile, sep='\t', header=1)
    primary = pm[pm['Class'] == 1]
    primary = zip(primary['i'], primary['j'])
    secondary = pm[pm['Class'] == 2]
    secondary = zip(secondary['i'], secondary['j'])
    profile = pd.read_csv(profilefile, sep='\t')

    for i, j in ct_pairs:
        arc = mpl.patches.Wedge(((i+j)/2, 0), (j-i)/2, 0, 180,
                                color='black', width=1)
        ax.add_patch(arc)
    for i, j in secondary:
        arc = mpl.patches.Wedge(((i+j+2)/2, -2), 1+(j-i)/2, 180, 0, width=3,
                                ec='none', color=(0, 0, 1, 0.2))
        ax.add_patch(arc)
    for i, j in primary:
        arc = mpl.patches.Wedge(((i+j+2)/2, -2), 1+(j-i)/2, 180, 0, width=3,
                                ec='none', color=(1, 0, 0, 0.8))
        ax.add_patch(arc)
    ax.set_aspect('equal')
    ax.set(ylim=[-40, 40],
           xlim=[0, 267])
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    pt.addSeqBar(ax, profile, yvalue=-1.5)

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
    # read in ct values
    ct = RNA.CT(ctfile)
    ct_pairs = ct.pairList()
    ct_helices = ct.extractHelices()

    # read in pairmap primary signals
    pm = pd.read_csv(pairfile, sep='\t', header=1)
    primary = pm[pm['Class'] == 1]
    pm_pairs = [(pair['i'], pair['j']+2) for index, pair in primary.iterrows()]

    # initialize true positives and false positives
    tp_ppv = 0.0
    total_ppv = len(primary)
    tp_sens = 0.0
    total_sens = len(ct_helices)

    # determine number pairmap signals that correspond to helices (tp_ppv)
    for i, j in pm_pairs:
        matches = 0
        for x in range(-1, 4):
            if (i+x, j-x) in ct_pairs:
                matches += 1
        if matches >= 3:
            tp_ppv += 1

    # determine number of helices predicted by a pairmap signal (tp_sens)
    for helix in ct_helices.values():
        matches = 0
        for pair in helix[0:-2]:
            if pair in pm_pairs:
                matches += 1
        if matches > 0:
            tp_sens += 1

    # calculate and return sensitivity and ppv
    if len(primary) == 0:
        sens = 0
        ppv = 0
    else:
        sens = tp_sens/total_sens
        ppv = tp_ppv/total_ppv
    return sens, ppv


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


def arcPlot(ct=False, fasta=False, refct=False, probability=False, ringz=False,
            ringsig=False, pairmap=False, compare_pairmap=False, ntshape=False,
            dmsprofile=False, bottom=False, title=False, showGrid=False,
            bound=False, filternc=False, profile=False):
    # ct, Base CT file to plot. By default, the first (0th)
    #     structure is plotted. Alternative structures can be selected by
    #     passing ,# after the CT file name (e.g. --ct ctfile.ct,3)
    # fasta, Fasta sequence file
    # refct, Reference CT to compare base ct against. By
    #        default the first (0th) structure is plotted. Alternative
    #        structures can be selected by passing ,# after the CT file name
    #        (e.g. --refct ctfile.ct,3)")
    # probability, Pairing probability file in dotplot format.
    #              By default, arcs are drawn for the following probability
    #              intervals: [0.03,0.1], [0.1,0.3], [0.3,0.8], [0.8,1.0].
    #              These intervals can be modified by passing thresholds as
    #              comma-separated values. For example, --prob
    #              file.dp,0.03,0.1,0.3,0.8,1.0 specifies the default
    #              intervals. At most 4 intervals can be plotted, but fewer
    #              intervals are allowed (e.g. --prob file.dp,0.1,0.3,0.8,1.0
    #              will plot 3 intervals).
    # ringz, Plot Z-scores from ringmapper correlation file.
    #        Default color thresholds are Z=1 and Z=5. Can be modified by
    #        passing thresholds as comma-separated values (e.g. --corrz
    #        corrs.txt,1,5)
    # ringsig, Plot statistical significance from ringmapper
    #          file. Default color thresholds are [20,500]. Can be modified by
    #          passing thresholds as comma-separated values (e.g. --corrsig
    #          corrs.txt,20,500)
    # pairmap, Plot pairmap signals from pairmap file. By
    #          default plot principal & minor correlations. Can plot all
    #          complementary correlations by passing ,all (e.g. --pairmap
    #          pairmap.txt,all)
    # compare_pairmap, Plot pairmap signals from second pairmap
    #                  file. By default plot principal & minor correlations.
    #                  Can plot all complementary correlations by passing ,all
    #                  (e.g. --compare_pairmap pairmap.txt,all)
    # ntshape, Color nucs by shape reactivty in provided shape
    #          /map file
    # profile, Plot reactivity profile on top from shape/map
    #          file
    # dmsprofile, type=str, help='Normalize and plot DMS reactivity from
    #             profile file
    # bottom, action='store_true', help="Plot arcs on the bottom
    # title, type=str, default='', help='Figure title
    # showGrid, action="store_true", help="plot a grid on axis ticks
    # bound, comma separated bounds of region to plot (e.g.
    #        --bound 511,796)
    # filternc, action="store_true", help="filter out non-canonical pairs in ct
    if refct and not ct:
        exit("refct is invalid without ct")
    numplots = int(ct is not None)

    # subparse the prob argument
    if probability:
        numplots += 1
        spl = probability.split(',')
        try:
            if len(spl) > 1:
                prob_bins = list(map(float, spl[1:]))
                probability = spl[0]
        except:
            raise TypeError('Incorrectly formatted --probability'
                            ' argument {}'.format(probability))

    def subparse3(arg, name):

        outarg = arg
        bins = None

        spl = arg.split(',')
        try:
            if len(spl) == 3:
                bins = list(map(float, spl[1:]))
                outarg = spl[0]
            elif len(spl) != 1:
                raise TypeError
        except:
            raise TypeError('Incorrectly formatted {0} argument'
                            ' {1}'.format(name, arg))

        return outarg, bins

    if ringz:
        numplots += 1
        ringz, ringz_bins = subparse3(ringz, '--ringz')

    if ringsig:
        numplots += 1
        ringsig, ringsig_bins = subparse3(ringsig, '--ringsig')

    pairmap_all = False
    if pairmap:
        numplots += 1
        spl = pairmap.split(',')
        if len(spl) == 1:
            pass
        elif len(spl) == 2 and spl[1] == 'all':
            pairmap_all = True
            pairmap = spl[0]
        else:
            raise TypeError(
                'Incorrectly formatted --pairmap argument {}'.format(pairmap))

    if numplots > 2:
        exit('Too many plots! Please select at maximum 2 of [--ct,'
             ' --probability, --ringz, --ringsig, --pairmap,'
             ' --compare_pairmap]')

    # subparse the bounds argument
    if bound:
        bound = list(map(int, bound.split(',')))

    if ct:
        spl = ct.split(',')
        if len(spl) == 1:
            ctstructnum = 0
        else:
            ctstructnum = int(spl[1])
            ct = spl[0]

    if refct:
        spl = refct.split(',')
        if len(spl) == 1:
            refctstructnum = 0
        else:
            refctstructnum = int(spl[1])
            refct = spl[0]

    aplot = ap.ArcPlot(title=title, fasta=fasta)

    panel = 1
    if bottom:
        panel = -1

    if ct:
        if refct:
            aplot.compareCTs(RNA.CT(refct,
                                         structNum=refctstructnum,
                                         filterNC=filternc),
                             RNA.CT(ct, structNum=ctstructnum,
                                         filterNC=filternc),
                             panel=panel)
        else:
            aplot.addCT(RNA.CT(ct, structNum=ctstructnum,
                                    filterNC=filternc),
                        panel=panel)

        panel *= -1

    if probability:
        aplot.addPairProb(RNA.DotPlot(probability),
                          panel=panel, bins=prob_bins)
        panel *= -1

    if ringz:
        aplot.addRings(ringz, panel=panel, metric='z', bins=ringz_bins)
        panel *= -1

    if ringsig:
        aplot.addRings(ringsig, panel=panel, metric='sig', bins=ringsig_bins)
        panel *= -1

    if pairmap:

        from pmanalysis import PairMap

        aplot.addPairMap(PairMap(pairmap), panel=panel, plotall=pairmap_all)
        panel *= -1

    if compare_pairmap:

        from pmanalysis import PairMap

        aplot.addPairMap(PairMap(pairmap2), panel=panel, plotall=pairmap_all)
        panel *= -1

    # if arg.intDistance:
    #    aplot.addInteractionDistance(arg.intDistance, arg.depthThresh, panel)

    if ntshape:
        aplot.colorSeqByMAP(ntshape)

    if profile:
        aplot.readProfile(profile)

    if dmsprofile:
        aplot.readDMSProfile(dmsprofile)

    if showGrid:
        aplot.grid = True

    aplot.writePlot('show')


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
