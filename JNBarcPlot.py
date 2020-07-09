import arcPlot as ap
import RNAtools2 as RNA
import numpy as np
import pandas as pd


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
    ct = pd.read_csv('RMRP.ct', sep='\s+', names=names,
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
