import arcPlot as ap
import RNAtools2 as RNAtools
import numpy as np


def corrsToBitmap(corrfile, ctfile, window=1, limit=100):
    size = ctLength(ctfile)
    bitmap = np.zeros((size, size))
    with open(corrfile) as f:
        f.readline()
        f.readline()
        for line in f.readlines():
            line = [item.strip() for item in line.split()]
            i = int(line[0])
            j = int(line[1])
            sig = float(line[2])*int(line[3])
            bitmap[i:i+window, j:j+window] = sig
    bitmap[bitmap > limit] = limit
    bitmap[bitmap < -limit] = -limit
    bitmap = bitmap/limit
    bitmap[0, 0] = 1
    bitmap[0, 1] = -1
    bitmap += 1
    bitmap /= 2
    return bitmap


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


def ctToBitmap(ctfile):
    size = ctLength(ctfile)
    bitmap = np.ones((size, size))
    with open(ctfile) as f:
        f.readline()
        for line in f.readlines():
            line = [item.strip() for item in line.split()]
            i = int(line[0])
            j = int(line[4])
            if j != 0 and i < j:
                bitmap[i, j] = 0
                nt = 0
                helix = True
                while helix:
                    if bitmap[i-nt-1, j+nt+1] == 0:
                        nt += 1
                    else:
                        helix = False
                bitmap[i-nt:i+1, j:j+nt+1] = 0
    return bitmap


def pairmapToBitmap(pairfile, ctfile):
    size = ctLength(ctfile)
    bitmap = np.zeros((size, size))+0.5
    classdict = {"0": 0.5,  # null (white with bwr)
                 "1": 0,  # primary (blue with bwr)
                 "2": 1}  # secondary (red with bwr)
    with open(pairfile) as f:
        f.readline()
        f.readline()
        for line in f.readlines():
            line = line.strip().split()
            i = int(line[0])
            j = int(line[1])
            Class = line[3]
            bitmap[i:i+3, j:j+3] = classdict[Class]
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
            aplot.compareCTs(RNAtools.CT(refct,
                                         structNum=refctstructnum,
                                         filterNC=filternc),
                             RNAtools.CT(ct, structNum=ctstructnum,
                                         filterNC=filternc),
                             panel=panel)
        else:
            aplot.addCT(RNAtools.CT(ct, structNum=ctstructnum,
                                    filterNC=filternc),
                        panel=panel)

        panel *= -1

    if probability:
        aplot.addPairProb(RNAtools.DotPlot(probability),
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
