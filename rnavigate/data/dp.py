#!/usr/bin/env python

import numpy as np
import sys
from rnavigate.data import CT


class DotPlot:
    """ Modified by Tony Mustoe to be more efficient at caclulations by
    rearranging data structure. This increases the memory footprint, but
    greatly reduces looping
    """

    def __init__(self, fIN=None):
        # if givin an input file construct the dotplot object automatically
        self.name = None
        self.length = None
        self.partfun = None
        self.sequence = None

        if fIN:
            self.name = fIN
            self.readDP(fIN)

    def __str__(self):

        entries = np.sum([len(self.partfun[x]['pair'])
                          for x in range(self.length)])/2
        a = '{{ Name= {0}, len(RNA)= {1}, entries(DotPlot)= {2} }}'.format(
            self.name, self.length, entries)
        return a

    def getMaxPair(self, nt, log10cut=1e10):
        """ uses 1-based coords """

        try:
            index = self.partfun[nt-1]['log10'].argmin()
            if self.partfun[nt-1]['log10'][index] < log10cut:
                return self.partfun[nt-1]['pair'][index]+1

        except ValueError:
            pass

        return 0

    def readDP(self, fIN):

        fhandle = open(fIN)

        # read the expected two lines of headers
        length = int(fhandle.readline().split()[0])

        if self.length is not None and length != self.length:
            sys.exit('Dotp file and dotp object are different lengths')

        if self.partfun is not None:
            # need to convert log10 and pair to lists
            # convert to numpy array
            partfun = []
            for i, v in enumerate(self.partfun):
                partfun.append({'nt': v['nt']})
                for j in ('log10', 'pair'):
                    partfun[i][j] = list(v[j])

        else:
            # make an array element for every nt
            partfun = [{'log10': [], 'pair':[], 'nt':x} for x in range(length)]

        # pop off the second header line
        fhandle.readline()

        # now read in data
        for line in fhandle:
            spl = line.split()
            i = int(spl[0])
            j = int(spl[1])
            logBP = float(spl[2])

            # i,j are anticipated to be 1-indexed
            partfun[i-1]['log10'].append(logBP)
            partfun[j-1]['log10'].append(logBP)
            partfun[i-1]['pair'].append(j-1)
            partfun[j-1]['pair'].append(i-1)

        fhandle.close()

        # convert to numpy array
        for i, v in enumerate(partfun):
            for j in ('log10', 'pair'):
                partfun[i][j] = np.array(v[j])

        self.length = length
        self.partfun = partfun

    def writeDP(self, filename):
        """write the DP object to file in DP format"""

        out = open(filename, 'w')

        out.write("%d\n" % self.length)
        out.write("i\tj\t-log10(Probability)\n")

        for nt in self.partfun:
            num = nt['nt']

            for i, pnum in enumerate(nt['pair']):

                if pnum > num:
                    # need to add one since internal indexing is 0-based
                    out.write("%d %d %f\n" % (num+1, pnum+1, nt['log10'][i]))

        out.close()

    def returnCT(self, probcut=0.5, skipConflicting=True, filterNC=True,
                 filterSingle=True):
        """return a CT object constructed from base pairs above a given
        pairing probability.
        """

        newDP = self.requireProb(probcut)

        newCT = CT()
        newCT.pair2CT(newDP.pairList(), ['n']*self.length,
                      skipConflicting=skipConflicting, filterNC=filterNC,
                      filterSingle=filterSingle)

        return newCT

    def pairList(self):
        # returns a list of base pairs i< j from the dotplot
        out = []
        for nt in self.partfun:

            mask = (nt['pair'] > nt['nt'])
            pairs = [(nt['nt']+1, x+1) for x in nt['pair'][mask]]
            if len(pairs) > 0:
                out.extend(pairs)

        return out

    def requireProb(self, minProb, maxProb=1.0):
        """return a new DP object where basepairs below a min prob have been
        trimed
        """

        maxlogBP = -np.log10(minProb)
        minlogBP = -np.log10(maxProb)

        out = DotPlot()
        out.length = self.length
        out.name = self.name

        partfun = []

        # select which nts are between a certain cutoff
        for nt in self.partfun:
            mask = (nt['log10'] <= maxlogBP) & (nt['log10'] > minlogBP)
            partfun.append({'log10': nt['log10'][mask],
                            'pair': nt['pair'][mask],
                            'nt': nt['nt']})

        out.partfun = partfun

        return out

    def trimEnds(self, trimSize, which='both'):

        which = which.lower()
        prime5 = 0
        prime3 = self.length

        if which == '5prime' or which == 'both':
            prime5 += trimSize
        if which == '3prime' or which == 'both':
            prime3 -= trimSize

        out = DotPlot()
        out.length = self.length
        out.name = self.name

        partfun = []

        for i, nt in partfun:
            if nt['nt'] < prime5 or nt['nt'] > prime3:
                partfun.append({'log10': np.array([]),
                                'pair': np.array([]),
                                'nt': nt['nt']})

            else:
                mask = (nt['pair'] >= prime5) & (nt['pair'] <= prime3)

                partfun.append({'log10': nt['log10'][mask],
                                'pair': nt['pair'][mask],
                                'nt': nt['nt']})

        out.partfun = partfun

        return out

    def calcShannon(self, printOut=False, toFile=None):

        partfun = self.partfun

        if toFile:
            outf = open(toFile, 'w')

        shannon = np.zeros(self.length)

        for i, nt in enumerate(partfun):

            nlogn = nt['log10']*10**(-nt['log10'])  # this is -log(p)*p
            ntsum = np.sum(nlogn)

            # catch rounding errors:
            if ntsum < 0:
                ntsum = 0

            shannon[i] = ntsum

            if printOut:
                print(nt['nt'], ntsum)

            if toFile:
                outf.write("%s\t%.3e\n" % (nt['nt'], ntsum))

        if toFile:
            outf.close()

        return np.array(shannon)

    def partfunDifference(self, comp, region=None):
        """ region is 1-based, inclusive; i.e. nt numbering"""

        if region is None:
            region = [1, self.length]

        sumdiff, numnts = 0.0, 0
        sarr = np.zeros(self.length)
        carr = np.zeros(self.length)

        for i in range(region[0]-1, region[1]):

            snt = self.partfun[i]
            cnt = comp.partfun[i]
            sarr[:] = 0
            carr[:] = 0

            for pindex, j in enumerate(snt['pair']):
                sarr[j] = 10**(-snt['log10'][pindex])

            for pindex, j in enumerate(cnt['pair']):
                carr[j] = 10**(-cnt['log10'][pindex])

            diff = sarr-carr
            d1 = np.sum(diff[diff > 0])
            diff = carr-sarr
            d2 = np.sum(diff[diff > 0])

            # get the geometeric average for each nt.
            sumdiff += np.sqrt(d1*d2)
            numnts += 1

        return sumdiff/numnts

    def hack_Gunfold(self, region):
        """This is a hacky (i.e. very unrigorous) but quick way to estimate
        unfolding free energy for a specified region of nts.
        Region is 1-based, inclusive; i.e. nt numbering
        """

        if self.sequence is None:
            AttributeError("Sequence is not defined")

        bpenergy = {'AT': 1, 'TA': 1, 'GC': 2, 'CG': 2, 'TG': 1, 'GT': 1}

        unfold = 0

        # shift each by -1 since want to convert to 0-based index,
        # and then +1 to region[1]
        for i in range(region[0]-1, region[1]):

            nt = self.partfun[i]
            seqi = self.sequence[i]

            for pindex, j in enumerate(nt['pair']):

                # make sure that each pair is only counted once.
                # If k < i, but outside the region,
                # we don't have to worry about double counting
                if region[0] <= j < i:
                    continue

                seqj = self.sequence[j]
                try:
                    eij = bpenergy[seqi+seqj]
                except KeyError:
                    print(seqi, seqj, i, j)
                    eij = 0

                # this is pij * Gij
                unfold += 10**(-nt['log10'][pindex]) * eij

        return unfold

    def calcPairProb(self):
        """ return the pairing probability of each nt"""
        partfun = self.partfun

        prob = np.zeros(self.length)

        # do sum the niave way since we aren't worried about small prob diffs
        # with this calculation
        # Note that it is likely numerically unstable for small probs though...
        for i, nt in enumerate(partfun):
            prob[i] = np.sum(10**(-nt['log10']))

        return np.array(prob)

    def maxProb(self):
        """return the max pairing prob of each nt"""

        partfun = self.partfun
        prob = np.zeros(self.length)

        for i, nt in enumerate(partfun):
            try:
                prob[i] = 10**(-1*np.min(nt['log10']))
            except ValueError:
                prob[i] = 0.0

        return np.array(prob)

    def averageSlippedBPs(self, struct=None, predictedOnly=True):
        """replaces PlusandMinus script. If a helix in a predicted structure is
        slipped +/-1 nt, we need to sum the predicted probabilities.Otherwise,
        the predicted Shannon entropy will be artificially high.

        turning off predicted only will go through all i<j basepairs and merge
        them in preference of liklihood. This is more compuationally intensive
        """

        dotPlot = self.dp
        # dotPlotCopy = {'logBP':copy.deepcopy(dotPlot['logBP'])}

        # this is the value in -log10(prob), 2 = a prob of 0.01
        slippedCutoff = 2
        slipped = []

        # if a reference structure is given, merge pairs to it first
        if struct:
            for pair in range(1, len(struct.ct)-1):
                # define the base pairs
                pair_i = pair+1
                pair_j = struct.ct[pair]

                # skip non-paired nucleotides
                if pair_j == 0:
                    continue
                # skip pairs i > j so we don't double count
                if pair_j < pair_i:
                    continue

                # create some search filters
                filter_i = dotPlot['i'] == pair_i
                filter_j = dotPlot['j'] == pair_j

                filter_i_before = dotPlot['i'] == pair_i - 1
                filter_i_after = dotPlot['i'] == pair_i + 1

                filter_j_before = dotPlot['j'] == pair_j - 1
                filter_j_after = dotPlot['j'] == pair_j + 1

                # find i,j union before, at, and after
                filterDP = {}
                filterDP['before_j'] = filter_j_before * filter_i
                filterDP['before_i'] = filter_j * filter_i_before

                filterDP['after_j'] = filter_j_after * filter_i
                filterDP['after_i'] = filter_j * filter_i_after

                # define current point
                at = filter_j * filter_i

                # handle slippage, first define base prob
                prob_at = 10**(-dotPlot['logBP'][at])

                # cycle through all filter combinations
                for filterPair in list(filterDP.keys()):

                    # shorthand variable for current filter
                    curr = filterDP[filterPair]

                    # if pair exists ...
                    if np.sum(curr) == 1:
                        # add it to predicted pair probability
                        prob_at += 10**(-dotPlot['logBP'][curr])

                        # add pair to slipped list if it meets slip criteria
                        if dotPlot['logBP'][curr] < slippedCutoff:
                            slipped.append((pair_i, pair_j))

                        # now set it to a very low probability for zero
                        dotPlot['logBP'][curr] = 50

                # return to -log10 format
                dotPlot['logBP'][at] = -np.log10(prob_at)

        if not predictedOnly:
            # go through all i<j basepair combinations and check to see if
            # there is a slipped base pair
            for i, j in self.pairList():

                # correct for python counting starting at 0
                pair_i, pair_j = int(i+0), int(j+0)

                # see if there exists a basepair for this combination
                filter_i = (dotPlot['i'] == pair_i)
                filter_j = (dotPlot['j'] == pair_j)
                filter_union = (filter_i * filter_j)
                paired = np.sum(filter_union) == 1
                union_logBP = dotPlot['logBP'][filter_union]
                # only continue if the pair is reasonably likely
                # (>1% chance of forming)
                if paired and union_logBP < 3:
                    filterList = {}

                    # define the various types of slippage
                    filter_ibefore = (dotPlot['i'] == pair_i-1)
                    filter_jbefore = (dotPlot['j'] == pair_j-1)
                    filter_iafter = (dotPlot['i'] == pair_i+1)
                    filter_jafter = (dotPlot['j'] == pair_j+1)

                    # index filters to a dict
                    filterList['before_i'] = filter_ibefore * filter_j
                    filterList['before_j'] = filter_i * filter_jbefore
                    filterList['after_i'] = filter_iafter * filter_j
                    filterList['after_j'] = filter_i * filter_jafter

                    # define the prob at in normal normal space
                    prob_at = 10**(-dotPlot['logBP'][filter_union])

                    # go through each of the filters
                    for pairFilter in list(filterList.keys()):
                        curr = filterList[pairFilter]

                        # if the current pair exists in the dotplot
                        if np.sum(curr) == 1:
                            # check if it's less likely than current pairs
                            if union_logBP < dotPlot['logBP'][curr]:
                                # and add it to the current pair if it is
                                prob_at += 10**(-dotPlot['logBP'][curr])

                                # set to a very low probabliity afterwards
                                dotPlot['logBP'][curr] = 50
                    dotPlot['logBP'][filter_union] = -np.log10(prob_at)

        return slipped

    def pairingProb(self):
        pass
