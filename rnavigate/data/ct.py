#!/usr/bin/env python

###############################################################################
#  CT object code comes from RNAtools2.py and Superfold
#      Main contributors:
#           Gregg Rice
#           Anthony Mustoe
#           Tom Christy
#           Patrick Irving
###############################################################################

import sys
import numpy as np
import xml.etree.ElementTree as xmlet
from .data import Data
import pandas as pd
import json


def get_ss_class(filepath, **kwargs):
    if filepath is None:
        return None
    extension = filepath.split('.')[-1]
    instantiator = {
        "varna": VARNA,
        "xrna": XRNA,
        "nsd": NSD,
        "cte": CTE,
        "ct": CT,
        "dbn": DotBracket,
        "bracket": DotBracket,
        "json": JSON,
    }[extension]
    return instantiator(filepath=filepath, **kwargs)


class CT(Data):

    def __init__(self, datatype="ct", dataframe=None, filepath=None, **kwargs):
        """
        if givin an input file .ct construct the ct object automatically
        """
        self.datatype = datatype
        if dataframe is not None:
            self.from_dataframe(dataframe)
            self.filepath = "dataframe"
        elif filepath is not None:
            self.filepath = filepath
            self.read_file()

    def __str__(self):
        """overide the default print statement for the object"""
        a = '{ Name= %s, len(CT)= %s }' % (self.filepath, str(len(self.ct)))
        return a

    def from_dataframe(self, dataframe):
        columns = ["Nucleotide", "Sequence", "Pair"]
        for col in columns:
            assert col in dataframe.columns, f"{col} column missing."
        seq = [nt.upper() for nt in dataframe["Sequence"]]
        if "T" in seq:
            print("T nucleotides have been recoded as U")
            seq = ["U" if nt == "T" else nt for nt in seq]
        self.sequence = "".join(seq)
        self.num = dataframe["Nucleotide"].values
        self.ct = dataframe["Pair"].values

    def read_file(self, structNum=0, filterNC=False, filterSingle=False):
        """Loads CT information from a given ct file. Requires a header!

        Args:
            structNum (int, optional): If ct file contains multiple structures,
                uses the given structure. Defaults to 0.
            filterNC (bool, optional): If True, will filter out non-canonical
                base pairs. Defaults to False.
            filterSingle (bool, optional): If True, will filter out any
                singleton base pairs. Defaults to False.
        """
        fIN = self.filepath
        num, seq, bp, mask = [], '', [], []

        try:
            with open(fIN) as f:

                # process the first line (pops off the first element)
                spl = f.readline().split()
                numnucs = int(spl[0])
                header = ' '.join(spl[1:])
                curstruct = 0

                for i, line in enumerate(f):
                    spl = line.split()
                    if (i+1) % (numnucs+1) == 0:
                        # catch the next header
                        curstruct += 1

                        if curstruct > structNum:
                            break

                        header = ' '.join(spl[1:])
                        continue

                    if curstruct == structNum:
                        num.append(int(spl[0]))
                        seq += str(spl[1])
                        bp.append(int(spl[4]))

                        # check if there is masking info
                        if len(spl) == 7 and spl[6] == '1':
                            mask.append(1)
                        else:
                            mask.append(0)

        except Exception:
            raise IOError(
                "{0} has invalid format or does not exist".format(fIN))

        if len(num) == 0:
            sys.exit("Structure %d was not found in the ct file" % structNum)

        if 'T' in seq:
            print("Note: T nucleotides have been recoded as U")
            seq = seq.replace('T', 'U')
            seq = seq.replace('t', 'u')

        # check consistency!
        for i in range(len(bp)):
            if bp[i] != 0:
                p1 = (i+1, bp[i])
                p2 = (bp[bp[i]-1], bp[i])
                if p1 != p2:
                    print("WARNING: Inconsistent pair "
                          f"{p1[0]}-{p1[1]} vs. {p2[0]}-{p2[1]}")

        self.header = header
        self.num = num
        self.sequence = seq
        self.ct = np.array(bp)
        self.mask = mask

        if filterNC:
            self.filterNC()
        if filterSingle:
            self.filterSingleton()

    def get_dbn(self):
        for pair in self.pairList():
            # check for pseudoknots
            for pair2 in self.pairList():
                if (pair[0] < pair2[0] < pair[1]) and (pair[1] < pair2[1]):
                    print("get_dbn doesn't support PKs, returning empty dbn")
                    return '.' * self.length
        dbn = ''
        for i, pair in enumerate(self.ct):
            if pair == 0:
                dbn += '.'
            elif pair > (i+1):
                dbn += '('
            elif pair < (i+1):
                dbn += ')'
        return dbn

    def filterNC(self):
        """
        Filter out non-canonical basepairs from the ct datastructure
        """

        for i, nt in enumerate(self.ct):

            if nt == 0:
                continue

            pi = self.num.index(nt)

            bp = self.sequence[pi] + self.sequence[i]
            if bp not in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG', '  '):
                self.ct[i] = 0
                self.ct[pi] = 0
                # print 'Deleted %d %s' % (self.num[i], self.sequence[i])
                continue

    def filterSingleton(self):
        """
        Filter out singleton basepairs from the ct datastructure
        """

        for i, nt in enumerate(self.ct):

            if nt == 0:
                continue

            pi = self.num.index(nt)

            neigh = False
            for j in (-1, 1):
                try:
                    neigh = (neigh or (nt-self.ct[i+j] == j))
                except IndexError:
                    pass

            if not neigh:
                self.ct[i] = 0
                self.ct[pi] = 0
                # print 'Deleted %d %s *' % (self.num[i], self.sequence[i])

    def writeCT(self, fOUT, writemask=False):
        """
        writes a ct file from the ct object
        writemask = True will write out any masking info
        """

        try:
            mask = (writemask and (len(self.mask) == len(self.ct)))
        except AttributeError:
            mask = False

        # handle empty ct object case
        if not self.ct:
            print("empty ct object. Nothing to write")
            return

        w = open(fOUT, 'w')
        w.write('{0:6d} {1}\n'.format(len(self.num), self.filepath))
        for i in range(len(self.num)):
            nt = self.num[i]
            prev = nt-1
            next = nt+1 % (len(self.ct)+1)  # last nt goes back to zero
            seq = self.sequence[i]
            ct = self.ct[i]
            if mask and self.mask[i]:
                line = f'{nt:5d} {seq} {prev:5d} {next:5d} {ct:5d} {nt:5d} 1\n'
            else:
                line = f'{nt:5d} {seq} {prev:5d} {next:5d} {ct:5d} {nt:5d}\n'
            w.write(line)
        w.close()

    def copy(self):
        """
        returns a deep copy of the ct object
        """
        out = CT()
        out.filepath = self.filepath[:]
        out.num = self.num[:]
        out.sequence = self.sequence[:]
        out.ct = self.ct[:]
        return out

    def pairList(self):
        """
        # returns a list of base pairs i<j as a array of tuples:
        [(19,50),(20,49)....]
        """
        out = []
        for nt in range(len(self.ct)):
            if self.ct[nt] != 0 and self.ct[nt] > self.num[nt]:
                out.append((self.num[nt], self.ct[nt]))
        return out

    def pairedResidueList(self, paired=True):
        """Returns a list of residues that are paired.

        Args:
            paired (bool, optional): It False, returns single-stranded list.
                Defaults to True.
        """

        out = []
        for i, nt in enumerate(self.ct):
            if (nt != 0) == paired:
                out.append(self.num[i])

        return out

    def junctionResidues(self, incGU=False):
        """Returns a list of residues at junctions.

        Args:
            incGU (bool, optional): If True, includes GU pairs.
                Defaults to False.
        """

        jun = []
        for i, nt in enumerate(self.ct):

            if nt == 0:  # if its not paired
                continue

            pi = self.num.index(nt)

            neighbors = [i-1, i+1, pi-1, pi+1]
            neighbors_paired = any(self.ct[neighbors] == 0)
            neighbors_not_in_chain = any(
                [x not in self.num for x in neighbors])
            if neighbors_paired or neighbors_not_in_chain:
                jun.append(self.num[i])

        return jun

    def addPairs(self, pairs):
        """Add base pairs to current ct file
        pairs should be a list of pairs (1-indexed)
        """

        for i, j in pairs:

            if self.ct[i-1] != 0:
                print(('Warning: nt {} is already paired!!!'.format(i)))
            if self.ct[j-1] != 0:
                print(('Warning: nt {} is already paired!!!'.format(j)))

            self.ct[i-1] = j
            self.ct[j-1] = i

    def pair2CT(self, pairs, seq=None, name=None, skipConflicting=True,
                filterNC=False, filterSingle=False):
        """
        constructs a ct object from a list of base pairs and a sequence

        pairs are an array of bp tuples ( i < j )
           i.e. [(4,26),(5,25),...]
        length is implied from the given sequence
        """

        # See if sequence has been defined
        if seq is None:
            assert hasattr(self, "sequence"), "Sequence is not defined"
            assert self.sequence is not None, "Sequence is not defined"
        else:
            self.sequence = seq

        length = len(self.sequence)
        self.num = list(range(1, length+1))

        # give it a name if it has one
        if name:
            self.filepath = name
        else:
            self.filepath = 'RNA_'+str(length)

        self.ct = []
        for i in range(length):
            self.ct.append(0)

        for i, j in pairs:
            if self.ct[i-1] != 0:
                print('Warning: conflicting pairs, (%s - %s) : (%s - %s)' %
                      (str(i), str(j), str(self.ct[i-1]), str(i)))
                if skipConflicting:
                    continue
            if self.ct[j-1] != 0:
                print('Warning: conflicting pairs, (%s - %s) : (%s - %s)' %
                      (str(i), str(j), str(j), str(self.ct[j-1])))
                if skipConflicting:
                    continue
            self.ct[i-1] = j
            self.ct[j-1] = i

        if filterNC:
            self.filterNC()

        if filterSingle:
            self.filterSingleton()

    def getNTslice(self, start=None, end=None):
        offset = self.num[0]

        try:
            start = start - offset
        except TypeError:
            assert start is None

        try:
            end = end+1-offset
        except TypeError:
            assert end is None

        sel = slice(start, end)

        return sel

    def maskCT(self, start, end, nocross=True, inverse=False):
        """return a copy of the original CT file, but with the specified region
        masked out (ie base pairs set to 0)
        if inverse = True, do the inverse -- ie return a CT with base pairs
        only involving the specified region
        """

        # in case numbering differs from indexing, go by indexes
        out = self.copy()
        out.filepath += '_mask_'+str(start)+'_'+str(end)

        sele = self.getNTslice(start, end)
        offset = self.num[0]

        if inverse:
            out.ct[:] = [0]*len(out.ct)
            out.ct[sele] = self.ct[sele]
            for nt in self.ct[sele]:
                out.ct[nt-offset] = self.ct[nt-offset]

            return out

        # below code will only execute if not inverse

        # zero out the masked region
        out.ct[sele] = [0]*len(out.ct[sele])
        for nt in self.ct[sele]:
            if nt != 0:
                out.ct[nt-offset] = 0

        if nocross:
            for i in range(sele.start):
                nt = out.ct[i]-offset
                if nt >= sele.stop:
                    out.ct[i] = 0
                    out.ct[nt] = 0

        return out

    def cutCT(self, start, end):
        """Returns a new ct file containing only base pairs within the
        specified region
        """

        sel = self.getNTslice(start, end)

        offset = self.num[0]
        numnts = end-start+1

        out = CT()
        out.sequence = self.sequence[sel]
        out.num = list(range(1, numnts+1))
        out.filepath = self.filepath + '_cut_'+str(start)+'_'+str(end)

        out.ct = []
        temp = self.ct[sel]
        # renumber from 1
        for nt in temp:
            nt_out = nt-(start-offset)
            # cut out pairings that lie outside the window
            if nt_out <= 0 or nt_out > numnts:
                nt_out = 0
            out.ct.append(nt_out)

        return out

    def stripCT(self):
        """
        returns an array the length of the ct object
        in a non-redundant form - e.g. only gives pairs i<j
        """
        pairs = self.pairList()
        halfPlexCT = np.zeros_like(self.ct)
        for i, j in pairs:
            halfPlexCT[i-1] = j
        return halfPlexCT

    def ctToArcList(self):
        """
        returns a 2D array containing arc list elements
        """

        pairs = self.stripCT()

        outarr = []
        for i, v in enumerate(pairs):
            if v != 0:
                temp = np.zeros_like(pairs)
                temp[i] = v
                outarr.append(temp)
        return outarr

    def get_distance_matrix(self, recalculate=False):
        """Based on Tom's contactDistance function above, but instead returns
        the all pairs shortest paths matrix, and stores it as an attribute. If
        the attribute has already been set, it returns the attribute. This is
        faster than calling contactDistance pairwise to fill the matrix.

        Args:
            recalculate (bool, optional): Set to true to recalculate the matrix
                even if the attribute is set. In case changes to the ct have
                been made.
        """
        if hasattr(self, "distance_matrix") and not recalculate:
            return self.distance_matrix
        # this method will be used later to make sure a nucleotide hasn't
        # been visited and is within the bounds of the RNA

        def viable(nt):
            not_empty = (nt != 0)
            valid_so_far = not_empty
            if valid_so_far:
                in_range = (min(self.num) <= nt <= max(self.num))
                valid_so_far = in_range
            if valid_so_far:
                not_visited = (level[self.num.index(nt)] == -1)
                valid_so_far = not_visited
            return valid_so_far

        # The default value of -1 means that a nt hasn't been visited.
        distance_matrix = np.full((len(self.num), len(self.num)), -1)
        for i, nt in enumerate(self.num):
            # Level keeps track of how far each other nt is from nt
            level = distance_matrix[i, :]
            # create the queue and add nt, set its level as 0
            queue = list()
            queue.append(nt)
            level[i] = 0
            # while there are nucleotides that haven't been visited
            while len(queue) > 0:
                current_nt = queue.pop(0)
                current_i = self.num.index(current_nt)

                # Find all the neighbors of the current NT
                NTBack = current_nt - 1
                NTUp = current_nt + 1
                NTPair = self.ct[current_i]
                # check all the neighbors to see if they should be added
                # and if so, record their level(ie their Contact Distance)
                for neighbor in [NTBack, NTUp, NTPair]:
                    if viable(neighbor):
                        queue.append(neighbor)
                        level[self.num.index(neighbor)] = level[current_i] + 1
            # set row of distance matrix for nt
            distance_matrix[i, :] = level
        # store the distance matrix from the search and return
        self.distance_matrix = distance_matrix
        return self.distance_matrix.copy()

    def contactDistance(self, i, j):
        """
        Caclulates the contact distance between pairs i,j in
        the RNA using the RNAtools CT Object. This is different
        than the old contact distance function because it utilizes a breadth
        first search to determine the distance. While BFS may take longer
        than contact distance, it will always find the shortest distance.

        Added by Tom Christy
        """
        if hasattr(self, "distance_matrix"):
            i_index = self.num.index(i)
            j_index = self.num.index(j)
            return self.distance_matrix[i_index, j_index]

        # this method will be used later to make sure a nucleotide hasn't
        # been visited and is within the bounds of the RNA
        def viable(nt, rna):
            addMe = True
            # don't add if nt is the pair of an unpaired nt (0 in the ct)
            if(nt == 0):
                addMe = False
            # don't add if nt is smaller than the lowest nt in the ct
            elif(nt < rna.num[0]):
                addMe = False
            # don't add if nt is larger than the highest nt in the ct
            elif(nt > rna.num[-1]):
                addMe = False
            # don't add if this nt has already been visited
            elif(level[nt] >= 0):
                addMe = False

            return addMe

        # create an array to keep track of how far each nt is from i
        # The default value of -1 also means that an nt hasn't been visited.
        level = np.zeros(len(self.num) + 1)
        level = level - 1
        # each index matches to the nt number,
        # so index 0 is just an empty place holder

        # create the queue and add nt i, set its level as 0
        queue = list()
        queue.append(i)
        level[i] = 0
        # while the queue has members,
        # keep popping and searching until the NT is found
        notFound = True
        contactDistance = 0
        while len(queue) > 0 and notFound:
            currentNT = queue.pop(0)
            # if the current nucleotide is the one we're searching for,
            # record it's level as the contact distance
            # and break out of the search
            if(currentNT == j):
                contactDistance = level[currentNT]
                notFound = False
                break

            # grab all the neighbors of the current NT
            # and determine if they are inbounds and unvisited
            # If so, add them to the queue and set their levels.
            NTBack = currentNT - 1
            NTUp = currentNT + 1
            # get the pair of NT -1 b/c of the weird way coded the ct file,
            # it's off by 1
            NTPair = self.ct[currentNT-1]

            # check all the neighbors to see if they should be added to the
            # search, and if so, record their level(ie their Contact Distance)
            if viable(NTBack, self):
                queue.append(NTBack)
                level[NTBack] = level[currentNT] + 1
            if viable(NTUp, self):
                queue.append(NTUp)
                level[NTUp] = level[currentNT] + 1
            if viable(NTPair, self):
                queue.append(NTPair)
                level[NTPair] = level[currentNT] + 1

        # return the contactDistance from the search
        return contactDistance

    def extractHelices(self, fillPairs=True, splitHelixByBulge=True):
        """Returns a list of helices in a given CT file and the nucleotides
        making up the list as a dict object. e.g. {0:[(1,50),(2,49),(3,48)}

        defaults to filling 1,1 and 2,2 mismatches

        Since April 20, 2018 an option has been added to no longer split
        helices by single nucleotide bulges. Prior to this option bulges on the
        5' of the helix split the helix but those on the 3' end did not.
        """
        # first step is to find all the helices
        rna = self.copy()
        if fillPairs:
            rna = rna.fillPairs()
        helices = {}
        nt = 0
        heNum = 0

        while nt < len(rna.ct):
            if rna.ct[nt] == 0:
                nt += 1
                continue
            else:
                # skip dups
                if nt > rna.ct[nt]:
                    nt += 1
                    continue

                tempPairs = []
                stillHelix = True

                previous = rna.ct[nt]

                # this is the default behavior which splits helices that are
                # separated by a single nucleotide bulge
                if splitHelixByBulge:
                    while stillHelix:
                        # see if the helix juts up against another one
                        if abs(rna.ct[nt] - previous) > 1:
                            break
                        # add the pairs
                        elif rna.ct[nt] != 0:
                            tempPairs.append((nt+1, rna.ct[nt]))
                            previous = rna.ct[nt]
                            nt += 1
                        else:
                            break

                else:
                    while stillHelix:
                        # Handle bulges on the 3' end of the helix
                        pair = rna.ct[nt]
                        if pair != 0 and not(abs(pair - previous) > 2):
                            tempPairs.append((nt+1, pair))
                            previous = pair
                            nt += 1
                        # Handle bulges on the 5' end of the helix
                        elif rna.ct[nt+1] == rna.ct[nt-1]-1:
                            nt += 1
                        else:
                            break

                # remove single bp helices
                if len(tempPairs) <= 1:
                    continue
                helices[heNum] = tempPairs
                heNum += 1
        return helices

    def fillPairs(self, mismatch=1):
        """
        fills 1,1 and/or 2,2 mismatches in an RNA structure
        mismatch = 1 will fill only 1,1 mismatches
        mismatch = 2 will fill both
        """

        rna = self.copy()
        # fill in 1,1 mismatch, 2,2 mismatch
        for i in range(len(rna.ct)-3):
            if rna.ct[i+1] == 0:
                if rna.ct[i] - rna.ct[i+2] == 2:
                    rna.ct[i+1] = rna.ct[i] - 1
            if mismatch == 2 and rna.ct[i+1] + rna.ct[i+2] == 0:
                if rna.ct[i] - rna.ct[i+3] == 3:
                    rna.ct[i+1] = rna.ct[i] - 1
                    rna.ct[i+2] = rna.ct[i] - 2
        return rna

    def extractPK(self, fillPairs=True):
        """
        returns the pk1 and pk2 pairs from a CT file. Ignores single base
        pairs. PK1 is defined as the helix crossing the most other bps.
        if there is a tie, the most 5' helix is called pk1

        Code from original RNAtools was updated to more robustly assign pk1/pk2

        returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
        """

        def checkOverlap(h1, h2):
            # only need to check one set of pairs from each
            # of the helices. Test is to see if they form
            # a cross hatching pattern
            h1_max = max(h1[0])
            h1_min = min(h1[0])
            h2_max = max(h2[0])
            h2_min = min(h2[0])
            if h1_max > h2_min > h1_min and h2_max > h1_max:
                return True
            if h2_max > h1_min > h2_min and h1_max > h2_max:
                return True

            return False

        # make a copy so we don't destroy the original object
        rna = self.copy()

        # get the helices by calling the extract helix function
        helices = rna.extractHelices(fillPairs=fillPairs)
        heNum = len(helices)

        # do the helices overlap? Check for a crosshatching pattern
        # append them to a list if they have it.
        overlaps = []  # stores the helix number

        for i in range(heNum):
            for j in range(i+1, heNum):
                if checkOverlap(helices[i], helices[j]):
                    overlaps.append((i, j))

        # if there are no overlapping bps, return none
        if len(overlaps) == 0:
            return None, None

        # iterate through all pairs of overlapping helices and count
        # the number of pairs in each to determine helix crossing most bps
        crossbps = {}
        for i, j in overlaps:
            if i not in crossbps:
                crossbps[i] = 0
            if j not in crossbps:
                crossbps[j] = 0
            crossbps[i] += len(helices[j])
            crossbps[j] += len(helices[i])

        # convert to list and sort
        crossbps = list(crossbps.items())

        # first sort list by helix number -- i.e. 5' to 3'
        crossbps.sort(key=lambda x: x[0])
        # now sort by which helix crosses the most BPs. This two-step procedure
        # ensures that in the case of ties, we always choose 5'-most helix
        crossbps.sort(key=lambda x: x[1], reverse=True)

        pk1 = helices[crossbps[0][0]]
        pk2 = []
        for k, v in crossbps[1:]:
            pk2.extend(helices[k])

        return pk1, pk2

    def padCT(self, referenceCT, giveAlignment=False):
        """Aligns the target CT to the reference CT and pads the referece
        CT file with 000s in order to input into CircleCompare"""
        out = CT()
        out.sequence = referenceCT.seq
        out.num = referenceCT.num

        # align target to reference
        seed = 200
        if len(self.seq) <= seed:
            seed = len(self.seq) - 1
        pos = 0
        maxScore = 0
        # print len(referenceCT.seq)
        for i in range(len(referenceCT.seq)-seed):
            a, b = referenceCT.seq[i:i+seed], self.seq[:seed]
            s = 0
            # s = # of identical nts across the alignment
            for k, l in zip(a, b):
                if k == l:
                    s += 1
            if s == seed:
                pos = i
                maxScore += 1
        # handle the exception when target and reference do not match
        if maxScore != 1:
            print('reference and target do not match <EXIT>')
            sys.exit()

        # create the renumbered ct to fit within the reference
        ct = []
        for i in range(len(referenceCT.seq)):
            # if the target falls within the range of the reference ct
            #     then change the numbers
            # else pad the files with 000's
            if i >= pos and i < pos+len(self.seq):
                val = self.ct[i-pos]
                if val > 0:
                    val += pos
                ct.append(val)
            else:
                ct.append(0)

        # set to the new ct file and return it
        out.ct = ct
        out.filepath = self.filepath+'_renum_'+str(pos)
        if giveAlignment:
            return out, pos
        else:
            return out

    def computePPVSens(self, compCT, exact=True, mask=False):
        """Compute the PPV and sensitivity between self and another CT.

        Args:
            compCT (CT object): The CT to compare to.
            exact (bool, optional): True requires BPs to be exactly correct.
                                    False allows +/-1 bp slippage.
                                    Defaults to True.
            mask (bool, optional): Exclude masked regions from calculation.
                                   Defaults to False.

        Returns:
            float, float, tuple: PPV, Sensitivity, (TP, TP+FP, TP+FN)
        """
        # check mask is properly set up if using
        if mask:
            assert hasattr(self, "mask"), "Mask has not been defined."
            assert len(self.mask) != len(self.ct), "Mask is incorrect length"

        assert len(self.ct) == len(compCT.ct), 'CT objects are different sizes'

        # compute totals
        totr, totc = 0, 0
        for i, v in enumerate(self.ct):

            if mask and self.mask[i]:
                continue
            if v != 0:
                totr += 1
            if compCT.ct[i] != 0:
                totc += 1

        totr /= 2
        totc /= 2

        sharedpairs = 0
        for i, v in enumerate(self.ct):

            if v == 0 or v < i or (mask and self.mask[i]):
                continue

            try:
                if v == compCT.ct[i]:
                    sharedpairs += 1
                elif not exact and (v in [compCT.ct[i]-1, compCT.ct[i]+1]):
                    sharedpairs += 1
                elif not exact and (v in [compCT.ct[i-1], compCT.ct[i+1]]):
                    sharedpairs += 1

            except IndexError:
                pass

        try:
            ppv = sharedpairs/float(totc)
            sens = sharedpairs/float(totr)
        except ZeroDivisionError:
            ppv, sens = 0.0, 0.0

        return sens, ppv, (sharedpairs, totc, totr)

    def writeSTO(self, outfile, name='seq'):
        """"write structure file out into Stockholm (STO) file format
        for use for infernal searches"""

        with open(outfile, 'w') as out:

            out.write('# STOCKHOLM 1.0\n\n')

            namelen = max(len(name), 12)

            out.write('{0} '.format(name.ljust(namelen)))
            for nt in self.sequence:
                out.write(nt)
            out.write('\n')

            out.write('{0} '.format('#=GC SS_cons'.ljust(namelen)))

            for nt, pair in enumerate(self.ct):
                if pair == 0:
                    out.write('.')
                elif pair > nt:
                    out.write('(')
                else:
                    out.write(')')

            out.write('\n//\n')

    def writeRNAstructureSeq(self, writename):
        with open(writename, 'w') as out:
            out.write(';\nSequence from {0}\n'.format(self.filepath))
            out.write('{0}1'.format(''.join(self.sequence)))

    def write_cte(self, outputPath):
        """writes the current structure out to CTE format for Structure Editor.

        Args:
            outputPath (string): path to output cte file to be created
        """
        # set scaling factors based on data source.
        xscale = {'xrna': 1.525*20, 'varna': 0.469*65,
                  'nsd': 30.5, 'cte': 30.5}[self.ss_type]
        yscale = {'xrna': -1.525*20, 'varna': 0.469*65,
                  'nsd': 30.5, 'cte': 30.5}[self.ss_type]
        w = open(outputPath, 'w')
        w.write('{0:6d} {1}\n'.format(self.length, self.filepath))
        line = ('{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d} ' +
                ';! X: {5:1f} Y: {6:1f}\n')
        for i in range(self.length):
            xcoord = round(xscale*self.xcoordinates[i])
            ycoord = round(yscale*self.ycoordinates[i])
            num, seq, cti = self.num[i], self.sequence[i], self.ct[i]
            # no nums after self.length, resets to zero
            nextnum = (num+1) % (self.length+1)
            cols = [num, seq, num-1, nextnum, cti, xcoord, ycoord]
            w.write(line.format(*cols))
        w.close()

    def fit_to(self, fit_to):
        am = self.get_alignment_map(fit_to=fit_to)
        ct = [0] * fit_to.length
        for i, j in self.pairList():
            new_i = am[i-1]+1
            new_j = am[j-1]+1
            if (new_i != 0) and (new_j != 0):
                ct[new_i] = new_j
                ct[new_j] = new_i
        data = pd.DataFrame({
            'Nucleotide': [i+1 for i in range(fit_to.length)],
            'Sequence': [nt for nt in fit_to.sequence],
            'Pair': ct
        })
        self.fitted = CT(dataframe=data)

    def get_ij_colors(self, compct=None, fitted=True):
        if fitted and hasattr(self, 'fitted'):
            ct1 = self.fitted
            if compct is not None and hasattr(compct, 'fitted'):
                ct2 = compct.fitted
            else:
                ct2 = compct
        else:
            ct1 = self
            ct2 = compct
        if ct2 is not None:
            assert ct1.length == ct2.length, "CT objects are not equal lengths"

        i_list, j_list, colors = [], [], []

        def add_ij_color(i, j, color):
            i_list.append(i)
            j_list.append(j)
            colors.append(color)

        if ct2 is None:
            ct_pairs = ct1.pairList()
            for i, j in ct_pairs:
                add_ij_color(i, j, (0.5, 0.5, 0.5, 0.7))
            return (i_list, j_list, colors)
        ct1 = set(ct1.pairList())
        ct2 = set(ct2.pairList())
        shared = ct1.intersection(ct2)
        ref = ct1.difference(ct2)
        comp = ct2.difference(ct1)
        sharedcolor = (150/255., 150/255., 150/255., 0.7)
        refcolor = (38/255., 202/255., 145/255., 0.7)
        compcolor = (153/255., 0.0, 1.0, 0.7)
        for i, j in comp:
            add_ij_color(i, j, compcolor)
        for i, j in ref:
            add_ij_color(i, j, refcolor)
        for i, j in shared:
            add_ij_color(i, j, sharedcolor)
        return (i_list, j_list, colors)

    def normalize_coordinates(self):
        # scale so that the median base pair is 1 unit distance
        bp_list = self.pairList()
        bp_distances = []
        for bp in bp_list:
            y_dist = self.ycoordinates[bp[0]-1] - self.ycoordinates[bp[1]-1]
            x_dist = self.xcoordinates[bp[0]-1] - self.xcoordinates[bp[1]-1]
            bp_distances.append((y_dist**2 + x_dist**2)**0.5)
        scale_factor = np.median(bp_distances)
        self.ycoordinates /= scale_factor
        self.xcoordinates /= scale_factor

        # shift so that center of plot is [0,0]
        x_center = (max(self.xcoordinates) + min(self.xcoordinates))/2
        self.xcoordinates -= x_center
        y_center = (max(self.ycoordinates) + min(self.ycoordinates))/2
        self.ycoordinates -= y_center


###############################################################################
# end of CT class
###############################################################################


class VARNA(CT):
    def __init__(self, filepath, datatype="ct", **kwargs):
        """
        if givin an input file .varna construct the ct object automatically
        """
        super().__init__(filepath=filepath, datatype=datatype, **kwargs)

    def read_file(self, **kwargs):
        """Generates CT object data from an ss file, including nucleotide x and
        y coordinates.

        Args:
            ss (file path): path to xrna, varna, nsd, or cte file
        """
        # get and check file extension, then read
        self.ss_type = "varna"
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        tree = xmlet.parse(self.filepath)
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
        # store attributes
        self.sequence = sequence.upper().replace("T", "U")
        self.num = np.arange(len(self.sequence))
        self.pair2CT(basepairs, **kwargs)
        self.xcoordinates = np.array(xcoords)
        self.ycoordinates = np.array(ycoords)
        self.normalize_coordinates()


class XRNA(CT):
    def __init__(self, datatype="ct", dataframe=None, filepath=None, **kwargs):
        super().__init__(datatype, dataframe, filepath, **kwargs)

    def read_file(self, **kwargs):
        """Generates CT object data from an ss file, including nucleotide x and
        y coordinates.

        Args:
            ss (file path): path to xrna, varna, nsd, or cte file
        """
        # get and check file extension, then read
        self.ss_type = "xrna"
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        tree = xmlet.parse(self.filepath)
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
        # store attributes
        coord_scale_factor = 1/20
        self.sequence = sequence.upper().replace("T", "U")
        self.num = np.arange(len(self.sequence))
        self.pair2CT(basepairs, **kwargs)
        self.xcoordinates = np.array(xcoords)
        self.ycoordinates = np.array(ycoords)
        self.normalize_coordinates()


class CTE(CT):
    def __init__(self, datatype="ct", dataframe=None, filepath=None, **kwargs):
        super().__init__(datatype, dataframe, filepath, **kwargs)

    def read_file(self, **kwargs):
        """Generates CT object data from an ss file, including nucleotide x and
        y coordinates.

        Args:
            ss (file path): path to xrna, varna, nsd, or cte file
        """
        # get and check file extension, then read
        self.ss_type = "cte"
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        names = ['nuc', 'seq', 'pair', 'xcoords', 'ycoords']
        ct = pd.read_csv(self.filepath, sep=r'\s+', usecols=[0, 1, 4, 8, 10],
                         names=names, header=0)
        sequence = ''.join(list(ct['seq']))
        xcoords = np.array([float(x) for x in ct.xcoords])
        ycoords = np.array([-float(y) for y in ct.ycoords])
        ct = ct[ct.nuc < ct.pair]
        basepairs = [[int(i), int(j)] for i, j in zip(ct.nuc, ct.pair)]
        # store attributes
        self.sequence = sequence.upper().replace("T", "U")
        self.num = np.arange(len(self.sequence))
        self.pair2CT(basepairs, **kwargs)
        self.xcoordinates = xcoords
        self.ycoordinates = ycoords
        self.normalize_coordinates()


class NSD(CT):
    def __init__(self, datatype="ct", dataframe=None, filepath=None, **kwargs):
        super().__init__(datatype, dataframe, filepath, **kwargs)

    def read_file(self):
        """Generates CT object data from an ss file, including nucleotide x and
        y coordinates.

        Args:
            ss (file path): path to xrna, varna, nsd, or cte file
        """
        # get and check file extension, then read
        self.ss_type = "nsd"
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], '', [], []
        with open(self.filepath, 'r') as file:
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
                    nt = {}
                    for field in line:
                        if ':' in field:
                            key, value = field.split(':')
                            nt[key] = value
                    sequence += nt["Base"]
                    xcoords.append(float(nt["X"]))
                    ycoords.append(-float(nt["Y"]))
                elif item == "pairs":
                    for field in line:
                        if field.startswith("Pair"):
                            field = field.strip('Pair:"').split(":")
                            basepairs.append([int(nuc) for nuc in field])
        # store attributes
        self.sequence = sequence.upper().replace("T", "U")
        self.num = np.arange(len(self.sequence))
        self.pair2CT(basepairs)
        self.xcoordinates = np.array(xcoords)
        self.ycoordinates = np.array(ycoords)
        self.normalize_coordinates()


class DotBracket(CT):
    def __init__(self, datatype="ct", dataframe=None, filepath=None, **kwargs):
        super().__init__(datatype, dataframe, filepath, **kwargs)

    def read_file(self, filterNC=False, filterSingle=False):
        """Generates CT object from a dot-bracket notation file. Lines
        starting with ">" or "#" are ignored. First non-ignored line is the
        sequence and the second is the structure.

        Args:
            filterNC (bool, optional): If True, will filter out non-canonical
                base pairs. Defaults to False.
            filterSingle (bool, optional): If True, will filter out any
                singleton base pairs. Defaults to False.
        """
        header, seq, bp_str, = "", "", ""
        with open(self.filepath) as f:
            for line in f:
                if line[0] in ">#":
                    header += line
                    continue
                if all(nt.upper() in "GUACT" for nt in line.strip()):
                    seq += line.strip()
                elif all(s in ".][)(}{><" for s in line.strip()):
                    bp_str += line.strip()
        assert len(bp_str) == len(seq), "db: seq and bp strings are mismatched"
        num = list(range(1, len(seq)+1))
        bp = [0 for _ in num]
        opens = {"[": [], "{": [], "(": [], "<": []}
        sym_pair = {"]": "[", "}": "{", ")": "(", ">": "<"}
        for i, sym in enumerate(bp_str):
            if sym == ".":
                continue
            elif sym in "[{(<":
                opens[sym].append(i)
            elif sym in "]})>":
                j = opens[sym_pair[sym]].pop()
                bp[i] = j + 1
                bp[j] = i + 1

        self.header = header
        self.num = num
        self.sequence = seq
        self.ct = bp
        self.mask = [0 for _ in seq]

        if filterNC:
            self.filterNC()
        if filterSingle:
            self.filterSingleton()


class JSON(CT):
    def __init__(self, datatype="ct", dataframe=None, filepath=None, **kwargs):
        super().__init__(datatype, dataframe, filepath, **kwargs)

    def read_file(self):
        """Generates CT object data from an ss file, including nucleotide x and
        y coordinates.

        Args:
            ss (file path): path to xrna, varna, nsd, or cte file
        """
        # get and check file extension, then read
        self.ss_type = "nsd"
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], '', [], []
        with open("Lys_RS-RF00168.colored.json", 'r') as file:
            file = json.load(file)
        for nt in file["rnaComplexes"][0]['rnaMolecules'][0]['sequence']:
            if nt["residueName"] not in ["5'", "3'"]:
                sequence += nt["residueName"]
                xcoords.append(nt["x"])
                ycoords.append(nt["y"])
        sequence = sequence.replace("5'", "").replace("3'", "")
        for bp in file["rnaComplexes"][0]['rnaMolecules'][0]['basePairs']:
            i = bp["residueIndex1"]
            j = bp["residueIndex2"]
            if [i, j] not in basepairs:
                basepairs.append([i, j])

        # store attributes
        self.sequence = sequence.upper().replace("T", "U")
        self.pair2CT(basepairs)
        self.xcoordinates = np.array(xcoords)
        self.ycoordinates = np.array(ycoords)
        self.normalize_coordinates()
