#!/usr/bin/env python
###############################################################################
# GPL statement:
# This file is part of SuperFold.
#
# SuperFold is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SuperFold is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SuperFold.  If not, see <http://www.gnu.org/licenses/>.
#
# 17 Nov 2014
# Copywrite 2014
# Greggory M Rice
# all rights reserved
#
# V2 version of RNAtools modified by Anthony Mustoe
# 2016
###############################################################################
########################################################
#                                                      #
#      RNAtools2.py CT and dotplot data structures     #
#               v 2.0       2016                       #
#                                                      #
#          author: Gregg Rice                          #
#                  Anthony Mustoe                      #
########################################################

# CHANGE-LOG
#   Feb 2020 Added separate filterNC, filter single functionality
#            -CT masking functionality,
#            -bp consistency checking
#            -writeRNAstructureSeq
#            -addPairs
#
#   12/6/18 Full merging of Tony's RNAtools
#
#   v0.8.1  contact distance fixed
#    Tom Christy: Contact Distance Fixed Again
#                 now uses a breadth first search algorithm
#
#   v0.8    numpy dotplot functions added
#
#   Added get_distance_matrix based on contactDistance, but faster for all pairs
#   made compatible with python 3
#   formatting changes in keeping with JNBTools project

import sys
import numpy as np


class CT(object):

    def __init__(self, fIN=None, **kwargs):
        """
        if givin an input file .ct construct the ct object automatically
        """
        if fIN:
            #self.name = fIN
            #self.num,self.seq,self.ct = self.readCT(fIN)
            self.readCT(fIN, **kwargs)

    def __str__(self):
        """overide the default print statement for the object"""
        a = '{ Name= %s, len(CT)= %s }' % (self.name, str(len(self.ct)))
        return a

    def readFasta(self, fastapath):
        """Assign sequence from fastafile"""

        with open(fastapath) as inp:

            seq = ''
            inp.readline()  # pop off the header
            for line in inp:
                if line[0] == '>':
                    print(
                        ("WARNING: Multiple sequences in {}. Using the first sequence.".format(fastapath)))
                    break
                seq += line.strip()

        if 'T' in seq or 't' in seq:
            print(("WARNING: replacing 'T' with 'U' from {}".format(fastapath)))
            seq = seq.replace('T', 'U')
            seq = seq.replace('t', 'u')

        self.seq = list(seq)

    def readCT(self, fIN, structNum=0, filterNC=False, filterSingle=False):
        """
        reads a ct file, !requires header!
        structNum allows selection of a given SS if multiple are stored in the same file
        filterNC will check the bps and filter out any between NC
        filterSingle will filter out any singleton bps
        """
        num, seq, bp, mask = [], [], [], []

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
                        seq.append(str(spl[1]))
                        bp.append(int(spl[4]))

                        # check if there is masking info
                        if len(spl) == 7 and spl[6] == '1':
                            mask.append(1)
                        else:
                            mask.append(0)

        except Exception as e:
            raise IOError(
                "{0} has invalid format or does not exist".format(fIN))

        if len(num) == 0:
            sys.exit("Structure %d was not found in the ct file" % structNum)

        if 'T' in seq:
            print("Note: T nucleotides have been recoded as U")
            seq = ['U' if x == 'T' else x for x in seq]

        # check consistency!
        for i in range(len(bp)):
            if bp[i] != 0:
                p1 = (i+1, bp[i])
                p2 = (bp[bp[i]-1], bp[i])
                if p1 != p2:
                    print(
                        ("WARNING: Inconsistent pair {0[0]}-{0[1]} vs. {1[0]}-{1[1]}".format(p1, p2)))

        self.name = fIN
        self.header = header
        self.num = num
        self.seq = seq
        self.ct = bp
        self.mask = mask

        if filterNC:
            self.filterNC()

        if filterSingle:
            self.filterSingleton()

    def filterNC(self):
        """
        Filter out non-canonical basepairs from the ct datastructure
        """

        for i, nt in enumerate(self.ct):

            if nt == 0:
                continue

            pi = self.num.index(nt)

            bp = self.seq[pi] + self.seq[i]
            if bp not in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG', '  '):
                self.ct[i] = 0
                self.ct[pi] = 0
                # print 'Deleted %d %s' % (self.num[i], self.seq[i])
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
                # print 'Deleted %d %s *' % (self.num[i], self.seq[i])

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
        w.write('{0:6d} {1}\n'.format(len(self.num), self.name))
        for i in range(len(self.num)-1):

            if mask and self.mask[i]:
                line = '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d} 1\n'.format(
                    self.num[i], self.seq[i], self.num[i]-1, self.num[i]+1, self.ct[i])
            else:
                line = '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(
                    self.num[i], self.seq[i], self.num[i]-1, self.num[i]+1, self.ct[i])
            w.write(line)

        # last line is different
        i = len(self.num)-1

        if mask and self.mask[i]:
            line = '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d} 1\n'.format(
                self.num[i], self.seq[i], self.num[i]-1, 0, self.ct[i])
        else:
            line = '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(
                self.num[i], self.seq[i], self.num[i]-1, 0, self.ct[i])

        w.write(line)
        w.close()

    def copy(self):
        """
        returns a deep copy of the ct object
        """
        out = CT()
        out.name = self.name[:]
        out.num = self.num[:]
        out.seq = self.seq[:]
        out.ct = self.ct[:]
        return out

    def pairList(self):
        """
        #returns a list of base pairs i<j as a array of tuples:
        [(19,50),(20,49)....]
        """
        out = []
        for nt in range(len(self.ct)):
            if self.ct[nt] != 0 and self.ct[nt] > self.num[nt]:
                out.append((self.num[nt], self.ct[nt]))
        return out

    def pairedResidueList(self, paired=True):
        """ return a list of residues that are paired (or single-stranded, if paired=False)"""

        out = []
        for i, nt in enumerate(self.ct):
            if (nt != 0) == paired:
                out.append(self.num[i])

        return out

    def junctionResidues(self, incGU=False):
        """ return a list of residues at junctions. Write in functionality to include opiton of GU pairs"""

        jun = []
        for i, nt in enumerate(self.ct):

            if nt == 0:  # if its not paired
                continue

            pi = self.num.index(nt)

            try:
                if self.ct[i-1] == 0 or self.ct[i+1] == 0 or \
                        self.ct[pi-1] == 0 or self.ct[pi+1] == 0:
                    jun.append(self.num[i])

            except IndexError:
                # this means we're at the beginning or end of the chain, hence must be junction
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

        # See if sequence has been defined either as self.seq or as seq keyword
        if seq is None:
            assert hasattr(self, "seq"), "Sequence is not defined"
            assert self.seq is not None, "Sequence is not defined"
        else:
            self.seq = seq

        length = len(self.seq)
        self.num = list(range(1, length+1))

        # give it a name if it has one
        if name:
            self.name = name
        else:
            self.name = 'RNA_'+str(length)

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
        if inverse = True, do the inverse -- ie return a CT with base pairs only
        involving the specified region"""

        # in case numbering differs from indexing, go by indexes
        out = self.copy()
        out.name += '_mask_'+str(start)+'_'+str(end)

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
        """
        returns a new ct file containing only base pairs within the specified region
        """

        sel = self.getNTslice(start, end)

        offset = self.num[0]
        numnts = end-start+1

        out = CT()
        out.seq = self.seq[sel]
        out.num = list(range(1, numnts+1))
        out.name = self.name + '_cut_'+str(start)+'_'+str(end)

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

    def get_distance_matrix(self, recalculate="False"):
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
            return distance_matrix
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
        return self.distance_matrix

    def contactDistance(self, i, j, matrix="False"):
        """
        Caclulates the contact distance between pairs i,j in
        the RNA using the RNAtools CT Object. This is different
        than the old contact distance function because it utilizes a breadth
        first search to determine the distance. While BFS may take longer
        than contact distance, it will always find the shortest distance.

        Added by Tom Christy
        """

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
        # each index matches to the nt number, so index 0 is just an empty place holder

        # create the queue and add nt i, set its level as 0
        queue = list()
        queue.append(i)
        level[i] = 0
        # while the queue has members, keep popping and searching until the NT is found
        notFound = True
        contactDistance = 0
        while len(queue) > 0 and notFound:
            currentNT = queue.pop(0)
            # if the current nucleotide is the one we're searching for, record it's level as the contact distance
            # and break out of the search
            if(currentNT == j):
                contactDistance = level[currentNT]
                notFound = False
                break

            # grab all the neighbors of the current NT and determine if they are inbounds and unvisited
            # If so, add them to the queue and set their levels.
            NTBack = currentNT - 1
            NTUp = currentNT + 1
            # get the pair of NT -1 b/c of the weird way coded the ct file, it's off by 1
            NTPair = self.ct[currentNT-1]

            # check all the neighbors to see if they should be added to the search, and if so, record their level(ie their Contact Distance)
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

    def GreggContactDistance(self, i, j):
        """
        calculates the contact distance between pairs i,j in
        in the RNA using the RNAtools CT object. Method for
        calculating the contact is described in Hajdin et al
        (2013).


        THIS FUNCTION IS BROKEN. It doesn't handle multi helix 
        junctions correctly and will often return a longer distance
        than the real answer.
        """

        print("WARNING: GreggContactDistance is broken and should not be used!!!")

        # correct for indexing offset
        i, j = i-1, j-1

        # error out if nucleotide out of range
        if max(i, j) > len(self.ct):
            print('Error!, nucleotide {0} out of range!'.format(max(i, j)+1))
            return

        # i must always be less than j, correct for this
        if j < i:
            i, j = j, i

        count = 0.0
        tempBack = 10000.0
        k = int(i)

        def backTrace(rna, j, k):
            bcount = 2
            k -= 1
            if k-1 == j:
                return bcount
            bcount += 1
            while k > j:
                if rna.ct[k] == 0:
                    k -= 1
                    bcount += 1
                # if not single stranded, exit the backtrace
                else:
                    return None
            return bcount
        # search forward through sequence
        while k < j:
            # debuging stuff, prevent infinite loop
            # if count > 200:
            #    print i,j
            #    break

            # nonpaired nucleotides are single beads
            if self.ct[k] == 0:
                k += 1
                #count += 6.5
                count += 1

            # branches through helices can't be skipped
            # treated as single beads
            elif self.ct[k] > j + 1:
                # try backtracing a few if it is close (within 10)
                if self.ct[k] - 10 < j:
                    back = backTrace(self, j, self.ct[k])
                    # if the backtracing is able to reach
                    # your nt, add its length to the count
                    # and break
                    if back:
                        # print 'backitup'
                        # store the backtracing count for later
                        # if it ends up being lower than the final
                        # we will use it instead
                        if count + back < tempBack:
                            tempBack = count + back
                k += 1
                #count += 6.5
                count += 1

            # simple stepping
            elif self.ct[k] < i + 1:
                k += 1
                #count += 6.5
                count += 1

            elif self.ct[k] < k+1:
                k += 1
                count += 1
                # print "Backtracking prevented, going forward to "+str(k+1)

            # handle branching, jumping the base of a helix is only 1
            else:
                # one count for the step to the next
                count += 1
                k = self.ct[k] - 1

                # if by jumping you land on your nt
                # stop counting
                if k - 1 == j:
                    break

                # one count for jumping the helix
                #count += 15.0
                #count += 1
            # print i,k,j
        finalCount = min(count, tempBack)
        return finalCount

    def extractHelices(self, fillPairs=True, splitHelixByBulge=True):
        """
        returns a list of helices in a given CT file and the nucleotides
        making up the list as a dict object. e.g. {0:[(1,50),(2,49),(3,48)}

        defaults to filling 1,1 and 2,2 mismatches

        Since April 20, 2018 an option has been added to no longer split helices by single
        nucleotide bulges. Prior to this option bulges on the 5' of the helix split the helix
        but those on the 3' end did not.
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

                # this is the default behavior which splits helices that are separated by a single nucleotide bulge
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
                        if rna.ct[nt] != 0 and not(abs(rna.ct[nt] - previous) > 2):
                            tempPairs.append((nt+1, rna.ct[nt]))
                            previous = rna.ct[nt]
                            nt += 1
                        # This if statement handles bulges on the 5' end of the helix by reading through them
                        # The prior if statement handles them on the 3' end of the helix
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
            if max(h1[0]) > min(h2[0]) > min(h1[0]) and max(h2[0]) > max(h1[0]):
                return True
            if max(h2[0]) > min(h1[0]) > min(h2[0]) and max(h1[0]) > max(h2[0]):
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
        """
        utilizes the global padCT method on this object
        """
        return padCT(self, referenceCT, giveAlignment)

    def readSHAPE(self, fIN):
        """
        utilizes the global readSHAPE method and appends a the data to the object as .shape
        """
        self.shape, seq = readSHAPE(fIN)
        if len(self.shape) < len(self.ct):
            print("warning! shape array is smaller than the CT range")

    def writeSHAPE(self, fOUT):
        """
        utilizes the global writeSHAPE method, and writes the .shape array attached to the object to a file
        """
        try:
            writeSHAPE(self.shape, fOUT)
        except:
            print("No SHAPE data present")
            return

    def computePPVSens(self, compCT, exact=True, mask=False):
        """
        compute the ppv and sensitivity between the current CT and the passed CT
        exact = True will require BPs to be exactly correct. 
                False allows +/-1 bp slippage (RNAstructure convention)
        mask = True will exclude masked regions from calculation
        """

        # check mask is properly set up if using
        if mask:
            try:
                if len(self.mask) != len(self.ct):
                    raise AttributeError()
            except:
                raise AttributeError('Mask is not properly initialized')

        if len(self.ct) != len(compCT.ct):
            raise IndexError('CT objects are different sizes')

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
                elif not exact and (v == compCT.ct[i]-1 or v == compCT.ct[i]+1):
                    sharedpairs += 1
                elif not exact and (v == compCT.ct[i-1] or v == compCT.ct[i+1]):
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
            for nt in self.seq:
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
            out.write(';\nSequence from {0}\n'.format(self.name))
            out.write('{0}1'.format(''.join(self.seq)))


#################################################################################
# end of CT class
#################################################################################


def padCT(targetCT, referenceCT, giveAlignment=False):
    """Aligns the target CT to the reference CT and pads the referece
    CT file with 000s in order to input into CircleCompare"""
    out = CT()
    out.seq = referenceCT.seq
    out.num = referenceCT.num

    # align target to reference
    seed = 200
    if len(targetCT.seq) <= seed:
        seed = len(targetCT.seq) - 1
    pos = 0
    maxScore = 0
    # print len(referenceCT.seq)
    for i in range(len(referenceCT.seq)-seed):
        a, b = referenceCT.seq[i:i+seed], targetCT.seq[:seed]
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
        if i >= pos and i < pos+len(targetCT.seq):
            val = targetCT.ct[i-pos]
            if val > 0:
                val += pos
            ct.append(val)
        else:
            ct.append(0)

    # set to the new ct file and return it
    out.ct = ct
    out.name = targetCT.name+'_renum_'+str(pos)
    if giveAlignment:
        return out, pos
    else:
        return out


def readSHAPE(fIN):
    """
    reads an RNA structure .shape or .map file. Returns an array of the SHAPE data
    """
    shape = []
    seq = ''

    with open(fIN, "rU") as inp:

        for line in inp:
            spl = line.split()
            shape.append(float(spl[1]))

            if len(spl) == 4:
                seq += spl[3][0]

    if len(seq) != len(shape):
        seq = ''

    return shape, seq


def writeSHAPE(shape, fOUT):
    """
    writes the data from a shape array into the file fOUT
    """
    w = open(fOUT, "w")

    for i in range(len(shape)):
        line = "{0}\t{1}\n".format(i+1, shape[i])
        w.write(line)
    w.close()


def readSeq(fIN, type='RNAstructure'):
    """
    reads an RNAstructure sequence file format and converts it to an
    array of nucleotides. e.g.: ['A','G','C','C'...]

    also returns the name of the sequence from the file
    """

    # strip the input file of comments
    seqRaw = []
    for i in open(fIN, "rU").read().split():
        if len(i) == 0:
            continue
        if i[0] == ";":
            continue
        seqRaw.append(i)

    name = seqRaw[0]
    seqJoin = ''.join(seqRaw[1:])
    seq = []
    for i in seqJoin:
        if i == '1':
            break
        seq.append(i)
    return processed, name


#################################################################################
# Start dotp class
#################################################################################


class DotPlot:
    """ Modified by Tony Mustoe to be more efficient at caclulations by rearranging data structure
    This increases the memory footprint, but greatly reduces looping
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

        except ValueError as e:
            pass

        return 0

    def readDP(self, fIN):

        fhandle = open(fIN)

        # read the expected two lines of headers
        length = int(fhandle.readline().split()[0])

        if self.length is not None and length != self.length:
            sys.exit('Dotp file is not the same length as existing dotp object')

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

    def returnCT(self, probcut=0.5, skipConflicting=True, filterNC=True, filterSingle=True):
        """
        return a CT object constructed from base pairs above a given pair. prob
        """

        newDP = self.requireProb(probcut)

        newCT = CT()
        newCT.pair2CT(newDP.pairList(), ['n']*self.length,
                      skipConflicting=skipConflicting, filterNC=filterNC, filterSingle=filterSingle)

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
        """
        return a new DP object where basepairs below a min prob have been trimed      
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
                partfun.append(
                    {'log10': np.array([]), 'pair': np.array([]), 'nt': nt['nt']})

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
        """This is a hacky (i.e. very unrigorous) but quick way to estimate unfolding
        free energy for a specified region of nts
        Region is 1-based, inclusive; i.e. nt numbering
        """

        if self.sequence is None:
            AttributeError("Sequence is not defined")

        bpenergy = {'AT': 1, 'TA': 1, 'GC': 2, 'CG': 2, 'TG': 1, 'GT': 1}

        unfold = 0

        # shift each by -1 since want to convert to 0-based index, and then +1 to region[1]
        for i in range(region[0]-1, region[1]):

            nt = self.partfun[i]
            seqi = self.sequence[i]

            for pindex, j in enumerate(nt['pair']):

                # make sure that each pair is only counted once. If k < i, but outside the
                # region, we don't have to worry about double counting
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

        # do sum the niave way since we aren't worried about small prob diffs with this calculation
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
        """
        replaces PlusandMinus script. If a helix in a predicted structure is slipped +/-1 nt
        we need to sum the predicted probabilities otherwise the predicted Shannon entropy
        will be artificially high.

        turning off predicted only will go through all i<j basepairs and merge them in preference
        of liklihood. This is more compuationally intensive
        """

        dotPlot = self.dp
        #dotPlotCopy = {'logBP':copy.deepcopy(dotPlot['logBP'])}

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
            # go through all i<j basepair combinations and check to see if there is a slipped base pair
            for i, j in self.pairList():

                # correct for python counting starting at 0
                pair_i, pair_j = int(i+0), int(j+0)

                # see if there exists a basepair for this combination
                filter_i = (dotPlot['i'] == pair_i)
                filter_j = (dotPlot['j'] == pair_j)

                filter_union = (filter_i * filter_j)

                # only continue if the pair is reasonably likely (>1% chance of forming)
                if np.sum(filter_union) == 1 and dotPlot['logBP'][filter_union] < 3:
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
                            # check to see if it's less likely than current pairs
                            if dotPlot['logBP'][filter_union] < dotPlot['logBP'][curr]:
                                # and add it to the current pair if it is
                                prob_at += 10**(-dotPlot['logBP'][curr])

                                # set to a very low probabliity afterwards
                                dotPlot['logBP'][curr] = 50
                    dotPlot['logBP'][filter_union] = -np.log10(prob_at)

        return slipped

    def pairingProb(self):
        pass


###############################################################################

def readShannonFile(fname):

    shan = []
    with open(fname) as inp:
        for line in inp:
            spl = line.split()
            shan.append(float(spl[1]))

    return np.array(shan)


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
