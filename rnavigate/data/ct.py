#!/usr/bin/env python

###############################################################################
#  SecondaryStructure code comes from CT object in RNAtools2.py and Superfold
#      Main contributors:
#           Gregg Rice
#           Anthony Mustoe
#           Tom Christy
#           Patrick Irving
###############################################################################

import sys
from os.path import isfile
import numpy as np
import xml.etree.ElementTree as xmlet
from rnavigate import data
import pandas as pd
import json


class SecondaryStructure(data.Sequence):
    """Base class for secondary structures.

    Args:
        Data (class): Parent class

    Attributes:
        input_data (str | pandas.DataFrame): structure file or dataframe
        sequence (str): sequence string
        num (numpy.array): nucleotide positions
        ct (numpy.array): paired nucleotide for each nucleotide position
        header (str): header information from CT file
        mask (list): list of 0's and 1's for each nucleotide position
        xcoordinates (numpy array): x-coordinate of each nucleotide
        ycoordinates (numpy array): y-coordinate of each nucleotide
    """

    def __init__(self, input_data, extension=None, **kwargs):
        """Creates a SecondaryStructure object from a given file or dataframe.

        Args:
            input_data (str | pandas.DataFrame):
                A dataframe or filepath containing a secondary structure
                DataFrame should contain these columns:
                    ["Nucleotide", "Sequence", "Pair"]
                "Pair" column must be redundant.
                Filepath parsing is determined by file extension:
                    varna, xrna, nsd, cte, ct, dbn, bracket, json (R2DT), forna
        """
        if isinstance(input_data, pd.DataFrame):
            self.data = input_data
            self.filepath = "dataframe"
        elif isfile(input_data):
            if extension is None:
                extension = input_data.split('.')[-1].lower()
            read_file = {
                "varna": self.read_varna,
                "xrna": self.read_xrna,
                "nsd": self.read_nsd,
                "cte": self.read_cte,
                "ct": self.read_ct,
                "dbn": self.read_dotbracket,
                "bracket": self.read_dotbracket,
                "json": self.read_r2dt,
                "forna": self.read_forna
            }[extension]
            self.filepath = input_data
            self.data = read_file(**kwargs)
        super().__init__(self.data)
        self.normalize_coordinates()

    # properties to provide backwards compatibility

    @property
    def num(self):
        return self.data['Nucleotide'].values

    @property
    def ct(self):
        return self.data["Pair"].values

    @property
    def ycoordinates(self):
        return self.data["Y_coordinate"].values
    
    @ycoordinates.setter
    def ycoordinates(self, values):
        self.data["Y_coordinate"] = values

    @property
    def xcoordinates(self):
        return self.data["X_coordinate"].values

    @xcoordinates.setter
    def xcoordinates(self, values):
        self.data["X_coordinate"] = values

    def __str__(self):
        """print the filepath and length of the RNA"""
        a = f'Name = {self.filepath}, length = {len(self.ct)}'
        return a

    # Loading Data

    def read_ct(self, structNum=0):
        """Loads secondary structure information from a given ct file.

        Requires a properly formatted header!

        Args:
            structNum (int, optional): If ct file contains multiple structures,
                uses the given structure. Defaults to 0.
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

        # check consistency!
        for i in range(len(bp)):
            if bp[i] != 0:
                p1 = (i+1, bp[i])
                p2 = (bp[bp[i]-1], bp[i])
                if p1 != p2:
                    print("WARNING: Inconsistent pair "
                          f"{p1[0]}-{p1[1]} vs. {p2[0]}-{p2[1]}")

        self.header = header
        return pd.DataFrame({
            'Nucleotide': num,
            'Sequence': list(seq),
            "Pair": bp,
            "Mask": mask})

    def read_varna(self):
        """Generates SecondaryStructure object data from a VARNA file,
        including nucleotide x and y coordinates.
        """
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
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i-1] = j
            pairs[j-1] = i
        return pd.DataFrame({
            'Nucleotide': np.arange(len(sequence))+1,
            'Sequence': list(sequence),
            'Pair': pairs,
            'X_coordinate': xcoords,
            'Y_coordinate': ycoords})

    def read_xrna(self):
        """Generates SecondaryStructure object data from an XRNA file,
        including nucleotide x and y coordinates.
        """
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
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i-1] = j
            pairs[j-1] = i
        return pd.DataFrame({
            'Nucleotide': np.arange(len(sequence))+1,
            'Sequence': list(sequence),
            'Pair': pairs,
            'X_coordinate': xcoords,
            'Y_coordinate': ycoords})

    def read_cte(self):
        """Generates SecondaryStructure object data from an ss file, including
        nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        return pd.read_table(
            self.filepath,
            sep=r'\s+',
            usecols=[0, 1, 4, 8, 10],
            names=[
                'Nucleotide',
                'Sequence',
                'Pair',
                'X_coordinate',
                'Y_coordinate'],
            header=0)

    def read_nsd(self):
        """Generates SecondaryStructure object data from an NSD file
        (format for RNAStructure StructureEditor), including nucleotide x and
        y coordinates.
        """
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
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i-1] = j
            pairs[j-1] = i
        return pd.DataFrame({
            'Nucleotide': np.arange(len(sequence))+1,
            'Sequence': list(sequence),
            'Pair': pairs,
            'X_coordinate': xcoords,
            'Y_coordinate': ycoords})

    def read_dotbracket(self):
        """Generates SecondaryStructure object data from a dot-bracket notation
        file, including nucleotide x and y coordinates.
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
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        sym_pair |= {c: c.upper() for c in alphabet}
        for i, sym in enumerate(bp_str):
            if sym == ".":
                continue
            elif sym in "[{(<" + alphabet.upper():
                opens[sym].append(i)
            elif sym in "]})>" + alphabet:
                j = opens[sym_pair[sym]].pop()
                bp[i] = j + 1
                bp[j] = i + 1
        # store attributes
        return pd.DataFrame({
            'Nucleotide': np.arange(len(seq))+1,
            'Sequence': list(seq),
            'Pair': bp})

    def read_r2dt(self):
        """Generates SecondaryStructure object data from an R2DT JSON file,
        including nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], '', [], []
        with open(self.filepath, 'r') as file:
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
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i-1] = j
            pairs[j-1] = i
        return pd.DataFrame({
            'Nucleotide': np.arange(len(sequence))+1,
            'Sequence': list(sequence),
            'Pair': pairs,
            'X_coordinate': xcoords,
            'Y_coordinate': ycoords})

    def read_forna(self):
        """Generates SecondaryStructure object data from a FORNA JSON file,
        including nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], [], [], []
        with open(self.filepath, 'r') as file:
            file = json.load(file)
        rna_key = list(file["rnas"].keys())[0]
        rna = file["rnas"][rna_key]
        rna_length = rna["rnaLength"]

        for nt in rna["nodes"][:rna_length]:
            sequence.append(nt["name"])
            xcoords.append(nt["x"])
            ycoords.append(nt["y"])
        for pair in rna["pseudoknotPairs"]:
            basepairs.append(pair)
        for i, j in enumerate(rna["pairtable"][-rna_length:]):
            if (j != 0) and ([j, i+1] not in basepairs):
                basepairs.append([i+1, j])
        # store attributes
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i-1] = j
            pairs[j-1] = i
        return pd.DataFrame({
            'Nucleotide': np.arange(len(sequence))+1,
            'Sequence': list(sequence),
            'Pair': pairs,
            'X_coordinate': xcoords,
            'Y_coordinate': ycoords})

    # writing files

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

    def get_dbn(self):
        """Returns a dotbracket notation representation of the secondary
        structure.

        Pseudoknot levels:
            1: ()
            2: []
            3: {}
            4: <>
            5: Aa
            6: Bb
            7: Cc
            etc...
        """
        dbn = ['.'] * self.length
        pair_list = self.pairList()
        pair_levels = np.zeros(len(pair_list), dtype=int)
        levels = ['()', '[]', '{}', '<>']

        # check for pseudoknots
        level = 0
        while level in pair_levels:
            print(level)
            indexes = np.where(pair_levels == level)[0]
            for i in indexes:
                for j in indexes[indexes > i]:
                    left1, right1 = pair_list[i]
                    left2, right2 = pair_list[j]
                    if pair_levels[i] != pair_levels[j]:
                        continue
                    elif (left1 < left2 < right1) and (right1 < right2):
                        pair_levels[j] += 1
            level += 1

        for pair, level in zip(pair_list, pair_levels):
            dbn[pair[0]-1] = levels[level][0]
            dbn[pair[1]-1] = levels[level][1]
        dbn = ''.join(dbn)
        return dbn

    def get_human_readable_dbn(self):
        """Returns dotbracket notation string representing SecondaryStructure
        object. This is an experimental format designed to be more human
        readable, i.e. no counting of brackets required.


        1)  Letters, instead of brackets, are used to denote nested base pairs.
        2)  Each helix is assigned a letter, which is incremented one letter
            alphabetically from the nearest enclosing stem.
        3)  Non-nested helices (pseudoknots) are assigned canonical brackets.

        From this canonical dbn string:
            how many bases are in the base stem?
            how many nested helices are there?
            ((((....(((.[[..)))))(((...(((..]].))))))))
        Same question, new format:
            AABB....CCC.[[..cccbbBBB...CCC..]].cccbbbaa

        Pseudoknot levels:
            1: Aa, Bb, Cc, etc.
            2: [], 3: {}, 4: <>
        """
        dbn = list(self.get_dbn())
        alphabet = 'abcdefghijklmnopqrstuvwxyz'
        closed = {n: 0 for n in alphabet}
        helices = self.extractHelices(
            splitHelixByBulge=False, keep_singles=True)
        for pairs in helices.values():
            left_most, right_most = pairs[0]
            if dbn[left_most-1] not in '()':
                continue
            for n in alphabet:
                if left_most > closed[n]:
                    letter = n
                    closed[n] = right_most
                    break
            for left, right in pairs:
                dbn[left-1] = letter.capitalize()
                dbn[right-1] = letter
        dbn = ''.join(dbn)
        return dbn

    def filterNC(self):
        """
        Removes non-canonical basepairs from the ct datastructure. Warning,
        this deletes information.
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
        Removes singleton basepairs from the ct datastructure. Warning, this
        deletes information.
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
        Writes a ct file from the ct object.

        Args:
            writemask (bool, optional)= True will write out any masking info
                Defaults to False.
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
        for i, nt in enumerate(self.num):
            prev = nt-1
            next = nt+1 % (self.length + 1)  # last nt goes back to zero
            seq = self.sequence[i]
            ct = self.ct[i]
            if mask and self.mask[i]:
                line = f'{nt:5d} {seq} {prev:5d} {next:5d} {ct:5d} {nt:5d} 1\n'
            else:
                line = f'{nt:5d} {seq} {prev:5d} {next:5d} {ct:5d} {nt:5d}\n'
            w.write(line)
        w.close()

    def copy(self):
        """Returns a deep copy of the ct object."""
        out = SecondaryStructure()
        out.filepath = self.filepath[:]
        out.num = self.num[:]
        out.sequence = self.sequence[:]
        out.ct = self.ct[:]
        return out

    def pairList(self):
        """
        Returns a non-redundant list of base pairs i < j as a array of tuples.
        e.g., [(19,50),(20,49)....]
        """
        out = []
        for left, right in zip(self.num, self.ct):
            if right != 0 and right > left:
                out.append((left, right))
        return out

    def pairedResidueList(self, paired=True):
        """Returns a list of residues that are paired.

        Args:
            paired (bool, optional): If False, returns unpaired residues.
                Defaults to True.
        """

        out = []
        for i, nt in enumerate(self.ct):
            if (nt != 0) == paired:
                out.append(self.num[i])

        return out

    def junctionResidues(self):
        """
        Returns a list of residues at junctions (paired, but adjacent to
        an unpaired residue or the end of a chain)
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
        """
        Add base pairs to current ct file.

        Args:
            pairs (list): 1-indexed list of paired residues.
                e.g. [(1, 20), (2, 19)]
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
        """Reconstructs this SecondaryStructure object from a list of base
        pairs and a sequence.

        Args:
            pairs (list): 1-indexed list of pairs. e.g. [(1, 20), (2, 19)]
            seq (str, optional): sequence string. e.g., "AUCGUGUCAUGCUA"
                Defaults to None.
            name (str, optional): name stored as self.filepath.
                Defaults to None.
            skipConflicting (bool, optional): whether to skip conflicting pairs
                Defaults to True.
            filterNC (bool, optional): whether to remove non-canonical pairs.
                Defaults to False.
            filterSingle (bool, optional): whether to remove singlet pairs.
                Defaults to False.
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
        """Returns a slicing object from start to end (1-indexed, inclusive)

        Args:
            start (int, optional): start position. Defaults to None.
            end (int, optional): end position. Defaults to None.

        Returns:
            slice: slicing object for the region from start to end
        """
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
        """Returns a deep copy of SecondaryStructure object, with pairs between
        start and end set to unpaired, (or the inverse).

        Args:
            start (int): start position
            end (int): end position
            nocross (bool, optional): _description_. Defaults to True.
            inverse (bool, optional): invert the behavior. Defaults to False.

        Returns:
            SecondaryStructure: new SecondaryStructure with pairs removed
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

        Args:
            start (int): start position
            end (int): end position

        Returns:
            SecondaryStructure: a new SecondaryStructure cut from self
        """

        sel = self.getNTslice(start, end)

        offset = self.num[0]
        numnts = end-start+1

        out = SecondaryStructure()
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
        """Returns the ct attribute in non-redundant form - only pairs i<j"""
        pairs = self.pairList()
        halfPlexCT = np.zeros_like(self.ct)
        for i, j in pairs:
            halfPlexCT[i-1] = j
        return halfPlexCT

    def ctToArcList(self):
        """Returns a 2D array containing arc list elements"""

        pairs = self.stripCT()

        outarr = []
        for i, v in enumerate(pairs):
            if v != 0:
                temp = np.zeros_like(pairs)
                temp[i] = v
                outarr.append(temp)
        return outarr

    def get_distance_matrix(self, recalculate=False):
        """Based on Tom's contactDistance function below, but instead returns
        the all pairs shortest paths matrix, and stores it as an attribute. If
        the attribute has already been set, it returns the attribute. This is
        faster than calling contactDistance pairwise to fill the matrix.

        Args:
            recalculate (bool, optional): Set to true to recalculate the matrix
                even if the attribute is set. In case changes to the structure
                have been made.
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
        """Contact distance is the shortest distance in a secondary structure
        (graph) between two nucleotides (nodes), where jumping across a
        base-pair or to an adjacent nucleotide counts as 1 distance
        (edge length = 1). This function utilizes a breadth first search (BFS)
        to determine the distance. BFS may take a long time on long RNAs (human
        LSU took 10 min on my laptop), but it always finds the shortest path.

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

    def extractHelices(self, fillPairs=True, splitHelixByBulge=True,
                       keep_singles=False):
        """Returns a dictionary of helices from the secondary structure.
        Keys are equivalent to list indices. Values are lists of paired
        nucleotides (1-indexed) in that helix. e.g. {0:[(1,50),(2,49),(3,48)}

        Args:
            fillPairs (bool, optional):
                Whether 1-1 and 2-2 bulges are replaced with base pairs.
                Defaults to True.
            splitHelixByBulge (bool, optional):
                Whether to split helices on bulges.
                Defaults to True.
            keep_singles (bool, optional):
                Whether to return helices that contain only 1 base-pair.
                Defaults to False.
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
                if len(tempPairs) <= 1 and not keep_singles:
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
        Returns the pk1 and pk2 pairs from the secondary structure. Ignores
        single base pairs. PK1 is defined as the helix crossing the most other
        bps. if there is a tie, the most 5' helix is called pk1

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
        out = SecondaryStructure()
        out.sequence = referenceCT.seq
        out.num = referenceCT.num

        # align target to reference
        seed = 200
        if len(self.sequence) <= seed:
            seed = len(self.sequence) - 1
        pos = 0
        maxScore = 0
        # print len(referenceCT.seq)
        for i in range(len(referenceCT.sequence)-seed):
            a, b = referenceCT.sequence[i:i+seed], self.sequence[:seed]
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
        for i in range(len(referenceCT.sequence)):
            # if the target falls within the range of the reference ct
            #     then change the numbers
            # else pad the files with 000's
            if i >= pos and i < pos+len(self.sequence):
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
        """Compute the PPV and sensitivity between self and another
        SecondaryStructure object.

        Args:
            compCT (SecondaryStructure): The SecondaryStructure to compare to.
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

    def get_aligned_data(self, alignment):
        """Returns a new SecondaryStructure matching the alignment target.

        Args:
            alignment (data.Alignment): an alignment object used to map values
        """
        df = alignment.map_nucleotide_dataframe(self.data)
        df['Pair'].fillna(0, inplace=True)
        mask = df['Pair'] != 0
        df.loc[mask, 'Pair'] = alignment.map_positions(
            df.loc[mask, "Pair"].values)
        return SecondaryStructure(input_data=df)

    def get_ij_colors(self, compct=None):
        """Gets i, j, and colors lists for plotting base pairs. i and j are the
        5' and 3' ends of each pair, and colors is the color to use for each
        pair (all grey for SecondaryStructure). If compct is provided, i and j
        are the union of pairs in self and compct, and colors represents
        whether the pair came from self, compct, or both.

        Args:
            compct (SecondaryStructure, optional):
                SecondaryStructure to compare to.
                Defaults to None.

        Returns:
            list, list, list: 5' and 3' ends of each pair, color for each pair
        """
        if compct is not None:
            alignment = data.SequenceAlignment(compct, self)
            compct = compct.get_aligned_data(alignment)
        ct1 = self
        ct2 = compct

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
        """Normalizes (x, y) coordinates of each nucleotide such that the
        structure is centered on (0, 0) and the median base-pair distance is 1.
        """
        # scale so that the median base pair is 1 unit distance
        if 'X_coordinate' not in self.data.columns:
            return
        df = self.data
        bp_list = self.pairList()
        bp_distances = []
        for bp in bp_list:
            y_dist = df.Y_coordinate[bp[0]-1] - df.Y_coordinate[bp[1]-1]
            x_dist = df.X_coordinate[bp[0]-1] - df.X_coordinate[bp[1]-1]
            bp_distances.append((y_dist**2 + x_dist**2)**0.5)
        scale_factor = np.median(bp_distances)
        df.Y_coordinate /= scale_factor
        df.X_coordinate /= scale_factor

        # shift so that center of plot is [0,0]
        x_center = (max(df.X_coordinate) + min(df.X_coordinate))/2
        df.X_coordinate -= x_center
        y_center = (max(df.Y_coordinate) + min(df.Y_coordinate))/2
        df.Y_coordinate -= y_center


###############################################################################
# end of SecondaryStructure class
###############################################################################
