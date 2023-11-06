"""Alignment objects map coordinates, vectors, and dataframes to a new sequence

Classes:
    SequenceAlignment (BaseAlignment)
        aligns one sequence another sequence
    RegionAlignment (BaseAlignment)
        cuts a sequence between a start and end position
    AlignmentChain (BaseAlignment)
        allows chaining of above alignments
"""

from abc import ABC, abstractmethod
from Bio.pairwise2 import align
import numpy as np
import pandas as pd
from rnavigate import data

# store sequence alignments
_alignments_cache = {}

# structure alignment parameters
# conversion of nt+pairing to pseudo amino acid sequence
def convert_sequence(aas, nts, dbn):
    nts_key = 'AAAAACCCCCUUUUUGGGGG-'
    dbn_key = '([.])([.])([.])([.])-'
    aas_key = 'ACDEFGHIKLMNPQRSTVWY-'
    to_aa = {nt+db: aa for nt, db, aa in zip(nts_key, dbn_key, aas_key)}
    to_nt = {aa: nt for nt, aa in zip(nts_key, aas_key)}
    to_db = {aa: db for db, aa in zip(dbn_key, aas_key)}
    if aas is True:
        # Make all nts uppercase and 'U' instead of 'T'
        nts = nts.upper().replace('T', 'U')
        # Make all pseudoknots level 1 (assumes they were assigned correctly)
        pk_level = {left: '[' for left in '{<ABCDEFGHIJKLMNOPQRSTUVWXYZ'}
        pk_level |= {right: ']' for right in '}>abcdefghijklmnopqrstuvwxyz'}
        pk_level |= {keep: keep for keep in '.()[]'}
        dbn = [pk_level[db] for db in dbn]
        aas = ''.join(to_aa[nt+db] for nt, db in zip(nts, dbn))
        return aas
    elif nts is True:
        nts = ''.join([to_nt[aa] for aa in aas])
        return nts
    elif dbn is True:
        dbn = ''.join([to_db[aa] for aa in aas])
        return dbn
    raise ValueError('Please specify which sequence to return by setting aas, '
                     'nts, or dbn to  True.')

# code used to generate structure_scoring_dict
# pair_scores = ({}
#     # both single-stranded = +6
#     | {'..': 6}
#     # pairs in the right orientation and the same nesting level = +5
#     | {pair: 5 for pair in ['((', '))', '[[', ']]']}
#     # pairs in the right orientation and different nesting level = +2
#     | {pair: 2 for pair in ['([', '])', '[(', ')]']}
#     # single stranded to double stranded = -8
#     | {bracket+'.': -8 for bracket in '([])'}
#     | {'.'+bracket: -8 for bracket in '([])'}
#     # pairs in the wrong orientation = -10
#     | {pair: -10 for pair in ['()', ')(', '[)', '(]', ')[', '](', '[]', '][']}
# )
# # matching nucleotides = +5, non-matching nucleotides = 0
# nt_scores = {n1+n2: 5  if n1==n2 else 0 for n1 in 'AUCG' for n2 in 'AUCG'}
# scoring_dict = {}
# for (nt1, pair1), aa1 in nt_dbn_to_aa_dict.items():
#     for (nt2, pair2), aa2 in nt_dbn_to_aa_dict.items():
#         nt_score = nt_scores[nt1+nt2]
#         pair_score = pair_scores[pair1+pair2]
#         scoring_dict[(aa1, aa2)] = nt_score + pair_score

structure_scoring_dict = {
    ('A', 'A'): 10, ('A', 'C'): 7, ('A', 'D'): -3, ('A', 'E'): -5, ('A', 'F'): -5, ('A', 'G'): 5, ('A', 'H'): 2, ('A', 'I'): -8, ('A', 'K'): -10, ('A', 'L'): -10, ('A', 'M'): 5, ('A', 'N'): 2, ('A', 'P'): -8, ('A', 'Q'): -10, ('A', 'R'): -10, ('A', 'S'): 5, ('A', 'T'): 2, ('A', 'V'): -8, ('A', 'W'): -10, ('A', 'Y'): -10,
    ('C', 'A'): 7, ('C', 'C'): 10, ('C', 'D'): -3, ('C', 'E'): -5, ('C', 'F'): -5, ('C', 'G'): 2, ('C', 'H'): 5, ('C', 'I'): -8, ('C', 'K'): -10, ('C', 'L'): -10, ('C', 'M'): 2, ('C', 'N'): 5, ('C', 'P'): -8, ('C', 'Q'): -10, ('C', 'R'): -10, ('C', 'S'): 2, ('C', 'T'): 5, ('C', 'V'): -8, ('C', 'W'): -10, ('C', 'Y'): -10,
    ('D', 'A'): -3, ('D', 'C'): -3, ('D', 'D'): 11, ('D', 'E'): -3, ('D', 'F'): -3, ('D', 'G'): -8, ('D', 'H'): -8, ('D', 'I'): 6, ('D', 'K'): -8, ('D', 'L'): -8, ('D', 'M'): -8, ('D', 'N'): -8, ('D', 'P'): 6, ('D', 'Q'): -8, ('D', 'R'): -8, ('D', 'S'): -8, ('D', 'T'): -8, ('D', 'V'): 6, ('D', 'W'): -8, ('D', 'Y'): -8,
    ('E', 'A'): -5, ('E', 'C'): -5, ('E', 'D'): -3, ('E', 'E'): 10, ('E', 'F'): 7, ('E', 'G'): -10, ('E', 'H'): -10, ('E', 'I'): -8, ('E', 'K'): 5, ('E', 'L'): 2, ('E', 'M'): -10, ('E', 'N'): -10, ('E', 'P'): -8, ('E', 'Q'): 5, ('E', 'R'): 2, ('E', 'S'): -10, ('E', 'T'): -10, ('E', 'V'): -8, ('E', 'W'): 5, ('E', 'Y'): 2,
    ('F', 'A'): -5, ('F', 'C'): -5, ('F', 'D'): -3, ('F', 'E'): 7, ('F', 'F'): 10, ('F', 'G'): -10, ('F', 'H'): -10, ('F', 'I'): -8, ('F', 'K'): 2, ('F', 'L'): 5, ('F', 'M'): -10, ('F', 'N'): -10, ('F', 'P'): -8, ('F', 'Q'): 2, ('F', 'R'): 5, ('F', 'S'): -10, ('F', 'T'): -10, ('F', 'V'): -8, ('F', 'W'): 2, ('F', 'Y'): 5,
    ('G', 'A'): 5, ('G', 'C'): 2, ('G', 'D'): -8, ('G', 'E'): -10, ('G', 'F'): -10, ('G', 'G'): 10, ('G', 'H'): 7, ('G', 'I'): -3, ('G', 'K'): -5, ('G', 'L'): -5, ('G', 'M'): 5, ('G', 'N'): 2, ('G', 'P'): -8, ('G', 'Q'): -10, ('G', 'R'): -10, ('G', 'S'): 5, ('G', 'T'): 2, ('G', 'V'): -8, ('G', 'W'): -10, ('G', 'Y'): -10,
    ('H', 'A'): 2, ('H', 'C'): 5, ('H', 'D'): -8, ('H', 'E'): -10, ('H', 'F'): -10, ('H', 'G'): 7, ('H', 'H'): 10, ('H', 'I'): -3, ('H', 'K'): -5, ('H', 'L'): -5, ('H', 'M'): 2, ('H', 'N'): 5, ('H', 'P'): -8, ('H', 'Q'): -10, ('H', 'R'): -10, ('H', 'S'): 2, ('H', 'T'): 5, ('H', 'V'): -8, ('H', 'W'): -10, ('H', 'Y'): -10,
    ('I', 'A'): -8, ('I', 'C'): -8, ('I', 'D'): 6, ('I', 'E'): -8, ('I', 'F'): -8, ('I', 'G'): -3, ('I', 'H'): -3, ('I', 'I'): 11, ('I', 'K'): -3, ('I', 'L'): -3, ('I', 'M'): -8, ('I', 'N'): -8, ('I', 'P'): 6, ('I', 'Q'): -8, ('I', 'R'): -8, ('I', 'S'): -8, ('I', 'T'): -8, ('I', 'V'): 6, ('I', 'W'): -8, ('I', 'Y'): -8,
    ('K', 'A'): -10, ('K', 'C'): -10, ('K', 'D'): -8, ('K', 'E'): 5, ('K', 'F'): 2, ('K', 'G'): -5, ('K', 'H'): -5, ('K', 'I'): -3, ('K', 'K'): 10, ('K', 'L'): 7, ('K', 'M'): -10, ('K', 'N'): -10, ('K', 'P'): -8, ('K', 'Q'): 5, ('K', 'R'): 2, ('K', 'S'): -10, ('K', 'T'): -10, ('K', 'V'): -8, ('K', 'W'): 5, ('K', 'Y'): 2,
    ('L', 'A'): -10, ('L', 'C'): -10, ('L', 'D'): -8, ('L', 'E'): 2, ('L', 'F'): 5, ('L', 'G'): -5, ('L', 'H'): -5, ('L', 'I'): -3, ('L', 'K'): 7, ('L', 'L'): 10, ('L', 'M'): -10, ('L', 'N'): -10, ('L', 'P'): -8, ('L', 'Q'): 2, ('L', 'R'): 5, ('L', 'S'): -10, ('L', 'T'): -10, ('L', 'V'): -8, ('L', 'W'): 2, ('L', 'Y'): 5,
    ('M', 'A'): 5, ('M', 'C'): 2, ('M', 'D'): -8, ('M', 'E'): -10, ('M', 'F'): -10, ('M', 'G'): 5, ('M', 'H'): 2, ('M', 'I'): -8, ('M', 'K'): -10, ('M', 'L'): -10, ('M', 'M'): 10, ('M', 'N'): 7, ('M', 'P'): -3, ('M', 'Q'): -5, ('M', 'R'): -5, ('M', 'S'): 5, ('M', 'T'): 2, ('M', 'V'): -8, ('M', 'W'): -10, ('M', 'Y'): -10,
    ('N', 'A'): 2, ('N', 'C'): 5, ('N', 'D'): -8, ('N', 'E'): -10, ('N', 'F'): -10, ('N', 'G'): 2, ('N', 'H'): 5, ('N', 'I'): -8, ('N', 'K'): -10, ('N', 'L'): -10, ('N', 'M'): 7, ('N', 'N'): 10, ('N', 'P'): -3, ('N', 'Q'): -5, ('N', 'R'): -5, ('N', 'S'): 2, ('N', 'T'): 5, ('N', 'V'): -8, ('N', 'W'): -10, ('N', 'Y'): -10,
    ('P', 'A'): -8, ('P', 'C'): -8, ('P', 'D'): 6, ('P', 'E'): -8, ('P', 'F'): -8, ('P', 'G'): -8, ('P', 'H'): -8, ('P', 'I'): 6, ('P', 'K'): -8, ('P', 'L'): -8, ('P', 'M'): -3, ('P', 'N'): -3, ('P', 'P'): 11, ('P', 'Q'): -3, ('P', 'R'): -3, ('P', 'S'): -8, ('P', 'T'): -8, ('P', 'V'): 6, ('P', 'W'): -8, ('P', 'Y'): -8,
    ('Q', 'A'): -10, ('Q', 'C'): -10, ('Q', 'D'): -8, ('Q', 'E'): 5, ('Q', 'F'): 2, ('Q', 'G'): -10, ('Q', 'H'): -10, ('Q', 'I'): -8, ('Q', 'K'): 5, ('Q', 'L'): 2, ('Q', 'M'): -5, ('Q', 'N'): -5, ('Q', 'P'): -3, ('Q', 'Q'): 10, ('Q', 'R'): 7, ('Q', 'S'): -10, ('Q', 'T'): -10, ('Q', 'V'): -8, ('Q', 'W'): 5, ('Q', 'Y'): 2,
    ('R', 'A'): -10, ('R', 'C'): -10, ('R', 'D'): -8, ('R', 'E'): 2, ('R', 'F'): 5, ('R', 'G'): -10, ('R', 'H'): -10, ('R', 'I'): -8, ('R', 'K'): 2, ('R', 'L'): 5, ('R', 'M'): -5, ('R', 'N'): -5, ('R', 'P'): -3, ('R', 'Q'): 7, ('R', 'R'): 10, ('R', 'S'): -10, ('R', 'T'): -10, ('R', 'V'): -8, ('R', 'W'): 2, ('R', 'Y'): 5,
    ('S', 'A'): 5, ('S', 'C'): 2, ('S', 'D'): -8, ('S', 'E'): -10, ('S', 'F'): -10, ('S', 'G'): 5, ('S', 'H'): 2, ('S', 'I'): -8, ('S', 'K'): -10, ('S', 'L'): -10, ('S', 'M'): 5, ('S', 'N'): 2, ('S', 'P'): -8, ('S', 'Q'): -10, ('S', 'R'): -10, ('S', 'S'): 10, ('S', 'T'): 7, ('S', 'V'): -3, ('S', 'W'): -5, ('S', 'Y'): -5,
    ('T', 'A'): 2, ('T', 'C'): 5, ('T', 'D'): -8, ('T', 'E'): -10, ('T', 'F'): -10, ('T', 'G'): 2, ('T', 'H'): 5, ('T', 'I'): -8, ('T', 'K'): -10, ('T', 'L'): -10, ('T', 'M'): 2, ('T', 'N'): 5, ('T', 'P'): -8, ('T', 'Q'): -10, ('T', 'R'): -10, ('T', 'S'): 7, ('T', 'T'): 10, ('T', 'V'): -3, ('T', 'W'): -5, ('T', 'Y'): -5,
    ('V', 'A'): -8, ('V', 'C'): -8, ('V', 'D'): 6, ('V', 'E'): -8, ('V', 'F'): -8, ('V', 'G'): -8, ('V', 'H'): -8, ('V', 'I'): 6, ('V', 'K'): -8, ('V', 'L'): -8, ('V', 'M'): -8, ('V', 'N'): -8, ('V', 'P'): 6, ('V', 'Q'): -8, ('V', 'R'): -8, ('V', 'S'): -3, ('V', 'T'): -3, ('V', 'V'): 11, ('V', 'W'): -3, ('V', 'Y'): -3,
    ('W', 'A'): -10, ('W', 'C'): -10, ('W', 'D'): -8, ('W', 'E'): 5, ('W', 'F'): 2, ('W', 'G'): -10, ('W', 'H'): -10, ('W', 'I'): -8, ('W', 'K'): 5, ('W', 'L'): 2, ('W', 'M'): -10, ('W', 'N'): -10, ('W', 'P'): -8, ('W', 'Q'): 5, ('W', 'R'): 2, ('W', 'S'): -5, ('W', 'T'): -5, ('W', 'V'): -3, ('W', 'W'): 10, ('W', 'Y'): 7,
    ('Y', 'A'): -10, ('Y', 'C'): -10, ('Y', 'D'): -8, ('Y', 'E'): 2, ('Y', 'F'): 5, ('Y', 'G'): -10, ('Y', 'H'): -10, ('Y', 'I'): -8, ('Y', 'K'): 2, ('Y', 'L'): 5, ('Y', 'M'): -10, ('Y', 'N'): -10, ('Y', 'P'): -8, ('Y', 'Q'): 2, ('Y', 'R'): 5, ('Y', 'S'): -5, ('Y', 'T'): -5, ('Y', 'V'): -3, ('Y', 'W'): 7, ('Y', 'Y'): 10
}

def set_alignment(sequence1, sequence2, alignment1, alignment2):
    """Add an alignment. When objects with these sequences are aligned for
    visualization, this alignment is used instead of an automated pairwise
    sequence alignment. Alignment 1 and 2 must have matching lengths.
    e.g.:
        sequence1 ="AAGCUUCGGUACAUGCAAGAUGUAC"
        sequence2 ="AUCGAUCGAGCUGCUGUGUACGUAC"
        alignment1="AAGCUUCG---------GUACAUGCAAGAUGUAC"
        alignment2="AUCGAUCGAGCUGCUGUGUAC---------GUAC"
                     |mm|   | indel |    | indel |

    Args:
        sequence1 (str): the first sequence
        sequence2 (str): the second sequence
        alignment1 (str): first sequence, plus dashes indicating indels
        alignment2 (str): second sequence, plus dashes indicating indels
    """
    # Normalize sequences
    sequence1 = sequence1.upper().replace("T", "U")
    sequence2 = sequence2.upper().replace("T", "U")
    hash1 = hash(sequence1)
    hash2 = hash(sequence2)
    alignment1 = alignment1.upper().replace("T", "U")
    alignment2 = alignment2.upper().replace("T", "U")

    # Check for consistency
    if len(alignment1) != len(alignment2):
        raise ValueError("Alignment 1 and 2 must having matching length.")
    if alignment1.replace("-", "") != sequence1:
        raise ValueError("Alignment 1 does not match sequence 1")
    if alignment2.replace("-", "") != sequence2:
        raise ValueError("Alignment 2 does not match sequence 2")
    # Store dictionaries in _alignments_cache
    if hash1 not in _alignments_cache and hash2 not in _alignments_cache:
        _alignments_cache[hash1] = {
            hash2: {
                "seqA": alignment1,
                "seqB": alignment2
                }
            }
    elif hash1 in _alignments_cache:
        _alignments_cache[hash1].update({
            hash2: {
                "seqA": alignment1,
                "seqB": alignment2
                }
            })
    elif hash2 in _alignments_cache:
        _alignments_cache[hash2].update({
            hash1: {
                "seqA": alignment2,
                "seqB": alignment1
                }
            })

def lookup_alignment(sequence1, sequence2):
    """look up a previously set alignment in the _alignments_cache
    
    Args:
        sequence1 (string)
            The first sequence to align
        sequence2 (string)
            The second sequence to be aligned to
    
    Returns:
        if 
    """
    # Normalize sequences
    sequence1 = sequence1.upper().replace("T", "U")
    sequence2 = sequence2.upper().replace("T", "U")
    # create hash keys
    hash1 = hash(sequence1)
    hash2 = hash(sequence2)
    # lookup and return alignment for sequence1 to sequence2
    try:
        return _alignments_cache[hash1][hash2]
    except KeyError:
        try:
            inverse_alignment = _alignments_cache[hash2][hash1]
            return {
                'seqA': inverse_alignment['seqB'],
                'seqB': inverse_alignment['seqA']
                }
        except KeyError:
            return None


class BaseAlignment(ABC):
    """Abstract base class for alignments

    Attributes:
        starting_sequence (str):
            the beginning sequence
        mapping (numpy.array): the alignment map array.
            index of starting_sequence is mapping[index] of target_sequence
        target_sequence (str):
            the portion of starting sequence that is mapped

    Methods:
        All map_functions map from starting sequence to target sequence.
        map_values: maps per-nucleotide values
        map_indices: maps a list of indices
        map_positions: maps a list of positions
        map_dataframe: maps a dataframe with multiple position columns
            (rows that cannot be mapped are dropped)
        map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
            (rows that cannot be mapped are dropped)
            (missing rows filled with NaN)
    """
    def __init__(self, starting_sequence, target_length):
        """Creates a BaseAlignment with starting and target sequences.

        Args:
            starting_sequence (str): the starting sequence
            target_sequence (str): the target sequence
        """
        self.starting_sequence = starting_sequence
        self.mapping = self.get_mapping()
        self.target_length = target_length
        self.target_sequence = self.get_target_sequence()

    @abstractmethod
    def get_mapping(self):
        """Alignments require a mapping from starting to target sequence"""
        return

    @abstractmethod
    def get_inverse_alignment(self):
        """Alignments require a method to get the inverted alignment"""
        return

    def get_target_sequence(self):
        """Gets the portion of starting sequence that fits the alignment"""
        return ''.join(self.map_values(list(self.starting_sequence), '-'))

    def map_values(self, values, fill=np.nan):
        """Takes an array of length equal to starting sequence and maps them to
        target sequence, unmapped positions in starting sequence are dropped
        and unmapped positions in target sequence are filled with fill value.

        Args:
            values (iterable): values to map to target sequence
            fill (any, optional): a value for unmapped positions in target.
            Defaults to numpy.nan.

        Returns:
            numpy.array: an array of values equal in length to target sequence
        """
        new_values = np.full(self.target_length, fill)
        for idx1, value in enumerate(values):
            idx2 = self.mapping[idx1]
            if idx2 == -1:
                continue
            new_values[idx2] = value
        return new_values

    def map_indices(self, indices, keep_minus_one=True):
        """Takes a list of indices (0-index) and maps them to target sequence

        Args:
            indices (int | list): a single or list of integer indices
            keep_minus_one (bool, optional): whether to keep unmapped starting
                sequence indices (-1) in the returned array. Defaults to True.

        Returns:
            numpy.array: the equivalent indices in target sequence
        """
        indices = np.array(indices, dtype=int)
        new_indices = self.mapping[indices]
        if not keep_minus_one:
            new_indices = new_indices[new_indices != -1]
        return new_indices

    def map_positions(self, positions, keep_zero=True):
        """Takes a list of positions (1-index) and maps them to target sequence

        Args:
            positions (int | list): a single or list of integer positions
            keep_zero (bool, optional): whether to keep unmapped starting
                sequence positions (0) in the returned array. Defaults to True.

        Returns:
            numpy.array: the equivalent positions in target sequence
        """
        positions = np.array(positions, dtype=int)
        new_positions = self.mapping[positions - 1] + 1
        if not keep_zero:
            new_positions = new_positions[new_positions != 0]
        return new_positions

    def map_dataframe(self, dataframe, position_columns):
        """Takes a dataframe and maps position columns to target sequence.
        Unmapped positions are dropped.

        Args:
            dataframe (pandas.DataFrame): a dataframe with position columns
            position_columns (list of str): a list of columns containing
                positions to map

        Returns:
            pandas.DataFrame: a new dataframe (copy) with position columns
                mapped or dropped
        """
        dataframe = dataframe.copy()
        for col in position_columns:
            dataframe[col] = self.map_positions(dataframe[col].values)
            dataframe = dataframe[dataframe[col] != 0]
        return dataframe.copy()

    def map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide',
                                 sequence_column='Sequence'):
        """Takes a dataframe which must have 1 row per nucleotide in starting
        sequence, with a position column and a sequence column. Dataframe is
        mapped to have the same format, but for target sequence.

        Args:
            dataframe (pandas.DataFrame): a per-nucleotide dataframe
            position_column (str, optional): name of the position column. Defaults to 'Nucleotide'.
            sequence_column (str, optional): name of the sequence column. Defaults to 'Sequence'.

        Returns:
            pandas.DataFrame: a new dataframe (copy) mapped to target sequence.
                Unmapped starting sequence positions are dropped and unmapped
                target sequence positions are filled.
        """
        dataframe = dataframe.copy()
        dataframe[position_column] = self.map_positions(dataframe[position_column])
        dataframe = dataframe[dataframe[position_column] != 0]
        new_dataframe = pd.DataFrame({
            position_column: np.arange(self.target_length)+1,
        })
        new_dataframe = new_dataframe.merge(dataframe, 'left', position_column)
        new_dataframe[sequence_column] = list(self.target_sequence)
        return new_dataframe


class SequenceAlignment(BaseAlignment):
    """The most useful feature of RNAvigate. Maps positions from one sequence
    to a totally different sequence using user-defined pairwise alignment or
    automatic pairwise alignment.

    Attributes:
        sequence1 (str): the sequence to be aligned
        sequence2 (str): the sequence to align to
        alignment1 (str): the alignment string matching sequence1 to sequence2
        alignment2 (str): the alignment string matching sequence2 to sequence1
        starting_sequence (str): sequence1
        target_sequence(str): sequence2 if full is False, else alignment2
        mapping (numpy.array): the alignment map array.
            index of starting_sequence is mapping[index] of target_sequence

    Methods:
        All map_functions map from starting sequence to target sequence.
        map_values: maps per-nucleotide values
        map_indices: maps a list of indices
        map_positions: maps a list of positions
        map_dataframe: maps a dataframe with multiple position columns
            (rows that cannot be mapped are dropped)
        map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
            (rows that cannot be mapped are dropped)
            (missing rows filled with NaN)
    """
    def __init__(
            self, sequence1, sequence2, align_kwargs=None, full=False,
            use_previous=True
            ):
        """Creates an alignment from sequence1 to sequence2.

        Args:
            sequence1 (str): the starting sequence
            sequence2 (str): the target sequence
            full (bool, optional): whether to keep unmapped starting sequence
                positions. Defaults to False.
        """
        if align_kwargs is None:
            self.align_kwargs = {
                "match": 1,
                "mismatch": 0,
                "open": -5,
                "extend": -0.1}
        if isinstance(sequence1, data.Sequence):
            sequence1 = sequence1.sequence
        if isinstance(sequence2, data.Sequence):
            sequence2 = sequence2.sequence
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.alignment1, self.alignment2 = self.get_alignment()
        self.full = full
        self.use_previous = use_previous
        if self.full:
            target_length = len(self.alignment2)
        else:
            target_length = len(self.sequence2)
        super().__init__(sequence1, target_length)

    def __repr__(self):
        """a nice text only representation of an alignment"""
        return f"""alignment:
        {self.alignment1}
        {''.join(['X '[n1==n2] for n1, n2 in zip(self.alignment1, self.alignment2)])}
        {self.alignment2}"""

    def get_inverse_alignment(self):
        return SequenceAlignment(self.sequence2, self.sequence1, self.full)

    def get_alignment(self):
        """Gets an alignment that has either been user-defined or previously
        calculated or produces a new pairwise alignment between two sequences.

        Returns:
            (tuple of 2 str): alignment1 and alignment2
        """
        # Normalize sequences
        seq1 = self.sequence1.upper().replace("T", "U")
        seq2 = self.sequence2.upper().replace("T", "U")
        # Check if sequences match
        if seq1 == seq2:
            return (seq1, seq2)
        else:
            # look in _alignments_cache, if not found, do a pairwise alignment
            alignments = lookup_alignment(seq1, seq2)
            if alignments is None:
                alignment = align.globalms(
                    seq1, seq2, penalize_end_gaps=False,
                    one_alignment_only=True,**self.align_kwargs
                    )
                set_alignment(
                    sequence1=seq1,
                    sequence2=seq2,
                    alignment1=alignment[0].seqA,
                    alignment2=alignment[0].seqB
                )
                alignments = lookup_alignment(seq1, seq2)
            align1, align2 = alignments['seqA'], alignments['seqB']
            return (align1, align2)

    def get_mapping(self):
        """Calculates a mapping from starting sequence to target sequence.

        Returns:
            numpy.array: an array of length of starting sequence that maps to
                an index of target sequence. Stored as self.mapping
                starting_sequence[idx] == target_sequence[self.mapping[idx]]
        """
        align1 = self.alignment1
        align2 = self.alignment2
        # get an index mapping from sequence 1 to the full alignment
        seq1_to_align = np.where([nt != '-' for nt in align1])[0]
        # if we want a mapping to a position in the full alignment, this is it.
        if self.full:
            return seq1_to_align
        # extra steps to get to sequence 2 positions
        # positions that are removed when plotting on sequence 2
        align_mask = np.array([nt != '-' for nt in align2])
        # an index mapping from the full alignment to position in sequence 2
        align_to_seq2 = np.full(len(align2), -1)
        align_to_seq2[align_mask] = np.arange(len(self.sequence2))
        # an index mapping from sequence 1 to sequence 2 positions
        seq1_to_seq2 = align_to_seq2[seq1_to_align]
        return seq1_to_seq2


class AlignmentChain(BaseAlignment):
    """Combines a list of alignments into one.

    Attributes:
        alignments (list): the constituent alignments
        starting_sequence (str): starting sequence of alignments[0]
        target_sequence (str): target sequence of alignments[-1]
        mapping (numpy.array): a vector that maps from starting to target
            index of starting_sequence is mapping[index] of target sequence

    Methods:
        All map_functions map from starting sequence to target sequence.
        map_values: maps per-nucleotide values
        map_indices: maps a list of indices
        map_positions: maps a list of positions
        map_dataframe: maps a dataframe with multiple position columns
            (rows that cannot be mapped are dropped)
        map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
            (rows that cannot be mapped are dropped)
            (missing rows filled with NaN)
    """
    def __init__(self, *alignments):
        """Creates a single alignment from multiple alignments

        Raises:
            ValueError: if the target sequence of one alignment doesn't match
                the starting sequence of the next.
        """
        next_sequence_len = len(alignments[0].starting_sequence)
        for alignment in alignments:
            if next_sequence_len == len(alignment.starting_sequence):
                next_sequence_len = len(alignment.target_sequence)
            else:
                raise ValueError("Alignments do not chain.")

        self.alignments = alignments
        super().__init__(self.alignments[0].starting_sequence,
                         len(self.alignments[-1].target_sequence))

    def get_mapping(self):
        """combines mappings from each alignment.

        Returns:
            numpy.array: a mapping from initial starting sequence to final
                target sequence
        """
        indices = self.alignments[0].mapping
        for alignment in self.alignments[1:]:
            valid = indices != -1
            indices[valid] = alignment.map_indices(indices[valid])
        return indices

    def get_inverse_alignment(self):
        alignments = [a.get_inverse_alignment() for a in self.alignments[::-1]]
        return AlignmentChain(*alignments)

class StructureAlignment(BaseAlignment):
    """Experimental secondary structure alignment based on RNAlign2D algorithm
    (https://doi.org/10.1186/s12859-021-04426-8)

    Attributes:
        sequence1 (str): the sequence to be aligned
        sequence2 (str): the sequence to align to
        structure1
        structure2
        aa_sequence1
        aa_sequence2
        alignment1 (str): the alignment string matching sequence1 to sequence2
        alignment2 (str): the alignment string matching sequence2 to sequence1
        starting_sequence (str): sequence1
        target_sequence(str): sequence2 if full is False, else alignment2
        mapping (numpy.array): the alignment map array.
            index of starting_sequence is mapping[index] of target_sequence

    Methods:
        All map_functions map from starting sequence to target sequence.
        map_values: maps per-nucleotide values
        map_indices: maps a list of indices
        map_positions: maps a list of positions
        map_dataframe: maps a dataframe with multiple position columns
            (rows that cannot be mapped are dropped)
        map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
            (rows that cannot be mapped are dropped)
            (missing rows filled with NaN)
    """
    def __init__(self, sequence1, sequence2, structure1=None, structure2=None,
                 full=False):
        """Creates an alignment from structure1 to structure2.

        Args:
            structure1 (str): the starting structure
            structure2 (str): the target structure
            full (bool, optional): whether to keep unmapped starting sequence
                positions. Defaults to False.
        """
        if structure1 is None:
            structure1 = sequence1
        if structure2 is None:
            structure2 = sequence2
        if isinstance(sequence1, data.Sequence):
            sequence1 = sequence1.sequence
        if isinstance(structure1, data.SecondaryStructure):
            structure1 = structure1.get_dotbracket()
        self.sequence1 = sequence1
        self.structure1 = structure1
        if isinstance(sequence2, data.Sequence):
            sequence2 = sequence2.sequence
        if isinstance(structure2, data.SecondaryStructure):
            structure2 = structure2.get_dotbracket()
        self.sequence2 = sequence2
        self.structure2 = structure2
        self.alignment1, self.alignment2 = self.get_alignment()
        self.full = full
        if self.full:
            target_length = len(self.alignment2)
        else:
            target_length = len(self.sequence2)
        super().__init__(self.sequence1, target_length)

    def get_alignment(self):
        """Aligns pseudo-amino-acid sequences according to RNAlign2D rules.

        Returns:
            (tuple of 2 str): alignment1 and alignment2
        """
        # Normalize sequences
        seq1 = convert_sequence(
            aas=True, nts=self.sequence1, dbn=self.structure1
            )
        seq2 = convert_sequence(
            aas=True, nts=self.sequence2, dbn=self.structure2
            )
        # Check if sequences match
        if seq1 == seq2:
            return (seq1, seq2)
        else:
            alignments = lookup_alignment(seq1, seq2)
            if alignments is None:
                alignment = align.globalds(
                    seq1, seq2, structure_scoring_dict, -12, -1,
                    penalize_end_gaps=False, one_alignment_only=True)
                alignment1 = alignment[0].seqA
                alignment2 = alignment[0].seqB
                set_alignment(
                    sequence1=seq1,
                    sequence2=seq2,
                    alignment1=alignment[0].seqA,
                    alignment2=alignment[0].seqB
                )
            alignment1, alignment2 = alignments['seqA'], alignments['seqB']
        # convert pseudo-amino acid alignments back into nucleotide alignments
        alignment1 = convert_sequence(aas=alignment1, nts=True, dbn=False)
        alignment2 = convert_sequence(aas=alignment2, nts=True, dbn=False)
        return alignment1, alignment2

    def get_mapping(self):
        """Calculates a mapping from starting sequence to target sequence.

        Returns:
            numpy.array: an array of length of starting sequence that maps to
                an index of target sequence. Stored as self.mapping
                starting_sequence[idx] == target_sequence[self.mapping[idx]]
        """
        align1 = self.alignment1
        align2 = self.alignment2
        # get an index mapping from sequence 1 to the full alignment
        seq1_to_align = np.where([nt != '-' for nt in align1])[0]
        # if we want a mapping to a position in the full alignment, this is it.
        if self.full:
            return seq1_to_align
        # extra steps to get to sequence 2 positions
        # positions that are removed when plotting on sequence 2
        align_mask = np.array([nt != '-' for nt in align2])
        # an index mapping from the full alignment to position in sequence 2
        align_to_seq2 = np.full(len(align2), -1)
        align_to_seq2[align_mask] = np.arange(len(self.sequence2))
        # an index mapping from sequence 1 to sequence 2 positions
        seq1_to_seq2 = align_to_seq2[seq1_to_align]
        return seq1_to_seq2

    def get_inverse_alignment(self):
        return StructureAlignment(
            self.sequence2, self.sequence1, self.structure2, self.structure1
            )

    def set_as_default_alignment(self):
        set_alignment(
            sequence1=self.sequence1, sequence2=self.sequence2,
            alignment1=self.alignment1, alignment2=self.alignment2
            )
