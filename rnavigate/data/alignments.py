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

_alignments_cache = {}
_globalms_params = {
    "match": 1,
    "mismatch": 0,
    "open": -5,
    "extend": -0.1}


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
    alignment1 = alignment1.upper().replace("T", "U")
    alignment2 = alignment2.upper().replace("T", "U")

    # Check for consistency
    if len(alignment1) != len(alignment2):
        raise ValueError("Alignment 1 and 2 must having matching length.")
    if alignment1.replace("-", "") != sequence1:
        raise ValueError("Alignment 1 does not match sequence 1")
    if alignment2.replace("-", "") != sequence2:
        raise ValueError("Alignment 2 does not match sequence 2")

    # Set up alignment dictionaries
    seq12 = {
        sequence2: {
            "seqA": alignment1,
            "seqB": alignment2}}
    seq21 = {
        sequence1: {
            "seqA": alignment2,
            "seqB": alignment1}}

    # Store dictionaries in _alignments_cache
    if sequence1 not in _alignments_cache:
        _alignments_cache[sequence1] = seq12
    else:
        _alignments_cache[sequence1].update(seq12)
    if sequence2 not in _alignments_cache:
        _alignments_cache[sequence2] = seq21
    else:
        _alignments_cache[sequence2].update(seq21)

class BaseAlignment(ABC):
    """Abstract base class for alignments

    Attributes:
        starting_sequence (str):
            the beginning sequence
        target_sequence (str):
            the sequence to map to
        mapping (numpy.array):
            a vector of size (len(starting_sequence)) which maps to an index in
            target_sequence, i.e.:
            starting_sequence[10] is aligned to target_sequence[mapping[10]]
            if mapping[10] == -1, starting_sequence[10] is unmapped (deleted)

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
    def __init__(self, starting_sequence, target_sequence):
        """Creates a BaseAlignment with starting and target sequences.

        Args:
            starting_sequence (str): the starting sequence
            target_sequence (str): the target sequence
        """
        self.starting_sequence = starting_sequence
        self.target_sequence = target_sequence
        self.mapping = self.get_mapping()

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
        new_values = np.full(len(self.target_sequence), fill)
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
        dataframe[position_column] = self.map_positions(
            dataframe[position_column])
        dataframe = dataframe[dataframe[position_column] != 0]
        new_dataframe = pd.DataFrame({
            position_column: np.arange(len(self.target_sequence))+1
        })
        new_dataframe = new_dataframe.merge(dataframe, 'left', 'Nucleotide')
        new_dataframe[sequence_column] = list(self.target_sequence)
        return dataframe.copy()

    @abstractmethod
    def get_mapping(self):
        """Alignments require a mapping from starting to target sequence"""
        pass


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
    def __init__(self, sequence1, sequence2, full=False):
        """Creates an alignment from sequence1 to sequence2.

        Args:
            sequence1 (str): the starting sequence
            sequence2 (str): the target sequence
            full (bool, optional): whether to keep unmapped starting sequence
                positions. Defaults to False.
        """
        if isinstance(sequence1, data.Sequence):
            sequence1 = sequence1.sequence
        if isinstance(sequence2, data.Sequence):
            sequence2 = sequence2.sequence
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.alignment1, self.alignment2 = self.get_alignment()
        self.full = full
        if full:
            target_sequence = self.alignment1
        else:
            target_sequence = self.sequence2
        super().__init__(sequence1, target_sequence)

    def __repr__(self):
        """a nice text only representation of an alignment"""
        return f"""alignment:
        {self.alignment1}
        {''.join(['X|'[n1==n2] for n1, n2 in zip(self.alignment1, self.alignment2)])}
        {self.alignment2}"""

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
            try:
                align1, align2 = _alignments_cache[seq1][seq2].values()
            except KeyError:
                alignment = align.globalms(seq1, seq2,
                                           penalize_end_gaps=False,
                                           **_globalms_params)
                set_alignment(
                    sequence1=seq1,
                    sequence2=seq2,
                    alignment1=alignment[0].seqA,
                    alignment2=alignment[0].seqB
                )
                align1, align2 = _alignments_cache[seq1][seq2].values()
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
            ValueError: if the target sequence of one alignment doesn't match the starting sequence of the next.
        """
        def normalize(sequence):
            return sequence.upper().replace('T', 'U')

        next_sequence = normalize(alignments[0].starting_sequence)
        for alignment in alignments:
            if next_sequence == normalize(alignment.starting_sequence):
                next_sequence = normalize(alignment.target_sequence)
            else:
                raise ValueError("Alignments do not chain.")

        self.alignments = alignments
        super().__init__(
            self.alignments[0].starting_sequence,
            self.alignments[-1].target_sequence)

    def get_mapping(self):
        """combines mappings from each alignment.

        Returns:
            numpy.array: a mapping from initial starting sequence to final
                target sequence
        """
        indices = self.alignments[0].mapping
        for alignment in self.alignments[1:]:
            indices = alignment.map_indices(indices)
        return indices

class RegionAlignment(BaseAlignment):
    """An alignment that drops values that are not within a region

    Attributes:
        starting_sequence (str): the starting sequence
        target_sequence (str): a subsequence of starting sequence
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
    def __init__(self, sequence, region):
        """Creates an alignment from a sequence to a substring of itself

        Args:
            sequence (str | rnavigate.data.Sequence): the starting sequence
            region (list of 2 int): the 1-indexed, inclusive starting and
                ending position of subsequence
        """
        if isinstance(sequence, data.Data):
            sequence = sequence.sequence
        start, end = region
        self.start = start
        self.end = end
        super().__init__(sequence, sequence[start-1: end])

    def get_mapping(self):
        """returns a mapping from sequence to subsequence

        Returns:
            numpy.array: the mapping array
        """
        mapping = np.full(len(self.starting_sequence), -1, dtype=int)
        mapping[self.start-1:self.end] = np.arange(self.end-self.start+1,
                                                   dtype=int)
        return mapping
