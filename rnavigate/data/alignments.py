from Bio.pairwise2 import align
import numpy as np
from rnavigate import data
from abc import ABC, abstractmethod
import pandas as pd

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

class Alignment(ABC):
    def __init__(self, target, mapping):
        self.target = target
        self.mapping = mapping

    def map_values(self, values, fill=np.nan):
        new_values = np.full(len(self.target), fill)
        for idx1, value in enumerate(values):
            idx2 = self.mapping[idx1]
            if idx2 == -1:
                continue
            new_values[idx2] = value
        return new_values

    def map_indices(self, indices, keep_minus_one=True):
        indices = np.array(indices, dtype=int)
        new_indices = self.mapping[indices]
        if not keep_minus_one:
            new_indices = new_indices[new_indices != -1]
        return new_indices

    def map_positions(self, positions, keep_zero=True):
        positions = np.array(positions, dtype=int)
        new_positions = self.mapping[positions - 1] + 1
        if not keep_zero:
            new_positions = new_positions[new_positions != 0]
        return new_positions

    def map_dataframe(self, dataframe, position_columns):
        dataframe = dataframe.copy()
        for col in position_columns:
            dataframe[col] = self.map_positions(dataframe[col].values)
            dataframe = dataframe[dataframe[col] != 0]
        return dataframe.copy()

    def map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide',
                                 sequence_column='Sequence'):
        dataframe = dataframe.copy()
        dataframe[position_column] = self.map_positions(
            dataframe[position_column])
        dataframe = dataframe[dataframe[position_column] != 0]
        new_dataframe = pd.DataFrame({
            position_column: np.arange(len(self.target))+1
        })
        new_dataframe = new_dataframe.merge(dataframe, 'left', 'Nucleotide')
        new_dataframe[sequence_column] = list(self.target)
        return dataframe.copy()

    @abstractmethod
    def get_mapping(self):
        pass


class SequenceAlignment(Alignment):
    def __init__(self, sequence1, sequence2, full=False):
        if isinstance(sequence1, data.Sequence):
            sequence1 = sequence1.sequence
        if isinstance(sequence2, data.Sequence):
            sequence2 = sequence2.sequence
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.alignment1, self.alignment2 = self.get_alignment()
        self.full = full
        if full:
            target = self.alignment1
        else:
            target = self.sequence2
        mapping = self.get_mapping()
        super().__init__(target=target, mapping=mapping)

    def __repr__(self):
        return f"""alignment:
        {self.alignment1}
        {''.join(['X|'[n1==n2] for n1, n2 in zip(self.alignment1, self.alignment2)])}
        {self.alignment2}"""

    def get_alignment(self):
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


class AlignmentChain(Alignment):
    def __init__(self, *alignments):
        self.alignments = alignments
        super().__init__(
            target=self.alignments[-1].target,
            mapping=self.get_mapping())

    def get_mapping(self):
        indices = self.alignments[0].mapping
        for alignment in self.alignments[1:]:
            indices = alignment.map_indices(indices)
        return indices

class RegionAlignment(Alignment):
    def __init__(self, sequence, start, end):
        if isinstance(sequence, data.Data):
            sequence = sequence.sequence
        self.sequence = sequence
        self.start = start
        self.end = end
        super().__init__(
            target=sequence[start-1, end-1],
            mapping=self.get_mapping())

    def get_mapping(self):
        mapping = np.full(len(self.sequence), -1, int)
        mapping[self.start-1:self.end] = np.arange(self.end-self.start+1, int)
        return mapping
