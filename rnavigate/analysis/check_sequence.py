"""This module contains the SequenceChecker used to inspect the sequences.

Given a list of samples, we can inspect which data keywords belong to the
samples, which sequences match up perfectly, and inspect the differences
between sequences.
"""
import pandas as pd
import numpy as np
from rnavigate import data


class SequenceChecker():
    """Check the sequences stored in a list of samples.

    Attributes:
        samples: the list of samples
        sequences: a list of all unique sequence strings stored in the list of
            samples. These are converted to an all uppercase RNA alphabet.
        keywords: a list of unique data keywords stored in the list of samples.
    """
    def __init__(self, samples):
        """Creates an instance of SequenceChecker given a list of samples

        Arguments:
            samples (list of rnav.Sample)
                samples for which to compare data keywords and sequences.
        """
        self.samples = samples
        self._keywords = []
        self._sequences = []

    @property
    def keywords(self):
        """A list of all unique data keywords across samples."""
        return list(set([dkw for s in self.samples for dkw in s.data]))

    @property
    def sequences(self):
        """A list of all unique sequences (uppercase RNA) across samples."""
        sequences = []
        for sample in self.samples:
            for dkw in self.keywords:
                if dkw not in sample.data:
                    continue
                seq = sample.get_data(dkw).sequence
                seq = seq.upper().replace("T", "U")
                if seq not in self.sequences:
                    sequences.append(seq)
        return sequences

    @property
    def which_sequences(self):
        """A DataFrame of sequence IDs (integers) for each data keyword."""
        sequences = self.sequences
        data_keywords = self.keywords
        df = {key: [] for key in ["Sample"] + data_keywords}
        for sample in self.samples:
            df["Sample"].append(sample.sample)
            for dkw in data_keywords:
                if dkw not in sample.data:
                    df[dkw].append(np.nan)
                else:
                    seq = sample.get_data(dkw).sequence
                    seq = seq.upper().replace("T", "U")
                    df[dkw].append(sequences.index(seq))
        self.which_sequences = pd.DataFrame(df)

    def print_which_sequences(self):
        """Print sequence ID (integer) for each data keyword and sample."""
        which_sequences = self.which_sequences
        for _, row in which_sequences.iterrows():
            for column in which_sequences.columns:
                if column == "Sample":
                    print(row[column])
                else:
                    if np.isnan(row[column]):
                        continue
                    print(f"    '{column}': Sequence {row[column]}")
            print()

    def print_alignments(self, print_format="long", which="all"):
        """Print alignments in the given format for sequence IDs provided.

        Optional arguments:
            print_format (string)
                What format to print the alignments in:
                "cigar" prints the cigar string
                "short" prints the numbers of mismatches and indels
                "long" prints the location and nucleotide identity of all
                    mismatches, insertions and deletions.
                Defaults to "long".
            which (pair of integers)
                two sequence IDs to compare.
                Defaults to every pairwise comparison.
        """
        kwargs = {"print_format": print_format}
        if which == "all":
            num = range(len(self.sequences))
            for s1 in num:
                for s2 in num:
                    if s1 < s2:
                        self.print_alignments(which=(s1, s2), **kwargs)
            return
        i, j = int(which[0]), int(which[1])
        a = self.sequences[i]
        b = self.sequences[j]
        alignment = data.SequenceAlignment(sequence1=a, sequence2=b, full=True)
        print(f"Sequence {i} to {j}:")
        alignment.print(print_format=print_format)

    def write_fasta(self, filename, which="all"):
        """Write all unique sequences to a fasta file.

        This is very useful for using external multiple sequence aligners such
        as ClustalOmega.

            1) go to https://www.ebi.ac.uk/Tools/msa/clustalo/
            2) upload new fasta file
            3) under STEP 2 output format, select Pearson/FASTA
            4) click 'Submit'
            5) wait for your alignment to finish
            6) download the alignment fasta file
            7) use rnav.data.set_alignments_from_pearsonfasta() to set the new
               default alignments

        Required arguments:
            filename (string)
                path to a new file to which fasta entries are written

        Optional arguments:
            which (list of integers)
                Sequence IDs to write to file.
        """
        sequences = self.sequences
        if which == "all":
            which = range(len(sequences))
        with open(filename, "w") as fasta:
            for i in which:
                fasta.write(f">Sequence_{i}\n")
                fasta.write(f"{sequences[i]}\n")
