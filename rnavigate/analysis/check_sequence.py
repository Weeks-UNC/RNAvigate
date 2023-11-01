import matplotlib.pyplot as plt
import seaborn as sns
from difflib import SequenceMatcher
from rnavigate.data import set_alignment

class SequenceChecker():
    def __init__(self, samples):
        self.samples = samples
        self.sequences = []
        self.data_keywords = list(set([dkw  for sample in samples for dkw in sample.data]))
        # 1 row per sequence, 1 column per data keyword, value==sequence # or N/A
        self.sample_dkw_seq = []
        # group each sample by sequence, print the sequence, sample names and data type
        for sample in samples:
            dkw_seq = []
            for dkw in self.data_keywords:
                if dkw not in sample.data:
                    which_seq = 'N/A'
                else:
                    seq = sample.get_data(dkw).sequence
                    seq = seq.upper().replace('T', 'U')
                    if seq not in self.sequences:
                        self.sequences.append(seq)
                    which_seq = self.sequences.index(seq)
                dkw_seq.append(which_seq)
            self.sample_dkw_seq.append(dkw_seq)
        num_seq = range(len(self.sequences))
        self.alignments = [[None for _ in num_seq] for _ in num_seq]

    def plot_which_sequence(self):
        colors = sns.color_palette("rainbow", len(self.sequences))
        colors = {i: color for i, color in enumerate(colors)} | {'N/A': 'w'}
        colors = [[colors[col] for col in row] for row in self.sample_dkw_seq]
        num_samples = len(self.samples)
        num_dkws = len(self.data_keywords)
        fig, ax = plt.subplots(1, figsize=(1.5*num_dkws, 0.75*num_samples))
        ax.table(
            cellText=self.sample_dkw_seq,
            cellColours=colors,
            cellLoc='center',
            bbox=[0, 0, 1, 1],
            rowLabels=[sample.sample for sample in self.samples],
            rowLoc='center',
            colLabels=[f'"{dkw}"' for dkw in self.data_keywords],
        )
        ax.set(yticks=[], xticks=[])
        return fig, ax

    def get_alignments(self, which='all', print_format=None):
        kwargs = {'print_format': print_format}
        if which == 'all':
            num = range(len(self.sequences))
            return [self.get_alignments(which=(s1, s2), **kwargs
                    ) for s1 in num for s2 in num if s1 < s2]
        i, j = which
        seq_match = self.alignments[i][j]
        a = self.sequences[i].upper().replace('T', 'U')
        b = self.sequences[j].upper().replace('T', 'U')
        if seq_match is None:
            seq_match = SequenceMatcher(a=a, b=b, autojunk=False)
            self.alignments[i][j] = seq_match
        if print_format is None:
            return seq_match
        elif print_format == 'cigar':
            tags = {
                'insert': 'I', 'delete': 'D', 'equal': '=', 'replace': 'X',
            }
            cigar = ''
            for tag, i1, i2, j1, j2 in seq_match.get_opcodes():
                if tag in ['equal', 'replace', 'delete']:
                    length = i2-i1
                elif tag == 'insert':
                    length = j2-j1
                cigar += f'{length}{tags[tag]}'
            print(f'Sequence {i} to {j}:')
            print(cigar)
        elif print_format == "long":
            print(f'Sequence {i} to {j}:')
            for tag, i1, i2, j1, j2 in seq_match.get_opcodes():
                s_a = a[i1:i2]
                s_b = b[j1:j2]
                if tag == 'equal':
                    continue
                elif tag == 'replace':
                    s = f'{"mismatch":<9} {i1+1:>6} {s_a} --> {s_b}'
                elif tag == 'delete':
                    s = f'{"delete":<9} {i1+1:>6} {s_a}'
                elif tag == 'insert':
                    s = f'{"insert":<9} {i1+1:>6} {s_b}'
                print(s)
            print()

    def set_as_default_alignments(self, which='all'):
        if which == 'all':
            for i in range(len(self.sequences)):
                for j in range(len(self.sequences)):
                    if i < j:
                        self.set_as_default_alignments(which=(i, j))
            return
        i, j = which
        seq_match = self.get_alignments(which=(i, j))
        alignment1, alignment2 = '', ''
        sequence1 = self.sequences[i]
        sequence2 = self.sequences[j]
        for tag, i1, i2, j1, j2 in seq_match.get_opcodes():
            s_a = sequence1[i1:i2]
            s_b = sequence2[j1:j2]
            if tag in ['equal', 'replace']:
                alignment1 += s_a
                alignment2 += s_b
            elif tag == 'delete':
                alignment1 += s_a
                alignment2 += '-' * len(s_a)
            elif tag == 'insert':
                alignment1 += '-' * len(s_b)
                alignment2 += s_b
        set_alignment(sequence1, sequence2, alignment1, alignment2)
