from Bio.pairwise2 import align
import Bio.SeqIO
import numpy as np
import matplotlib.pyplot as plt


def get_pairs_sens_PPV(self, ct="ct"):
    "Returns sensitivity and PPV for pair data to the ct structure"
    import pmanalysis as pma
    pm = pma.PairMap(self.paths["pairs"])
    ct = getattr(self, ct).copy()
    ct.filterNC()
    ct.filterSingleton()
    p, s = pm.ppvsens_duplex(ct, ptype=1, exact=False)
    return p, s


class Data():
    def __init__(self, sequence=None, fasta=None):
        if sequence is not None:
            self.sequence = sequence
        elif fasta is not None:
            self.read_fasta(fasta)

    def read_fasta(self, fasta):
        fasta = list(Bio.SeqIO.parse(open(fasta), 'fasta'))
        self.sequence = str(fasta[0].seq).upper().replace("T", "U")
        self.gene = fasta[0].id

    @property
    def length(self):
        return len(self.sequence)

    def get_alignment_map(self, fit_to, full=False):
        alignment = align.globalxs(self.sequence, fit_to.sequence, -1, -0.1,
                                   penalize_end_gaps=False)
        # get an index map from this sequence to that.
        alignment_map = []
        i = 0
        for nt1, nt2 in zip(alignment[0].seqA, alignment[0].seqB):
            #  012-34567 index 1
            #  AUC-UGGCU sequence 1
            #  AUCGUG-CU sequence 2
            #  012345-67 index 2
            #  012345678 index "full"
            #  desired: alignmen_map[index 1] == index 2
            #  012 45-67 not full
            # desired: alignment_map[index 1] == index "full"
            #  012 45678 full
            if nt1 == '-':
                i += 1
            elif nt2 == '-':
                if full:
                    alignment_map.append(i)
                    i += 1
                else:
                    alignment_map.append(-1)
            else:
                alignment_map.append(i)
                i += 1
        return alignment_map

    def get_colorby_sequence(self, colors='new'):
        seq = self.sequence
        colors = np.array([get_nt_color(nt.upper(), colors) for nt in seq])
        return colors

    def get_colorby_position(self, cmap='rainbow'):
        cmap = plt.get_cmap(cmap)
        colors = np.array([cmap(n/self.length) for n in range(self.length)])
        return colors
