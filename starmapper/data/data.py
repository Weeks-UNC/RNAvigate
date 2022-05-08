from Bio.pairwise2 import align
import Bio.SeqIO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like


def get_nt_color(nt, colors="new"):
    nt_color = {"old": {"A": "#f20000",  # red
                        "U": "#f28f00",  # yellow
                        "G": "#00509d",  # blue
                        "C": "#00c200"},  # green
                "new": {"A": "#366ef0",  # blue
                        "U": "#9bb9ff",  # light blue
                        "G": "#f04c4c",  # red
                        "C": "#ffa77c"}  # light red
                }[colors][nt]
    return nt_color


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
        if full:
            return alignment_map, i
        return np.array(alignment_map)

    def get_colors(self, source, nt_colors='new', pos_cmap='turbo',
                   profile=None, ct=None):
        if source == "sequence":
            seq = self.sequence
            colors = np.array([get_nt_color(nt.upper(), nt_colors)
                               for nt in seq])
            return colors
        elif source == "position":
            cmap = plt.get_cmap(pos_cmap)
            colors = np.array([cmap(n/self.length)
                               for n in range(self.length)])
            return colors
        elif source == "profile":
            assert type(profile).__name__ == "Profile", "Invalid profile"
            if profile.datatype == 'RNP':
                cmap = np.array(["silver", "limegreen"])
                prof_colors = cmap[profile.data["RNPsite"]]
            elif profile.datatype in ['profile', 'dance']:
                cmap = np.array(['gray', 'black', 'orange', 'red'])
                bins = np.array([0, 0.4, 0.85])
                with np.errstate(invalid='ignore'):  # always false for nans
                    prof_colors = cmap[[sum(p > bins)
                                        for p in profile.data.Norm_profile]]
            colors = np.full(self.length, 'gray', dtype='<U16')
            am = profile.get_alignment_map(self)
            for i, i2 in enumerate(am):
                if i2 != -1:
                    colors[i2] = prof_colors[i]
            return colors
        elif source == "structure":
            assert type(ct).__name__ == "CT", "Invalid ct"
            cmap = np.array(['C0', 'C1'])
            ct_colors = cmap[[int(nt == 0) for nt in ct.ct]]
            colors = np.full(self.length, 'gray', dtype='<U8')
            am = ct.get_alignment_map(self)
            for i, i2 in enumerate(am):
                if i2 != -1:
                    colors[i2] = ct_colors[i]
            return colors
        elif (isinstance(source, list) and (len(source) == self.length)
              and all(is_color_like(c) for c in source)):
            return np.array(source)
        elif is_color_like(source):
            return np.full(self.length, source, dtype="<U16")
        else:
            print("Invalid colors: choices = profile, sequence, position, " +
                  "a list of mpl colors, or a single mpl color. " +
                  "Defaulting to sequence.")
            return self.get_colors("sequence")
