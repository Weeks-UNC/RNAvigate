from Bio.pairwise2 import align
import Bio.SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpc

_alignments_cache = {}


def set_alignment(sequence1, sequence2, alignment1, alignment2):
    seq12 = {
        sequence2.upper().replace("T", "U"): {
            "seqA": alignment1,
            "seqB": alignment2}}
    seq21 = {
        sequence1.upper().replace("T", "U"): {
            "seqA": alignment2,
            "seqB": alignment1
        }
    }
    if sequence1 not in _alignments_cache.keys():
        _alignments_cache[sequence1] = seq12
    else:
        _alignments_cache[sequence1].update(seq12)
    if sequence2 not in _alignments_cache.keys():
        _alignments_cache[sequence2] = seq21
    else:
        _alignments_cache[sequence2].update(seq21)


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
    def __init__(self, sequence=None, filepath=None, dataframe=None):
        if sequence is not None:
            self.sequence = sequence
        elif filepath is not None:
            self.read_fasta(filepath)
        elif dataframe is not None:
            self.get_seq_from_data(dataframe)
        else:
            print(f"{self.datatype} initialized without sequence.")
        if not hasattr(self, "datatype"):
            self.datatype = "data"

    def read_fasta(self, fasta):
        fasta = list(Bio.SeqIO.parse(open(fasta), 'fasta'))
        self.sequence = str(fasta[0].seq).upper().replace("T", "U")
        self.gene = fasta[0].id

    def get_seq_from_data(self, dataframe):
        sequence = ''.join(dataframe["Sequence"].values)
        self.sequence = sequence.upper().replace("T", "U")

    @property
    def length(self):
        return len(self.sequence)

    def get_cmap(self, cmap):
        if mpc.is_color_like(cmap):
            cmap = mpc.ListedColormap([cmap])
        elif isinstance(cmap, list) and all(mpc.is_color_like(c) for c in cmap):
            cmap = mpc.ListedColormap(cmap)
        else:
            assert cmap in plt.colormaps(), ("cmap must be one of: valid mpl "
                                             "color, list of mpl colors, or "
                                             "mpl colormap."
                                             + str(cmap))
        cmap = plt.get_cmap(cmap)
        return cmap

    def set_alignment(self, fit_to, seqA, seqB):
        fit_to = fit_to.upper()
        self.alignments[fit_to] = {
            "seqA": seqA,
            "seqB": seqB,
        }

    def get_alignment_map(self, fit_to, full=False, print_sequences=False):
        seq1 = self.sequence.upper().replace("T", "U")
        seq2 = fit_to.sequence.upper().replace("T", "U")
        if seq1 == seq2:
            if print_sequences:
                print(seq1, seq2, sep="\n")
            return np.arange(len(self.sequence))
        else:
            try:
                alignment = _alignments_cache[seq1][seq2]
            except KeyError:
                alignment = align.globalxs(self.sequence, fit_to.sequence,
                                           -1, -0.1, penalize_end_gaps=False)
                set_alignment(
                    sequence1=seq1,
                    sequence2=seq2,
                    alignment1=alignment[0].seqA,
                    alignment2=alignment[0].seqB
                )
                alignment = _alignments_cache[seq1][seq2]
        if print_sequences:
            print(alignment["seqA"], alignment["seqB"], sep="\n")
        # get an index map from this sequence to that.
        alignment_map = []
        i = 0
        for nt1, nt2 in zip(alignment["seqA"], alignment["seqB"]):
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

    def get_colors(self, source, nt_colors='new', pos_cmap='rainbow',
                   profile=None, ct=None):
        if isinstance(source, str) and (source == "sequence"):
            seq = self.sequence
            colors = np.array([get_nt_color(nt.upper(), nt_colors)
                               for nt in seq])
            return colors
        elif isinstance(source, str) and (source == "position"):
            cmap = plt.get_cmap(pos_cmap)
            colors = np.array([cmap(n/self.length)
                               for n in range(self.length)])
            return colors
        elif isinstance(source, str) and (source == "profile"):
            prof_colors = profile.colors
            colors = np.full(self.length, 'gray', dtype='<U16')
            am = profile.get_alignment_map(self)
            for i, i2 in enumerate(am):
                if i2 != -1:
                    colors[i2] = mpc.to_hex(prof_colors[i])
            return colors
        elif isinstance(source, str) and (source == "structure"):
            # TODO: this should be implemented in CT object for reusability
            cmap = np.array(['C0', 'C1'])
            ct_colors = cmap[[int(nt == 0) for nt in ct.ct]]
            colors = np.full(self.length, 'gray', dtype='<U8')
            alignment_map = ct.get_alignment_map(self)
            for i, i2 in enumerate(alignment_map):
                if i2 != -1:
                    colors[i2] = ct_colors[i]
            return colors
        elif mpc.is_color_like(source):
            return np.full(self.length, source, dtype="<U16")
        elif ((len(source) == self.length)
              and all(mpc.is_color_like(c) for c in source)):
            return np.array(list(source))
        else:
            print("Invalid colors: choices = profile, sequence, position, " +
                  "a list of mpl colors, or a single mpl color. " +
                  "Defaulting to sequence.")
            return self.get_colors("sequence")
