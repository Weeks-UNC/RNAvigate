from Bio.pairwise2 import align
import Bio.SeqIO


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
            #  AUC-UGGCU
            #  AUCGUG-CU
            #  012-3467 if not full
            #  012-34567 if full
            if nt1 == '-':
                alignment_map.append(-1)  # these will be masked when plotting
                i += 1
            elif nt2 == '-':
                if full:
                    alignment_map.append(i)
                    i += 1
            else:
                alignment_map.append(i)
                i += 1
        return alignment_map

    def get_colorby_sequence(self, colors='new'):
        """Returns list of mpl colors representing sequence.

        Args:
            sequence (str): string matching attribute with sequence to color
            colors (str, optional): Options are "new" or "old".
                    "new": A=blue, U=light blue, G=red, C=light red
                    "old": A=red, U=yellow, G=blue, C=green
                    Defaults to 'new'.

        Returns:
            list of mpl colors: color list representing sequence
        """
        seq = self.sequence
        colors = np.array([get_nt_color(nt.upper(), colors) for nt in seq])
        return colors

    def get_colorby_position(self, cmap='rainbow'):
        """Returns list of mpl colors that spans the rainbow. Fits length of
        given sequence.

        Args:
            sequence (str): string matching sequence to fit colorlist to
            cmap (str, optional): mpl colormap to use. Defaults to 'rainbow'.

        Returns:
            list of mpl colors: spectrum of colors with same length as sequence
        """
        cmap = plt.get_cmap(cmap)
        colors = np.array([cmap(n/self.length) for n in range(self.length)])
        return colors
