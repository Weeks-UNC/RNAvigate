from Bio.pairwise2 import align
import Bio.SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
from ..styles import get_nt_color

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
        sequence1="AUCGAUCGGUACAUGUGAUGUAC"
        sequence2="AUCGAUCGAGCGUCAUGACGUCGAUGUACGUAC"
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
        sequence2.upper().replace("T", "U"): {
            "seqA": alignment1,
            "seqB": alignment2}}
    seq21 = {
        sequence1.upper().replace("T", "U"): {
            "seqA": alignment2,
            "seqB": alignment1}}

    # Store dictionaries in _alignments_cache
    if sequence1 not in _alignments_cache.keys():
        _alignments_cache[sequence1] = seq12
    else:
        _alignments_cache[sequence1].update(seq12)
    if sequence2 not in _alignments_cache.keys():
        _alignments_cache[sequence2] = seq21
    else:
        _alignments_cache[sequence2].update(seq21)


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
        """Constructs a Data object given a sequence string, fasta file, or
        dataframe containing a "Sequence" column.

        Args:
            sequence (str, optional): sequence string. Defaults to None.
            filepath (str, optional): path to fasta file. Defaults to None.
            dataframe (pandas DataFrame, optional): must contain a "Sequence"
                column. Defaults to None.
        """
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
        """Parse a fasta file for the first sequence. Store the sequence name
        as self.gene and the sequence string as self.sequence.

        Args:
            fasta (str): path to fasta file
        """
        fasta = list(Bio.SeqIO.parse(open(fasta), 'fasta'))
        self.sequence = str(fasta[0].seq).upper().replace("T", "U")
        self.gene = fasta[0].id

    def get_seq_from_data(self, dataframe):
        """Parse a dataframe for the sequence string, store as self.sequence.

        Args:
            dataframe (pandas DataFrame): must contain a "Sequence" column
        """
        sequence = ''.join(dataframe["Sequence"].values)
        self.sequence = sequence.upper().replace("T", "U")

    @property
    def length(self):
        """Get the length of the sequence

        Returns:
            int: the length of self.sequence
        """
        return len(self.sequence)

    def get_cmap(self, cmap):
        """Given a matplotlib color, list of colors, or colormap name, return
        a colormap object

        Args:
            cmap (str | tuple | float | list): A valid mpl color, list of valid
                colors or a valid colormap name

        Returns:
            matplotlib colormap: listed colormap matching the input
        """
        if mpc.is_color_like(cmap):
            cmap = mpc.ListedColormap([cmap])
        elif (isinstance(cmap, list)
              and all(mpc.is_color_like(c) for c in cmap)):
            cmap = mpc.ListedColormap(cmap)
        else:
            assert cmap in plt.colormaps(), ("cmap must be one of: valid mpl "
                                             "color, list of mpl colors, or "
                                             "mpl colormap."
                                             + str(cmap))
        cmap = plt.get_cmap(cmap)
        return cmap

    def get_alignment_map(self, fit_to, full=False, return_alignment=False):
        """Get an array to assist in repositioning data to the fit_to
        sequence. old_index (1-based) is mapped to new_index (1-based) thusly:
            new_index = alignment_map[old_index - 1] -1
        If full is False, new_index is a position within fit_to sequence
            deleted positions are dropped (new_idx == 0)
        If full is True, new_index is a position within the alignment of these
            two sequences. No positions are dropped.

        Args:
            fit_to (Data or subclass): Data object to fit new positions to
            full (bool, optional): Whether to drop deletions from the new index
                Defaults to False.
            return_alignment (bool, optional): whether to return the alignment
                instead of the alignment map. Defaults to False.

        Returns:
            numpy array: a mapping of old positions to new positions
        """
        # Normalize sequences
        seq1 = self.sequence.upper().replace("T", "U")
        seq2 = fit_to.sequence.upper().replace("T", "U")
        seq2_to_final_map = np.where([nt != '.' for nt in seq2])[0]
        seq2 = seq2.replace('.', '')
        # Check if sequences match
        if seq1 == seq2:
            if return_alignment:
                return (seq1, seq2, seq1, seq2)
            else:
                return seq2_to_final_map
        else:
            # look in _alignments_cache, if not found, do alignment
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
        # skip creating alignment_map if return_alignment
        if return_alignment:
            return (seq1, seq2, align1, align2)
        # get an index map from this sequence to that.
        seq1_to_align = np.where([nt != '-' for nt in align1])[0]
        if full:
            i = len(align1)
            return seq1_to_align, i
        else:
            align_mask = np.where([nt != '-' for nt in align2])[0]
            align_to_seq2 = np.full(len(align2), -1)
            align_to_seq2[~align_mask] = np.arange(len(seq2))
            seq1_to_seq2 = align_to_seq2[seq1_to_align]
            keepers = seq1_to_seq2[seq1_to_seq2 != -1]
            seq1_to_seq2[seq1_to_seq2 != -1] = seq2_to_final_map[keepers]
            return seq1_to_seq2

    def get_colors(self, source, nt_colors='new', pos_cmap='rainbow',
                   profile=None, ct=None, annotations=None):
        """Get a numpy array of colors that fits the current sequence.

        Args:
            source (str | array of color-like): One of the following:
                "position": colors represent position in sequence
                "sequence": colors represent nucleotide identity
                "annotations": colors represent sequence annotations
                "profile": colors represent per-nucleotide data
                "structure": colors represent base-pairing status
                matplotlib color-like: all colors are this color
                array of color like: must match length of sequence
            nt_colors (str, optional): 'new' or 'old' as defined in
                rnavigate.style. Defaults to 'new'.
            pos_cmap (str, optional): cmap used if source="position".
                Defaults to 'rainbow'.
            profile (Profile or subclass, optional): Data object containing
                per-nucleotide information. Defaults to None.
            ct (CT of subclass, optional): Data object containing secondary
                structure information. Defaults to None.
            annotations (list of Annotations or subclass, optional): list of
                Data objects containing annotations. Defaults to None.

        Returns:
            numpy array: one matplotlib color-like value for each nucleotide in
                self.sequence
        """
        if isinstance(source, str) and (source == "sequence"):
            seq = self.sequence
            colors = np.array([get_nt_color(nt, nt_colors) for nt in seq])
            return colors
        elif isinstance(source, str) and (source == "position"):
            cmap = plt.get_cmap(pos_cmap)
            cmap_values = np.arange(self.length)/self.length
            return cmap(cmap_values)
        elif isinstance(source, str) and (source == "profile"):
            prof_colors = profile.colors
            colors = np.full(self.length, 'gray', dtype='<U16')
            am = profile.get_alignment_map(self)
            for i, i2 in enumerate(am):
                if i2 != -1:
                    colors[i2] = mpc.to_hex(prof_colors[i])
            return colors
        elif isinstance(source, str) and (source == "annotations"):
            colors = np.full(self.length, 'gray', dtype='<U16')
            for annotation in annotations:
                for site, color in zip(*annotation.get_sites_colors(self)):
                    colors[site-1] = color
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
            print('Invalid colors:\n\tchoices are "profile", "sequence", '
                  '"position", "structure", "annotations", a list of mpl '
                  'colors, or a single mpl color.\nDefaulting to sequence.')
            return self.get_colors("sequence")
