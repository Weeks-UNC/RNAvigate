import re
from .data import Data
import pandas as pd


class Annotation(Data):
    def __init__(self,
                 name=None,
                 datatype="annotation",
                 filepath="",
                 fasta=None,
                 sequence=None,
                 annotation_type=None,
                 site_list=None,
                 span_list=None,
                 groups=None,
                 primer_list=None,
                 annotations=None,
                 color="blue"):
        """Base annotation class to store 1D features of an RNA. This can
        include groups of separted nucleotides (e.g. binding pocket), spans of
        nucleotides (e.g. coding sequence, Alu elements), list of sites (e.g.
        m6A locations) or primer binding sites.

        Args:
            name (str, optional): Name of annotation.
                Defaults to None.
            datatype (str, optional): annotation datatype.
                Defaults to "annotation".
            filepath (str, optional): path to file containing annotations.
                Defaults to "".
            fasta (str, optional): path to fasta file containing one sequence.
                Defaults to None.
            sequence (str, optional): Nucleotide sequence.
                Defaults to None.
            annotation_type (str, optional): "groups", "sites", "spans", or
                "primers". Must match the type used below.
                Defaults to None.
            site_list (list of int, optional): 1-indexed location of sites of
                interest within sequence or fasta.
                Defaults to None.
            span_list (list of pairs, optional): 1-indexed locations of spans
                of interest within sequence or fasta.
                e.g. [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30.
                Defaults to None.
            groups (list of lists, optional): 1-indexed locations of groups of
                sites of interest.
                Defaults to None.
            primer_list (list of pairs, optional): Similar to span_list above,
                but reverse primers are reverse ordered.
                e.g. [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20.
                Defaults to None.
            annotations (list, optional): catchall for the above list formats.
                List will be treated according to annotation_type argument.
                Defaults to None.
            color (matplotlib color-like, optional): Color to be used for
                displaying this annotation on plots.
                Defaults to "blue".
        """
        self.name = name
        self.color = color
        self.datatype = datatype
        super().__init__(filepath=fasta, sequence=sequence)
        self.annotation_type = annotation_type
        # TODO: make sure whichever list is stored matches expected format
        for anno in [annotations, site_list, span_list, groups, primer_list]:
            if anno is not None:
                self._list = anno

    # TODO: create a method for reading annotations from a file.
    def read_sites(self, filepath, sep, read_csv_kw):
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)

    # TODO: implement fit_to for Motif subclass
    def fit_to(self, fit_to):
        """Creates a new Annotation, stored as self.fitted, which maps the
        indices to a new sequence.

        Args:
            fit_to (rnavigate.data.Data): A data object containing a sequence.
        """
        am = self.get_alignment_map(fit_to=fit_to)

        def recursive_fit_to(indices):
            # given nested list of indices, returns same shape with indices
            # mapped to fit_to sequence
            new_list = []
            for idx in indices:
                if isinstance(idx, list):
                    new_list.append(recursive_fit_to(idx))
                else:
                    new_list.append(am[idx-1]+1)
            return new_list

        self.fitted = Annotation(
            name=self.name,
            color=self.color,
            sequence=fit_to,
            annotation_type=self.annotation_type,
            annotations=recursive_fit_to(self._list)
        )

    def __getitem__(self, i):
        return self._list[i]

    def __len__(self):
        return len(self._list)


class Motif(Annotation):
    def __init__(self,
                 name=None,
                 filepath="",
                 fasta=None,
                 sequence=None,
                 motif=None,
                 color="blue"):
        """Creates a Motif annotation, which acts like a span Annotation, for
        highlighting a sequence motif of interest, given with conventional
        nucleotide codes. e.g. "DRACH"

        Args:
            name (str, optional): name of this annotation.
                Defaults to None.
            filepath(str, optional): Defaults to "".
            fasta (str, optional): path to fasta file containing sequence.
                Defaults to None.
            sequence (str, optional): sequence to be searched.
                Defaults to None.
            motif (str, optional): sequence motif to be searched for.
                Defaults to None.
            color (str, optional): color used to display these motif locations.
                Defaults to "blue".
        """
        span_list = self.get_spans_from_motif(sequence, motif)
        super().__init__(name=name, fasta=fasta, sequence=sequence,
                         annotation_type='spans',
                         span_list=span_list, color=color)

    def get_spans_from_motif(self, sequence, motif):
        """Returns a list of spans [[start, end], [start, end]] for each
        location of motif found within sequence, using conventional nucleotide
        codes.

        Args:
            sequence (str): sequence to be searched
            motif (str): sequence motif to be searched for.

        Returns:
            _type_: _description_
        """
        nuc_codes = {"A": "A", "T": "T", "U": "U", "G": "G", "C": "C",
                     "B": "[CGTU]", "D": "[ATUG]", "H": "[ATUC]", "V": "[ACG]",
                     "W": "[ATU]", "S": "[CG]",  # strong and weak
                     "M": "[AC]", "K": "[GTU]",  # amino and ketone
                     "R": "[AG]", "Y": "[CTU]",  # purine and pyrimidine
                     "N": "[ATUGC]"}  # any nuc
        re_pattern = ''.join([nuc_codes[n] for n in motif])
        spans = []
        for match in re.finditer(re_pattern, sequence):
            start, end = match.span()
            spans.append([start+1, end])
        return spans


class ORFs(Annotation):
    def __init__(self,
                 name=None,
                 filepath=None,
                 fasta=None,
                 sequence=None,
                 color="blue"):
        span_list = self.get_spans_from_orf(sequence)
        super().__init__(name=name, fasta=fasta, sequence=sequence,
                         annotation_type='spans',
                         span_list=span_list, color=color)

    def get_spans_from_orf(self, sequence):
        spans = []
        stop_codons = "UAA|UAG|UGA"
        stop_sites = []
        for match in re.finditer(stop_codons, sequence):
            stop_sites.append(match.span()[1])
        start_codon = "AUG"
        start_sites = []
        for match in re.finditer(start_codon, sequence):
            start_sites.append(match.span()[0]+1)
        for start in start_sites:
            for stop in stop_sites:
                if ((stop-start) % 3 == 2) and (start < stop):
                    spans.append([start, stop])
        return spans
