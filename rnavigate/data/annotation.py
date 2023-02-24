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
                 annotations=None,
                 site_list=None,
                 span_list=None,
                 groups=None,
                 primer_list=None,
                 color="blue"):
        self.name = name
        self.color = color
        self.datatype = datatype
        super().__init__(filepath=fasta, sequence=sequence)
        self.annotation_type = annotation_type
        for anno in [annotations, site_list, span_list, groups, primer_list]:
            if anno is not None:
                self._list = anno

    def read_sites(self, filepath, sep, read_csv_kw):
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)

    def fit_to(self, fit_to):
        am = self.get_alignment_map(fit_to=fit_to)

        def recursive_fit_to(indices):
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
                 filepath=None,
                 fasta=None,
                 sequence=None,
                 motif=None,
                 color="blue"):
        span_list = self.get_spans_from_motif(sequence, motif)
        super().__init__(name=name, fasta=fasta, sequence=sequence,
                         annotation_type='spans',
                         span_list=span_list, color=color)

    def get_spans_from_motif(self, sequence, motif):
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
