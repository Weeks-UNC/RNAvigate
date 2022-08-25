import re
from .data import Data
import pandas as pd


class Annotation(Data):
    def __init__(self, name=None, datatype="annotation",
                 dataframe=None, filepath=None, read_csv_kw={},
                 fasta=None, sequence=None,
                 site_list=None, span_list=None, groups=None,
                 motif=None, orf=False,
                 color="blue"):
        super().__init__(filepath=fasta, sequence=sequence)
        self.name = name
        self.color = color
        self.datatype = datatype
        if dataframe is not None:
            self.data = dataframe
        elif filepath is not None:
            self.read_sites(filepath=filepath, read_csv_kw=read_csv_kw)
        super().__init__(filepath=fasta, sequence=sequence,
                         dataframe=dataframe)

        self.sites = []
        self.spans = []
        self.groups = []
        if site_list is not None:
            self.annotation_type = "sites"
            self.sites = site_list
        elif span_list is not None:
            self.annotation_type = "spans"
            self.spans = span_list
        elif groups is not None:
            self.annotation_type = "groups"
            self.groups = groups
        elif motif is not None:
            self.annotation_type = "spans"
            self.get_spans_from_motif(motif)
        elif orf:
            self.annotation_type = "spans"
            self.get_spans_from_orf()

    def read_sites(self, filepath, sep, read_csv_kw):
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)

    def get_spans_from_motif(self, motif):
        nuc_codes = {"A": "A", "T": "T", "U": "U", "G": "G", "C": "C",
                     "B": "[CGTU]", "D": "[ATUG]", "H": "[ATUC]", "V": "[ACG]",
                     "W": "[ATU]", "S": "[CG]",  # strong and weak
                     "M": "[AC]", "K": "[GTU]",  # amino and ketone
                     "R": "[AG]", "Y": "[CTU]",  # purine and pyrimidine
                     "N": "[ATUGC]"}  # any nuc
        re_pattern = ''.join([nuc_codes[n] for n in motif])
        self.spans = []
        for match in re.finditer(re_pattern, self.sequence):
            start, end = match.span()
            self.spans.append([start+1, end])

    def get_spans_from_orf(self):
        spans = []
        stop_codons = "UAA|UAG|UGA"
        stop_sites = []
        for match in re.finditer(stop_codons, self.sequence):
            stop_sites.append(match.span()[1])
        start_codon = "AUG"
        start_sites = []
        for match in re.finditer(start_codon, self.sequence):
            start_sites.append(match.span()[0]+1)
        for start in start_sites:
            for stop in stop_sites:
                if (stop-start) % 3 == 2:
                    spans.append([start, stop])
        self.spans = spans
