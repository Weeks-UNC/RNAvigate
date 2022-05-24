import re
import numpy as np
from .data import Data
import pandas as pd


class Annotation(Data):
    def __init__(self, name, filepath=None, fasta=None, sequence=None,
                 sites=None, motif=None, color="blue"):
        self.name = name
        self.color = color
        if filepath is not None:
            self.read_sites(filepath)
        elif fasta is not None:
            self.read_fasta(fasta)
        elif sequence is not None:
            self.sequence = sequence
        if sites is not None:
            assert all(x in [0, 1] for x in sites), "sites should be 0s and 1s"
            assert len(self.sequence) == len(
                sites), "Sequence and sites lengths differ."
            self.sites = np.array(sites, dtype=bool)
        elif motif is not None:
            self.get_sites_from_motif(motif)

    def read_sites(self, filepath):
        data = pd.read_csv(filepath, sep='\t')
        self.sequence = ''.join(data["Sequence"].values)
        self.sites = np.array(data["sites"].values, dtype=bool)

    def get_sites_from_motif(self, motif):
        nuc_codes = {"A": "A", "T": "T", "U": "U", "G": "G", "C": "C",
                     "B": "[CGTU]", "D": "[ATUG]", "H": "[ATUC]", "V": "[ACG]",
                     "W": "[ATU]", "S": "[CG]",  # strong and weak
                     "M": "[AC]", "K": "[GTU]",  # amino and ketone
                     "R": "[AG]", "Y": "[CTU]",  # purine and pyrimidine
                     "N": "[ATUGC]"}  # any nuc
        self.sites = np.zeros(len(self.sequence), dtype=bool)
        re_pattern = ''.join([nuc_codes[n] for n in motif])
        for match in re.finditer(re_pattern, self.sequence):
            start, end = match.span()
            for i in range(start, end):
                self.sites[i] = True
