import pandas as pd
from .data import Data
import numpy as np


class Profile(Data):
    def __init__(self, filepath, datatype="profile", component=None):
        self.datatype = datatype
        if datatype == "profile":
            self.read_profile(filepath, '\t')
        elif datatype == "RNP":
            self.read_profile(filepath, ',')
        elif datatype == "dance":
            assert component is not None, "Must pass a dance component number"
            self.read_dance_reactivities(filepath, component)

    def read_profile(self, profile, sep):
        self.data = pd.read_csv(profile, sep=sep)
        sequence = ''.join(self.data["Sequence"].values)
        self.sequence = sequence.upper().replace("T", "U")

    def read_dance_reactivities(self, filepath, component):
        # parse header
        with open(filepath) as f:
            header1 = f.readline()
            header2 = f.readline()
        self.header = header1 + header2
        self.components = int(header1.strip().split()[0])
        self.percents = [float(x) for x in header2.strip().split()[1:]]
        self.component = component
        self.percent = self.percents[component]
        # parse datatable
        read_kwargs = {}
        read_kwargs['names'] = ["Nucleotide", "Sequence", "Norm_profile",
                                "Modified_rate", "Untreated_rate"]
        col_offset = 3 * component
        bg_col = 3 * self.components + 2
        read_kwargs["usecols"] = [0, 1, 2+col_offset, 3+col_offset, bg_col]
        self.data = pd.read_csv(filepath, sep='\t', header=2, **read_kwargs)
        stripped = []
        for x in self.data["Untreated_rate"]:
            if type(x) is float:
                stripped.append(x)
            else:
                stripped.append(float(x.rstrip(' i')))
        self.data["Untreated_rate"] = stripped
        self.data["Reactivity_profile"] = (self.data["Modified_rate"] -
                                           self.data["Untreated_rate"])
        sequence = ''.join(self.data["Sequence"].values)
        self.sequence = sequence.upper().replace("T", "U")
