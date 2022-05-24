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

    def set_dms_profile(self):
        """Perform normalization of data based on DMS reactivities (AC seperate
        from GU), and replace data columns "Norm_profile" and "Norm_stderr"
        """
        profile = self.data["HQ_profile"]
        error = self.data["HQ_stderr"]
        # initialize the profile and error array
        dms_profile = np.array(profile)
        with np.errstate(invalid='ignore'):
            norm_error = np.zeros(error.shape)
            mask = profile > 0
            norm_error[mask] = (error[mask]/profile[mask])**2

        def get_norm_factors(data):
            finite_data = data[np.isfinite(data)]
            lower, upper = np.percentile(finite_data, [90., 99.])  # 99
            mask = (finite_data >= lower) & (finite_data <= upper)
            norm_set = finite_data[mask]
            average = np.mean(norm_set)
            std = np.std(norm_set)
            return average, std/np.sqrt(len(norm_set))

        ac_mask = np.isin(self.data["Sequence"], ['A', 'C'])
        ac_norm_factor, ac_norm_error = get_norm_factors(dms_profile[ac_mask])
        gu_mask = np.isin(self.data["Sequence"], ['G', 'U'])
        gu_norm_factor, gu_norm_error = get_norm_factors(dms_profile[gu_mask])
        dms_profile[ac_mask] /= ac_norm_factor
        dms_profile[gu_mask] /= gu_norm_factor
        if norm_error is not None:
            norm_error[ac_mask] += (ac_norm_error/ac_norm_factor)**2
            norm_error[gu_mask] += (gu_norm_error/gu_norm_factor)**2
            for mask in [ac_mask, gu_mask]:
                norm_error[mask] = np.sqrt(norm_error[mask])
                norm_error[mask] *= np.abs(dms_profile[mask])

        norm_factors = {'G': gu_norm_factor,
                        'U': gu_norm_factor,
                        'A': ac_norm_factor,
                        'C': ac_norm_factor}

        self.data["Norm_profile"] = dms_profile
        self.data["Norm_stderr"] = norm_error
