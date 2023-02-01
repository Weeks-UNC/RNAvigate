import pandas as pd
from .data import Data
import numpy as np
import matplotlib.colors as mpc
import matplotlib.pyplot as plt


class Profile(Data):
    def __init__(self, datatype="profile",
                 column=None, ap_scale_factor=1,
                 dataframe=None, sequence=None, fasta=None,
                 filepath=None, sep=None, read_csv_kw=None,
                 cmap=None, norm_method=None, norm_values=None,
                 color_column=None, colors=None):
        if sep is None:
            sep = '\t'
        if read_csv_kw is None:
            read_csv_kw = {}
        self.datatype = datatype
        self.default_column = column
        self.ap_scale_factor = ap_scale_factor
        self.set_color_defaults(column=color_column, cmap=cmap,
                                norm_method=norm_method,
                                norm_values=norm_values)
        # assign data
        if dataframe is not None:
            self.data = dataframe
        elif filepath is not None:
            self.read_file(filepath=filepath, sep=sep, read_csv_kw=read_csv_kw)
        else:
            print(f"{self.datatype} initialized without data.")
        # assign sequence
        super().__init__(filepath=fasta, sequence=sequence,
                         dataframe=self.data)
        # assign colors
        if colors is not None:
            self.colors = colors
        elif "Colors" in self.data.columns:
            self.colors = self.data["Colors"]
        else:
            self.set_profile_colors(start_from="defaults")

    def read_file(self, filepath, sep, read_csv_kw):
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)
        self.filepath = filepath

    def set_color_defaults(self, column, cmap, norm_method, norm_values):
        self._color_defaults = {"column": column,
                                "cmap": self.get_cmap(cmap),
                                "norm_method": norm_method,
                                "norm_values": norm_values}

    def set_color_values(self, start_from, **kwargs):
        if start_from == "defaults":
            self.color_values = self._color_defaults
        elif start_from == "current":
            pass
        for key in kwargs:
            if key not in self.color_values:
                print(f"{key} is not a valid color argument.")
            elif key == "cmap":
                self.color_values[key] = self.get_cmap(kwargs[key])
            elif kwargs[key] is not None:
                self.color_values[key] = kwargs[key]

    def set_profile_colors(self, start_from="current", **kwargs):
        self.set_color_values(start_from=start_from, **kwargs)
        cv = self.color_values
        if cv["norm_method"] == "bins":
            norm = mpc.BoundaryNorm(cv["norm_values"],
                                    cv["cmap"].N-1, extend="both")
        elif cv["norm_method"] == "min_max":
            norm = plt.Normalize(cv["norm_values"][0], cv["norm_values"][1])
        elif cv["norm_method"] == "0_1":
            norm = plt.Normalize()
        elif cv["norm_method"] is None:
            def norm(x): return x  # does nothing to values
        values = self.data[cv["column"]].values
        values = norm(values)
        self.colors = cv["cmap"](values)


class SHAPEMaP(Profile):
    def __init__(self, filepath, dms=False, datatype="shapemap",
                 read_csv_kw=None, **kwargs):
        if filepath.endswith(".map"):
            read_csv_kw = {"names": ["Nucleotide", "Norm_profile",
                                     "Norm_stderr", "Sequence"],
                           "na_values": "-999"}
        super().__init__(filepath=filepath,
                         datatype=datatype,
                         sep="\t",
                         read_csv_kw=read_csv_kw,
                         column="Norm_profile",
                         ap_scale_factor=5,
                         color_column="Norm_profile",
                         cmap=["grey", "black", "orange", "grey", "red"],
                         norm_method="bins",
                         norm_values=[-0.4, 0.4, 0.85],
                         **kwargs)
        if dms:
            self.ap_scale_factor = 10
            self.set_dms_profile()
            self.set_profile_colors(start_from="defaults")

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

        self.norm_factors = {'G': gu_norm_factor,
                             'U': gu_norm_factor,
                             'A': ac_norm_factor,
                             'C': ac_norm_factor}

        self.data["Norm_profile"] = dms_profile
        self.data["Norm_stderr"] = norm_error


class DanceMaP(SHAPEMaP):
    def __init__(self, filepath, component, datatype="dancemap", **kwargs):
        self.component = component
        super().__init__(filepath=filepath, datatype=datatype,
                         read_csv_kw={}, **kwargs)
        self.ap_scale_factor = 10

    def read_file(self, filepath, sep='\t', read_csv_kw={}):
        # parse header
        self.filepath = filepath
        with open(self.filepath) as f:
            header1 = f.readline()
            header2 = f.readline()
        self.header = header1 + header2
        self.components = int(header1.strip().split()[0])
        self.percents = [float(x) for x in header2.strip().split()[1:]]
        self.percent = self.percents[self.component]
        # parse datatable
        read_csv_kw['names'] = ["Nucleotide", "Sequence", "Norm_profile",
                                "Modified_rate", "Untreated_rate"]
        col_offset = 3 * self.component
        bg_col = 3 * self.components + 2
        read_csv_kw["usecols"] = [0, 1, 2+col_offset, 3+col_offset, bg_col]
        self.data = pd.read_csv(self.filepath, sep=sep, header=2,
                                **read_csv_kw)
        # some rows have an "i" added to the final column
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


class RNPMaP(Profile):
    def __init__(self, filepath, datatype="rnpmap", **kwargs):
        super().__init__(filepath=filepath,
                         datatype=datatype,
                         sep=",",
                         column="NormedP",
                         color_column="RNPsite",
                         cmap=["silver", "limegreen"],
                         **kwargs)


class DeltaProfile(Profile):
    def __init__(self, profile1, profile2, column=None, norm_method="min_max",
                 norm_values=[-0.8, 0.8], cmap="coolwarm"):
        if column is None:
            column = profile1.default_column
        columns = ["Nucleotide", "Sequence", column]
        new_data = profile1.data[columns].copy()
        new_data = pd.merge(profile1.data[columns], profile2.data[columns],
                            how="outer", on=["Nucleotide", "Sequence"],
                            suffixes=["_1", "_2"])
        new_data.eval(f"Delta_profile = {column}_1 - {column}_2",
                      inplace=True)
        super().__init__(datatype="deltaprofile",
                         column="Delta_profile",
                         ap_scale_factor=5,
                         dataframe=new_data,
                         cmap=cmap,
                         norm_method=norm_method,
                         norm_values=norm_values,
                         color_column="Delta_profile")
