import pandas as pd
from .data import Data
from .fragmap import Fragmapper
import numpy as np
import matplotlib.colors as mpc
import matplotlib.pyplot as plt


class Profile(Data):
    def __init__(self, datatype="profile",
                 column=None, err_column=None, ap_scale_factor=1,
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
        self.default_err_column = err_column
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
        if (colors is None) and ("Colors" in self.data.columns):
            colors = self.data["Colors"]
        else:
            self.set_color_values(start_from="defaults", colors=colors)

    def read_file(self, filepath, sep, read_csv_kw):
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)
        self.filepath = filepath

    def set_color_defaults(self, column, cmap, norm_method, norm_values):
        self._color_defaults = {"column": column,
                                "cmap": cmap,
                                "norm_method": norm_method,
                                "norm_values": norm_values}

    def set_color_values(self, start_from="current", column=None, cmap=None,
                         norm_method=None, norm_values=None, colors=None):
        if start_from == "defaults":
            self.color_values = self._color_defaults
        elif start_from == "current":
            pass
        if column is not None:
            self.color_values["column"] = column
        if cmap is not None:
            self.color_values["cmap"] = cmap
        if norm_method is not None:
            self.color_values["norm_method"] = norm_method
        if norm_values is not None:
            self.color_values["norm_values"] = norm_values
        if colors is not None:
            self._colors = colors
        elif colors is None:
            self._colors = None

    @property
    def colors(self):
        if self._colors is not None:
            return self._colors
        cv = self.color_values
        cmap = self.get_cmap(cv["cmap"])
        if cv["norm_method"] == "bins":
            norm = mpc.BoundaryNorm(cv["norm_values"], cmap.N, extend="both")
        elif cv["norm_method"] == "min_max":
            norm = plt.Normalize(cv["norm_values"][0], cv["norm_values"][1])
        elif cv["norm_method"] == "0_1":
            norm = plt.Normalize()
        elif cv["norm_method"] == "none":
            def norm(x): return x  # does nothing to values
        values = self.data[cv["column"]].values
        colors = np.full(len(values), "#888888")
        mask = ~np.isnan(values)
        values = norm(values[mask])
        colors[mask] = np.array([mpc.to_hex(color) for color in cmap(values)])
        return colors

    def fit_to(self, fit_to):
        am = self.get_alignment_map(fit_to=fit_to)
        self.data["nt_offset"] = [am[nt-1]+1 for nt in self.data["Nucleotide"]]
        self.data["mask"] = self.data["nt_offset"] != 0
        self.fit_to_length = fit_to.length

    def get_plotting_dataframe(self, all_columns=False, column=None,
                               err_column=None):
        if all_columns:
            columns = self.data.columns.values
            columns = np.delete(columns, columns == "Nucleotide")
            data = self.data.loc[self.data["mask"], columns].copy()
            data = data.rename(columns={"nt_offset": "Nucleotide"})
        else:
            columns = ["nt_offset", "mask"]
            column_names = ["Nucleotide", "mask"]
            if column is not None:
                columns.append(column)
            else:
                columns.append(self.default_column)
            column_names.append("Values")
            if err_column is not None:
                columns.append(err_column)
                column_names.append("Errors")
            elif self.default_err_column is not None:
                columns.append(self.default_err_column)
                column_names.append("Errors")
            data = self.data.loc[self.data["mask"], columns].copy()
            data.columns = column_names
        data["Colors"] = self.colors[self.data["mask"]]
        plotting_df = pd.DataFrame(
            {"Nucleotide": np.arange(self.fit_to_length)+1})
        plotting_df = plotting_df.merge(data, how='outer', on="Nucleotide")
        plotting_df.fillna({"Colors": mpc.to_hex("grey")}, inplace=True)
        return plotting_df


class SHAPEMaP(Profile):
    def __init__(self, filepath, dms=False, datatype="shapemap",
                 read_csv_kw=None, err_column="Norm_stderr", **kwargs):
        if filepath.endswith(".map"):
            read_csv_kw = {"names": ["Nucleotide", "Norm_profile",
                                     "Norm_stderr", "Sequence"],
                           "na_values": "-999"}
        super().__init__(filepath=filepath,
                         datatype=datatype,
                         sep="\t",
                         read_csv_kw=read_csv_kw,
                         column="Norm_profile",
                         err_column=err_column,
                         ap_scale_factor=5,
                         color_column="Norm_profile",
                         cmap=["grey", "black", "orange", "red"],
                         norm_method="bins",
                         norm_values=[-0.4, 0.4, 0.85],
                         **kwargs)
        if dms:
            self.ap_scale_factor = 10
            self.set_dms_profile()

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
    def __init__(self, filepath, component, datatype="dancemap",
                 err_column=None, **kwargs):
        self.component = component
        super().__init__(filepath=filepath, datatype=datatype,
                         read_csv_kw={}, err_column=err_column, **kwargs)
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
                         norm_method="none",
                         **kwargs)


class DeltaProfile(Profile):
    def __init__(self, profile1, profile2, column=None, norm_method="min_max",
                 norm_values=[-0.8, 0.8], cmap="coolwarm",
                 ap_scale_factor=None):
        if column is None:
            column = profile1.default_column
        if ap_scale_factor is None:
            ap_scale_factor = profile1.ap_scale_factor * 2
        columns = ["Nucleotide", "Sequence", column]
        profile2.fit_to(profile1)
        profile2 = profile2.get_plotting_dataframe(all_columns=True)
        new_data = profile1.data[columns].merge(
            profile2[columns], how="left", on=["Nucleotide"],
            suffixes=["_1", "_2"])
        new_data.eval(f"Delta_profile = {column}_1 - {column}_2",
                      inplace=True)
        super().__init__(datatype="deltaprofile",
                         column="Delta_profile",
                         sequence=profile1.sequence,
                         ap_scale_factor=ap_scale_factor,
                         dataframe=new_data,
                         cmap=cmap,
                         norm_method=norm_method,
                         norm_values=norm_values,
                         color_column="Delta_profile")


class FragmapProfile(Profile):
    def __init__(self, sample1, sample2, column=None, norm_method="min_max",
                 norm_values=[-0.8, 0.8], cmap="coolwarm",
                 ap_scale_factor=None, filepath=None, sequence=None, **kwargs):

        self.fragmap_data = Fragmapper(
            sample1=sample1, sample2=sample2, **kwargs)

        super().__init__(datatype="fragmap",
                         column="Fragmap_profile",
                         sequence=sample1.data['shapemap'].sequence,
                         ap_scale_factor=ap_scale_factor,
                         dataframe=self.fragmap_data.data,
                         cmap=cmap,
                         norm_method=norm_method,
                         norm_values=norm_values,
                         color_column="Fragmap_profile")
