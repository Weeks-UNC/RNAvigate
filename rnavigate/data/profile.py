from os.path import isfile
import pandas as pd
import numpy as np
from rnavigate import data


class Profile(data.Data):
    def __init__(self, input_data, metric='default', metric_defaults=None,
                 read_table_kw=None, sequence=None):
        if metric_defaults is None:
            metric_defaults = {}
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw)

    @property
    def recreation_kwargs(self):
        return {}

    def get_aligned_data(self, alignment):
        dataframe = alignment.map_nucleotide_dataframe(self.data)
        return self.__class__(
            input_data=dataframe,
            metric=self._metric,
            metric_defaults=self.metric_defaults,
            sequence=alignment.target_sequence,
            **self.recreation_kwargs)

    def get_plotting_dataframe(self):
        new_names = ["Nucleotide"]
        old_names = ["Nucleotide"]
        old_names.append(self.metric)
        new_names.append("Values")
        if self.error_column is not None:
            old_names.append(self.error_column)
            new_names.append("Errors")
        plotting_dataframe = self.data[old_names].copy()
        plotting_dataframe.columns = new_names
        plotting_dataframe["Colors"] = self.colors
        return plotting_dataframe


class SHAPEMaP(Profile):
    def __init__(self, input_data, dms=False, read_table_kw=None,
                 sequence=None, metric='Norm_profile', metric_defaults=None):

        if metric_defaults is None:
            metric_defaults = {}

        metric_defaults = {
            'Norm_profile': {
                'metric_column': 'Norm_profile',
                'error_column': 'Norm_stderr',
                'cmap': ["grey", "black", "orange", "red", "red"],
                'normalization': "bins",
                'values': [-0.4, 0.4, 0.85, 2],
                'title': 'SHAPE Reactivity',
                'extend': 'both'}
            } | metric_defaults
        if (isinstance(input_data, str)
                and input_data.endswith(".map")
                and read_table_kw is None):
            read_table_kw = {
                "names": ["Nucleotide", "Norm_profile",
                          "Norm_stderr", "Sequence"],
                "na_values": "-999"}

        super().__init__(input_data=input_data,
                         read_table_kw=read_table_kw,
                         sequence=sequence,
                         metric=metric,
                         metric_defaults=metric_defaults)
        if dms:
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

        def get_norm_factors(profile):
            finite_data = profile[np.isfinite(profile)]
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
    def __init__(self, input_data, component, read_table_kw=None,
                 sequence=None, metric='Norm_profile', metric_defaults=None):
        self.component = component
        super().__init__(input_data=input_data,
                         read_table_kw=read_table_kw,
                         sequence=sequence,
                         metric=metric,
                         metric_defaults=metric_defaults)

    @property
    def recreation_kwargs(self):
        return {'component': self.component}


    def read_file(self, input_data, read_table_kw={}):
        # parse header
        self.filepath = input_data
        with open(self.filepath) as f:
            header1 = f.readline()
            header2 = f.readline()
        self.header = header1 + header2
        self.components = int(header1.strip().split()[0])
        self.percents = [float(x) for x in header2.strip().split()[1:]]
        self.percent = self.percents[self.component]
        # parse datatable
        read_table_kw['names'] = ["Nucleotide", "Sequence", "Norm_profile",
                                  "Modified_rate", "Untreated_rate"]
        col_offset = 3 * self.component
        bg_col = 3 * self.components + 2
        read_table_kw["usecols"] = [0, 1, 2+col_offset, 3+col_offset, bg_col]
        df = pd.read_table(self.filepath, header=2, **read_table_kw)
        # some rows have an "i" added to the final column
        stripped = []
        for x in df["Untreated_rate"]:
            if type(x) is float:
                stripped.append(x)
            else:
                stripped.append(float(x.rstrip(' i')))
        df["Untreated_rate"] = stripped
        df = df.eval("Reactivity_profile = Modified_rate - Untreated_rate")
        return df


class RNPMaP(Profile):
    def __init__(self, input_data, read_table_kw=None, sequence=None,
                 metric="NormedP", metric_defaults=None):
        if metric_defaults is None:
            metric_defaults = {}
        if read_table_kw is None:
            read_table_kw = {}
        metric_defaults = {
            'NormedP': {
                'metric_column': 'NormedP',
                'color_column': 'RNPsite',
                'cmap': ["silver", "limegreen"],
                'normalization': "none",
                'values': None}
            } | metric_defaults
        read_table_kw = {
            'sep': ','
            } | read_table_kw
        super().__init__(input_data=input_data,
                         read_table_kw=read_table_kw,
                         sequence=sequence,
                         metric=metric,
                         metric_defaults=metric_defaults)


class DeltaProfile(Profile):
    def __init__(self, profile1, profile2, metric=None, metric_defaults=None):
        if metric is None:
            metric = profile1.metric

        columns = ["Nucleotide", "Sequence", metric]
        alignment = data.SequenceAlignment(profile2, profile1)
        profile2 = profile2.get_aligned_data(alignment)
        new_data = pd.merge(
            profile1.data[columns],
            profile2.data[columns],
            how="left", on=["Nucleotide"], suffixes=["_1", "_2"])
        new_data.eval(f"Delta_profile = {metric}_1 - {metric}_2", inplace=True)
        metric_defaults = {
            'Delta_profile': {
                'metric_column': 'Delta_profile',
                'error_column': None,
                'color_column': None,
                'cmap': 'coolwarm',
                'normalization': "min_max",
                'values': [-0.8, 0.8]}
            } | metric_defaults
        super().__init__(
            input_data=new_data,
            sequence=profile1.sequence,
            metric="Delta_profile",
            metric_defaults=metric_defaults)
