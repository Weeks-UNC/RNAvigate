"""LogCompare compares reactivity profiles for significant differences.

This analysis requires replicates.
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from rnavigate import Sample, data


# TODO: refactor as a subclass of rnavigate.Sample
class LogCompare(Sample):
    """Compares 2 experimental samples, given replicates of each sample.

    Algorithm
    ---------
    1. Calculate the ln(modified/untreated) rate for each replicate.
    2. Scale these values to minimize the median of the absolute difference
    between samples.
    3. Calculate the standard error in these values for each replicate.
    4. Calculate the difference between samples.
    5. Calculate z-scores between samples.
    6. Plot the results in two panels:
    (1) the scaled log10(modified/untreated) rate for each sample with error bars, and
    (2) the difference between samples, colored by z-score.

    Methods
    -------
    __init__: computes log10(modified/untreated) rates, rescales the data, then calls make_plot()
    get_profile_sequence: gets log10(m/u) rate and sequence from sample
    rescale: rescales a profile to minimize difference to another profile
    load_replicates: calculates average and standard error of replicates
    make_plots: displays the two panels described above.

    Attributes
    ----------
    data : str
        a key of sample.data to retrieve per-nucleotide data
    groups : dict
        a dictionary with keys 1 and 2, each containing: self.data (averaged scaled log10(m/u)),
        "stderr" (standard errors), "stacked" (2d array of scaled log10(m/u) per replicate),
        "seq" (the sequence string)
    """

    def __init__(
        self,
        samples1,
        samples2,
        name1,
        name2,
        profile_kw,
        sequence=None,
        inherit=None,
    ):
        """Takes replicates of two samples for comparison. Replicates are
        required. Calculates the log division profile
        (log10(modified/untreated)) and minimizes the median of the absolute
        difference between these. Finally, creates a plot.

        Args:
            samples1 (list of Sample objects): Replicates of the first sample
            samples2 (list of Sample objects): Replicates of the second sample
            name1 (string): name of first sample
            name2 (string): name of second sample
            profile (str, optional): Datatype to compare. Defaults to "profile".
        """
        super().__init__(
            sample=f"{name1} vs {name2} - logcompare",
            inherit=inherit,
        )
        if sequence is None:
            sequence = samples1[0].get_data(profile_kw)
        self.set_data(data_keyword="sequence", inputs=sequence)
        sequence = self.get_data("sequence")
        profiles1 = []
        for i, sample in enumerate(samples1):
            profile = sample.get_data(profile_kw)
            profile = profile.get_aligned_data(
                data.SequenceAlignment(profile, sequence)
            )
            self.set_data(data_keyword=f"Sample1_{i + 1}", inputs=profile)
            profiles1.append(profile)
        profiles2 = []
        for i, sample in enumerate(samples2):
            profile = sample.get_data(profile_kw).get_aligned_data(
                data.SequenceAlignment(profile, sequence)
            )
            self.set_data(data_keyword=f"Sample2_{i + 1}", inputs=profile)
            profiles2.append(profile)
        self.set_data(
            data_keyword="log_compare",
            inputs=LogProfile(input_data=[profiles1, profiles2]),
        )


class LogProfile(data.Profile):
    """A class for log10(Modified_rate/Untreated_rate) profiles."""

    def __init__(
        self,
        input_data,
        metric="mean_diff",
        metric_defaults=None,
        sequence=None,
        **kwargs,
    ):
        """Create a new LogProfile object.

        Parameters
        ----------
        input_data (list of rnav.data.Profile)
            list of profiles
        """
        if not isinstance(input_data, pd.DataFrame):
            dfs = []
            for sample in input_data:
                df = sample[0].data[["Nucleotide", "Sequence"]].copy()
                cols = [f"profile{i + 1}" for i in range(len(sample))] + [
                    "avg",
                    "stderr",
                ]
                df[cols] = np.vstack(self.load_replicates(sample)).T
                dfs.append(df)
            input_data = pd.merge(
                left=dfs[0],
                right=dfs[1],
                how="outer",
                on=["Nucleotide", "Sequence"],
                suffixes=("_1", "_2"),
            )
            stack1 = input_data["avg_1"]
            stack2 = input_data["avg_2"]
            # calculate mean difference, standard error, and z-scores
            stack1 = stack1 + 10
            stack1 = self.rescale(stack1, stack2)
            meandiff = stack2 - stack1
            std_err = input_data["stderr_1"] + input_data["stderr_2"]
            z_scores = meandiff / std_err
            input_data["avg_1"] = stack1
            input_data["mean_diff"] = meandiff
            input_data["std_err"] = std_err
            input_data["z_scores"] = z_scores

        if metric_defaults is None:
            metric_defaults = {
                f"avg_{i}": {
                    "metric_column": f"avg_{i}",
                    "error_column": f"stderr_{i}",
                    "color_column": None,
                    "cmap": "viridis",
                    "normalization": "min_max",
                    "values": [-1, 3],
                    "extend": "neither",
                    "title": "ln(Modified/Untreated)",
                    "alpha": 0.7,
                }
                for i in [1, 2]
            } | {
                "mean_diff": {
                    "metric_column": "mean_diff",
                    "error_column": "std_err",
                    "color_column": None,
                    "cmap": "bwr",
                    "normalization": "min_max",
                    "values": [-10, 10],
                    "extend": "neither",
                    "title": "Standard score of difference",
                    "alpha": 1.0,
                },
            }
        super().__init__(
            input_data=input_data,
            metric=metric,
            metric_defaults=metric_defaults,
        )

    def calc_profile(self, profile):
        """Calculate log10(Modified_rate/Untreated_rate) for the given sample/profile.

        Args:
            sample (rnavigate.Sample): an rnavigate sample

        Returns:
            np.array: log profile
        """
        df = profile.data
        sequence = profile.sequence
        plus = df.Modified_rate.values.copy()
        minus = df.Untreated_rate.values.copy()
        nans = np.full(len(sequence), np.nan)
        profile = np.divide(plus, minus, out=nans, where=(0 < minus) & (minus < 0.05))
        profile[profile == 0] = np.nan
        profile[minus > 0.05] = np.nan
        profile[sequence.islower()] = np.nan
        profile = np.log(profile)
        return profile

    def rescale(self, profile, target_profile):
        """scales profile to minimize difference to target_profile.

        Args:
            profile (np.array): log10 profile to scale
            target_profile (np.array): 2nd log10 profile

        Returns:
            np.array: scaled profile
        """

        def f(offset):
            return np.nanmedian(np.abs(profile + offset - target_profile))

        result = minimize_scalar(f, bounds=[-20, 20])
        offset = result.x
        return profile + offset

    def load_replicates(self, profiles):
        """calculates log profiles, avg and sterr for a group of replicates.

        Args:
            *profiles (list of rnavigate.Sample): replicates to load
        """
        profiles = [self.calc_profile(profile) for profile in profiles]
        rescaled_profiles = []
        for profile in profiles[1:]:
            rescaled_profiles.append(self.rescale(profile, profiles[0]))
        profiles = [profiles[0]] + rescaled_profiles

        stacked = np.vstack(profiles)
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "Mean of empty slice")
            avgprofile = np.nanmean(stacked, axis=0)
        stderr = np.std(stacked, axis=0)
        return profiles + [avgprofile, stderr]
