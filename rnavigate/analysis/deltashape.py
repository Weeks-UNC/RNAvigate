#!/usr/bin/env python
"""
deltaSHAPE software for detecting meaningful changes in SHAPE reactivity
between two samples. Parameters are optimized for detecting in cell vs. cell
free protein protections and enhancements, but useful for identifying any
useful differences.

Copyright Matthew J. Smola 2015
Largely rewritten for RNAvigate by Patrick Irving 2023
"""

###########################################################################
# GPL statement:                                                          #
#                                                                         #
# This program is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
###########################################################################

from rnavigate import data, plots, Sample
from scipy.stats import zscore



class DeltaSHAPE(Sample):
    def __init__(
            self, sample1, sample2, profile='shapemap', smoothing_window=3,
            zf_coeff=1.96, ss_thresh=1, site_window=3, site_nts=2,
            ):
        """Detects meaningful differences in chemical probing reactivity according
        to the deltaSHAPE algorithm (doi:10.1021/acs.biochem.5b00977)

        Args:
            sample1 (rnavigate.Sample): First sample to compare
            sample2 (rnavigate.Sample): Second sample to compare
            smoothing_window (int, optional): Size of windows for data smoothing
                Defaults to 1
            zf_coeff (float, optional): Ajust the Z-factor stringency by changing the equation coefficient.
                Defaults to 1.96
            ss_thresh (int, optional): Set the cutoff threshold of standard score filtering
                Defaults to 1
            find_site (tuple, optional): tuple of integers indicating the window pad size and number of required hits when finding binding sites
                Defaults to (2, 3)
        """

        # STEP ONE
        # extract relevant information from sample 1 and 2 profiles
        # sequence, normalized profile, normalized standard error
        self.parent = sample1
        self.sample = f"{sample1.sample} vs. {sample2.sample}"
        self.data = {}
        columns = ["Nucleotide", "Sequence", "Norm_profile", "Norm_stderr"]
        profile_1 = sample1.get_data(profile)
        profile_1 = profile_1.get_aligned_data(
            data.SequenceAlignment(profile_1, profile_1)
            )
        profile_2 = sample2.get_data(profile)
        profile_2 = profile_2.get_aligned_data(
            data.SequenceAlignment(profile_2, profile_1)
            )
        self.data['profile_1'] = profile_1
        self.data['profile_2'] = profile_2
        df_1 = profile_1.data[columns].copy()
        df_2 = profile_2.data[columns].copy()
        self.data['deltashape'] = DeltaSHAPEProfile(
            input_data=df_1.merge(
                df_2, how="left", on=["Nucleotide", "Sequence"],
                suffixes=("_1", "_2")
                )
            )
        self.calculate_deltaSHAPE(
            smoothing_window=smoothing_window, zf_coeff=zf_coeff,
            ss_thresh=ss_thresh, site_window=site_window, site_nts=site_nts
            )

    def calculate_deltaSHAPE(
            self, smoothing_window=3, zf_coeff=1.96, ss_thresh=1,
            site_window=2, site_nts=3
            ):
        deltashape = self.data['deltashape']
        deltashape.calculate_deltaSHAPE(
            smoothing_window=smoothing_window, zf_coeff=zf_coeff,
            ss_thresh=ss_thresh, site_window=site_window, site_nts=site_nts
            )
        self.data['protections'] = deltashape.get_protections_annotation(
            site_window
            )
        self.data['enhancements'] = deltashape.get_enhancements_annotation(
            site_window
            )

    def plot(self, region='all'):
        plot = plots.Profile(1, self.data['deltashape'].length, region=region)
        plot.plot_data(
            profile=self.data['deltashape'],
            annotations=self.get_data(['protections', 'enhancements']),
            domains=None, plot_error=True, label=self.sample)
        plot.set_figure_size()
        return plot



class DeltaSHAPEProfile(data.Profile):
    def __init__(
            self, input_data, metric='Smooth_diff', metric_defaults=None,
            read_table_kw=None, sequence=None,
            ):
        metric_defaults = {
            'Smooth_diff': {
                'metric_column': 'Smooth_diff',
                'error_column': 'Smooth_diff_stderr',
                'color_column': 'Class',
                'cmap': ['lightgrey', '#7F3B95', '#3EB452'],
                'normalization': 'none',
                'values': None,
                'extend': 'neither',
                'title': 'deltaSHAPE protections\nand enhancements',
                'ticks': [0, 1, 2],
                'tick_labels': ['other', 'protection', 'enhancement'],
                'alpha': 0.7},
        }
        super().__init__(
            input_data, metric, metric_defaults, read_table_kw, sequence
            )

    def calculate_deltaSHAPE(
            self, smoothing_window=3, zf_coeff=1.96, ss_thresh=1,
            site_window=2, site_nts=3):
        # STEP TWO
        # smooth data and errors
        def propagate_errors(errors):
            errors = errors[~errors.isna()]
            return (errors ** 2).sum() ** 0.5 / len(errors)

        self.calculate_windows(
            column='Norm_profile_1', window=smoothing_window,
            new_name='Smooth_profile_1', method='mean', minimum_points=1,
            mask_na=True
            )
        self.calculate_windows(
            column='Norm_stderr_1', window=smoothing_window,
            new_name='Smooth_stderr_1', method=propagate_errors,
            minimum_points=1, mask_na=True
            )
        self.calculate_windows(
            column='Norm_profile_2', window=smoothing_window,
            new_name='Smooth_profile_2', method='mean', minimum_points=1,
            mask_na=True
            )
        self.calculate_windows(
            column='Norm_stderr_2', window=smoothing_window,
            new_name='Smooth_stderr_2', method=propagate_errors,
            minimum_points=1, mask_na=True
            )

        # STEP THREE
        # subtract raw and smoothed data
        self.data['Diff_profile'] = self.data.eval(
            "Norm_profile_1 - Norm_profile_2"
            )
        self.data['Diff_stderr'] = self.data.eval(
            "Norm_stderr_1 + Norm_stderr_2"
            )
        self.data['Smooth_diff'] = self.data.eval(
            "Smooth_profile_1 - Smooth_profile_2"
            )
        self.data['Smooth_diff_stderr'] = self.data.eval(
            'Smooth_stderr_1 + Smooth_stderr_2'
            )
        self.data['Positive'] = self.data["Smooth_diff"] > 0

        # STEP FOUR
        # calculate Z-factors from smoothed data and smoothed errs
        confidence_interval = self.data.eval(
            f'Smooth_stderr_1 + Smooth_stderr_2 * {zf_coeff}')
        difference = self.data['Smooth_diff'].abs()
        self.data['Z_factor'] = 1 - (confidence_interval/difference)

        # STEP FIVE
        # calculate Z-scores from difference of smoothed data
        self.data["Z_score"] = abs(zscore(
            self.data["Smooth_diff"], nan_policy='omit'))

        # STEP SIX
        # identify site_window nt windows where site_nts are significant
        def determine_significance(series):
            df = self.data.loc[series.index, ['Z_factor', 'Z_score']]
            significant = df['Z_factor'] > 0
            significant &= df['Z_score'] > ss_thresh
            return sum(significant) >= site_nts

        windows = self.data['Nucleotide'].rolling(site_window, center=True)
        self.data['Significant'] = windows.apply(determine_significance)
        self.data['Significant'] = self.data['Significant'].fillna(0).astype(bool)
        self.data['Class'] = 0
        self.data.loc[self.data.eval('Significant & ~ Positive'), 'Class'] = 1
        self.data.loc[self.data.eval('Significant & Positive'), 'Class'] = 2

    def get_protections_annotation(self, site_window):
        is_protected = self.data.eval('Significant & ~ Positive')
        return data.Annotation.from_boolean_array(
            values=is_protected, window=site_window, annotation_type='spans',
            name="Protections", sequence=self.sequence, color='#7F3B95',
            )

    def get_enhancements_annotation(self, site_window):
        is_enhanced = self.data.eval('Significant & Positive')
        return data.Annotation.from_boolean_array(
            values=is_enhanced, window=site_window, annotation_type='spans',
            name="Enhancements", sequence=self.sequence, color='#3EB452',
            )
