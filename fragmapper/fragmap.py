#!/usr/bin/env python3

"""Fragmapper analysis tools.

Description:
Fraggmapper is designed to compare reactivity profile differences between
RNAvigate sample objects. The intended application of Fragmapper is to
detect fragment or ligand crosslinking sites in RNA through mutational profiling.
"""

from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from .annotation import Annotation

__version__ = '0.0.3'
__author__ = 'Seth D. Veenbaas'
__maintainer__ = 'Seth D. Veenbaas'
__email__ = 'sethv@live.unc.edu'


@dataclass()
class Fragmapper:
    sample1: object = field(repr=True)
    sample2: object = field(repr=True)
    mutation_rate_thresh: float = field(repr=False, default=0.025)
    depth_thresh: int = field(repr=False, default=5000)
    seq_source: str = field(repr=False, default='shapemap')

    p_significant: float = field(repr=False, default=0.01)
    ss_thresh: float = field(repr=False, default=0.05)
    correction: bool = field(repr=False, default=False)
    correction_method: str = field(repr=False, default='fdr_bh')

    by_seq: bool = field(repr=False, default=True)

    _run_post_init: bool = field(default=True, repr=False, init=False)

    def __post_init__(self) -> None:

        Fragmapper._get_data(self)

        Fragmapper._apply_filters(self)

        Fragmapper._calc_profile(self)

        Fragmapper._get_FDR_correction(self)

        Fragmapper._calc_significance(self)

        Fragmapper._calc_delta_rate(self)

        Fragmapper._get_fragmap_sites(self)

        Fragmapper._set_fragmap_sites_annotation(self)

        Fragmapper._set_fragmap_sites_data(self)

    def _get_data(self):
        self.columns = ["Nucleotide", "Sequence", "Modified_mutations",
                        "Modified_effective_depth", "Modified_rate"]
        profile_1_data = self.sample1.data['shapemap'].data[self.columns].copy(
        )
        profile_1_data.eval(
            'Stderr = (Modified_rate*(1-Modified_rate)/Modified_effective_depth)**0.5', inplace=True)
        profile_2_data = self.sample2.data['shapemap'].data[self.columns].copy(
        )
        profile_2_data.eval(
            'Stderr = (Modified_rate*(1-Modified_rate)/Modified_effective_depth)**0.5', inplace=True)

        assert all(profile_1_data["Sequence"] == profile_2_data["Sequence"]
                   ), "Sequences don't match"
        self.data = profile_1_data.merge(profile_2_data, how="left",
                                         on=["Nucleotide", "Sequence"],
                                         suffixes=("_1", "_2"))

    def _apply_filters(self):
        # Effective read depths > depth_thresh
        self.data.loc[(self.data['Modified_effective_depth_1'] <
                       self.depth_thresh), ['Modified_rate_1', 'Modified_rate_2']] = np.nan
        self.data.loc[(self.data['Modified_effective_depth_2'] <
                       self.depth_thresh), ['Modified_rate_1', 'Modified_rate_2']] = np.nan

        # Sample_2 (control) mutation rate < mutation rate threshold
        self.data.loc[(self.data['Modified_rate_2'] >
                       self.mutation_rate_thresh), ['Modified_rate_1', 'Modified_rate_2']] = np.nan

    def _calc_profile(self):

        if self.by_seq and len(self.sample1.data['shapemap'].data['Sequence']) > 200:
            self.data['zscore_1'] = np.nan
            self.data.loc[self.data['Sequence'].isin(['A']), 'zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'].loc[self.data['Sequence'].isin(['A'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['C']), 'zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'].loc[self.data['Sequence'].isin(['C'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['U']), 'zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'].loc[self.data['Sequence'].isin(['U'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['G']), 'zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'].loc[self.data['Sequence'].isin(['G'])], nan_policy='omit')

            self.data['zscore_2'] = np.nan
            self.data.loc[self.data['Sequence'].isin(['A']), 'zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'].loc[self.data['Sequence'].isin(['A'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['C']), 'zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'].loc[self.data['Sequence'].isin(['C'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['U']), 'zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'].loc[self.data['Sequence'].isin(['U'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['G']), 'zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'].loc[self.data['Sequence'].isin(['G'])], nan_policy='omit')

        elif self.by_seq:
            self.data['zscore_1'] = np.nan
            self.data.loc[self.data['Sequence'].isin(['A', 'C']), 'zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'].loc[self.data['Sequence'].isin(['A', 'C'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['U', 'G']), 'zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'].loc[self.data['Sequence'].isin(['U', 'G'])], nan_policy='omit')
            self.data['zscore_2'] = np.nan
            self.data.loc[self.data['Sequence'].isin(['A', 'C']), 'zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'].loc[self.data['Sequence'].isin(['A', 'C'])], nan_policy='omit')
            self.data.loc[self.data['Sequence'].isin(['U', 'G']), 'zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'].loc[self.data['Sequence'].isin(['U', 'G'])], nan_policy='omit')

        else:
            self.data['zscore_1'] = stats.zscore(
                self.data['Modified_rate_1'], nan_policy='omit')
            self.data['zscore_2'] = stats.zscore(
                self.data['Modified_rate_2'], nan_policy='omit')

        # self.data.loc[self.data['Sequence'].isin(['U', 'G']), 'Fragmap_profile'] = self.data['zscore_1'] - self.data['zscore_2']
        self.data['Fragmap_profile'] = self.data['zscore_1'] - \
            self.data['zscore_2']

        self.data['Fragmap_profile_no_nan'] = self.data['Fragmap_profile'].fillna(
            0)
        self.data['Fragmap_pvalue'] = stats.norm.sf(
            self.data['Fragmap_profile'])

    def _get_FDR_correction(self, column_name: str = 'Fragmap_pvalue', corrected_column_name: str = None):
        """Perform False Discovery Rate (FDR) correction on a specified column of a pandas DataFrame.

        Args:
            df (pandas.DataFrame): Input DataFrame.
            column_name (str): Name of the column containing the p-values.
            alpha (float): The desired FDR threshold (default: 0.05).
            method (str): The FDR correction method to use.
                        Options: 'benjamini-hochberg' (default), 'bonferroni', 'sidak', and more.
            corrected_column_name (str): Name for the column to store the corrected p-values.
                                        If None, the original column will be overwritten (default: None).

        Returns:
            None
        """
        if self.correction:
            p_values = self.data[column_name].values
            reject, pvals_corrected, _, _ = multipletests(
                p_values, alpha=self.ss_thresh, method=self.correction_method)
            if corrected_column_name is None:
                corrected_column_name = 'Fragmap_pvalue_corrected'
            self.data['Reject'] = reject
            self.data[corrected_column_name] = pvals_corrected

        else:
            self.data['Reject'] = False
            self.data.loc[self.data['Fragmap_pvalue']
                          < self.ss_thresh, 'Reject'] = True
            self.data['Fragmap_pvalue_corrected'] = np.nan

    def _calc_significance(self):
        self.data['tstat'], self.data['pvalue'] = stats.ttest_ind_from_stats(mean1=self.data["Modified_rate_1"], std1=self.data["Stderr_1"], nobs1=self.data[
                                                                             "Modified_effective_depth_1"], mean2=self.data["Modified_rate_2"], std2=self.data["Stderr_2"], nobs2=self.data["Modified_effective_depth_2"], equal_var=False)
        # Bonferoni correction
        self.p_significant = self.p_significant/len(self.data)
        self.data['Significant'] = False
        self.data.loc[self.data['pvalue'] <
                      self.p_significant, 'Significant'] = True

    def _calc_delta_rate(self):
        self.data.eval(
            'Delta_rate = Modified_rate_1 - Modified_rate_2', inplace=True)
        self.data['Delta_rate_zscore'] = stats.zscore(
            self.data['Delta_rate'], nan_policy='omit')

    def _get_fragmap_sites(self):
        # self.zscore_thresh = stats.norm.ppf(1 - self.ss_thresh)
        # self.data['Site'] = False
        # self.data.loc[(self.data['Fragmap_profile'] > self.zscore_thresh) &
        #               (self.data['Significant'] == True), 'Site'] = True

        self.data['Site'] = False
        self.data.loc[(self.data['Reject'] == True) &
                      (self.data['Significant'] == True), 'Site'] = True

        self.fragmap_sites_data = self.data.loc[self.data['Site'] == True].copy(
        )
        fragment_name = self.sample1.sample.split('_')[0]
        self.fragmap_sites_data['Fragment'] = fragment_name

        self.fragmap_sites_data = self.fragmap_sites_data[[
            'Fragment', 'Nucleotide', 'Sequence',
            'Fragmap_profile', 'Fragmap_pvalue', 'Reject', 'Fragmap_pvalue_corrected',
            'tstat', 'pvalue',
            'Delta_rate', 'Delta_rate_zscore',
            'zscore_1', 'zscore_2',
            'Modified_mutations_1',
            'Modified_effective_depth_1',
            'Modified_rate_1', 'Stderr_1',
            'Modified_mutations_2',
            'Modified_effective_depth_2',
            'Modified_rate_2', 'Stderr_2'
        ]]

        self.fragmap_sites = self.fragmap_sites_data['Nucleotide'].to_list()

        self.pymol_sites = '+'.join([str(nt) for nt in self.fragmap_sites])

    def _set_fragmap_sites_annotation(self):

        self.sample1.set_data(
            name="fragmap_sites",
            seq_source=self.seq_source,
            instantiator=Annotation,
            sites=self.fragmap_sites,
            filepath=None,
            color='green')

    def _set_fragmap_sites_data(self):

        # self.sample1.data['fragmap_sites'] = lambda: None
        setattr(
            self.sample1.data['fragmap_sites'], 'data', self.fragmap_sites_data)
        setattr(
            self.sample1.data['fragmap_sites'], 'list', self.fragmap_sites)
        setattr(
            self.sample1.data['fragmap_sites'], 'pymol_string', self.pymol_sites)
        setattr(
            self.sample1.data['fragmap_sites'], 'self', self)

    def plot_scatter(self, column='Modified_rate', xmax=None, ymax=None, xmin=None, ymin=None) -> 'scatter_plot':
        """Generates scatter plots useful for fragmapper quality control.

        Args:
            column (str, optional): Dataframe column containing data to plot 
                                    (must be avalible for the sample and control).
                                    Defaults to 'Modified_rate'.
            x_max (_type_, optional): Maximun values for the x-axis. Defaults to None.
            y_max (_type_, optional): Maximum value for the y-axis. Defaults to None.

        Returns:
            scatter_plot: Scatter plot with control values on the x-axis,
                          sample values on the y-axis, and each point
                          representing a nucleotide not filtered out in the
                          fragmapper pipeline.
        """
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        fig, ax = plt.subplots(1, 1, figsize=(7, 7))

        scatter_data = self.data.copy()

        scatter_data.loc[(self.data['Modified_effective_depth_1'] > self.depth_thresh) & (
            scatter_data['Modified_effective_depth_2'] > self.depth_thresh), 'Region'] = True

        scatter_data.loc[scatter_data['Site'] == True, 'Positive_mask'] = True

        # scatter_data.loc[(scatter_data['Fragmap_profile'] < -4)
        #                  & (scatter_data['Significant'] == True), 'Negative_mask'] = True

        pos_sites = scatter_data.loc[scatter_data['Positive_mask']
                                     == True, 'Nucleotide'].tolist()
        # neg_sites = scatter_data.loc[scatter_data['Negative_mask']
        #                              == True, 'Nucleotide'].tolist()

        # pos_sites = positive_mask.loc[positive_mask == True].index.tolist()
        # neg_sites = negative_mask.loc[negative_mask == True].index.tolist()
        # sites = pos_sites + neg_sites
        sites = pos_sites

        x_max = self.data[f'{column}_2'].max()*1.1
        x_min = self.data[f'{column}_2'].min()*1.1
        
        y_max = self.data[f'{column}_1'].max()*1.1
        y_min = self.data[f'{column}_1'].min()*1.1

        for site in sites:
            x = scatter_data.loc[scatter_data['Nucleotide']
                                 == site, f'{column}_2'].values
            y = scatter_data.loc[scatter_data['Nucleotide']
                                 == site, f'{column}_1'].values
            # x = x[0]
            # y = y[0]

            ax.text(x, y, site, fontsize=12)

            # if x*1.1 > x_max:
            #     x_max = x*1.1
            # if y*1.1 > y_max:
            #     y_max = y*1.1

        ax.scatter(scatter_data.loc[scatter_data['Region'] == True, f'{column}_2'].values,
                   scatter_data.loc[scatter_data['Region'] == True, f'{column}_1'].values, s=5)
        # ax.scatter(scatter_data.loc[scatter_data['Negative_mask'] == True, f'{column}_2'].values,
        #            scatter_data.loc[scatter_data['Negative_mask'] == True, f'{column}_1'].values, s=15)
        ax.scatter(scatter_data.loc[scatter_data['Positive_mask'] == True, f'{column}_2'].values,
                   scatter_data.loc[scatter_data['Positive_mask'] == True, f'{column}_1'].values, s=25)

        # square_size = max([x_max, y_max])

        if xmax is not None:
            x_max = xmax
        if xmax is not None:
            x_min = xmin
        if ymax is not None:
            y_max = ymax
        if ymax is not None:
            y_min = ymin

        ax.set(xlim=[x_min, x_max],
               ylim=[y_min,  y_max])

        sample_name = self.sample1.sample.split('_')[0]
        control_name = self.sample2.sample.split('_')[0]

        ax.set_xlabel(f'{control_name} {column}', fontsize=14)
        ax.set_ylabel(f'{sample_name} {column}', fontsize=14)
        # ax.set_aspect('equal')

        # plt.legend(['Uncalled Nucleotides', 'High Methyl Sites',
        #             'Frag-MaP Sites'], fontsize=10)
        plt.legend(['Uncalled Nucleotides',
                    'Frag-MaP Sites'], fontsize=12)
