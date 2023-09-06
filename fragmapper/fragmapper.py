#!/usr/bin/env python3

"""Fragmapper analysis tools.

Description:
Fraggmapper is designed to compare reactivity profile differences between
RNAvigate sample objects. The intended application of Fragmapper is to
detect fragment or ligand crosslinking sites in RNA through mutational profiling.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from rnavigate.data import Profile, Annotation
from rnavigate.rnavigate import Sample

__version__ = '0.1.0'
__author__ = 'Seth D. Veenbaas'
__maintainer__ = 'Seth D. Veenbaas'
__email__ = 'sethv@live.unc.edu'



class FragMaP(Profile):
    def __init__(self, profile1, profile2, parameters,
                 column='Fragmap_profile',
                 norm_method='min_max', norm_values=[-3, 3],
                 cmap='coolwarm', ap_scale_factor=20):
        dataframe = self.get_dataframe(profile1, profile2, **parameters)
        super().__init__(datatype='fragmap',
                         column=column,
                         sequence=profile1.sequence,
                         ap_scale_factor=ap_scale_factor,
                         dataframe=dataframe,
                         cmap=cmap,
                         norm_method=norm_method,
                         norm_values=norm_values,
                         color_column=column)

    def get_dataframe(self, profile1, profile2, mutation_rate_threshold,
                      depth_threshold, p_significant, ss_threshold,
                      correction_method, by_nucleotide):
        columns = ['Nucleotide', 'Sequence', 'Modified_mutations',
                        'Modified_effective_depth', 'Modified_rate', 'Std_err']
        dataframe = pd.merge(profile1.data[columns], profile2.data[columns],
            how='left', on=['Nucleotide', 'Sequence'], suffixes=('_1', '_2'))

        # Filter data
        dataframe['filter'] = True
        # Effective read depths must be greater than threshold
        dataframe['filter'] &= dataframe['Modified_effective_depth_1'] > depth_threshold
        dataframe['filter'] &= dataframe['Modified_effective_depth_2'] > depth_threshold
        # Control mutation rate (profile 2) must be less than threshold
        dataframe['filter'] &= dataframe['Modified_rate_2'] < mutation_rate_threshold
        valid = dataframe['filter']

        # Calculate Z-scores for profile 1 and profile 2
        dataframe[['zscore_1', 'zscore_2']] = np.nan
        # Z-scores per nucleotide
        if by_nucleotide and len(dataframe) > 200:
            for nt in 'AUCG':
                valid_nt = valid  & (dataframe['Sequence'] == nt)
                dataframe.loc[valid_nt, 'zscore_1'] = stats.zscore(
                    dataframe.loc[valid_nt, 'Modified_rate_1'])
                dataframe.loc[valid_nt, 'zscore_2'] = stats.zscore(
                    dataframe.loc[valid_nt, 'Modified_rate_2'])
        # Z-scores per nucleotide group: A/C and U/G
        elif by_nucleotide:
            for nts in [['A', 'C'], ['U', 'G']]:
                valid_nts = valid & dataframe['sequence'].isin(nts)
                dataframe.loc[valid_nts, 'zscore_1'] = stats.zscore(
                    dataframe.loc[valid_nts, 'Modified_rate_1'])
                dataframe.loc[valid_nts, 'zscore_2'] = stats.zscore(
                    dataframe.loc[valid_nts, 'Modified_rate_2'])
        # Z-scores for all nucleotides
        else:
            dataframe.loc[valid, 'zscore_1'] = stats.zscore(
                dataframe.loc[valid, 'Modified_rate_1'])
            dataframe.loc[valid, 'zscore_2'] = stats.zscore(
                dataframe.loc[valid, 'Modified_rate_2'])

        # Frag-MaP profile is the difference in Z-scores
        dataframe.eval('Fragmap_profile = zscore_1 - zscore_2', inplace=True)
        dataframe['Fragmap_pvalue'] = stats.norm.sf(
            dataframe['Fragmap_profile'])

        # Calculate T-statistic and P-value
        dataframe['tstat'], dataframe['pvalue'] = stats.ttest_ind_from_stats(
            mean1=dataframe['Modified_rate_1'],
            std1=dataframe['Std_err_1'],
            nobs1=dataframe['Modified_effective_depth_1'],
            mean2=dataframe['Modified_rate_2'],
            std2=dataframe['Std_err_2'],
            nobs2=dataframe['Modified_effective_depth_2'],
            equal_var=False)

        # Bonferoni p-value correction
        p_significant = p_significant/len(dataframe)
        dataframe['Significant'] = dataframe['pvalue'] < p_significant

        # Calculate Delta rate
        dataframe.eval('Delta_rate = Modified_rate_1 - Modified_rate_2',
                       inplace=True)
        dataframe['Delta_rate_zscore'] = stats.zscore(dataframe['Delta_rate'])

        # Perform False Discovery Rate (FDR) correction
        if correction_method is not None:
            p_values = dataframe['Fragmap_pvalue'].values
            reject, pvals_corrected, _, _ = multipletests(
                p_values, alpha=ss_threshold, method=correction_method)
            dataframe['Reject'] = reject
            dataframe['Fragmap_pvalue_corrected'] = pvals_corrected
        else:
            dataframe['Reject'] = dataframe['Fragmap_pvalue'] < ss_threshold
            dataframe['Fragmap_pvalue_corrected'] = np.nan

        # Called Frag-MaP sites
        dataframe['Site'] = (dataframe['Reject']) & (dataframe['Significant'])

        return dataframe

    def get_annotation(self):
        return Annotation(
            'Frag-MaP sites',
            sequence=self.sequence,
            sites=self.data.loc[self.data['Site'], "Nucleotide"].to_list(),
            color="green")


class Fragmapper(Sample):
    def __init__(self, sample1, sample2, parameters=None, profile='shapemap'):
        if sample1.data[profile].sequence != sample2.data[profile].sequence:
            raise ValueError('Profiles must have the same sequence')
        super().__init__(f'FragMaP: {sample1.sample}/{sample2.sample}')
        self.parent = sample1
        self.sample1 = sample1
        self.sample2 = sample2
        if parameters is None:
            parameters = {}
        self.parameters = {
            'mutation_rate_threshold': 0.025,
            'depth_threshold': 5000,
            'p_significant': 0.01,
            'ss_threshold': 0.05,
            'correction_method': 'fdr_bh',
            'by_nucleotide': True
        }
        self.parameters |= parameters

        self.data['fragmap'] = FragMaP(sample1.get_data(profile),
                                       sample2.get_data(profile),
                                       parameters=self.parameters)
        self.data['fragmap_sites'] = self.data['fragmap'].get_annotation()


    def plot_scatter(self, column='Modified_rate'):
        """Generates scatter plots useful for fragmapper quality control.

        Args:
            column (str, optional): Dataframe column containing data to plot 
                                    (must be avalible for the sample and control).
                                    Defaults to 'Modified_rate'.

        Returns:
            scatter_plot: Scatter plot with control values on the x-axis,
                          sample values on the y-axis, and each point
                          representing a nucleotide not filtered out in the
                          fragmapper pipeline.
        """
        fig, ax = plt.subplots(1, 1, figsize=(7, 7))

        columns = [f'{column}_1', f'{column}_2']
        scatter_data = self.data["fragmap"].data
        scatter_data = scatter_data[scatter_data['filter']].copy()
        sites = scatter_data[scatter_data['Site']]
        non_sites = scatter_data[~scatter_data['Site']]

        for _, (nt, y, x) in sites[['Nucleotide']+columns].iterrows():
            ax.text(x, y, int(nt), fontsize=12)

        ax.scatter(non_sites[columns[1]], non_sites[columns[0]], s=5)
        ax.scatter(sites[columns[1]], sites[columns[0]], s=25)

        ax.set_xlabel(f'{self.sample2.sample} {column}', fontsize=14)
        ax.set_ylabel(f'{self.sample1.sample} {column}', fontsize=14)

        plt.legend(['Uncalled Nucleotides', 'Frag-MaP Sites'], fontsize=12)

        return fig, ax
