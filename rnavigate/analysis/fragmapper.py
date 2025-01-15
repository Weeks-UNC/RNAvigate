#!/usr/bin/env python3

"""Fragmapper analysis tools.

Description:
FragMapper compares reactivity profile differences between SHAPE-MaP profiles.
The intended application of Fragmapper is to detect fragment or ligand
crosslinking sites in RNA.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import rnavigate as rnav
from rnavigate import data, Sample
from textwrap import wrap

__version__ = "0.2.0"
__author__ = "Seth D. Veenbaas"
__maintainer__ = "Seth D. Veenbaas"
__email__ = "sethv@live.unc.edu"


class FragMaP(data.Profile):
    def __init__(
        self,
        input_data,
        parameters,
        metric="Delta_zscore",
        metric_defaults=None,
        read_table_kw=None,
        sequence=None,
        name=None,
    ):
        self.parameters = parameters
        if isinstance(input_data, (list, tuple)) and len(input_data) == 2:
            profile1, profile2 = input_data
            input_data = self.get_dataframe(profile1, profile2, **parameters)
        elif not isinstance(input_data, pd.DataFrame):
            raise ValueError(
                "Input data must be a list of 2 shapemap profiles or a "
                "pre-initialized dataframe"
            )
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "Delta_zscore": {
                "metric_column": "Delta_zscore",
                "error_column": "Delta_zscore_err",
                "color_column": "Site",
                # "cmap": ["C0", "C1"],
                "cmap": ["#9FB3B3", "#E86500"],
                "normalization": "none",
                "values": None,
                "extend": "neither",
                "alpha": 1.0,
                "title": "Frag-MaP sites",
                "ticks": [0, 1],
                "tick_labels": ["other", "sites"],
            }
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            sequence=sequence,
            name=name,
        )

    @property
    def recreation_kwargs(self):
        return {"parameters": self.parameters}
    
    def get_dataframe(
        self,
        profile1,
        profile2,
        mutation_rate_threshold,
        depth_threshold,
        delta_rate_threshold,
        zscore_threshold,
        zscore_min_threshold
    ):
        columns = [
            "Nucleotide",
            "Sequence",
            "Modified_mutations",
            "Modified_effective_depth",
            "Modified_rate",
            "Std_err",
        ]
        dataframe = pd.merge(
            profile1.data[columns],
            profile2.data[columns],
            how="left",
            on=["Nucleotide", "Sequence"],
            suffixes=("_1", "_2"),
        )

        # Filter data
        dataframe["Valid"] = dataframe.eval(
            # Effective read depths must be greater than threshold
            f"Modified_effective_depth_1 > {depth_threshold} "
            f"& Modified_effective_depth_2 > {depth_threshold} "
            # Control mutation rate (profile 2) must be less than threshold
            f"& Modified_rate_2 < {mutation_rate_threshold}"
        ).astype(bool)
        valid = dataframe["Valid"]

        # Calculate Delta rate
        dataframe["Delta_rate"] = dataframe.eval("Modified_rate_1 - Modified_rate_2")
        dataframe["Delta_std_err"] = dataframe.eval("Std_err_1 + Std_err_2")
        dataframe["Delta_rate_min"] = dataframe.eval("Delta_rate - Delta_std_err")
                     
        # Calculate rolling Delta rate modified zscore
        rolling_median = dataframe["Delta_rate"].rolling(window=50, min_periods=25, center=True).median()
        rolling_mad = dataframe["Delta_rate"].rolling(window=50, min_periods=25, center=True).apply(lambda x: np.median(np.abs(x - np.median(x))), raw=False)
        
        # Calculate the rolling modified z-score
        dataframe["Delta_zscore"] = 0.6745 * (dataframe["Delta_rate"] - rolling_median) / rolling_mad

        # Calculate the minimum rolling modified z-score
        dataframe["Delta_zscore_min"] = 0.6745 * (dataframe["Delta_rate_min"] - rolling_median) / rolling_mad

        # Calculate the rolling modified z-score error
        dataframe["Delta_zscore_err"] = (dataframe["Delta_zscore"] - dataframe["Delta_zscore_min"])
        
        # Define peaks based on z-score amplitude
        dataframe["Site"] = (
            (dataframe["Delta_rate"] > delta_rate_threshold) &
            (dataframe["Delta_zscore"] > zscore_threshold) &
            (dataframe["Delta_zscore_min"] > zscore_min_threshold) &
            (dataframe["Valid"] == True)
        )

        return dataframe


    def get_annotation(self):
        return data.Annotation(
            input_data=self.data.loc[
                self.data["Site"] & self.data["Valid"], "Nucleotide"
            ].to_list(),
            name="Frag-MaP sites",
            sequence=self.sequence,
            annotation_type="sites",
            color="#E86500",
        )


class Fragmapper(Sample):
    def __init__(self, sample1, sample2, parameters=None, profile="shapemap"):
        if sample1.data[profile].sequence != sample2.data[profile].sequence:
            raise ValueError("Profiles must have the same sequence")
        if parameters is None:
            parameters = {}
        self.parameters = {
            "mutation_rate_threshold": 0.025,
            "depth_threshold": 50000,
            "delta_rate_threshold": 0.015,
            "zscore_threshold": 30,
            "zscore_min_threshold": 5,
        }
        self.parameters |= parameters
        self.fragmap = FragMaP(
            input_data=[sample1.get_data(profile), sample2.get_data(profile)],
            parameters=self.parameters,
        )
        super().__init__(
            sample=f"FragMaP: {sample1.sample}/{sample2.sample}",
            inherit=[sample1, sample2],
            fragmap=self.fragmap,
            fragmap_sites=self.fragmap.get_annotation(),
        )
        self.sample1 = sample1
        self.sample2 = sample2
        
    
    def update_annotation(self):
        self.data['fragmap_sites']=self.fragmap.get_annotation()
        

    def plot_scatter(self, column="Modified_rate"):
        """Generates scatter plots useful for fragmapper quality control.

        Args:
            column (str, optional):
                Dataframe column containing data to plot (must be avalible for
                the sample and control).
                Defaults to "Modified_rate".

        Returns:
            (matplotlib figure, matplotlib axis)
                Scatter plot with control values on the x-axis, sample values
                on the y-axis, and each point representing a nucleotide not
                filtered out in the fragmapper pipeline.
        """
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

        columns = [f"{column}_1", f"{column}_2"]
        scatter_data = self.data["fragmap"].data
        scatter_data = scatter_data[scatter_data["Valid"]].copy()
        sites = scatter_data[scatter_data["Site"]]
        non_sites = scatter_data[~scatter_data["Site"]]

        for _, (nt, y, x) in sites[["Nucleotide"] + columns].iterrows():
            ax.text(x, y + 0.002, int(nt), size=6)

        ax.scatter(x=non_sites[columns[1]], y=non_sites[columns[0]], s=5, alpha=0.25)
        ax.scatter(x=sites[columns[1]], y=sites[columns[0]], s=25, alpha=1)

        ax.set_xlabel(f"{self.sample2.sample}: {column}")
        ax.set_ylabel(f"{self.sample1.sample}: {column}")

        plt.legend(["Uncalled Nucleotides", "Frag-MaP Sites"])

        return fig, ax


class FragmapperReplicates(Sample):
    def __init__(
        self,
        samples_1 : list,
        samples_2 : list,
        parameters=None,
        profile="shapemap",
    ):

        # if len(samples_1) != len(samples_2):
        #     raise ValueError("The samples lists must be the same length.")

        for idx_1 in range(len(samples_1)):
            sequence_1 = samples_1[idx_1].data[profile].sequence
            for idx_2 in range(len(samples_2)):
                sequence_2 = samples_2[idx_2].data[profile].sequence
                if sequence_1 != sequence_2:
                    raise ValueError(
                        "Profiles must have the same sequence. "
                        f"{samples_1[idx_1]} sequence != {samples_2[idx_2]} sequence"
                    )

        if parameters is None:
            parameters = {}
        self.parameters = {
            "mutation_rate_threshold": 0.025,
            "depth_threshold": 50000,
            "delta_rate_threshold": 0.015,
            "zscore_threshold": 30,
            "zscore_min_threshold": 5,
        }
        self.parameters |= parameters

        # Process samples_1
        self.merged_1 = self.merge_samples(samples_1, profile)
        avg_merged_1 = self.average_columns(self.merged_1)
        fragmap_rep1 = rnav.Sample(
            sample=f'fragmap_rep1',
            shapemap={'shapemap': avg_merged_1}
        )

        # Process samples_2
        self.merged_2 = self.merge_samples(samples_2, profile)
        avg_merged_2 = self.average_columns(self.merged_2)
        fragmap_rep2 = rnav.Sample(
            sample=f'fragmap_rep2',
            shapemap={'shapemap': avg_merged_2}
        )

        self.fragmap = FragMaP(
            input_data=[
                fragmap_rep1.get_data(profile),
                fragmap_rep2.get_data(profile),
            ],
            parameters=self.parameters,
        )

        super().__init__(
            sample=f"FragMaP: {samples_1[0].sample}/{samples_2[0].sample}",
            inherit=samples_1 + samples_2,
            fragmap=self.fragmap,
            fragmap_sites=self.fragmap.get_annotation(),
        )

        self.sample1 = samples_1[0]
        self.sample2 = samples_2[0]


    def update_annotation(self):
        self.data['fragmap_sites']=self.fragmap.get_annotation()


    def merge_samples(
        self,
        samples : list,
        profile : str = 'shapemap',
        suffix : str = "rep",
        columns : list = [
            "Nucleotide",
            "Sequence",
            "Modified_mutations",
            "Modified_effective_depth",
            "Modified_rate", 
        ],
        exceptions : list = [
            'Nucleotide',
            'Sequence',
        ],
    ):

        for index, sample in enumerate(samples):
            df = sample.get_data(profile).data[columns]
            suffix_idx = f"_{suffix}{index+1}"
            if index == 0:
                merged_df = df.copy()
                merged_df = merged_df.rename(
                    columns={
                        col: col + suffix_idx for col in df.columns if col not in exceptions
                    }
                )

            else:
                df = df.rename(
                    columns={
                        col: col + suffix_idx for col in df.columns if col not in exceptions
                    }
                )
                merged_df = pd.merge(
                    merged_df,
                    df,
                    how="left",
                    on=["Nucleotide", "Sequence"],
                )


        return merged_df

    
    def average_columns(
        self,
        df : pd.DataFrame,
        avg_columns : list[str] = [
            "Modified_mutations",
            "Modified_effective_depth",
            "Modified_rate",
        ],
        sem_column : list[str] = [
            "Modified_rate"
        ],
    ):
        avg_df = df[["Nucleotide", "Sequence"]].copy()

        for pattern in avg_columns:
            matching_cols = [col for col in df.columns if pattern in col]
            if matching_cols:
                avg_df[f'{pattern}'] = df[matching_cols].mean(axis=1)

        for pattern in sem_column:
            matching_cols = [col for col in df.columns if pattern in col]
            if matching_cols:
                avg_df['Std_err'] = df[matching_cols].sem(axis=1)

        return avg_df


    def plot_scatter(
        self,
        column:str="Modified_rate",
        error:str="Std_err",
        label_size:int=None,
        ylabel:str=None,
        xlabel:str=None,
    ):
        """Generates scatter plots useful for fragmapper quality control.

        Args:
            column (str, optional):
                Dataframe column containing data to plot (must be avalible for
                the sample and control).
                Defaults to "Modified_rate".

        Returns:
            (matplotlib figure, matplotlib axis)
                Scatter plot with control values on the x-axis, sample values
                on the y-axis, and each point representing a nucleotide not
                filtered out in the fragmapper pipeline.
        """
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

        columns = [f"{column}_1", f"{column}_2"]
        scatter_data = self.data["fragmap"].data
        scatter_data = scatter_data[scatter_data["Valid"]].copy()
        sites = scatter_data[scatter_data["Site"]]
        non_sites = scatter_data[~scatter_data["Site"]]
        
        if error is not None:
            errors = [f"{error}_1", f"{error}_2"]
            (_, caps, _) = ax.errorbar(
                x=sites[columns[1]],
                y=sites[columns[0]],
                xerr=sites[errors[1]],
                yerr=sites[errors[0]],
                zorder=1,
                fmt='none',
                ecolor="#9FB3B3",
                elinewidth=1,
                alpha=0.4,
                capsize=1.25,
            )
            for cap in caps:
                cap.set_markeredgewidth(1)
                cap.set_alpha(1)
                cap.set_color("#9FB3B3")

        ax.scatter(
            x=non_sites[columns[1]],
            y=non_sites[columns[0]],
            zorder=2,
            s=5,
            alpha=0.25,
            c="#9FB3B3",
        )
        ax.scatter(
            x=sites[columns[1]],
            y=sites[columns[0]],
            zorder=3,
            s=25,
            alpha=1,
            color="#E86500",
        )

        if label_size is not None:
            for _, (nt, y, x) in sites[["Nucleotide"] + columns].iterrows():
                ax.text(x, y+0.001, int(nt), zorder=4, size=label_size)

        if ylabel is None:
            ylabel = str(self.sample1.sample)
        if xlabel is None:
            xlabel = str(self.sample2.sample)

        ax.set_ylabel(f"{ylabel}: {column}")
        ax.set_xlabel(f"{xlabel}: {column}")

        # legend = plt.legend(["Uncalled Nucleotides", "Frag-MaP Sites"])

        # Add a legend
        labels = [ '\n'.join(wrap(l, 11)) for l in ["Uncalled Nucleotides", "Frag-MaP Sites"]]

        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * 0.9, pos.height])
        ax.legend(
            labels,
            loc='center right',
            bbox_to_anchor=(1.5, 0.5),
            fontsize="small",
        )
 
        return fig, ax