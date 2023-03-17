import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from .plots import Plot


class QC(Plot):

    def __init__(self, num_samples):
        super().__init__(num_samples)
        if self.length == 1:
            self.ax_muts_unt = self.axes[0, 0]
            self.ax_muts_mod = self.axes[0, 0]
            self.ax_read_unt = self.axes[0, 1]
            self.ax_read_mod = self.axes[0, 1]
            self.ax_boxplot = self.axes[0, 2]
        elif self.length > 1:
            self.ax_muts_unt = self.axes[0, 0]
            self.ax_muts_mod = self.axes[0, 1]
            self.ax_muts_mod.get_shared_x_axes().join(self.ax_muts_mod,
                                                      self.ax_muts_unt)
            self.ax_muts_mod.get_shared_y_axes().join(self.ax_muts_mod,
                                                      self.ax_muts_unt)
            self.ax_muts_mod.set_yticklabels([])
            self.ax_read_unt = self.axes[0, 2]
            self.ax_read_mod = self.axes[0, 3]
            self.ax_muts_mod.get_shared_x_axes().join(self.ax_read_mod,
                                                      self.ax_read_unt)
            self.ax_muts_mod.get_shared_y_axes().join(self.ax_read_mod,
                                                      self.ax_read_unt)
            self.ax_read_mod.set_yticklabels([])
            gs = self.axes[1, 0].get_gridspec()
            for ax in self.axes[1, :]:
                ax.remove()
            self.ax_boxplot = self.fig.add_subplot(gs[1, :])
        self.ax_muts_unt.set(xlabel="Mutations per molecule",
                             ylabel="Percentage of Reads",
                             title="Untreated")
        self.ax_muts_mod.set(xlabel="Mutations per molecule",
                             title="Modified")
        self.ax_read_unt.set(xticks=range(12),
                             xlabel='Read Length',
                             ylabel='Percentage of Reads',
                             title='Untreated')
        self.ax_read_mod.set(xticks=range(12),
                             xlabel='Read Length',
                             title='Modified')
        if self.length == 1:
            self.ax_muts_mod.set_title("Mutations per Molecule Distribution")
            self.ax_read_mod.set_title("Read Length Distribution")
            self.ax_boxplot.set_title("Mutation rates")
        xticklabels = [f"{x*50}" for x in range(12)]
        self.ax_read_unt.set_xticklabels(xticklabels, rotation=45)
        self.ax_read_mod.set_xticklabels(xticklabels, rotation=45)
        self.profiles = []
        plt.tight_layout()

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=None,
                        width_ax_in=7, height_ax_in=7,
                        height_gap_in=1, width_gap_in=0.5,
                        top_in=1, bottom_in=0.5,
                        left_in=0.5, right_in=0.5):
        super().set_figure_size(fig=fig, ax=ax, rows=rows, cols=cols,
                                height_ax_rel=height_ax_rel,
                                width_ax_rel=width_ax_rel,
                                width_ax_in=width_ax_in,
                                height_ax_in=height_ax_in,
                                height_gap_in=height_gap_in,
                                width_gap_in=width_gap_in, top_in=top_in,
                                bottom_in=bottom_in, left_in=left_in,
                                right_in=right_in)

    def get_rows_columns(self, rows=None, cols=None):
        if self.length == 1:
            return (1, 3)
        else:
            return (2, 4)

    def get_figsize(self):
        if self.length == 1:
            return (30, 10)
        else:
            return (40, 20)

    def plot_data(self, log, profile, label):
        self.plot_MutsPerMol(log, label)
        self.plot_ReadLength(log, label)
        self.profiles.append(profile)
        self.i += 1
        if self.i == self.length:
            handles, labels = self.ax_muts_unt.get_legend_handles_labels()
            if self.length == 1:
                labels = ["Modified", "Untreated"]
            self.ax_muts_unt.legend(handles, labels, title="Samples", loc=1)
            if self.length == 1:
                labels = [label]
            self.make_boxplot(labels)

    def plot_MutsPerMol(self, log, label, upper_limit=12):
        x = log.data.loc[:upper_limit, 'Mutation_count']
        y1 = log.data.loc[:upper_limit, 'Modified_mutations_per_molecule']
        self.ax_muts_mod.plot(x, y1, label=label)
        y2 = log.data.loc[:upper_limit, 'Untreated_mutations_per_molecule']
        self.ax_muts_unt.plot(x, y2, label=label)

    def plot_ReadLength(self, log, label, upper_limit=12):
        if self.length == 1:
            width = 0.4
            x1 = np.arange(upper_limit) - 0.6 + width
            x2 = x1 + width
        else:
            width = 0.8/self.length
            x1 = np.arange(upper_limit) - 0.4 - (width/2) + (width*self.i)
            x2 = x1
        y1 = log.data.loc[:upper_limit-1, 'Modified_read_length']
        self.ax_read_mod.bar(x1, y1, width, label=label)
        y2 = log.data.loc[:upper_limit-1, 'Untreated_read_length']
        self.ax_read_unt.bar(x2, y2, width, label=label)

    def make_boxplot(self, labels):
        cols = ["Modified_rate", "Untreated_rate"]
        data = []
        for i, profile in enumerate(self.profiles):
            data.append(profile.data[cols].copy().assign(Sample=i+1))
        data = pd.concat(data)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        ax = sns.violinplot(x="Sample", y="value", hue='Rate', data=data,
                            ax=self.ax_boxplot, split=True, scale='width')
        ax.set(
            #    yscale='log',
            #    ylim=(0.00003, 0.3),
            ylabel="Mutation Rate",
            xticklabels=labels)
