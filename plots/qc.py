import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from .plots import Plot


class QC(Plot):

    def __init__(self, num_samples):
        self.num_samples = num_samples
        self.fig = plt.figure(figsize=self.get_figsize(),
                              constrained_layout=True)
        gs = self.fig.add_gridspec(2, 4)
        # Muts per molecule plots
        self.ax1 = self.fig.add_subplot(gs[0, 0])
        self.ax2 = self.fig.add_subplot(gs[0, 1],
                                        sharex=self.ax1, sharey=self.ax1)
        self.ax1.set(xlabel="Mutations per molecule",
                     ylabel="Percentage of Reads",
                     title="Untreated")
        self.ax2.set(xlabel="Mutations per molecule",
                     title="Modified")
        # Read length distribution
        self.ax3 = self.fig.add_subplot(gs[0, 2])
        self.ax4 = self.fig.add_subplot(gs[0, 3],
                                        sharex=self.ax3, sharey=self.ax3)
        self.ax3.set(xticks=range(12),
                     xlabel='Read Length',
                     ylabel='Percentage of Reads',
                     title='Untreated')
        self.ax4.set(xticks=range(12),
                     xlabel='Read Length',
                     title='Modified')
        xticklabels = [f"{x*50}" for x in range(12)]
        self.ax3.set_xticklabels(xticklabels, rotation=45)
        self.ax4.set_xticklabels(xticklabels, rotation=45)
        # Reactivities boxplot
        self.ax6 = self.fig.add_subplot(gs[1, :])
        self.i = 0
        self.colors = sns.color_palette("Paired")
        self.profiles = []
        self.pass_through = []

    def get_figsize(self):
        return (20, 10)

    def plot_data(self, log, profile, label):
        self.plot_MutsPerMol(log, label)
        self.plot_ReadLength(log, label)
        self.profiles.append(profile)
        self.i += 1
        if self.i == self.num_samples:
            handles, labels = self.ax1.get_legend_handles_labels()
            self.add_legend(handles, labels)
            self.make_boxplot(labels)

    def plot_MutsPerMol(self, log, label, upper_limit=12):
        x = log.data.loc[:upper_limit, 'Mutation_count']
        y1 = log.data.loc[:upper_limit, 'Modified_mutations_per_molecule']
        self.ax2.plot(x, y1, label=label)
        y2 = log.data.loc[:upper_limit, 'Untreated_mutations_per_molecule']
        self.ax1.plot(x, y2, label=label)

    def plot_ReadLength(self, log, label, upper_limit=12):
        width = 0.8/self.num_samples
        x = np.arange(upper_limit) - 0.4 - (width/2) + (width*self.i)
        y1 = log.data.loc[:upper_limit-1, 'Modified_read_length']
        self.ax4.bar(x, y1, width, label=label)
        y2 = log.data.loc[:upper_limit-1, 'Untreated_read_length']
        self.ax3.bar(x, y2, width, label=label)

    def add_legend(self, handles, labels):
        self.ax1.legend(handles, labels, title="Samples", loc=1)

    def make_boxplot(self, labels):
        cols = ["Modified_rate", "Untreated_rate"]
        data = []
        for i, profile in enumerate(self.profiles):
            data.append(profile.data[cols].copy().assign(Sample=i+1))
        data = pd.concat(data)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v', ax=self.ax6)
        ax.set(yscale='log',
               ylim=(0.0005, 0.5),
               ylabel="Mutation Rate",
               xticklabels=labels)
