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
        xticklabels = [f"{x*50}-\n{x*50+50}" for x in range(12)]
        self.ax3.set_xticklabels(xticklabels)
        self.ax4.set_xticklabels(xticklabels)
        # Legend
        self.ax5 = self.fig.add_subplot(gs[1, 0])
        self.ax5.set(title="Samples", axis="off", frame="on")
        # Reactivities boxplot
        self.ax6 = self.fig.add_subplot(gs[1, 1:])
        self.i = 0
        self.colors = sns.color_palette("Paired")

    def get_figsize(self):
        return (20, 10)

    def plot_data(self, log, profile, label):
        self.plot_MutsPerMol(log)
        self.make_boxplotprofiles.append(profile)
        self.labels.append(label)

    def plot_MutsPerMol(self, log, upper_limit=10):
        x = log.data.iloc[:upper_limit, 'Mutation_count']
        y1 = log.data.iloc[:upper_limit, 'Modified_mutations_per_molecule']
        self.ax2.plot(x, y1, color=self.colors[self.i])
        y2 = log.data.iloc[:upper_limit, 'Untreated_mutations_per_molecule']
        self.ax1.plot(x, y2, color=self.colors[self.i])

    def plot_ReadLength(self, log, upper_limit=10):
        width = 0.8/self.num_samples
        x = np.arange(upper_limit) - 0.4 - (width/2) + (width*self.i)
        y1 = log.data.iloc[:upper_limit, 'Modified_read_length']
        self.ax4.bar(x, y1, width)
        y2 = log.data.iloc[:upper_limit, 'Untreated_read_length']
        self.ax3.bar(x, y2, width)

    def add_to_legend(self, label):
        pass

    def make_boxplot(self, ax):
        cols = ["Modified_rate", "Untreated_rate"]
        xticklabels = [label for label in self.labels]
        data = []
        for i, profile in enumerate(self.profiles):
            data.append(profile.data[cols].copy().assign(Sample=i+1))
        data = pd.concat(data)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v', ax=ax)
        ax.set(yscale='log',
               ylim=(0.0005, 0.5),
               ylabel="Mutation Rate",
               xticklabels=xticklabels)
