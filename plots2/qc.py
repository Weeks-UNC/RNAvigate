import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class QC():
    def __init__(self, logs=None, profiles=[], labels=[]):
        self.logs = logs
        self.profiles = profiles
        self.labels = labels

    def add_data(self, log, profile, label):
        self.logs.append(log)
        self.profiles.append(profile)
        self.labels.append(label)

##############################################################################
# Mutations per molecule plotting
#   set, plot, make
##############################################################################

    def set_MutsPerMol(self, ax):
        ax.legend(title="Samples")
        ax.set(xlabel="Mutations per molecule",
               ylabel="Percentage of Reads",
               title='Mutations per molecule distribution')

    def plot_MutsPerMol(self, mod_ax, unt_ax, upper_limit=10):
        for log, label in zip(self.logs, self.labels):
            x = log.data['Mutation_count'][:upper_limit]
            y = log.data['Modified_mutations_per_molecule'][:upper_limit]
            mod_ax.plot(x, y, label=label+": Modified")
            y = log.data['Untreated_mutations_per_molecule'][:upper_limit]
            unt_ax.plot(x, y, label=label+": Untreated")

    def make_MutsPerMol(self, mod_ax, unt_ax):
        self.plot_MutsPerMol(mod_ax, unt_ax)
        self.set_MutsPerMol(mod_ax)
        self.set_MutsPerMol(unt_ax)

##############################################################################
# Read Length Distribution plotting
#   set, plot, make
##############################################################################

    def set_ReadLength(self, ax, upper_limit=10):
        ax.legend(title="Samples")
        ax.set(xticks=range(upper_limit),
               xlabel='Read Length',
               ylabel='Percentage of Reads',
               title='Read length distribution')
        ax.set_xticklabels(self.logs[0].data.loc[:upper_limit-1, "Read_length"],
                           rotation=45, ha='right')

    def plot_ReadLength(self, mod_ax, unt_ax, upper_limit=10):
        width = 0.8/len(self.logs)
        for i, (log, label) in enumerate(zip(self.logs, self.labels)):
            x = np.arange(upper_limit) - 0.4 - (width/2) + (width*i)
            y = log.data['Modified_read_length'][:upper_limit]
            mod_ax.bar(x, y, width, label=label+": Modified")
            y = log.data['Untreated_read_length'][:upper_limit]
            unt_ax.bar(x, y, width, label=label+": Untreated")

    def make_ReadLength(self, mod_ax, unt_ax):
        self.plot_ReadLength(mod_ax, unt_ax)
        self.set_ReadLength(mod_ax)
        self.set_ReadLength(unt_ax)

##############################################################################
# Reactivity boxplot plotting
##############################################################################

    def make_boxplot(self, ax):
        cols = ["Modified_rate", "Untreated_rate"]
        xticklabels = [label for label in self.labels]
        data = []
        for i, profile in enumerate(self.profiles):
            data.append(profile.data[cols].copy().assign(Sample=i+1))
        data = pd.concat(data)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v')
        ax.set(yscale='log',
               ylim=(0.0005, 0.5),
               ylabel="Mutation Rate",
               xticklabels=xticklabels)

##############################################################################
# Combination plotting
##############################################################################

    def make_plot(self):
        equal_data = len(self.profiles) == len(self.logs) == len(self.labels)
        assert equal_data, "Uneven number of profiles, logs, and labels."
        if len(self.profiles) > 1:
            fig = plt.figure(figsize=(20, 10))
            mpm_mod = fig.add_subplot(241)
            mpm_unt = fig.add_subplot(242)
            rl_mod = fig.add_subplot(243)
            rl_unt = fig.add_subplot(244)
            boxplot = fig.add_subplot(212)
        else:
            _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21, 7))
            mpm_mod, mpm_unt = ax1, ax1
            rl_mod, rl_unt = ax2, ax2
            boxplot = ax3
        self.make_MutsPerMol(mpm_mod, mpm_unt)
        self.make_ReadLength(rl_mod, rl_unt)
        self.make_boxplot(boxplot)
        plt.tight_layout()
