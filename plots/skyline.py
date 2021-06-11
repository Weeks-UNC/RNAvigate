import matplotlib.pyplot as plt
import matplotlib as mp
from plots.plots import get_rows_columns, same_lengths, add_sequence
from plots import *


class Skyline():

    def __init__(self, profiles=[], labels=[]):
        self.profiles = profiles
        self.labels = labels

    def add_profile(self, profile, label):
        self.profiles.append(profile)
        self.labels.append(label)

    def set_figsize(self):
        left_inches = 0.9
        right_inches = 0.4
        ax_width = self.profiles[0].length * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        self.figsize = (fig_width, fig_height)

    def set_axis(self, ax=None, xlim=None, xticks=20, xticks_minor=5):
        if ax is None:
            _, self.ax = plt.subplots(1, figsize=self.figsize)
        else:
            self.ax = ax
        # x axis appearance
        length = self.profiles[0].length
        if xlim is None:
            xlim = [0, length]
        self.ax.set_xlim()
        self.ax.set_xticks(range(xlim[0], xlim[1], xticks))
        self.ax.set_xticks(range(xlim[0], xlim[1], xticks_minor), minor=True)

    def set_labels(self, axis_title="Raw Reactivity Profile",
                   legend_title="Samples", xlabel="Nucleotide",
                   ylabel="Profile"):
        self.ax.set_title(axis_title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.legend(title=legend_title)

    def plot_profiles(self, column='Reactivity_profile'):
        for profile, label in zip(self.profiles, self.labels):
            x = [0.5]
            y = [0]
            # converts standard plot to skyline plot.
            for n, r in zip(profile.data['Nucleotide'], profile.data[column]):
                x.extend([n - 0.5, n + 0.5])
                y.extend([r, r])
            x.append(x[-1])
            y.append(0)
            self.ax.plot(x, y, label=label)

    def make_plot(self, ax=None, column="Reactivity_profile", **kwargs):
        self.set_figsize()
        self.set_axis(ax)
        self.plot_profiles(column)
        add_sequence(self.ax, self.profiles[0].sequence)
        self.set_labels(**kwargs)
