import matplotlib as mp
import matplotlib.pyplot as plt
import math
import numpy as np
from abc import ABC, abstractmethod, abstractproperty


class Plot(ABC):
    def __init__(self, num_samples, rows=None, cols=None):
        self.length = num_samples
        self.rows, self.columns = self.get_rows_columns(rows, cols)
        figsize = self.get_figsize()
        self.fig, self.axes = plt.subplots(self.rows, self.columns,
                                           figsize=figsize, squeeze=False)
        self.i = 0
        self.pass_through = []

    def get_ax(self, i=None):
        if i is None:
            i = self.i
        row = i // self.columns
        col = i % self.columns
        return self.axes[row, col]

    def add_sample(self, sample, **kwargs):
        if isinstance(sample, list):
            for s in sample:
                self.add_sample(s, **kwargs)
            return
        for key in kwargs.keys():
            if key not in self.pass_through:
                kwargs[key] = sample.get_data(kwargs[key])
        self.plot_data(**kwargs)

    @classmethod
    def view_colormap(self, ax=None, ij=None, metric=None, ticks=None, values=None,
                      title=None, cmap=None):
        if ij == "ct_compare":
            metric = "Pairing"
            ticks = [10/6, 30/6, 50/6]
            values = ["Shared", "Structure 1", "Structure 2"]
            title = "Base-pairing comparison"
            cmap = mp.colors.ListedColormap([(0.6, 0.6, 0.6, 0.7),
                                             (0.15, 0.8, 0.6, 0.7),
                                             (0.6, 0.0, 1.0, 0.7)])
        elif ij is not None:
            metric = ij.metric
        if ticks is None:
            if metric == "Class":
                ticks = [10/6, 30/6, 50/6]
            else:
                ticks = [0, 2, 4, 6, 8, 10]
        if values is None:
            if metric == "Class":
                values = ['Complementary', 'Primary', 'Secondary']
            else:
                mn, mx = ij.min_max
                values = [f"{mn + ((mx-mn)/5)*i:.1f}" for i in range(6)]
        if title is None:
            title = f"{ij.datatype.capitalize()}: {metric.lower()}"
        if cmap is None:
            cmap = ij.cmap
        else:
            cmap = plt.get_cmap(cmap)
        colors = cmap(np.arange(cmap.N))

        if ax is None:
            _, ax = plt.subplots(1, figsize=(6, 2))
        ax.imshow([colors], extent=[0, 10, 0, 1])
        ax.set_title(title)
        ax.set_xticks(ticks)
        ax.set_xticklabels(values)
        ax.set_yticks([])

    def get_rows_columns(self, rows=None, cols=None):
        has_rows = isinstance(rows, int)
        has_cols = isinstance(cols, int)
        if has_rows and has_cols:
            return rows, cols
        elif has_rows:
            cols = math.ceil(self.length / rows)
        elif has_cols:
            rows = math.ceil(self.length / cols)
        elif self.length < 10:
            rows, cols = [(0, 0), (1, 1), (1, 2), (1, 3), (2, 2),
                          (2, 3), (2, 3), (3, 3), (3, 3), (3, 3)
                          ][self.length]
        else:
            cols = 4
            rows = math.ceil(self.length / cols)
        return rows, cols

    @classmethod
    def add_sequence(self, ax, sequence, yvalue=0.005):
        # set font style and colors for each nucleotide
        font_prop = mp.font_manager.FontProperties(
            family="monospace", style="normal", weight="bold", size="12")
        color_dict = {"A": "#f20000", "U": "#f28f00",
                      "G": "#00509d", "C": "#00c200"}
        # transform yvalue to a y-axis data value
        ymin, ymax = ax.get_ylim()
        yvalue = (ymax-ymin)*yvalue + ymin
        for i, seq in enumerate(sequence):
            col = color_dict[seq.upper()]
            ax.annotate(seq, xy=(i + 1, yvalue), xycoords='data',
                        fontproperties=font_prop,
                        color=col, horizontalalignment="center")

    @abstractmethod
    def get_figsize(self):
        pass

    @abstractmethod
    def plot_data(self):
        pass
