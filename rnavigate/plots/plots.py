import matplotlib as mp
import matplotlib.pyplot as plt
import math
import numpy as np
from abc import ABC, abstractmethod
from ..data import Data


class Plot(ABC):
    def __init__(self, num_samples, rows=None, cols=None, figsize=None,
                 **kwargs):
        self.length = num_samples
        self.rows, self.columns = self.get_rows_columns(rows, cols)
        if figsize is None:
            figsize = self.get_figsize()
        self.fig, self.axes = plt.subplots(self.rows, self.columns,
                                           figsize=figsize, squeeze=False,
                                           **kwargs)
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
                if isinstance(kwargs[key], Data):
                    continue
                kwargs[key] = sample.get_data_list(kwargs[key])
        self.plot_data(**kwargs)

    @classmethod
    def view_colormap(self, ax=None, interactions=None, metric=None, ticks=None,
                      values=None, title=None, cmap=None):
        if ((interactions is None or interactions.datatype == "ct") and
                (None in [ticks, values, title, cmap])):
            ax.remove()
            return
        elif interactions == "ct_compare":
            metric = "Pairing"
            ticks = [10/6, 30/6, 50/6]
            values = ["Shared", "Structure 1", "Structure 2"]
            title = "Base-pairing comparison"
            cmap = mp.colors.ListedColormap([(0.6, 0.6, 0.6, 0.7),
                                             (0.15, 0.8, 0.6, 0.7),
                                             (0.6, 0.0, 1.0, 0.7)])
        elif interactions is not None:
            metric = interactions.metric
        if ticks is None:
            if metric == "Class":
                ticks = [10/6, 30/6, 50/6]
            else:
                ticks = [0, 2, 4, 6, 8, 10]
        if values is None:
            if metric == "Class":
                values = ['Complementary', 'Primary', 'Secondary']
            else:
                mn, mx = interactions.min_max
                values = [f"{mn + ((mx-mn)/5)*i:.1f}" for i in range(6)]
        if title is None:
            title = f"{interactions.datatype.capitalize()}: {metric.lower()}"
        if cmap is None:
            cmap = interactions.cmap
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

    def add_sequence(self, ax, sequence, yvalue=0, ytrans="axes"):
        # set font style and colors for each nucleotide
        font_prop = mp.font_manager.FontProperties(
            family="monospace", style="normal", weight="bold", size="12")
        color_dict = {"A": "#f20000", "U": "#f28f00",
                      "G": "#00509d", "C": "#00c200"}
        # transform yvalue to a y-axis data value
        if ytrans == "axes":
            trans = ax.get_xaxis_transform()
        elif ytrans == "data":
            trans = ax.transData
        sequence = sequence[self.region[0]-1: self.region[1]]
        for i, seq in enumerate(sequence):
            col = color_dict[seq.upper()]
            ax.text(i + self.region[0], yvalue, seq, fontproperties=font_prop,
                    transform=trans, color=col, horizontalalignment="center",
                    verticalalignment="center")

    @abstractmethod
    def get_figsize(self):
        pass

    @abstractmethod
    def plot_data(self):
        pass

    def save(self, filename):
        self.fig.savefig(filename)


def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
        else:
            spine.set_color('none')  # don't draw spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([])


def clip_spines(ax, spines):
    for spine in spines:
        if spine in ['left', 'right']:
            ticks = ax.get_yticks()
        if spine in ['top', 'bottom']:
            ticks = ax.get_xticks()
        ax.spines[spine].set_bounds((min(ticks), max(ticks)))
