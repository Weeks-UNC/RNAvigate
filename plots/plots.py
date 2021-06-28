import matplotlib as mp
import matplotlib.pyplot as plt
import math
import numpy as np


def view_colormap(ij=None, metric=None, ticks=None, values=None,
                  title=None, cmap=None):
    """Given an ij_data (ij data) will display a colorbar for the default
    values (metric, cmap, min_max).

    Args:
        ij_data (str, optional): string matching an ij data type.
            Options are "rings", "pairs" or "deletions".
            Defaults to None.
        metric (str, optional): string matching column name of ij data.
            Default determined by get_default_metric.
        ticks (list, optional): locations to add ticks. scale is 0-10.
            Defaults to [0.5, 0.95], or [10/6, 30/6/, 50/6] for "Class" metric.
        title (str, optional): string for title of colorbar.
            Defaults to "{ij_data}: {metric}"
        cmap (str, optional): string matching a valid matplotlib colormap.
            Default determined by get_default_cmap.
    """
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

    _, ax = plt.subplots(1, figsize=(6, 2))
    ax.imshow([colors], extent=[0, 10, 0, 1])
    ax.set_title(title)
    ax.set_xticks(ticks)
    ax.set_xticklabels(values)
    ax.set_yticks([])


def get_rows_columns(number_of_samples, rows=None, cols=None):
    if isinstance(rows, int) and cols is None:
        cols = math.ceil(number_of_samples / rows)
    elif isinstance(cols, int) and rows is None:
        rows = math.ceil(number_of_samples / cols)
    elif number_of_samples < 10:
        rows, cols = [(0, 0), (1, 1), (1, 2), (1, 3), (2, 2),  # 0-4 samples
                      (2, 3), (2, 3), (3, 3), (3, 3), (3, 3)  # 5-9 samples
                      ][number_of_samples]
    else:
        cols = 4
        rows = math.ceil(number_of_samples / cols)
    return rows, cols


def same_lengths(*lists):
    it = iter(lists)
    the_len = len(next(it))
    return all(len(l) == the_len for l in it)


def add_sequence(ax, sequence, yvalue=0.005):
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
