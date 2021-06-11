import matplotlib as mp
import math


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
