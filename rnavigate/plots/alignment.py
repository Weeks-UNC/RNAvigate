from rnavigate import plots
import numpy as np


class Alignment(plots.Plot):
    def __init__(self, num_samples, rows=None, cols=1, figsize=None, **kwargs):
        super().__init__(num_samples, rows, cols, figsize, **kwargs)

    def get_figsize(self):
        return (50, 1.2)

    def plot_data(self, alignment, label, axis=None):
        if axis is None:
            axis = self.get_ax()
        alignment1 = alignment.alignment1
        self.region = (1, len(alignment1))
        plot_alignment(self, axis, alignment, label,
                       spines_positions={"top": 1, "bottom": -1})
        offset = 0.78
        axis.set(xlim=(0.5, len(alignment1)+0.5),
               ylim=(-1, 1),
               yticks=[0-offset, 0+offset],
               yticklabels=label[::-1]
               )
        for spine in ['top', 'bottom', 'left', 'right']:
            axis.spines[spine].set_color(None)

    def set_figure_size(self, fig=None, axis=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=0.1,
                        width_ax_in=None, height_ax_in=1.2,
                        height_gap_in=1, width_gap_in=0.5,
                        top_in=1, bottom_in=0.5,
                        left_in=0.5, right_in=0.5):
        super().set_figure_size(fig=fig, axis=axis, rows=rows, cols=cols,
                                height_ax_rel=height_ax_rel,
                                width_ax_rel=width_ax_rel,
                                width_ax_in=width_ax_in,
                                height_ax_in=height_ax_in,
                                height_gap_in=height_gap_in,
                                width_gap_in=width_gap_in, top_in=top_in,
                                bottom_in=bottom_in, left_in=left_in,
                                right_in=right_in)


def plot_alignment(plot, axis, alignment, label, center=-0.04, offset=0.78,
                   spines_positions=None):
    if spines_positions is None:
        spines_positions = {"top": center+offset*1.1,
                            "bootom": center-offset*1.1}
    seq1 = alignment.sequence1
    seq2 = alignment.sequence2
    al1 = alignment.alignment1
    al2 = alignment.alignment2
    al1 = al1.replace(".", "-")
    plot.add_sequence(axis, sequence=al1, yvalue=center+offset, ytrans="data")
    plot.add_sequence(axis, sequence=al2, yvalue=center-offset, ytrans="data")

    am2 = np.array([i for i, nt in enumerate(al2) if nt != "-"])
    xtick_labels2 = np.arange(10, len(seq2)+1, 10)
    xticks2 = am2[xtick_labels2 - 1] + 1
    axis.set(xticks=xticks2, xticklabels=xtick_labels2)
    axis.spines["bottom"].set(position=("data", spines_positions["bottom"]),
                            visible=False)
    for label in axis.get_xticklabels():
        label.set_bbox({"facecolor": "white",
                        "edgecolor": "None",
                        "alpha": 0.5,
                        "boxstyle": "round,pad=0.1,rounding_size=0.2"})

    am1 = np.array([i for i, nt in enumerate(al1) if nt != "-"])
    xtick_labels1 = np.arange(10, len(seq1)+1, 10)
    xticks1 = am1[xtick_labels1 - 1] + 1
    xlims = [plot.region[0]-0.5, plot.region[1]+0.5]
    ax2 = axis.twiny()
    ax2.set(xticks=xticks1, xticklabels=xtick_labels1, xlim=xlims)
    ax2.spines["top"].set(position=("data", spines_positions["top"]),
                          visible=False)
    for label in ax2.get_xticklabels():
        label.set_bbox({"facecolor": "white",
                        "edgecolor": "None",
                        "alpha": 0.5,
                        "boxstyle": "round,pad=0.1,rounding_size=0.2"})

    for spine in ["top", "bottom", "left", "right"]:
        ax2.spines[spine].set_color(None)
    for label in ax2.get_xticklabels():
        label.set_bbox({"facecolor": "white",
                        "edgecolor": "None",
                        "alpha": 0.5,
                        "boxstyle": "round,pad=0.1,rounding_size=0.2"})
    for idx, (nt1, nt2) in enumerate(zip(al1, al2)):
        axis.fill_between(
            x=[idx+0.5, idx+1.5],
            y1=[center+(0.6*offset)*(nt1 != "-")]*2,
            y2=[center-(0.6*offset)*(nt2 != "-")]*2,
            color="grey", ec="none")
        if nt1 != nt2 and "-" not in [nt1, nt2]:
            axis.fill_between(
                x=[idx+0.5, idx+1.5],
                y1=[center-0.6*offset]*2,
                y2=[center+0.6*offset]*2,
                color="red", ec="none")
