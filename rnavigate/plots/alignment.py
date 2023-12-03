from rnavigate import plots


class Alignment(plots.Plot):
    def __init__(self, num_samples, rows=None, cols=1, **kwargs):
        super().__init__(num_samples, rows, cols, **kwargs)

    def plot_data(self, alignment, label, ax=None):
        if ax is None:
            ax = self.get_ax()
        alignment1 = alignment.alignment1
        self.region = (1, len(alignment1))
        ax.set(xlim=(0, len(alignment1) + 1), ylim=(-5, 5), yticks=[])
        plots.plot_sequence_alignment(ax, alignment, label)
        for spine in ["top", "bottom", "left", "right"]:
            ax.spines[spine].set_color(None)

    def set_figure_size(
        self,
        fig=None,
        ax=None,
        rows=None,
        cols=None,
        height_ax_rel=0.03,
        width_ax_rel=0.03,
        width_ax_in=None,
        height_ax_in=None,
        height_gap_in=1,
        width_gap_in=0.5,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        super().set_figure_size(
            fig=fig,
            ax=ax,
            rows=rows,
            cols=cols,
            height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel,
            width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in,
            width_gap_in=width_gap_in,
            top_in=top_in,
            bottom_in=bottom_in,
            left_in=left_in,
            right_in=right_in,
        )
