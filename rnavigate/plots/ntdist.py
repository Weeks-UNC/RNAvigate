import numpy as np
import seaborn as sns
from rnavigate import plots, styles


class NucleotideDistribution(plots.Plot):
    def __init__(self, num_samples, **plot_kwargs):
        super().__init__(num_samples, sharex=True, cols=1, **plot_kwargs)

    def set_figure_size(
            self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None,
            width_ax_rel=None, width_ax_in=2, height_ax_in=2,
            height_gap_in=0.2, width_gap_in=0.4, top_in=1, bottom_in=1,
            left_in=1, right_in=1
            ):
        super().set_figure_size(
            fig=fig, ax=ax, rows=rows, cols=cols, height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel, width_ax_in=width_ax_in,
            height_ax_in=height_ax_in, height_gap_in=height_gap_in,
            width_gap_in=width_gap_in, top_in=top_in, bottom_in=bottom_in,
            left_in=left_in, right_in=right_in
            )

    def plot_data(self, profile, label, column=None, normalize=None, ax=None):
        profile = profile.copy()
        if ax is None:
            ax = self.get_ax()
            self.i += 1
        if column is None:
            column = profile.metric
        if normalize is not None:
            profile.normalize(**normalize)
        data = profile.data
        nt_idx = data[column] > 0
        sns.kdeplot(
            ax=ax,
            data=data.loc[nt_idx],
            x=np.log10(data.loc[nt_idx, column]),
            hue="Sequence",
            hue_order=["A", "U", "C", "G"],
            palette={nt: styles.get_nt_color(nt) for nt in "AUGC"},
            common_norm=False
            )
        ax.set(
            xlim=(-2.5, 1.5),
            xticks=[-2, -1, 0, 1],
            xticklabels=[0.01, 0.1, 1, 10],
            xlabel="Normalized profile (log scale)",
            yticks=[],
            ylabel=label,
            )
        ax.axvspan(-1, 0, color="lightgrey", alpha=0.4, ec="none")
