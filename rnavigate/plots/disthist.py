import numpy as np
from rnavigate import plots


class DistHist(plots.Plot):
    def __init__(self, num_samples, **plot_kwargs):
        super().__init__(num_samples, sharey=True, **plot_kwargs)
        self.axes2 = {}
        base_ax2 = None
        for row in self.axes:
            for ax in row:
                ax2 = ax.twinx()
                self.axes2[ax] = ax2
                if base_ax2 is None:
                    base_ax2 = ax2
                else:
                    ax2.sharey(base_ax2)

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=None,
                        width_ax_in=2, height_ax_in=2,
                        height_gap_in=1, width_gap_in=0.4,
                        top_in=1, bottom_in=1,
                        left_in=1, right_in=1):
        super().set_figure_size(fig=fig, ax=ax, rows=rows, cols=cols,
                                height_ax_rel=height_ax_rel,
                                width_ax_rel=width_ax_rel,
                                width_ax_in=width_ax_in,
                                height_ax_in=height_ax_in,
                                height_gap_in=height_gap_in,
                                width_gap_in=width_gap_in, top_in=top_in,
                                bottom_in=bottom_in, left_in=left_in,
                                right_in=right_in)

    def plot_data(self, structure, interactions, bg_interactions, label,
                  atom="O2'", ax=None):
        if ax is None:
            ax = self.get_ax()
        ax2 = self.axes2[ax]
        if bg_interactions is not None:
            self.plot_experimental_distances(ax=ax2, structure=structure,
                                             interactions=bg_interactions,
                                             atom=atom, histtype='step')
        else:
            self.plot_structure_distances(
                ax=ax2, structure=structure, atom=atom)
        self.plot_experimental_distances(ax=ax, structure=structure,
                                         interactions=interactions, atom=atom)
        ax.set(title=label)
        self.i += 1
        if self.i == self.length:
            for leftmost_axis in self.axes[:, 0]:
                leftmost_axis.set(ylabel="Experimental")
            for rightmost_axis in self.axes[:, -1]:
                self.axes2[rightmost_axis].set(ylabel="Pairwise")
            for top_axis in self.axes[-1, :]:
                top_axis.set(xlabel="3D distance")
            for row in self.axes:
                for not_rightmost_axis in row[:-1]:
                    self.axes2[not_rightmost_axis].yaxis.set_tick_params(
                        labelright=False)

    def plot_structure_distances(self, ax, structure, atom):
        matrix = structure.get_distance_matrix(atom=atom)
        dists = []
        for i in range(len(matrix)-6):
            for j in range(i+6, len(matrix)):
                if not np.isnan(matrix[i, j]):
                    dists.append(matrix.item(i, j))
        ax.hist(dists, bins=range(0, int(max(dists))+5, 5),
                histtype="step", color="0.5", label="All distances")

    def plot_experimental_distances(self, ax, structure, interactions, atom,
                                    histtype='bar'):
        interactions.set_3d_distances(structure, atom)
        ij_dists = interactions.data.loc[interactions.data['mask'], "Distance"]
        if (len(ij_dists) > 0) and (histtype == 'bar'):
            ax.hist(ij_dists, bins=range(0, int(max(ij_dists))+5, 5),
                      width=5, ec='none')
        elif (len(ij_dists) > 0) and (histtype == 'step'):
            ax.hist(ij_dists, bins=range(0, int(max(ij_dists))+5, 5),
                      histtype='step', color='0.5', label='All distances')
