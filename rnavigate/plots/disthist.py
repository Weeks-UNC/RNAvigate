from .plots import Plot


class DistHist(Plot):
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
        self.pass_through = ["atom", "ax"]

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
            for ax in self.axes[:, 0]:
                ax.set(ylabel="Experimental")
            for ax in self.axes[:, -1]:
                self.axes2[ax].set(ylabel="Pairwise")
            for ax in self.axes[-1, :]:
                ax.set(xlabel="3D distance")
            for row in self.axes:
                for ax in row[:-1]:
                    self.axes2[ax].yaxis.set_tick_params(labelright=False)

    def get_figsize(self):
        return (10*self.columns, 10*self.rows)

    def plot_structure_distances(self, ax, structure, atom):
        matrix = structure.get_distance_matrix(atom=atom)
        dists = []
        for i in range(len(matrix)-6):
            for j in range(i+6, len(matrix)):
                if matrix[i, j] != 1000:
                    dists.append(matrix.item(i, j))
        ax.hist(dists, bins=range(0, int(max(dists))+5, 5),
                histtype="step", color="0.5", label="All distances")

    def plot_experimental_distances(self, ax, structure, interactions, atom,
                                    histtype='bar'):
        interactions.set_3d_distances(structure, atom)
        ij_dists = interactions.data.loc[interactions.data["mask"], "Distance"]
        if (len(ij_dists) > 0) and (histtype == 'bar'):
            ax.hist(ij_dists, bins=range(0, int(max(ij_dists))+5, 5), width=5,
                    ec='none')
        elif (len(ij_dists) > 0) and (histtype == 'step'):
            ax.hist(ij_dists, bins=range(0, int(max(ij_dists))+5, 5),
                    histtype='step', color='0.5', label='All distances')
