from .plots import Plot


class DistHist(Plot):
    def __init__(self, num_samples):
        super().__init__(num_samples)
        self.pass_through = ["atom", "ax"]

    def plot_data(self, structure, interactions, label, atom="O2'", ax=None):
        if ax is None:
            ax = self.get_ax()
        self.plot_all_distances(ax, structure, atom)
        self.plot_experimental_distances(ax, structure, interactions, atom)
        ax.set(title=label,
               ylabel="Experimental",
               xlabel="3D distance")
        self.i += 1

    def get_figsize(self):
        return (10*self.columns, 10*self.rows)

    def plot_all_distances(self, ax, structure, atom):
        matrix = structure.get_distance_matrix(atom=atom)
        pdb_dists = []
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if i+6 > j and matrix[i, j] != 1000:
                    pdb_dists.append(matrix.item(i, j))
        ax2 = ax.twinx()
        ax2.set(ylabel="Pairwise")
        ax2.hist(pdb_dists, bins=range(0, int(max(pdb_dists)), 5),
                 histtype="step", color="0.5", label="All distances")

    def plot_experimental_distances(self, ax, structure, interactions, atom):
        interactions.set_3d_distances(structure, atom)
        ij_dists = interactions.data.loc[interactions.data["mask"], "Distance"]
        ax.hist(ij_dists, bins=range(0, int(max(ij_dists)), 5), width=5,
                ec='none')
