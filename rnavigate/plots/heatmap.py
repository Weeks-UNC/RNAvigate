from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from .plots import Plot


class Heatmap(Plot):
    def __init__(self, num_samples, structure):
        super().__init__(num_samples, sharey=True)
        self.structure = structure
        for i in range(num_samples):
            ax = self.get_ax(i)
            ax.set_aspect("equal")
            mn_mx = (0.5, structure.length + 0.5)
            ax.set(ylim=mn_mx, xlim=mn_mx,
                   xticks=list(range(50, structure.length+1, 50)))
        self.pass_through = ["levels", "regions", "interpolation", "atom",
                             "plot_type"]

    def plot_data(self, ij, label, levels=None, regions=None,
                  interpolation='none', atom="O2'", plot_type="heatmap"):
        ax = self.get_ax()
        if regions is not None:
            self.plot_contour_regions(ax, ij, regions)
        else:
            self.plot_contour_distances(ax, levels, atom)
        assert plot_type in ["heatmap",
                             "kde"], "plot_type must be heatmap or kde"
        if plot_type == "heatmap":
            self.plot_heatmap_data(ax, ij, interpolation)
        elif plot_type == "kde":
            self.plot_kde_data(ax, ij)
        ax.set_title(label)
        self.i += 1

    def get_figsize(self):
        return (7*self.columns, 7*self.rows)

    def plot_contour_regions(self, ax, ij, regions):
        matrix = np.full([ij.length, ij.length], 0)
        for (mn1, mx1), (mn2, mx2) in regions:
            matrix[mn1-1:mx1-1, mn2-1:mx2-1] += 1
        levels = [0.5]
        cmap = LinearSegmentedColormap.from_list('contours', ['black', 'gray'])
        x_y = list(range(1, ij.length+1))
        ax.contour(x_y, x_y, matrix, levels=levels, cmap=cmap, linewidths=1)

    def plot_contour_distances(self, ax, levels, atom):
        structure = self.structure
        distances = structure.get_distance_matrix(atom)
        for i in range(structure.length):
            for j in range(i, structure.length):
                distances[i, j] = 0
        if levels is None:
            levels = {"ct": [5],
                      "ss": [5],
                      "pdb": [20]
                      }[structure.datatype]
        cmap = LinearSegmentedColormap.from_list('contours', ['black', 'gray'])
        x_y = list(range(1, structure.length+1))
        ax.contour(x_y, x_y, distances, levels=levels, cmap=cmap,
                   linewidths=1)

    def plot_heatmap_data(self, ax, ij, interpolation):
        structure = self.structure
        data = ij.data.copy()
        metric = ij.metric
        columns = ["i", "j", metric]
        if ij.datatype == "rings":
            data = data[columns+["Sign"]]
            data[metric] = data[metric]*data["Sign"]
        data = data[columns]
        am = ij.get_alignment_map(structure)
        length = structure.length
        data_im = np.full([length, length], ij.fill)
        window = ij.window
        for _, i, j, value in data.itertuples():
            i = am[i-1]
            j = am[j-1]
            data_im[j:j+window, i:i+window] = value
        min_max = ij.min_max
        if metric == "Percentile":
            min_max = [0, 1]
            cmap = plt.get_cmap("jet")
            cmap = cmap(np.linspace(0, 1, 256))
            cmap[:, -1] = 0.6
        elif metric == "Class":
            data_im = data_im + 1
            min_max = [0, 3]
            cmap = ij.cmap
            cmap = cmap(np.arange(cmap.N))
            no_data = np.array([1, 1, 1, 0])
            cmap = np.vstack((no_data, cmap))
        else:
            cmap = ij.cmap
            cmap = cmap(np.arange(cmap.N))
        if ij.datatype != "rings":
            cmap[0, :] = [1, 1, 1, 0]
        cmap = ListedColormap(cmap)
        ax.imshow(data_im, cmap=cmap, vmin=min_max[0], vmax=min_max[1],
                  interpolation=interpolation)

    def plot_kde_data(self, ax, ij, **kwargs):
        data = ij.data.loc[ij.data["mask"]]
        sns.kdeplot(ax=ax, data=data, x="i_offset", y="j_offset",
                    fill=True, levels=5, bw_adjust=0.2, cmap=ij.cmap,
                    common_norm=True, ** kwargs)
