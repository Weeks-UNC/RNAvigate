
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
import numpy as np
from .plots import Plot


class Heatmap(Plot):
    def plot_data(self, structure, ij, levels=None):
        ax = self.get_ax()
        self.plot_contour_distances(ax, structure, ij, levels)
        self.plot_heatmap_data(ax, structure, ij)
        self.i += 1

    def get_figsize(self):
        return (7*self.columns, 7*self.rows)

    def plot_contour_distances(self, ax, structure, ij, levels):
        am, length = structure.get_alignment_map(ij, full=True)
        fill = 1000
        distances = np.full([length, length], fill)
        structure_distances = structure.get_distance_matrix()
        for i1, i2 in enumerate(am):
            for j1, j2 in enumerate(am):
                distances[i2, j2] = structure_distances[i1, j1]
        if levels is None:
            levels = {"ct": [5], "pdb": [20, 500]}[structure.datatype]
        cmap = LinearSegmentedColormap.from_list('contours', ['black', 'gray'])
        x_y = list(range(1, length+1))
        ax.contour(x_y, x_y, distances, levels=levels, cmap=cmap,
                   linewidths=0.5)

    def plot_heatmap_data(self, ax, structure, ij):
        data = ij.data.copy()
        metric = ij.metric
        columns = ["i", "j", metric]
        if ij.datatype == "rings":
            data = data[columns+["Sign"]]
            data[metric] = data[metric]*data["Sign"]
        data = data[columns]
        am, length = ij.get_alignment_map(structure, full=True)
        data_im = np.full([length, length], ij.fill)
        window = ij.window
        for _, i, j, value in data.itertuples():
            i = am[i-1]
            j = am[j-1]
            data_im[i:i+window, j:j+window] = value
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
                  interpolation='none')
