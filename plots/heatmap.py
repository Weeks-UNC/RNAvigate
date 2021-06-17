
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
import numpy as np


class Heatmap():
    def __init__(self, heatmap, contour):
        self.heatmap = heatmap
        self.contour = contour

    def plot_contour_distances(self, ax, levels):
        alignment_map = self.contour.get_alignment_map(self.heatmap, full=True)
        length = len(alignment_map)
        fill = 1000
        distances = np.full([length, length], fill)
        contour_distance_matrix = self.contour.get_distance_matrix()
        for i1, i2 in enumerate(alignment_map):
            for j1, j2 in enumerate(alignment_map):
                if i2 != -1 and j2 != -1:
                    distances[i1, j1] = contour_distance_matrix[i2, j2]
        if levels is None:
            datatype = self.contour.datatype
            if datatype == "ct":
                levels = [5]
            elif datatype == "pdb":
                levels = [20, 500]
        cmap = LinearSegmentedColormap.from_list('contours', ['black', 'gray'])
        x_y = list(range(1, length+1))
        ax.contour(x_y, x_y, distances, levels=levels, cmap=cmap,
                   linewidths=0.5)

    def plot_heatmap_data(self, ax):
        data = self.heatmap.data.copy()
        metric = self.heatmap.metric
        columns = ["i", "j", metric]
        datatype = self.heatmap.datatype
        if datatype == "rings":
            data = data[columns+["Sign"]]
            data[metric] = data[metric]*data["Sign"]
        data = data[columns]
        alignment_map = self.heatmap.get_alignment_map(self.contour, full=True)
        length = len(alignment_map)
        fill = self.heatmap.fill
        data_im = np.full([length, length], fill)
        window = self.heatmap.window
        for _, i, j, value in data.itertuples():
            i = alignment_map[i-1]
            j = alignment_map[j-1]
            data_im[i:i+window, j:j+window] = value
            data_im[j:j+window, i:i+window] = value
        min_max = self.heatmap.min_max
        if metric == "Percentile":
            min_max = [0, 1]
            cmap = plt.get_cmap("jet")
            cmap = cmap(np.linspace(0, 1, 256))
            cmap[:, -1] = 0.6
        elif metric == "Class":
            data_im = data_im + 1
            min_max = [0, 3]
            cmap = self.heatmap.cmap
            cmap = cmap(np.arange(cmap.N))
            no_data = np.array([1, 1, 1, 0])
            cmap = np.vstack((no_data, cmap))
        else:
            cmap = self.heatmap.cmap
            cmap = cmap(np.arange(cmap.N))
        if self.heatmap.datatype != "rings":
            cmap[0, :] = [1, 1, 1, 0]
        cmap = ListedColormap(cmap)
        ax.imshow(data_im, cmap=cmap, vmin=min_max[0], vmax=min_max[1],
                  interpolation='none')

    def make_plot(self, ax=None, levels=None):
        if ax is None:
            _, ax = plt.subplots(1, figsize=(10, 10))
        # self.set_heatmap(ax)
        self.plot_heatmap_data(ax)
        self.plot_contour_distances(ax, levels)
