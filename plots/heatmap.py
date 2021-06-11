
def get_distance_matrix(self, structure):
    if not hasattr(self, "distances"):
        self.distances = {}
    if structure in self.distances.keys():
        return self.distances[structure]
    if structure == "ct":
        matrix = self.ct.get_distance_matrix()
        self.distances["ct"] = matrix
    elif structure == "pdb":
        length = self.length["pdb"]
        fill = get_default_fill('Distance')
        matrix = np.full([length, length], fill)
        for i in range(length):
            for j in range(length):
                if matrix[i, j] == fill:
                    matrix[i, j] = self.get_3d_distance(i+1, j+1)
                    matrix[j, i] = matrix[i, j]
        self.distances["pdb"] = matrix
    return matrix


def plot_contour_distances(self, ax, structure, ij_data, levels):
    clip_pad = self.get_clip_pad(structure, ij_data)
    length = self.length[structure]+sum(clip_pad[1])
    fill = get_default_fill('Distance')
    distances = np.full([length, length], fill)
    start = clip_pad[1][0]
    end = length - clip_pad[1][1]
    distances[start:end, start:end] = self.get_distance_matrix(structure)
    if levels is None:
        levels = {"ct": [5], "pdb": [20, 500]}[structure]
    cmap = mp.colors.LinearSegmentedColormap.from_list('contours',
                                                       ['black', 'gray'])
    x_y = list(range(1, length+1))
    ax.contour(x_y, x_y, distances, levels=levels, cmap=cmap,
               linewidths=0.5)


def plot_heatmap_data(self, ax, ij_data, structure, metric=None):
    clip_pad = self.get_clip_pad(ij_data, structure)
    data = self.ij_data[ij_data].copy()
    if metric is None:
        metric = get_default_metric(ij_data)
    columns = ["i", "j", metric]
    if ij_data == "rings":
        data = data[columns+["+/-"]]
        data[metric] = data[metric]*data["+/-"]
        data = data[columns]
    else:
        data = data[columns]
    data[["i", "j"]] += clip_pad[1][0]
    length = self.length[ij_data] + sum(clip_pad[1])
    fill = get_default_fill(metric)
    data_im = np.full([length, length], fill)
    window = self.window[ij_data]
    for _, i, j, value in data.itertuples():
        data_im[i-1:i-1+window, j-1:j-1+window] = value
        data_im[j-1:j-1+window, i-1:i-1+window] = value
    min_max = get_default_min_max(metric)
    if metric == "Percentile":
        min_max = [0, 1]
        cmap = plt.get_cmap("jet")
        cmap = cmap(np.linspace(0, 1, 256))
        cmap[:, -1] = 0.6
    elif metric == "Class":
        data_im = data_im + 1
        min_max = [0, 3]
        cmap = get_default_cmap(metric)
        cmap = cmap(np.arange(cmap.N))
        no_data = np.array([1, 1, 1, 0])
        cmap = np.vstack((no_data, cmap))
    else:
        cmap = get_default_cmap(metric)
        cmap = cmap(np.arange(cmap.N))
    if ij_data != "rings":
        cmap[0, :] = [1, 1, 1, 0]
    cmap = mp.colors.ListedColormap(cmap)
    ax.imshow(data_im, cmap=cmap, vmin=min_max[0], vmax=min_max[1],
              interpolation='none')


def make_heatmap(self, ij_data, structure, metric=None, ax=None,
                 levels=None):
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(10, 10))
    # self.set_heatmap(ax)
    self.plot_heatmap_data(ax, ij_data, structure, metric)
    self.plot_contour_distances(ax, structure, ij_data, levels)
