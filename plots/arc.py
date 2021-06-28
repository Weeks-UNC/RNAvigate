from plots import *
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt


class AP():
    def __init__(self, structures=[], ijs=[], profiles=[], labels=[]):
        self.structures = []
        self.ijs = []
        self.profiles = []
        self.labels = []
        self.ij_windows = []
        for sample in zip(structures, ijs, profiles, labels):
            self.add_sample(sample)

    def add_sample(self, structure, ij, profile, label):
        self.structures.append(structure)
        self.profiles.append(profile)
        self.labels.append(label)
        self.ijs.append(ij.get_ij_colors())
        self.ij_windows.append(ij.window)
        view_colormap(ij)

    def add_arc(self, ax, i, j, window, color, panel):
        if panel == "top":
            center = ((i+j)/2., 0)
            theta1 = 0
            theta2 = 180
        if panel == "bottom":
            center = ((i+j+window-1)/2., -2)
            theta1 = 180
            theta2 = 360
        radius = 0.5+(j+window-1-i)/2.
        arc = mp.patches.Wedge(center, radius, theta1, theta2, color=color,
                               width=window, ec='none')
        ax.add_patch(arc)

    def get_figsize(self, rows, cols, sequence="ct"):
        dim = self.profiles[0].length * 0.1 + 1
        return (dim*cols, dim*rows)

    def set_plot(self, ax, i):
        ax.set_aspect('equal')
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        width = self.profiles[i].length
        height = min(300, width/2)
        ax.set(xlim=(0, width),
               ylim=(-height-5, height+1))

    def plot_ct(self, ax, ct, panel):
        ct_pairs = ct.pairList()
        for i, j in ct_pairs:
            self.add_arc(ax, i, j, 1, (0.1, 0.1, 0.1, 0.7), panel)

    def plot_ctcompare(self, ax, ct1, ct2, panel):
        ct1 = set(ct1.pairList())
        ct2 = set(ct2.pairList())
        shared = ct1.union(ct2)
        ref = ct1.difference(ct2)
        comp = ct2.difference(ct1)
        sharedcolor = (150/255., 150/255., 150/255., 0.7)
        refcolor = (38/255., 202/255., 145/255., 0.7)
        compcolor = (153/255., 0.0, 1.0, 0.7)
        for i, j in comp:
            self.add_arc(ax, i, j, 1, compcolor, panel)
        for i, j in ref:
            self.add_arc(ax, i, j, 1, refcolor, panel)
        for i, j in shared:
            self.add_arc(ax, i, j, 1, sharedcolor, panel)
        handles = [mp.patches.Patch(color=sharedcolor, alpha=0.7),
                   mp.patches.Patch(color=refcolor, alpha=0.7),
                   mp.patches.Patch(color=compcolor, alpha=0.7)]
        labels = ["Shared pairs", "Reference only", "Compared only"]
        ax.legend(title="top panel", handles=handles, labels=labels,
                  loc="upper right")

    def plot_profile(self, ax, profile):
        near_black = (0, 0, 1 / 255.0)
        orange_thresh = 0.4
        red_thresh = 0.85
        values = profile.data['Norm_profile']
        cindex = np.zeros(len(values), dtype=int)
        # where values are not NaNs, add 1 to color index array
        cindex[np.array(np.logical_not(np.isnan(values)), dtype=bool)] += 1
        # where values excede orange threshold, add 1 to color index array
        cindex[np.array(values > orange_thresh, dtype=bool)] += 1
        # where values excede red threshold (0.85), add 1 to color index array
        cindex[np.array(values > red_thresh, dtype=bool)] += 1
        # create color map array based on cindex
        colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
        ax.bar(profile.data['Nucleotide'], values*5, align="center",
               width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
               yerr=profile.data['Norm_stderr'], ecolor=near_black, capsize=1)

    def plot_ij(self, ax, ij, panel, metric=None):
        ij_colors = ij.get_ij_colors(metric)
        window = ij.window
        for i, j, color in zip(*ij_colors):
            self.add_arc(ax, i, j, window, color, panel)

    def plot_data(self, ax, index=0):
        def cls_name(data):
            return type(data).__name__

        def plot_func(ax, data, panel):
            f = {'CT': self.plot_ct, 'IJ': self.plot_ij}
            if isinstance(data, list):
                if cls_name(data[0]) == cls_name(data[1]) == "CT":
                    self.plot_ctcompare(ax, data[0], data[1], panel)
                else:
                    for datum in data:
                        cls = cls_name(datum)
                        f[cls](ax, datum, panel)
            else:
                cls = cls_name(data)
                f[cls](ax, data, panel)

        plot_func(ax, self.structures[index], "top")
        window = self.ij_windows[index]
        for i, j, color in zip(*self.ijs[index]):
            self.add_arc(ax, i, j, window, color, "bottom")
        self.plot_profile(ax, self.profiles[index])

    def make_plot(self, ax=None):
        length = len(self.labels)
        rows, columns = get_rows_columns(length)
        if ax is None:
            figsize = self.get_figsize(rows, columns)
            _, axes = plt.subplots(rows, columns, figsize=figsize,
                                   squeeze=False)
        for i in range(length):
            row = i // columns
            col = i % columns
            ax = axes[row, col]
            self.set_plot(ax, i)
            self.plot_data(ax, i)
            add_sequence(ax, self.profiles[i].sequence, yvalue=0.5)
            ax.annotate(self.labels[i], xy=(0.1, 0.9),
                        xycoords="axes fraction", fontsize=60, ha='left')
