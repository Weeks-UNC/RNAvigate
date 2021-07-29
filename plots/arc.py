from .plots import Plot
import numpy as np
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection


class AP(Plot):
    def __init__(self, num_samples, nt_length):
        self.nt_length = nt_length
        super().__init__(num_samples)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.set_aspect('equal')
            ax.yaxis.set_visible(False)
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
            ax.spines['bottom'].set_position('zero')
            ax.spines['top'].set_color('none')
            width = self.nt_length
            height = min(300, width/2)
            ax.set(xlim=(0, width),
                   ylim=(-height-5, height+1))

    def plot_data(self, ct, comp, ij, ij2, profile, label):
        ax = self.get_ax()
        if ij is not None:
            ax_ins1 = ax.inset_axes([0.05, 0.2, 0.3, 0.03])
            self.view_colormap(ax_ins1, ij)
        if ij2 is not None:
            ax_ins2 = ax.inset_axes([0.05, 0.3, 0.3, 0.03])
            self.view_colormap(ax_ins2, ij2)
        if comp is not None:
            ax_ins3 = ax.inset_axes([0.05, 0.8, 0.3, 0.03])
            self.view_colormap(ax_ins3, "ct_compare")
        self.add_patches(ax, ct, "top", comp)
        self.add_patches(ax, ij, "bottom")
        self.add_patches(ax, ij2, "bottom")
        self.add_sequence(ax, ct.sequence, yvalue=0.5)
        self.plot_profile(ax, profile, ct)
        self.add_title(ax, label)
        self.i += 1

    def add_patches(self, ax, data, panel, comp=None):
        if comp is not None:
            ij_colors = data.get_ij_colors(comp)
        elif data is not None:
            ij_colors = data.get_ij_colors()
        else:
            return
        patches = []
        for i, j, color in zip(*ij_colors):
            if panel == "top":
                center = ((i+j)/2., 0)
                theta1 = 0
                theta2 = 180
            elif panel == "bottom":
                center = ((i+j)/2., -2)
                theta1 = 180
                theta2 = 360
            radius = 0.5+(j-i)/2.
            patches.append(Wedge(center, radius, theta1, theta2, color=color,
                                 width=1, ec='none'))
        ax.add_collection(PatchCollection(patches, match_original=True))

    def add_title(self, ax, label):
        ax.annotate(label, xy=(0.1, 0.9), xycoords="axes fraction",
                    fontsize=60, ha='left')

    def get_figsize(self):
        dim = self.nt_length * 0.1 + 1
        return (dim*self.columns, dim*self.rows)

    def plot_profile(self, ax, profile, ct):
        columns = ["Norm_profile", "Reactivity_profile", "Modified_rate"]
        for column in columns:
            if column in profile.data.columns:
                break
        near_black = (0, 0, 1 / 255.0)
        orange_thresh = 0.4
        red_thresh = 0.85
        am = profile.get_alignment_map(ct)
        values = np.full(ct.length, np.nan)
        nts = np.arange(ct.length)+1
        yerr = np.full(ct.length, np.nan)
        for i1, i2 in enumerate(am):
            if i2 != -1:
                values[i2] = profile.data.loc[i1, column]
                if 'Norm_stderr' in profile.data.columns:
                    yerr[i2] = profile.data.loc[i1, 'Norm_stderr']
        cindex = np.zeros(len(values), dtype=int)
        # where values are not NaNs, add 1 to color index array
        not_nan_values = values[~np.isnan(values)]
        not_nans = np.ones(len(not_nan_values), dtype=int)
        # where values excede orange threshold, add 1 to color index array
        not_nans[np.array(not_nan_values > orange_thresh, dtype=bool)] += 1
        # where values excede red threshold (0.85), add 1 to color index array
        not_nans[np.array(not_nan_values > red_thresh, dtype=bool)] += 1
        cindex[~np.isnan(values)] += not_nans
        # create color map array based on cindex
        colormap = np.array(["0.80", "black", "orange", "red"])[cindex]
        ax = self.get_ax()
        ax.bar(nts, values*5, align="center",
               width=1.05, color=colormap, edgecolor=colormap, linewidth=0.0,
               yerr=yerr, ecolor=near_black, capsize=1)
