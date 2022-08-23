from .plots import Plot
import numpy as np
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection


class AP(Plot):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1]-region[0]+1
            self.region = region
        super().__init__(num_samples, **kwargs)
        for i in range(self.length):
            self.set_axis(self.get_ax(i))
        self.pass_through = ["ax", "colorbar", "seqbar", "title", "ij_panel",
                             "ij2_panel", "ct_panel"]

    def set_axis(self, ax, xticks=20, xticks_minor=5):
        ax.set_aspect('equal')
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        height = min(300, self.nt_length/2)
        mn, mx = self.region
        ax.set(xlim=(mn - 0.5, mx + 0.5),
               ylim=(-height-5, height+1),
               xticks=list(range(xticks, mx, xticks)))

    def plot_data(self, ct, comp, ij, ij2, profile, label, ax=None,
                  colorbar=True, seqbar=True, title=True, ij_panel="bottom",
                  ij2_panel="bottom", ct_panel="top"):
        ax = self.get_ax(ax)
        if colorbar:
            ax_ins1 = ax.inset_axes([0.05, 0.2, 0.3, 0.03])
            self.view_colormap(ax_ins1, ij)
            ax_ins2 = ax.inset_axes([0.05, 0.3, 0.3, 0.03])
            self.view_colormap(ax_ins2, ij2)
            ax_ins3 = ax.inset_axes([0.05, 0.8, 0.3, 0.03])
            if comp is not None:
                self.view_colormap(ax_ins3, "ct_compare")
            else:
                self.view_colormap(ax_ins3, ct)
        self.add_patches(ax, ct, ct_panel, comp)
        self.add_patches(ax, ij, ij_panel)
        self.add_patches(ax, ij2, ij2_panel)
        self.plot_profile(ax, profile, ct)
        if seqbar:
            self.add_sequence(ax, ct.sequence, yvalue=0.5)
        if title:
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
        mn, mx = self.region
        for i, j, color in zip(*ij_colors):
            if j > i:  # flip the order
                i, j = j, i
            if ((i < mn) and (j < mn)) or ((i > mx) and (j > mx)):
                continue
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
        width = self.nt_length * 0.1 + 1
        height = min(self.nt_length, 602) * 0.1 + 1
        return (width*self.columns, height*self.rows)

    def plot_profile(self, ax, profile, ct):
        if profile is None:
            return
        column = {"RNP": "normedP",
                  "profile": "Norm_profile",
                  "dance": "Norm_profile"}[profile.datatype]
        factor = {"RNP": 1,
                  "profile": 5,
                  "dance": 5}[profile.datatype]
        am = profile.get_alignment_map(ct)
        values = np.full(ct.length, np.nan)
        colormap = ct.get_colors("profile", profile=profile)
        nts = np.arange(ct.length)+1
        yerr = np.full(ct.length, np.nan)
        for i1, i2 in enumerate(am):
            if i2 != -1:
                values[i2] = profile.data.loc[i1, column]
                if 'Norm_stderr' in profile.data.columns:
                    yerr[i2] = profile.data.loc[i1, 'Norm_stderr']
        mn, mx = self.region
        ax.bar(nts[mn-1:mx], values[mn-1:mx]*factor, align="center",
               width=1.05, color=colormap[mn-1:mx],
               edgecolor=colormap[mn-1:mx], linewidth=0.0,
               yerr=yerr[mn-1:mx], ecolor=(0, 0, 1 / 255.0), capsize=1)
