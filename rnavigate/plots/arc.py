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
        self.pass_through = ["ax", "colorbar", "seqbar", "title", "ij_panel",
                             "ij2_panel", "ct_panel", "annotations"]

    def set_axis(self, ax, annotations=0, xticks=20, xticks_minor=10):
        def get_ticks(x, mn, mx):
            return [tick for tick in range(x, mx+1, x) if mn <= tick <= mx]

        ax.set_aspect('equal')
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set(position=('data', -2*annotations),
                                visible=False)
        ax.spines['top'].set_color('none')
        height = min(300, self.nt_length/2)
        mn, mx = self.region
        ax.set(xlim=(mn - 0.5, mx + 0.5),
               ylim=(-height-1-(2*annotations), height+1),
               xticks=get_ticks(xticks, mn, mx),
               axisbelow=False)
        ax.set_xticks(get_ticks(xticks_minor, mn, mx), minor=True)
        for label in ax.get_xticklabels():
            label.set_bbox({"facecolor": "white",
                            "edgecolor": "None",
                            "alpha": 0.5,
                            "boxstyle": "round,pad=0.1,rounding_size=0.2"})

    def plot_data(self, ct, comp, ij, ij2, profile, label, ax=None,
                  colorbar=True, seqbar=True, title=True, ij_panel="bottom",
                  ij2_panel="bottom", ct_panel="top", annotations=[]):
        ax = self.get_ax(ax)
        annotations = len(annotations) + seqbar
        self.set_axis(ax=ax, annotations=annotations)
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
        self.add_patches(ax=ax, data=ct, panel=ct_panel,
                         annotations=annotations, comp=comp)
        self.add_patches(ax=ax, data=ij, panel=ij_panel,
                         annotations=annotations)
        self.add_patches(ax=ax, data=ij2, panel=ij2_panel,
                         annotations=annotations)
        self.plot_profile(ax, profile, ct)
        if seqbar:
            self.add_sequence(ax, ct.sequence, yvalue=1-(2*annotations),
                              ytrans="data")
        if title:
            self.add_title(ax, label)
        self.i += 1

    def add_patches(self, ax, data, panel, annotations, comp=None):
        if comp is not None:
            ij_colors = data.get_ij_colors(comp)
        elif data is not None:
            ij_colors = data.get_ij_colors()
        else:
            return
        patches = []
        mn, mx = self.region
        for i, j, color in zip(*ij_colors):
            if j < i:  # flip the order
                i, j = j, i
            if ((i < mn) and (j < mn)) or ((i > mx) and (j > mx)):
                continue
            if panel == "top":
                center = ((i+j)/2., 0)
                theta1 = 0
                theta2 = 180
            elif panel == "bottom":
                center = ((i+j)/2., -2*annotations)
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
        column = profile.default_column
        factor = profile.ap_scale_factor
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
