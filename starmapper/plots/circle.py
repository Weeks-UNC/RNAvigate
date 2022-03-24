from matplotlib.collections import PatchCollection
from matplotlib.patches import Path, PathPatch
from math import sin, cos, pi
from .plots import Plot


class Circle(Plot):

    def __init__(self, num_samples, nt_length):
        self.nt_length = nt_length
        super().__init__(num_samples)
        self.x, self.y = [], []
        diameter = self.nt_length / pi
        for i in range(nt_length):
            theta = 2*pi * (i+4)/(nt_length+8)
            self.x.append(sin(theta)*10)
            self.y.append(cos(theta)*10)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.set_aspect('equal')
            ax.axis('off')
            ax.set(xlim=(-10, 10),
                   ylim=(-10, 10))

    def get_figsize(self):
        dim = self.nt_length / pi / 4
        return (dim * self.columns, dim * self.rows)

    def plot_data(self, ct, comp, ij, ij2, profile, label):
        ax = self.get_ax()
        if ij is not None:
            ax_ins1 = ax.inset_axes(
                [-5, -0.4, 10, 0.8], transform=ax.transData)
            self.view_colormap(ax_ins1, ij)
        if ij2 is not None:
            ax_ins2 = ax.inset_axes([10, -100, 100, 8], transform=ax.transData)
            self.view_colormap(ax_ins2, ij2)
        if comp is not None:
            ax_ins3 = ax.inset_axes([10, 80, 100, 8], transform=ax.transData)
            self.view_colormap(ax_ins3, "ct_compare")
        self.add_patches(ax, ct, comp)
        self.add_patches(ax, ij)
        self.add_patches(ax, ij2)
        # self.add_sequence(ax, ct.sequence, yvalue=0.5)
        # self.plot_profile(ax, profile, ct)
        ax.set_title(label)
        nuc_colors = profile.get_colors(profile)
        ax.scatter(self.x, self.y, marker='o', c=nuc_colors)
        self.i += 1

    def add_patches(self, ax, data, comp=None):
        if comp is not None:
            ij_colors = data.get_ij_colors(comp)
        elif data is not None:
            ij_colors = data.get_ij_colors()
        else:
            return
        patches = []
        for i, j, color in zip(*ij_colors):
            x_i = self.x[i-1]
            y_i = self.y[i-1]
            x_j = self.x[j-1]
            y_j = self.y[j-1]
            x_center = (x_i+x_j)/2
            y_center = (y_i+y_j)/2
            # scaling the center point towards zero depending on angle(i,j)
            # I should test taking the root of the center point
            f = (1 - ((j - i) / (self.nt_length + 8))) ** 4
            verts = [[x_i, y_i], [x_center*f, y_center*f], [x_j, y_j]]
            codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
            patches.append(PathPatch(Path(verts, codes), fc="none", ec=color))
        ax.add_collection(PatchCollection(patches, match_original=True))
