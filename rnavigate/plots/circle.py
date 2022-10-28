from matplotlib.collections import PatchCollection
from matplotlib.patches import Path, PathPatch
from matplotlib.colors import to_rgb
from math import sin, cos, pi
from .plots import Plot
import numpy as np


class Circle(Plot):

    def __init__(self, num_samples, seq_source):
        self.sequence = seq_source
        length = self.sequence.length
        super().__init__(num_samples)
        self.x, self.y = np.zeros(length), np.zeros(length)
        self.diameter = length / pi
        self.theta = np.array([2*pi * (i+4)/(length+8) for i in range(length)])
        self.x = np.sin(self.theta)*self.diameter
        self.y = np.cos(self.theta)*self.diameter
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.set_aspect('equal')
            ax.axis('off')
            # ax.set(xlim=(-10.1, 10.1),
            #        ylim=(-10.1, 10.1))
        self.pass_through = ["colors", "apply_color_to", "sequence", "colorbar",
                             "title", "positions"]
        self.zorder = {"annotations": 0,
                       "data": 5,
                       "nucleotide": 10,
                       "sequence": 15,
                       "position": 20}

    def get_figsize(self):
        dim = self.sequence.length / pi / 4
        return (dim * self.columns, dim * self.rows)

    def plot_data(self, ct, comp, interactions, interactions2, profile,
                  annotations, label,
                  colors="sequence", apply_color_to="sequence", colorbar=True,
                  title=True, positions=True):
        ax = self.get_ax()
        if interactions is not None and colorbar:
            ax_ins1 = ax.inset_axes(
                [-5, -0.4, 10, 0.8], transform=ax.transData)
            self.view_colormap(ax_ins1, interactions)
        if interactions2 is not None and colorbar:
            ax_ins2 = ax.inset_axes([10, -100, 100, 8], transform=ax.transData)
            self.view_colormap(ax_ins2, interactions2)
        if comp is not None and colorbar:
            ax_ins3 = ax.inset_axes([10, 80, 100, 8], transform=ax.transData)
            self.view_colormap(ax_ins3, "ct_compare")
        self.add_patches(ax, ct, comp)
        self.add_patches(ax, interactions)
        self.add_patches(ax, interactions2)
        # self.add_sequence(ax, ct.sequence, yvalue=0.5)
        # self.plot_profile(ax, profile, ct)
        self.plot_sequence(ax, profile, colors, apply_color_to)
        for annotation in annotations:
            self.plot_annotation(ax, annotation)
        if title:
            ax.set_title(label)
        if positions:
            self.plot_positions(ax)
        self.i += 1

    def plot_sequence(self, ax, profile, colors, apply_color_to):
        sequence = self.sequence
        nuc_z = self.zorder["nucleotide"]
        seq_z = self.zorder["sequence"]
        valid_apply = ["background", "sequence", None]
        message = f"invalid apply_color_to, must be in {valid_apply}"
        assert apply_color_to in valid_apply, message
        if colors is None or apply_color_to is None:
            colors = sequence.get_colors("black")
            apply_color_to = "sequence"
        if apply_color_to == "background":
            bg_color = sequence.get_colors(colors, profile=profile)
        if apply_color_to == "sequence":
            nt_color = sequence.get_colors(colors, profile=profile)
            bg_color = sequence.get_colors("white")
        elif sequence:
            nt_color = ['k'] * len(bg_color)
            for i, color in enumerate(bg_color):
                r, g, b = to_rgb(color)
                if (r*0.299 + g*0.587 + b*0.114) < 175/256:
                    nt_color[i] = 'w'
            nt_color = np.array(nt_color)
        ax.scatter(self.x, self.y, marker="o", c=bg_color, s=256, zorder=nuc_z)
        for nuc in "GUACguac":
            mask = [nt == nuc for nt in sequence.sequence]
            xcoords = self.x[mask]
            ycoords = self.y[mask]
            marker = "$\mathsf{"+nuc+"}$"
            ax.scatter(xcoords, ycoords, marker=marker, s=100,
                       c=nt_color[mask], lw=1, zorder=seq_z)

    def plot_positions(self, ax, interval=20):
        for i in np.arange(interval, self.sequence.length, interval):
            x = self.x[i]
            y = self.y[i]
            theta = self.theta[i] * -180/np.pi
            ax.plot([x, x*1.06], [y, y*1.06], c='k')
            ha = ['left', 'right'][x < 0]
            va = ['top', 'bottom'][y > 0]
            ax.text(x*1.08, y*1.08, i, ha='center',
                    va='center', rotation=theta)

    def add_patches(self, ax, data, comp=None):
        if comp is not None:
            ij_colors = data.get_ij_colors(comp)
        elif data is not None:
            ij_colors = data.get_ij_colors()
        else:
            return
        patches = []
        for i, j, color in zip(*ij_colors):
            if j < i:  # flip the order
                i, j = j, i
            x_i = self.x[i-1]
            y_i = self.y[i-1]
            x_j = self.x[j-1]
            y_j = self.y[j-1]
            x_center = (x_i+x_j)/2
            y_center = (y_i+y_j)/2
            # scaling the center point towards zero depending on angle(i,j)
            # I should test taking the root of the center point
            f = (1 - ((j - i) / (self.sequence.length + 8))) ** 4
            verts = [[x_i, y_i], [x_center*f, y_center*f], [x_j, y_j]]
            codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
            patches.append(PathPatch(Path(verts, codes), fc="none", ec=color))
        ax.add_collection(PatchCollection(patches, match_original=True))

    def plot_annotation(self, ax, annotation):
        color = annotation.color
        zorder = self.zorder["annotations"]
        if annotation.annotation_type == "spans":
            for start, end in annotation.spans:
                x = self.x[start-1:end]*1.03
                y = self.y[start-1:end]*1.03
                ax.plot(x, y, color=color, alpha=0.8, lw=10, zorder=zorder)
        elif annotation.annotation_type == "sites":
            x = self.x[annotation.sites]
            y = self.y[annotation.sites]
            ax.scatter(x, y, color=color, marker='o', ec="none", alpha=0.7,
                       s=30**2, zorder=zorder)
        elif annotation.annotation_type == "groups":
            for group in annotation.groups:
                x = self.structure.xcoordinates[group["sites"]]
                y = self.structure.ycoordinates[group["sites"]]
                color = group["color"]
                ax.plot(x, y, color=color, alpha=0.2, lw=30, zorder=zorder)
                ax.scatter(x, y, color=color, marker='o', ec="none", alpha=0.4,
                           s=30**2, zorder=zorder)
        elif annotation.annotation_type == "primers":
            for start, end in annotation.primers:
                if start < end:
                    index = np.arange(start-1, end)
                elif start > end:
                    index = np.arange(start-1, end-2, -1)
                index = np.append(index, index[-2])
                x = self.x[index] * 1.03
                y = self.y[index] * 1.03
                x[-1] *= 1.03
                y[-1] *= 1.03
                ax.plot(x, y, color=color, alpha=0.8, zorder=zorder)
