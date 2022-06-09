from turtle import color
from .plots import Plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.path import Path
from matplotlib.collections import LineCollection


class SS(Plot):
    def __init__(self, num_samples, structure, **kwargs):
        self.structure = structure
        xmin = min(self.structure.xcoordinates)
        xmax = max(self.structure.xcoordinates)
        xbuffer = 1
        ymin = min(self.structure.ycoordinates)
        ymax = max(self.structure.ycoordinates)
        ybuffer = 1
        super().__init__(num_samples, **kwargs)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set(xlim=[xmin-xbuffer, xmax+xbuffer],
                   ylim=[ymin-3*ybuffer, ymax+ybuffer])
        self.pass_through = ["colors", "sequence", "apply_color_to",
                             "colorbar"]

    def plot_data(self, ij, ij2, profile, label, colors="sequence",
                  sequence=False, apply_color_to="background", colorbar=True):
        ax = self.get_ax()
        self.plot_sequence(ax, profile, colors, sequence, apply_color_to)
        self.plot_ij(ax, ij, colorbar, 0)
        self.plot_ij(ax, ij2, colorbar, 1)
        ax.set_title(label)
        self.i += 1
        if self.i == self.length:
            plt.tight_layout()

    def get_figsize(self):
        ss = self.structure
        xmin = min(ss.xcoordinates)
        xmax = max(ss.xcoordinates)
        ymin = min(ss.ycoordinates)
        ymax = max(ss.ycoordinates)
        scale = 0.55
        width = (xmax-xmin)*scale
        height = (ymax-ymin)*scale
        return (width*self.columns, height*self.rows)

    def plot_structure(self, ax, struct_color):
        ss = self.structure

        x = ss.xcoordinates
        y = ss.ycoordinates

        path = Path(np.column_stack([x, y]))
        verts = path.interpolated(steps=2).vertices
        x, y = verts[:, 0], verts[:, 1]
        colors = [x for c in struct_color for x in [c, c]][1:-1]
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, colors=colors, zorder=0)
        ax.add_collection(lc)

        for pair in ss.pairList():
            xcoords = [ss.xcoordinates[pair[0]-1], ss.xcoordinates[pair[1]-1]]
            ycoords = [ss.ycoordinates[pair[0]-1], ss.ycoordinates[pair[1]-1]]
            ax.plot(xcoords, ycoords, color="grey",
                    linestyle=(0, (1, 1)), zorder=0)

    def plot_sequence(self, ax, profile, colors, sequence, apply_color_to):
        ss = self.structure
        valid_apply = ["background", "sequence", "structure", None]
        message = f"invalid apply_color_to, must be in {valid_apply}"
        assert apply_color_to in valid_apply, message
        if colors is None or apply_color_to is None:
            return
        elif apply_color_to == "structure":
            struct_color = ss.get_colors(colors, profile=profile,
                                         ct=self.structure)
            self.plot_structure(ax, struct_color)
            ax.scatter(ss.xcoordinates, ss.ycoordinates, marker=".",
                       c=struct_color, zorder=1)
            return
        elif apply_color_to == "background":
            bg_color = ss.get_colors(colors, profile=profile,
                                     ct=self.structure)
            self.plot_structure(ax, ss.get_colors("grey"))
            if sequence:
                nt_color = np.full(bg_color.shape, 'k')
                for i, color in enumerate(nt_color):
                    r, g, b = to_rgb(color)
                    if (r*0.299 + g*0.587 + b*0.114) < 200/256:
                        nt_color[i] = 'w'
        elif apply_color_to == "sequence":
            sequence = True
            nt_color = ss.get_colors(colors, profile=profile,
                                     ct=self.structure)
            bg_color = "white"
            self.plot_structure(ax, ss.get_colors('grey'))
        ax.scatter(ss.xcoordinates, ss.ycoordinates, marker="o",
                   c=bg_color, s=256, zorder=1)
        if sequence:
            for nuc in "GUAC":
                mask = [nt == nuc for nt in ss.sequence]
                xcoords = ss.xcoordinates[mask]
                ycoords = ss.ycoordinates[mask]
                marker = "$\mathsf{"+nuc+"}$"
                ax.scatter(xcoords, ycoords, marker=marker, s=100,
                           c=nt_color[mask], lw=1, zorder=2)

    def add_lines(self, ax, i, j, color, linewidth=1.5):
        ss = self.structure
        x = [ss.xcoordinates[i-1], ss.xcoordinates[j-1]]
        y = [ss.ycoordinates[i-1], ss.ycoordinates[j-1]]
        ax.plot(x, y, color=color, linewidth=linewidth)

    def plot_ij(self, ax, ij, colorbar, cmap_pos):
        if ij is None:
            return
        ij_colors = ij.get_ij_colors()
        window = ij.window
        if window == 1:
            lw = 3
        else:
            lw = 6
        for i, j, color in zip(*ij_colors):
            self.add_lines(ax, i, j, color, linewidth=lw)
        if colorbar:
            x, width = [(0, 0.49), (0.51, 0.49)][cmap_pos]
            ax_ins1 = ax.inset_axes([x, 0, width, 0.05])
            self.view_colormap(ax_ins1, ij)

    def plot_positions(self, ax, spacing=20):
        ss = self.structure
        for i in range(0, ss.length, spacing):
            ax.annotate(i+1, xy=(ss.xcoordinates[i], ss.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords="data", fontsize=16,
                        bbox=dict(fc="white", ec="none", pad=0.3))
