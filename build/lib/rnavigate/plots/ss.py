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
                             "colorbar", "title", "positions", "bp_style"]
        self.zorder = {"structure": 0,
                       "backbone": 0,
                       "annotations": 5,
                       "nucleotide": 10,
                       "data": 15,
                       "sequence": 20,
                       "position": 25}

    def plot_data(self, interactions, interactions2, profile, annotations,
                  label,
                  colors="sequence",
                  sequence=False,
                  apply_color_to="background",
                  colorbar=True,
                  title=True,
                  positions=False,
                  bp_style="dotted"):
        ax = self.get_ax()
        self.plot_sequence(ax=ax, profile=profile, colors=colors,
                           sequence=sequence, apply_color_to=apply_color_to,
                           positions=positions, bp_style=bp_style)
        self.plot_interactions(
            ax=ax, interactions=interactions, colorbar=colorbar, cmap_pos=0)
        self.plot_interactions(
            ax=ax, interactions=interactions2, colorbar=colorbar, cmap_pos=1)
        for annotation in annotations:
            self.plot_annotation(ax=ax, annotation=annotation)
        if title:
            ax.set_title(label)
        self.i += 1
        if self.i == self.length:
            plt.tight_layout()

    def get_figsize(self):
        ss = self.structure
        xmin = min(ss.xcoordinates)-1
        xmax = max(ss.xcoordinates)+1
        ymin = min(ss.ycoordinates)-1
        ymax = max(ss.ycoordinates)+3
        scale = 0.55
        width = (xmax-xmin)*scale
        height = (ymax-ymin)*scale
        return (width*self.columns, height*self.rows)

    def plot_structure(self, ax, struct_color, bp_style):
        bp_styles = ["conventional", "dotted", "line"]
        assert bp_style in bp_styles, f"bp_style must be one of {bp_styles}"

        ss = self.structure

        x = ss.xcoordinates
        y = ss.ycoordinates

        path = Path(np.column_stack([x, y]))
        verts = path.interpolated(steps=2).vertices
        x, y = verts[:, 0], verts[:, 1]
        colors = [x for c in struct_color for x in [c, c]][1:-1]
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        zorder = self.zorder["structure"]
        lc = LineCollection(segments, colors=colors, zorder=zorder)
        ax.add_collection(lc)

        for pair in ss.pairList():
            x = ss.xcoordinates[[p-1 for p in pair]]
            y = ss.ycoordinates[[p-1 for p in pair]]
            xdist = x[1]-x[0]
            ydist = y[1]-y[0]
            if bp_style == "dotted":
                x_dots = [x[0] + i * xdist / 6 for i in [2, 3, 4]]
                y_dots = [y[0] + i * ydist / 6 for i in [2, 3, 4]]
                ax.scatter(x_dots, y_dots, c="grey", marker='.', s=4,
                           zorder=zorder)
            if bp_style == "line":
                x_caps = [x[0] + i * xdist / 3 for i in [1, 2]]
                y_caps = [y[0] + i * ydist / 3 for i in [1, 2]]
                ax.plot(x_caps, y_caps, color="grey", zorder=zorder)
            if bp_style == "conventional":
                nts = ''.join([ss.sequence[p-1] for p in pair]).upper()
                x_caps = [x[0] + i * xdist / 7 for i in [2, 5]]
                y_caps = [y[0] + i * ydist / 7 for i in [2, 5]]
                if nts in ["UA", "AU", "GU", "UG"]:
                    ax.plot(x_caps, y_caps, color="grey", zorder=zorder,
                            solid_capstyle='butt')
                if nts in ["GU", "UG", "AG", "GA"]:
                    x_center = x[0] + xdist/2
                    y_center = y[0] + ydist/2
                    if nts in ["GU", "UG"]:
                        fc = 'white'
                        ec = 'grey'
                    elif nts in ["AG", "GA"]:
                        fc = 'grey'
                        ec = 'none'
                    ax.scatter(x_center, y_center, zorder=zorder+1, s=36,
                               color="grey", marker="o", fc=fc, ec=ec)
                if nts in ["GC", "CG"]:
                    ax.plot(x_caps, y_caps, color="grey", zorder=zorder,
                            linewidth=6, solid_capstyle='butt')
                    ax.plot(x_caps, y_caps, color="white", zorder=zorder,
                            linewidth=2)

    def plot_sequence(self, ax, profile, colors, sequence, apply_color_to,
                      positions, bp_style):
        ss = self.structure
        nuc_z = self.zorder["nucleotide"]
        seq_z = self.zorder["sequence"]
        valid_apply = ["background", "sequence", "structure", None]
        message = f"invalid apply_color_to, must be in {valid_apply}"
        assert apply_color_to in valid_apply, message
        if colors is None or apply_color_to is None:
            colors = ss.get_colors("gray")
            apply_color_to = "structure"
        if apply_color_to == "structure":
            struct_color = ss.get_colors(colors, profile=profile,
                                         ct=self.structure)
            self.plot_structure(ax, struct_color, bp_style=bp_style)
            ax.scatter(ss.xcoordinates, ss.ycoordinates, marker=".",
                       c=struct_color, zorder=seq_z)
            return
        if apply_color_to == "background":
            bg_color = ss.get_colors(colors, profile=profile,
                                     ct=self.structure)
            self.plot_structure(ax, ss.get_colors("grey"), bp_style=bp_style)
        if apply_color_to == "sequence":
            sequence = True
            nt_color = ss.get_colors(colors, profile=profile,
                                     ct=self.structure)
            bg_color = ss.get_colors("white")
            self.plot_structure(ax, ss.get_colors('grey'), bp_style=bp_style)
        elif sequence:
            nt_color = ['k'] * len(bg_color)
            for i, color in enumerate(bg_color):
                r, g, b = to_rgb(color)
                if (r*0.299 + g*0.587 + b*0.114) < 175/256:
                    nt_color[i] = 'w'
            nt_color = np.array(nt_color)

        if positions:
            self.plot_positions(ax, text_color=nt_color, bbox_color=bg_color)
        ax.scatter(ss.xcoordinates, ss.ycoordinates, marker="o",
                   c=bg_color, s=256, zorder=nuc_z)
        if sequence:
            for nuc in "GUACguac":
                mask = [nt == nuc for nt in ss.sequence]
                xcoords = ss.xcoordinates[mask]
                ycoords = ss.ycoordinates[mask]
                marker = "$\mathsf{"+nuc+"}$"
                ax.scatter(xcoords, ycoords, marker=marker, s=100,
                           c=nt_color[mask], lw=1, zorder=seq_z)

    def add_lines(self, ax, i, j, color, linewidth=1.5):
        ss = self.structure
        zorder = self.zorder["data"]
        x = [ss.xcoordinates[i-1], ss.xcoordinates[j-1]]
        y = [ss.ycoordinates[i-1], ss.ycoordinates[j-1]]
        ax.plot(x, y, color=color, linewidth=linewidth, zorder=zorder)

    def plot_interactions(self, ax, interactions, colorbar, cmap_pos):
        if interactions is None:
            return
        ij_colors = interactions.get_ij_colors()
        window = interactions.window
        if window == 1:
            lw = 3
        else:
            lw = 6
        for i, j, color in zip(*ij_colors):
            self.add_lines(ax, i, j, color, linewidth=lw)
        if colorbar:
            x, width = [(0, 0.49), (0.51, 0.49)][cmap_pos]
            ax_ins1 = ax.inset_axes([x, 0, width, 0.05])
            self.view_colormap(ax_ins1, interactions)

    def plot_positions(self, ax, text_color, bbox_color, spacing=20):
        ss = self.structure
        zorder = self.zorder["position"]
        for i in range(spacing, ss.length, spacing):
            ax.annotate(f"{ss.sequence[i-1]}\n{i}",
                        xy=(ss.xcoordinates[i-1], ss.ycoordinates[i-1]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        color=text_color[i-1],
                        xycoords="data", fontsize=12, zorder=zorder,
                        bbox=dict(boxstyle="Circle", pad=0.1, ec="none",
                                  fc=bbox_color[i-1]))

    def plot_annotation(self, ax, annotation):
        color = annotation.color
        zorder = self.zorder["annotations"]
        if annotation.annotation_type == "spans":
            for start, end in annotation.spans:
                x = self.structure.xcoordinates[start-1:end]
                y = self.structure.ycoordinates[start-1:end]
                ax.plot(x, y, color=color, alpha=0.2, lw=30, zorder=zorder)
        elif annotation.annotation_type == "sites":
            x = self.structure.xcoordinates[annotation.sites]
            y = self.structure.ycoordinates[annotation.sites]
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
