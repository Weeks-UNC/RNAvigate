from .plots import Plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.path import Path
from matplotlib.collections import LineCollection


class SS(Plot):
    def __init__(self, num_samples, structures, **kwargs):
        self.plot_params = {
            "structure_lw": 3,
            "structure_s": 10**2,

            "data_lw": 1.5,
            "data_a": None,
            "data_z": 15,

            "annotations_lw": 30,
            "annotations_s": 30**2,
            "annotations_a": 0.4,

            "basepair_z": 0,
            "structure_z": 2,
            "annotations_z": 5,
            "nucleotide_z": 10,
            "sequence_z": 20,
            "position_z": 25,
        }

        for kw in list(kwargs.keys()):
            if kw in self.plot_params:
                self.plot_params[kw] = kwargs.pop(kw)
        # ensure that all structures fit within axes limits
        self.xmin = min([min(ss.xcoordinates) for ss in structures])
        self.xmax = max([max(ss.xcoordinates) for ss in structures])
        self.ymin = min([min(ss.ycoordinates) for ss in structures])
        self.ymax = max([max(ss.ycoordinates) for ss in structures])
        super().__init__(num_samples, **kwargs)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set(xlim=[self.xmin-1, self.xmax+1],
                   ylim=[self.ymin-1, self.ymax+1])
        self.pass_through = ["colors", "sequence", "apply_color_to",
                             "colorbar", "title", "positions", "bp_style"]

    def plot_data(self, structure, interactions, interactions2, profile,
                  annotations, label,
                  colors="sequence",
                  sequence=False,
                  apply_color_to="background",
                  colorbar=True,
                  title=True,
                  positions=False,
                  bp_style="dotted"):
        profile = profile.fitted
        annotations = [annotation.fitted for annotation in annotations]
        ax = self.get_ax()
        self.plot_sequence(ax=ax, ss=structure, profile=profile, colors=colors,
                           sequence=sequence, apply_color_to=apply_color_to,
                           positions=positions, bp_style=bp_style)
        self.plot_interactions(ax=ax, ss=structure, interactions=interactions,
                               colorbar=colorbar, cmap_pos=0)
        self.plot_interactions(ax=ax, ss=structure, interactions=interactions2,
                               colorbar=colorbar, cmap_pos=1)
        for annotation in annotations:
            self.plot_annotation(ax=ax, ss=structure, annotation=annotation)
        if title:
            ax.set_title(label)
        self.i += 1
        if self.i == self.length:
            plt.tight_layout()

    def get_figsize(self):
        scale = 0.55
        width = (self.xmax-self.xmin)*scale
        height = (self.ymax-self.ymin)*scale
        return (width*self.columns, height*self.rows)

    def plot_structure(self, ax, ss, struct_color, bp_style):
        bp_styles = ["conventional", "dotted", "line"]
        assert bp_style in bp_styles, f"bp_style must be one of {bp_styles}"

        x = ss.xcoordinates
        y = ss.ycoordinates

        path = Path(np.column_stack([x, y]))
        verts = path.interpolated(steps=2).vertices
        x, y = verts[:, 0], verts[:, 1]
        colors = [x for c in struct_color for x in [c, c]][1:-1]
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        current_color = None
        new_segments = []
        new_colors = []
        for segment, color in zip(segments, colors):
            if current_color is None:
                current_color = color
                current_segment = list(segment)
            elif current_color == color:
                current_segment.append(segment[-1])
            else:
                new_segments.append(current_segment)
                new_colors.append(current_color)
                current_color = color
                current_segment = list(segment)
        new_segments.append(current_segment)
        new_colors.append(current_color)

        zorder = self.plot_params["structure_z"]
        lc = LineCollection(new_segments, colors=new_colors, zorder=zorder,
                            lw=self.plot_params["structure_lw"])
        ax.add_collection(lc)

        zorder = self.plot_params["basepair_z"]
        for pair in ss.pairList():
            x = ss.xcoordinates[[p-1 for p in pair]]
            y = ss.ycoordinates[[p-1 for p in pair]]
            xdist = x[1]-x[0]
            ydist = y[1]-y[0]
            if xdist != 0:
                angle_xy = np.arctan(ydist/xdist)
            elif ydist > 0:
                angle_xy = 1/2 * np.pi
            elif ydist < 0:
                angle_xy = 3/2 * np.pi
            if (xdist < 0):
                angle_xy += np.pi
            x_offset = np.cos(angle_xy) / 3
            y_offset = np.sin(angle_xy) / 3
            x_caps = [x[0] + x_offset, x[1] - x_offset]
            y_caps = [y[0] + y_offset, y[1] - y_offset]
            if bp_style == "dotted":
                x_caps_dist = x_caps[0]-x_caps[1]
                y_caps_dist = y_caps[0]-y_caps[1]
                caps_dist = (x_caps_dist**2 + y_caps_dist**2)**0.5
                segs = int(max([1/3, caps_dist]) * 3) * 2
                x_dots = [x_caps[0]-x_caps_dist*i/segs for i in range(segs+1)]
                y_dots = [y_caps[0]-y_caps_dist*i/segs for i in range(segs+1)]
                ax.scatter(x_dots, y_dots, c="grey", marker='.', s=4,
                           zorder=zorder)
            if bp_style == "line":
                ax.plot(x_caps, y_caps, color="grey", zorder=zorder)
            if bp_style == "conventional":
                nts = ''.join([ss.sequence[p-1] for p in pair]).upper()
                # x_caps = [x[0] + i * xdist / 7 for i in [2, 5]]
                # y_caps = [y[0] + i * ydist / 7 for i in [2, 5]]
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

    def plot_sequence(self, ax, ss, profile, colors, sequence,
                      apply_color_to, positions, bp_style):
        nuc_z = self.plot_params["nucleotide_z"]
        seq_z = self.plot_params["sequence_z"]
        valid_apply = ["background", "sequence", "structure", None]
        message = f"invalid apply_color_to, must be in {valid_apply}"
        assert apply_color_to in valid_apply, message
        if colors is None or apply_color_to is None:
            colors = ss.get_colors("gray")
            apply_color_to = "structure"
        if apply_color_to == "structure":
            struct_color = ss.get_colors(colors, profile=profile,
                                         ct=ss)
            self.plot_structure(ax=ax, ss=ss, struct_color=struct_color,
                                bp_style=bp_style)
            ax.scatter(ss.xcoordinates, ss.ycoordinates, marker=".",
                       c=struct_color, zorder=seq_z,
                       s=self.plot_params["structure_s"])
            return
        if apply_color_to == "background":
            bg_color = ss.get_colors(colors, profile=profile,
                                     ct=ss)
            self.plot_structure(ax=ax, ss=ss,
                                struct_color=ss.get_colors("grey"),
                                bp_style=bp_style)
        if apply_color_to == "sequence":
            sequence = True
            nt_color = ss.get_colors(colors, profile=profile,
                                     ct=ss)
            bg_color = ss.get_colors("white")
            self.plot_structure(ax=ax, ss=ss,
                                struct_color=ss.get_colors('grey'),
                                bp_style=bp_style)
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

    def add_lines(self, ax, ss, i, j, color):
        zorder = self.plot_params["data_z"]
        alpha = self.plot_params["data_a"]
        x = [ss.xcoordinates[i-1], ss.xcoordinates[j-1]]
        y = [ss.ycoordinates[i-1], ss.ycoordinates[j-1]]
        ax.plot(x, y, color=color, lw=self.plot_params[
            "data_lw"], zorder=zorder, alpha=alpha)

    def plot_interactions(self, ax, ss, interactions, colorbar,
                          cmap_pos):
        if interactions is None:
            return
        ij_colors = interactions.get_ij_colors()
        for i, j, color in zip(*ij_colors):
            self.add_lines(ax=ax, ss=ss, i=i, j=j, color=color)
        if colorbar:
            x, width = [(0, 0.49), (0.51, 0.49)][cmap_pos]
            ax_ins1 = ax.inset_axes([x, 0, width, 0.05])
            self.view_colormap(ax_ins1, interactions)

    def plot_positions(self, ax, ss, text_color, bbox_color, spacing=20):
        zorder = self.plot_params["position_z"]
        for i in range(spacing, ss.length, spacing):
            ax.annotate(f"{ss.sequence[i-1]}\n{i}",
                        xy=(ss.xcoordinates[i-1], ss.ycoordinates[i-1]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        color=text_color[i-1],
                        xycoords="data", fontsize=12, zorder=zorder,
                        bbox=dict(boxstyle="Circle", pad=0.1, ec="none",
                                  fc=bbox_color[i-1]))

    def plot_annotation(self, ax, ss, annotation):
        annotation = annotation.fitted
        color = annotation.color
        alpha = self.plot_params["annotations_a"]
        size = self.plot_params["annotations_s"]
        linewidth = self.plot_params["annotations_lw"]
        zorder = self.plot_params["annotations_z"]
        if annotation.annotation_type in ["spans", "primers"]:
            for start, end in annotation[:]:
                if start > end:
                    start, end = end, start
                x = ss.xcoordinates[start-1:end]
                y = ss.ycoordinates[start-1:end]
                ax.plot(x, y, color=color, alpha=alpha, lw=linewidth,
                        zorder=zorder)
        elif annotation.annotation_type == "sites":
            sites = np.array(annotation[:])-1
            x = ss.xcoordinates[sites]
            y = ss.ycoordinates[sites]
            ax.scatter(x, y, color=color, marker='o', ec="none", alpha=alpha,
                       s=size, zorder=zorder)
        elif annotation.annotation_type == "groups":
            for group in annotation[:]:
                x = ss.xcoordinates[group["sites"]]
                y = ss.ycoordinates[group["sites"]]
                color = group["color"]
                ax.plot(x, y, color=color, a=alpha,
                        lw=linewidth, zorder=zorder)
                ax.scatter(x, y, color=color, marker='o', ec="none", a=alpha,
                           s=size, zorder=zorder)
