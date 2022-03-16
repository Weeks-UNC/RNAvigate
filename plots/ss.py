from .plots import Plot
import numpy as np


class SS(Plot):
    def __init__(self, num_samples, structure):
        self.structure = structure
        xmin = min(self.structure.xcoordinates)
        xmax = max(self.structure.xcoordinates)
        xbuffer = (xmax - xmin) * 0.05
        ymin = min(self.structure.ycoordinates)
        ymax = max(self.structure.ycoordinates)
        ybuffer = (ymax-ymin) * 0.05
        super().__init__(num_samples)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set(xlim=[xmin-xbuffer, xmax+xbuffer],
                   ylim=[ymin-3*ybuffer, ymax+ybuffer])
        self.pass_through = ["nt_color"]

    def plot_data(self, ij, profile, label, nt_color="sequence"):
        ax = self.get_ax()
        self.plot_structure(ax)
        self.plot_sequence(ax, profile, nt_color)
        if ij is not None:
            self.plot_ij(ax, ij)
            ax_ins1 = ax.inset_axes([0.15, 0.05, 0.7, 0.05])
            self.view_colormap(ax_ins1, ij)
        ax.set_title(label, fontsize=30)
        self.i += 1

    def get_figsize(self):
        ss = self.structure
        xmin = min(ss.xcoordinates)
        xmax = max(ss.xcoordinates)
        ymin = min(ss.ycoordinates)
        ymax = max(ss.ycoordinates)
        scale = 0.5
        width = (xmax-xmin)*scale
        height = (ymax-ymin)*scale*1.05
        return (width*self.columns, height*self.rows)

    def plot_structure(self, ax):
        ss = self.structure
        for pair in ss.pairList():
            xcoords = [ss.xcoordinates[pair[0]-1], ss.xcoordinates[pair[1]-1]]
            ycoords = [ss.ycoordinates[pair[0]-1], ss.ycoordinates[pair[1]-1]]
            ax.plot(xcoords, ycoords, color="grey",
                    linestyle=(0, (1, 1)), zorder=0)
        ax.plot(ss.xcoordinates, ss.ycoordinates, color="grey", zorder=0)

    def plot_sequence(self, ax, profile, nt_color, markers="o"):
        ss = self.structure
        if isinstance(nt_color, list) and len(nt_color) == self.length["ss"]:
            self.colors = np.array(nt_color)
        elif nt_color == "sequence":
            nt_color = ss.get_colorby_sequence()
        elif nt_color == "position":
            nt_color = ss.get_colorby_position()
        elif nt_color == "profile":
            nt_color = profile.get_colors(ss)
        else:
            print("Invalid colors: choices are profile, sequence, position " +
                  "or a list of length equal to structure sequence. " +
                  "Defaulting to sequence.")
            nt_color = ss.get_colorby_sequence()
        if markers == "sequence":
            for nuc in "GUAC":
                mask = [nt == nuc for nt in ss.sequence]
                xcoords = ss.xcoordinates[mask]
                ycoords = ss.ycoordinates[mask]
                marker = "$"+nuc+"}$"
                marker = markers
                ax.scatter(xcoords, ycoords, marker=marker,
                           c=nt_color[mask])
        else:
            ax.scatter(ss.xcoordinates, ss.ycoordinates, marker=markers,
                       c=nt_color)

    def add_lines(self, ax, i, j, color, linewidth=1.5):
        ss = self.structure
        x = [ss.xcoordinates[i-1], ss.xcoordinates[j-1]]
        y = [ss.ycoordinates[i-1], ss.ycoordinates[j-1]]
        ax.plot(x, y, color=color, linewidth=linewidth)

    def plot_ij(self, ax, ij):
        ij_colors = ij.get_ij_colors()
        window = ij.window
        if window == 1:
            lw = 1.5
        else:
            lw = 6
        for i, j, color in zip(*ij_colors):
            for w in range(window):
                self.add_lines(ax, i+w, j+window-1-w, color, linewidth=lw)

    def plot_positions(self, ax, spacing=20):
        ss = self.structure
        for i in range(0, ss.length, spacing):
            ax.annotate(i+1, xy=(ss.xcoordinates[i], ss.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords="data", fontsize=16,
                        bbox=dict(fc="white", ec="none", pad=0.3))
