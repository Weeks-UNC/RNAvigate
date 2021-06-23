import matplotlib.pyplot as plt
from plots.plots import get_rows_columns


class SS():
    def __init__(self, structures, labels, ijs=None, profiles=None):
        self.structures = structures
        self.ijs = ijs
        self.profiles = profiles
        self.labels = labels

    def get_figsize(self, rows, columns):
        ss = self.structures[0]
        xmin = min(ss.xcoordinates)
        xmax = max(ss.xcoordinates)
        ymin = min(ss.ycoordinates)
        ymax = max(ss.ycoordinates)
        scale = 0.5
        width = (xmax-xmin)*scale
        height = (ymax-ymin)*scale
        return (width*columns, height*rows)

    def set_plot(self, ax, index):
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(self.labels[index], fontsize=30)

    def plot_structure(self, ax, index):
        ss = self.structures[index]
        for pair in ss.pairList():
            xcoords = [ss.xcoordinates[pair[0]-1], ss.xcoordinates[pair[1]-1]]
            ycoords = [ss.ycoordinates[pair[0]-1], ss.ycoordinates[pair[1]-1]]
            ax.plot(xcoords, ycoords, color='grey',
                    linestyle=(0, (1, 1)), zorder=0)
        ax.plot(ss.xcoordinates, ss.ycoordinates, color='grey', zorder=0)

    def plot_sequence(self, ax, index, colors, markers='o'):
        ss = self.structures[index]
        if isinstance(colors, list) and len(colors) == self.length["ss"]:
            self.colors = np.array(colors)
        elif colors == "sequence":
            colors = ss.get_colorby_sequence()
        elif colors == "position":
            colors = ss.get_colorby_position()
        elif colors == "profile":
            colors = self.profiles[index].get_colors(ss)
        else:
            print("Invalid colors: choices are profile, sequence, position " +
                  "or a list of length equal to structure sequence. " +
                  "Defaulting to sequence.")
            colors = ss.get_colorby_sequence()
        if markers == "sequence":
            for nuc in "GUAC":
                mask = [nt == nuc for nt in ss.sequence]
                xcoords = ss.xcoordinates[mask]
                ycoords = ss.ycoordinates[mask]
                marker = "$"+nuc+"}$"
                marker = markers
                ax.scatter(xcoords, ycoords, marker=marker,
                           c=colors[mask])
        else:
            ax.scatter(ss.xcoordinates, ss.ycoordinates, marker=markers,
                       c=colors)

    def add_lines(self, ax, i, j, color, index, linewidth=1.5):
        ss = self.structures[index]
        x = [ss.xcoordinates[i-1], ss.xcoordinates[j-1]]
        y = [ss.ycoordinates[i-1], ss.ycoordinates[j-1]]
        ax.plot(x, y, color=color, linewidth=linewidth)

    def plot_ij(self, ax, index):
        ij = self.ijs[index]
        ij_colors = ij.get_ij_colors()
        window = ij.window
        for i, j, color in zip(*ij_colors):
            if window == 1:
                self.add_lines(ax, i, j, color, index)
            else:
                for w in range(window):
                    self.add_lines(ax, i+w, j+window-1 -
                                   w, color, index, linewidth=6)

    def plot_positions(self, ax, index, spacing=20):
        ss = self.structures[index]
        for i in range(0, ss.length, spacing):
            ax.annotate(i+1, xy=(ss.xcoordinates[i], ss.ycoordinates[i]),
                        horizontalalignment="center",
                        verticalalignment="center",
                        xycoords='data', fontsize=16,
                        bbox=dict(fc="white", ec='none', pad=0.3))

    def make_plot(self, ax=None, positions=True, colors=None, markers='o'):
        length = len(self.labels)
        rows, cols = get_rows_columns(length)
        if ax is None:
            _, axes = plt.subplots(rows, cols, squeeze=False,
                                   figsize=self.get_figsize(rows, cols))
        for i in range(length):
            row = i // cols
            col = i % cols
            ax = axes[row, col]
            self.set_plot(ax, i)
            self.plot_structure(ax, i)
            if colors is None:
                if self.profiles is not None:
                    colors = "profile"
                else:
                    colors = "sequence"
            self.plot_sequence(ax, i, colors=colors, markers=markers)
            if positions is True:
                self.plot_positions(ax, i)
            if self.ijs is not None:
                self.plot_ij(ax, i)
