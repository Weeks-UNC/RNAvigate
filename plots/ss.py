
def get_ss_figsize(self, rows, columns):
    """Returns a tuple of the width and height for a secondary structure
    graph, based on the number of rows and columns. Pass this function call
    to the figsize argument of plt.subplots.

    Args:
        rows (int): number of rows in figure
        columns (int): number of columns in plot

    Returns:
        [type]: [description]
    """
    xmin = min(self.xcoordinates)
    xmax = max(self.xcoordinates)
    ymin = min(self.ycoordinates)
    ymax = max(self.ycoordinates)
    scale = 0.5
    width = (xmax-xmin)*scale
    height = (ymax-ymin)*scale
    return (width*columns, height*rows)


def set_ss(self, ax):
    """Turn off axes and axis border and set aspect ratio to equal

    Args:
        ax (pyplot axis): axis to set
    """
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(self.sample, fontsize=30)


def plot_ss_structure(self, ax):
    """Plots the skeleton of the secondary structure, with base pairs as
    dotted lines

    Args:
        ax (pyplot axis): axis on which to plot
    """
    for pair in self.basepairs:
        xcoords = [self.xcoordinates[pair[0]-1],
                   self.xcoordinates[pair[1]-1]]
        ycoords = [self.ycoordinates[pair[0]-1],
                   self.ycoordinates[pair[1]-1]]
        ax.plot(xcoords, ycoords, color='grey',
                linestyle=(0, (1, 1)), zorder=0)
    ax.plot(self.xcoordinates, self.ycoordinates, color='grey', zorder=0)


def get_ss_seqcolors(self, colors):
    """creates an array of white and black that contrasts well with the
    given color array
    """
    def whiteOrBlack(color):
        """excepts any valid mpl color and returns white or black,
        which ever has higher contrast"""
        try:
            c = mc.cnames[color]
        except KeyError:
            c = color
        if mc.rgb_to_hsv(mc.to_rgb(c))[2] < 0.179:
            return 'white'
        else:
            return 'black'
    letter_colors = np.apply_along_axis(whiteOrBlack, 0, colors)
    return letter_colors


def plot_ss_sequence(self, ax, colorby='sequence', markers='o'):
    """plots the location of nucleotides on top of the secondary structure
    skeleton.

    Args:
        ax (pyplot axis): axis on which to plot
        colorby (str or list, optional): Options are 'sequence', 'profile',
                'position' or a lis of valid matplotlib colors.
                Defaults to None.
        markers (str, optional): Options are a matplotlib marker type or
                'sequence', which uses the sequence letters as markers.
                Defaults to 'o' (a filled circle).
    """
    if isinstance(colorby, list) and len(colorby) == self.length["ss"]:
        self.colors = np.array(colorby)
    try:
        colors = getattr(self, f"get_colorby_{colorby}")("ss")
    except AttributeError:
        print("Invalid colorby: choices are profile, sequence, position" +
              " or a list of length equal to structure sequence, using" +
              " sequence.")
        self.get_colorby_sequence("ss")
    if markers == "sequence":
        for nuc in "GUAC":
            mask = [nt == nuc for nt in self.sequence["ss"]]
            xcoords = self.xcoordinates[mask]
            ycoords = self.ycoordinates[mask]
            marker = "$"+nuc+"}$"
            marker = markers
            ax.scatter(xcoords, ycoords, marker=marker,
                       c=colors[mask])
    else:
        ax.scatter(self.xcoordinates, self.ycoordinates, marker=markers,
                   c=colors)


def add_ss_lines(self, ax, i, j, color, linewidth=1.5):
    xi = self.xcoordinates[i-1]
    yi = self.ycoordinates[i-1]
    xj = self.xcoordinates[j-1]
    yj = self.ycoordinates[j-1]
    ax.plot([xi, xj], [yi, yj], color=color, linewidth=linewidth)


def plot_ss_data(self, ax, ij_data, metric=None, all_pairs=False,
                 **kwargs):
    self.filter_ij_data(ij_data, "ss", all_pairs=all_pairs, **kwargs)
    ij_colors = self.get_ij_colors(ij_data, metric)
    window = self.window[ij_data]
    for i, j, color in zip(*ij_colors):
        if window == 1:
            self.add_ss_lines(ax, i, j, color)
        else:
            for w in range(window):
                self.add_ss_lines(ax, i+w, j+window-1 -
                                  w, color, linewidth=6)


def plot_ss_positions(self, ax, spacing=20):
    """adds position labels to the secondary structure graph

    Args:
        ax (pyplot axis): axis on which to add labels
        spacing (int, optional): labels every nth nucleotide.
                Defaults to 20.
    """
    for i in range(0, self.length["ss"], spacing):
        ax.annotate(i+1, xy=(self.xcoordinates[i], self.ycoordinates[i]),
                    horizontalalignment="center",
                    verticalalignment="center",
                    xycoords='data', fontsize=16,
                    bbox=dict(fc="white", ec='none', pad=0.3))


def make_ss(self, ax=None, ij_data=None, metric=None, positions=True,
            colorby=None, markers='o', all_pairs=False, **kwargs):
    """Creates a full secondary structure graph with data plotted

    Args:
        ax (pyplot axis, option): axis to plot
        ij_data (str, optional): string matching the ij_data to plot.
                Options: rings, pairs, deletions.
        metric (str, optional): string matching a column name of the
                ij_data dataframe. Defaults defined at top of file.
        positions (bool, optional): whether to label positions.
                Defaults to True.
        colorby (str, optional): method used to color nucleotides.
                Options: sequence, profile, positions
                Defaults to 'profile',
                If profile is missing, sequence is used.
        markers (str, optional): marker type to use for nucleotides.
                Defaults to 'o'. (dots)
    """
    if ax is None:
        _, ax = plt.subplots(1, figsize=self.get_ss_figsize(1, 1))
    self.set_ss(ax)
    self.plot_ss_structure(ax)
    if colorby is None:
        if hasattr(self, "profile"):
            colorby = "profile"
        else:
            colorby = "sequence"
    self.plot_ss_sequence(ax, colorby=colorby, markers=markers)
    if positions is True:
        self.plot_ss_positions(ax)
    if metric == "Distance":
        self.set_3d_distances(ij_data)
    if ij_data is not None:
        self.plot_ss_data(ax, ij_data, metric, all_pairs, **kwargs)
    return ax
