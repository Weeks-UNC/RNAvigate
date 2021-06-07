

class skyline():

    def __init__(self):
        self.profiles = []
        self.labels = []

    def add_profile(self, profile, label):
        self.profiles.append(profile)
        self.labels.append(label)

    def set_figsize(self, profile):
        left_inches = 0.9
        right_inches = 0.4
        ax_width = profile.length * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        self.figsize = (fig_width, fig_height)

    def set_axis(self, profile, ax=None):
        if ax is None:
            _, self.ax = plt.subplots(1, figsize=self.figsize)
        else:
            self.ax = ax
        # x axis appearance
        ax.set_xlim([0, profile.length])
        ax.set_xticks(range(0, profile.length, 20))
        ax.set_xticks(range(0, profile.length, 5), minor=True)

    def set_labels(self, axis_title="Raw Reactivity Profile",
                   legend_title="Samples", xlabel="Nucleotide",
                   ylabel="Profile"):
        ax.set_title(axis_title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=legend_title)

    def plot_profile(self, profile, label, column='Reactivity_profile'):
        x = [0.5]
        y = [0]
        # converts standard plot to skyline plot.
        for n, r in zip(profile.data['Nucleotide'], profile.data[column]):
            x.extend([n - 0.5, n + 0.5])
            y.extend([r, r])
        x.append(x[-1])
        y.append(0)
        ax.plot(x, y, label=label)

    def add_sequence(self, sequence, yvalue=0.005):
        # set font style and colors for each nucleotide
        font_prop = mp.font_manager.FontProperties(
            family="monospace", style="normal", weight="bold", size="12")
        color_dict = {"A": "#f20000", "U": "#f28f00",
                      "G": "#00509d", "C": "#00c200"}
        # transform yvalue to a y-axis data value
        ymin, ymax = axis.get_ylim()
        yvalue = (ymax-ymin)*yvalue + ymin
        for i, seq in enumerate(sequence):
            col = color_dict[seq.upper()]
            axis.annotate(seq, xy=(i + 1, yvalue), xycoords='data',
                          fontproperties=font_prop,
                          color=col, horizontalalignment="center")

    def make_plot(self, ax=None, column, **kwargs):
        self.set_figsize(self.profile[0])
        self.set_axis(self.profile[0], ax)
        self.add_sequence(self.profile[0].sequence)
        for profile, label in zip(self.profiles, self.labels):
            self.plot_profile(profile, label, column)
        self.set_labels(**kwargs)
