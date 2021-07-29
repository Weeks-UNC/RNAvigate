from .plots import Plot


class Skyline(Plot):
    def __init__(self, nt_length):
        self.nt_length = nt_length
        super().__init__(1)
        self.set_axis()
        self.set_labels()
        self.ax = self.axes[0, 0]

    def plot_data(self, profile, label, column="Reactivity_profile"):
        self.plot_profile(profile, label, column)
        self.i += 1
        if self.i == self.length:
            self.add_sequence(self.ax, profile.sequence)

    def get_figsize(self):
        left_inches = 0.9
        right_inches = 0.4
        ax_width = self.nt_length * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        return (fig_width, fig_height)

    def set_axis(self, xlim=None, xticks=20, xticks_minor=5):
        if xlim is None:
            xlim = [0, self.nt_length]
        self.ax.set_xlim(xlim)
        self.ax.set_xticks(range(xlim[0], xlim[1], xticks))
        self.ax.set_xticks(range(xlim[0], xlim[1], xticks_minor), minor=True)

    def set_labels(self, axis_title="Raw Reactivity Profile",
                   legend_title="Samples", xlabel="Nucleotide",
                   ylabel="Profile"):
        self.ax.set_title(axis_title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.legend(title=legend_title)

    def plot_profile(self, profile, label, column):
        x = [0.5]
        y = [0]
        # converts standard plot to skyline plot.
        for n, r in zip(profile.data['Nucleotide'], profile.data[column]):
            x.extend([n - 0.5, n + 0.5])
            y.extend([r, r])
        x.append(x[-1])
        y.append(0)
        self.ax.plot(x, y, label=label)
