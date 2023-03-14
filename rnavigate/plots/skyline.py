from .plots import Plot, adjust_spines


class Skyline(Plot):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples, **kwargs)
        self.ax = self.axes[0, 0]
        self.pass_through = ["columns", "seqbar", "errorbars"]

    def get_rows_columns(self, number_of_samples=None, rows=None, cols=None):
        return (1, 1)

    def plot_data(self, profile, annotations, label,
                  columns="Reactivity_profile", seqbar=True, errorbars=None):
        profile = profile.fitted
        annotations = [annotation.fitted for annotation in annotations]
        self.plot_profile(profile, label, columns, errorbars)
        self.i += 1
        if self.i == self.length:
            if seqbar:
                self.add_sequence(self.ax, profile.sequence)
            if isinstance(columns, list):
                ylabel = [column.replace("_", " ") for column in columns]
                ylabel = ', '.join(ylabel)
            else:
                ylabel = columns.replace("_", " ")
            self.set_labels(ax=self.ax, ylabel=ylabel, axis_title=None,
                            legend_title=None)
            for annotation in annotations:
                self.plot_annotation(ax=self.ax, annotation=annotation)
            self.set_axis(self.ax)

    def get_figsize(self):
        left_inches = 0.9
        right_inches = 0.4
        ax_width = self.nt_length * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        return (fig_width, fig_height)

    def set_axis(self, ax, xticks=20, xticks_minor=5):
        xlim = self.region
        ax.set_xlim([xlim[0] - 0.5, xlim[1] + 0.5])
        xrange = range(xlim[0], xlim[1]+1)
        ax.set_xticks([x for x in xrange if (x % xticks) == 0])
        ax.set_xticks([x for x in xrange if (x % xticks_minor) == 0],
                      minor=True)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(axis='y', visible=True)
        ax.tick_params(axis='y', length=0, grid_alpha=0.4)
        adjust_spines(ax, ["left", "bottom"])

    def set_labels(self, ax, axis_title="Raw Reactivity Profile",
                   legend_title="Samples", xlabel="Nucleotide Position",
                   ylabel="Profile"):
        ax.set_title(axis_title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=legend_title, labelcolor="linecolor", frameon=False,
                  handlelength=0)

    def plot_profile(self, profile, label, columns, errorbars):
        if isinstance(columns, list):
            for column in columns:
                self.ax.plot(profile.data["Nucleotide"], profile.data[column],
                             label=f"{label}: {column.replace('_', ' ')}",
                             drawstyle="steps-mid")
        elif isinstance(columns, str):
            x = profile.data["Nucleotide"]
            profile_values = profile.data[columns]
            self.ax.plot(x, profile_values, label=label, drawstyle="steps-mid")
            if errorbars is not None:
                stderr = profile.data[errorbars]
                self.ax.fill_between(x, profile_values-stderr,
                                     profile_values+stderr, step='mid',
                                     color='C'+str(self.i), alpha=0.25, lw=0)

    def plot_annotation(self, ax, annotation, mode="vbar"):
        color = annotation.color
        a_type = annotation.annotation_type
        if a_type in ["primers", "spans"]:
            for start, end in annotation[:]:
                if start > end:
                    start, end = end, start
                ax.axvspan(start-0.5, end+0.5,
                           fc=color, ec="none", alpha=0.1)
        if a_type == "sites":
            for site in annotation[:]:
                ax.axvline(site, color=color, ls=":")
