from rnavigate import plots
import seaborn as sns


class Skyline(plots.Plot):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples, **kwargs)
        self.axis = self.axes[0, 0]
        self.pass_through = ["columns", "seqbar", "errorbars"]

    def set_figure_size(self, fig=None, axis=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=0.1,
                        width_ax_in=None, height_ax_in=6,
                        height_gap_in=1, width_gap_in=0.5,
                        top_in=1, bottom_in=0.5,
                        left_in=0.5, right_in=0.5):
        super().set_figure_size(fig=fig, axis=axis, rows=rows, cols=cols,
                                height_ax_rel=height_ax_rel,
                                width_ax_rel=width_ax_rel,
                                width_ax_in=width_ax_in,
                                height_ax_in=height_ax_in,
                                height_gap_in=height_gap_in,
                                width_gap_in=width_gap_in, top_in=top_in,
                                bottom_in=bottom_in, left_in=left_in,
                                right_in=right_in)

    def get_rows_columns(self, rows=None, cols=None):
        return (1, 1)

    def get_ax(self, i=None):
        return self.axis

    def plot_data(self, profile, annotations, label,
                  columns="Reactivity_profile", seqbar=True, plot_error=False):
        axis = self.get_ax()
        self.plot_profile(axis, profile, label, columns, plot_error)
        if seqbar and (self.i == 0):
            self.add_sequence(axis, profile.sequence)
        self.i += 1
        if self.i == self.length:
            if isinstance(columns, list):
                ylabel = [column.replace("_", " ") for column in columns]
                ylabel = ', '.join(ylabel)
            else:
                ylabel = columns.replace("_", " ")
            self.set_labels(axis=axis, ylabel=ylabel, axis_title=None,
                            legend_title=None)
            for annotation in annotations:
                self.plot_annotation(axis=axis, annotation=annotation)
            self.set_axis(axis)

    def get_figsize(self):
        left_inches = 0.9
        right_inches = 0.4
        ax_width = self.nt_length * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        return (fig_width, fig_height)

    def set_axis(self, axis, xticks=20, xticks_minor=5):
        xlim = self.region
        axis.set_xlim([xlim[0] - 0.5, xlim[1] + 0.5])
        xrange = range(xlim[0], xlim[1]+1)
        axis.set_xticks([x for x in xrange if (x % xticks) == 0])
        axis.set_xticks([x for x in xrange if (x % xticks_minor) == 0],
                      minor=True)
        axis.spines['bottom'].set_visible(False)
        axis.spines['left'].set_visible(False)
        axis.grid(axis='y', visible=True)
        axis.tick_params(axis='y', length=0, grid_alpha=0.4)
        sns.despine(ax=axis, bottom=True, left=True)

    def set_labels(self, axis, axis_title="Raw Reactivity Profile",
                   legend_title="Samples", xlabel="Nucleotide Position",
                   ylabel="Profile"):
        axis.set_title(axis_title, loc="left")
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.legend(title=legend_title, labelcolor="linecolor", frameon=False,
                  handlelength=0)

    def plot_profile(self, axis, profile, label, columns, plot_error):
        data = profile.data
        if isinstance(columns, list):
            for column in columns:
                axis.plot(data["Nucleotide"], data[column],
                        label=f"{label}: {column.replace('_', ' ')}",
                        drawstyle="steps-mid")
        elif isinstance(columns, str):
            x = data["Nucleotide"]
            profile.metric = columns
            profile_values = data[columns]
            axis.plot(x, profile_values, label=label, drawstyle="steps-mid")
            if plot_error and profile.error_column is not None:
                stderr = data[profile.error_column]
                axis.fill_between(x, profile_values-stderr,
                                profile_values+stderr, step='mid',
                                color='C'+str(self.i), alpha=0.25, lw=0)

    def plot_annotation(self, axis, annotation, mode="vbar"):
        color = annotation.color
        a_type = annotation.annotation_type
        if a_type in ["primers", "spans"]:
            for start, end in annotation:
                if start > end:
                    start, end = end, start
                axis.axvspan(start-0.5, end+0.5,
                           fc=color, ec="none", alpha=0.1)
        if a_type == "sites":
            for site in annotation:
                axis.axvline(site, color=color, ls=":")


class Profile(Skyline):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        super().__init__(num_samples, nt_length, region, sharey=True, **kwargs)
        self.pass_through = ["plot_errors", "column", "seqbar", "error_column"]
        self.axes[0, 0].set_ylim(-0.5, 4.5)

    def get_ax(self, i=None):
        return super(Skyline, self).get_ax(i)

    def get_rows_columns(self, rows=None, cols=None):
        return super(Skyline, self).get_rows_columns(rows, cols=1)

    def set_labels(self, axis, axis_title="Reactivity Profile",
                   xlabel="Nucleotide Position", ylabel="Reactivity"):
        axis.set_title(axis_title, loc="left")
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)

    def plot_data(self, seq, profile, annotations, label, plot_errors=True,
                  column=None, seqbar=True):
        axis = self.get_ax()
        if column is not None:
            profile.metric = column
        column = profile.metric
        self.plot_profile(axis=axis, profile=profile, plot_errors=plot_errors)
        if seqbar:
            self.add_sequence(axis, seq.sequence)
        self.i += 1
        ylabel = column.replace("_", " ")
        self.set_labels(axis=axis, ylabel=ylabel, axis_title=label)
        for annotation in annotations:
            self.plot_annotation(axis=axis, annotation=annotation)
        self.set_axis(axis)

    def plot_profile(self, axis, profile, plot_errors):
        data = profile.get_plotting_dataframe()
        values = data["Values"]
        colormap = data["Colors"]
        nts = data["Nucleotide"]
        mn, mx = self.region
        if plot_errors and ("Errors" in data.columns):
            yerr = data["Errors"]
            axis.bar(nts[mn-1:mx], values[mn-1:mx], align="center",
                   width=1, color=colormap[mn-1:mx],
                   edgecolor=colormap[mn-1:mx], linewidth=0.0,
                   yerr=yerr[mn-1:mx], ecolor=(0, 0, 1 / 255.0), capsize=1)
        else:
            axis.bar(nts[mn-1:mx], values[mn-1:mx], align="center",
                   width=1, color=colormap[mn-1:mx],
                   edgecolor=colormap[mn-1:mx], linewidth=0.0)

    def set_figure_size(self, fig=None, axis=None, rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=0.1, width_ax_in=None,
                        height_ax_in=3, height_gap_in=1.5, width_gap_in=0.5,
                        top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5):
        return super().set_figure_size(fig, axis, rows, cols, height_ax_rel,
                                       width_ax_rel, width_ax_in, height_ax_in,
                                       height_gap_in, width_gap_in, top_in,
                                       bottom_in, left_in, right_in)
