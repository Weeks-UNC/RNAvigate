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

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=0.1,
                        width_ax_in=None, height_ax_in=6,
                        height_gap_in=1, width_gap_in=0.5,
                        top_in=1, bottom_in=0.5,
                        left_in=0.5, right_in=0.5):
        super().set_figure_size(fig=fig, ax=ax, rows=rows, cols=cols,
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
        return self.ax

    def plot_data(self, profile, annotations, label,
                  columns="Reactivity_profile", seqbar=True, errorbars=None):
        ax = self.get_ax()
        annotations = [annotation.fitted for annotation in annotations]
        self.plot_profile(ax, profile, label, columns, errorbars)
        if seqbar and (self.i == 0):
            self.add_sequence(ax, profile.sequence)
        self.i += 1
        if self.i == self.length:
            if isinstance(columns, list):
                ylabel = [column.replace("_", " ") for column in columns]
                ylabel = ', '.join(ylabel)
            else:
                ylabel = columns.replace("_", " ")
            self.set_labels(ax=ax, ylabel=ylabel, axis_title=None,
                            legend_title=None)
            for annotation in annotations:
                self.plot_annotation(ax=ax, annotation=annotation)
            self.set_axis(ax)

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
        ax.set_title(axis_title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=legend_title, labelcolor="linecolor", frameon=False,
                  handlelength=0)

    def plot_profile(self, ax, profile, label, columns, errorbars):
        data = profile.get_plotting_dataframe(all_columns=True)
        if isinstance(columns, list):
            for column in columns:
                ax.plot(data["Nucleotide"], data[column],
                        label=f"{label}: {column.replace('_', ' ')}",
                        drawstyle="steps-mid")
        elif isinstance(columns, str):
            x = data["Nucleotide"]
            profile_values = data[columns]
            ax.plot(x, profile_values, label=label, drawstyle="steps-mid")
            if errorbars is not None:
                stderr = data[errorbars]
                ax.fill_between(x, profile_values-stderr,
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


class Profile(Skyline):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        super().__init__(num_samples, nt_length, region, sharey=True, **kwargs)
        self.pass_through = ["plot_errors", "column", "seqbar", "error_column"]
        self.axes[0, 0].set_ylim(-0.5, 4.5)

    def get_ax(self, i=None):
        return super(Skyline, self).get_ax(i)

    def get_rows_columns(self, rows=None, cols=None):
        return super(Skyline, self).get_rows_columns(rows, cols=1)

    def set_labels(self, ax, axis_title="Reactivity Profile",
                   xlabel="Nucleotide Position", ylabel="Reactivity"):
        ax.set_title(axis_title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    def plot_data(self, seq, profile, annotations, label, plot_errors=True,
                  column=None, seqbar=True, error_column=None):
        ax = self.get_ax()
        if column is None:
            column = profile.default_column
        annotations = [annotation.fitted for annotation in annotations]
        self.plot_profile(ax=ax, profile=profile, column=column,
                          error_column=error_column, plot_errors=plot_errors)
        if seqbar:
            self.add_sequence(ax, seq.sequence)
        self.i += 1
        ylabel = column.replace("_", " ")
        self.set_labels(ax=ax, ylabel=ylabel, axis_title=label)
        for annotation in annotations:
            self.plot_annotation(ax=ax, annotation=annotation)
        self.set_axis(ax)

    def plot_profile(self, ax, profile, column, error_column, plot_errors):
        data = profile.get_plotting_dataframe(column=column,
                                              err_column=error_column)
        values = data["Values"]
        colormap = data["Colors"]
        nts = data["Nucleotide"]
        mn, mx = self.region
        if plot_errors and ("Errors" in data.columns):
            yerr = data["Errors"]
            ax.bar(nts[mn-1:mx], values[mn-1:mx], align="center",
                   width=1, color=colormap[mn-1:mx],
                   edgecolor=colormap[mn-1:mx], linewidth=0.0,
                   yerr=yerr[mn-1:mx], ecolor=(0, 0, 1 / 255.0), capsize=1)
        else:
            ax.bar(nts[mn-1:mx], values[mn-1:mx], align="center",
                   width=1, color=colormap[mn-1:mx],
                   edgecolor=colormap[mn-1:mx], linewidth=0.0)

    def set_figure_size(self, fig=None, ax=None, rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=0.1, width_ax_in=None,
                        height_ax_in=3, height_gap_in=1.5, width_gap_in=0.5,
                        top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5):
        return super().set_figure_size(fig, ax, rows, cols, height_ax_rel,
                                       width_ax_rel, width_ax_in, height_ax_in,
                                       height_gap_in, width_gap_in, top_in,
                                       bottom_in, left_in, right_in)
