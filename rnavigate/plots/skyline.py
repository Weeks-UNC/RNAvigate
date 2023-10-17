from rnavigate import plots, styles
import seaborn as sns


class Skyline(plots.Plot):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples=num_samples, **kwargs)
        self.ax = self.axes[0, 0]
        self.track_height = 0

    def set_figure_size(
            self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None,
            width_ax_rel=0.03, width_ax_in=None, height_ax_in=2,
            height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=1,
            left_in=1, right_in=1
            ):
        super().set_figure_size(
            fig=fig, ax=ax, rows=rows, cols=cols, height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel, width_ax_in=width_ax_in,
            height_ax_in=height_ax_in, height_gap_in=height_gap_in,
            width_gap_in=width_gap_in, top_in=top_in,
            bottom_in=bottom_in+2*self.track_height,
            left_in=left_in, right_in=right_in
            )

    def get_rows_columns(self, rows=None, cols=None):
        return (1, 1)

    def get_ax(self, i=None):
        return self.ax

    def plot_data(
            self, profile, annotations=None, domains=None, label=None,
            columns=None, seqbar=True, errors=None,
            annotations_mode="track"
            ):
        if columns is None:
            columns = profile.metric
        track_unit = 0.03
        annotations_track = 2 * track_unit * len(annotations)
        annotations_track *= annotations_mode == 'track'
        domains_track = 3 * track_unit * (domains is not None)
        sequence_track = track_unit * seqbar
        self.track_height = sequence_track + annotations_track + domains_track

        ax = self.get_ax()
        plots.plot_profile_skyline(ax, profile, label, columns, errors)
        if self.i == 0:
            if seqbar:
                plots.plot_sequence_track(
                    ax=ax, sequence=profile.sequence,
                    yvalue=self.track_height * -1, height=sequence_track,
                    region=self.region, ytrans='axes')
                self.add_colorbar_args(styles.get_nt_cmap())
            if domains is not None:
                for domain in domains:
                    plots.plot_domain_track(
                        ax=ax, spans=domain,
                        yvalue=sequence_track - self.track_height,
                        height=domains_track, region=self.region,
                        ytrans='axes')
            for i, annotation in enumerate(annotations):
                plots.plot_annotation_track(
                    ax=ax, annotation=annotation, mode=annotations_mode,
                    yvalue=(sequence_track + domains_track
                            - self.track_height + (i+0.5) * (2 * track_unit)),
                    height=2*track_unit, ytrans='axes', region=self.region,
                    )

        self.i += 1
        if self.i == self.length:
            if isinstance(columns, list):
                ylabel = [column.replace("_", " ") for column in columns]
                ylabel = ', '.join(ylabel)
            else:
                ylabel = columns.replace("_", " ")
            self.set_labels(ax=ax, ylabel=ylabel, axis_title=None,
                            legend_title=None)
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
        ax.spines['bottom'].set(position=('axes', self.track_height * -1),
                                visible=False)
        ax.spines['left'].set_visible(False)
        ax.grid(axis='y', visible=True)
        ax.tick_params(axis='y', length=0, grid_alpha=0.4)
        sns.despine(ax=ax, bottom=True, left=True)

    def set_labels(self, ax, axis_title="Raw Reactivity Profile",
                   legend_title="Samples", xlabel="Nucleotide Position",
                   ylabel="Profile"):
        ax.set_title(axis_title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=legend_title, labelcolor="linecolor", frameon=False,
                  handlelength=0)
