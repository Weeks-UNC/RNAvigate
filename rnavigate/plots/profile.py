import seaborn as sns
from rnavigate import plots, styles

class Profile(plots.Plot):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples, sharey=True, **kwargs)
        self.track_height = 0

    def get_rows_columns(self, rows=None, cols=None):
        return super().get_rows_columns(rows=rows, cols=1)

    def set_labels(self, ax, axis_title="Reactivity Profile",
                   xlabel="Nucleotide Position", ylabel="Reactivity"):
        ax.set_title(axis_title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    def plot_data(self, profile, annotations, domains, label, plot_error=True,
                  column=None, seqbar=True, annotations_mode='track'):
        ax = self.get_ax()
        if column is not None:
            profile.metric = column
        column = profile.metric
        plots.plot_profile_bars(
            ax=ax, profile=profile, plot_error=plot_error, region=self.region
            )
        self.add_colorbar_args(profile.cmap)

        track_unit = 0.015
        annotations_track = 4 * track_unit * len(annotations)
        annotations_track *= annotations_mode == 'track'
        domains_track = 6 * track_unit * (domains is not None)
        sequence_track = 2 * track_unit * seqbar
        self.track_height = sequence_track + annotations_track + domains_track
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
                        - self.track_height + (i+0.5) * (4 * track_unit)),
                ytrans='axes', region=self.region,
                )

        self.i += 1
        ylabel = column.replace("_", " ")
        self.set_labels(ax=ax, ylabel=ylabel, axis_title=label)
        self.set_axis(ax)

    def set_figure_size(
            self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None,
            width_ax_rel=0.03, width_ax_in=None, height_ax_in=2,
            height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=1,
            left_in=1, right_in=1
            ):
        super().set_figure_size(
            fig=fig, ax=ax, rows=rows, cols=cols, height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel, width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in+2*self.track_height,
            width_gap_in=width_gap_in, top_in=top_in,
            bottom_in=bottom_in+2*self.track_height,
            left_in=left_in, right_in=right_in
            )

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
