from rnavigate import plots, styles


class AP(plots.Plot):
    def __init__(
            self, num_samples, nt_length, region="all", track_labels=True,
            **kwargs
            ):
        self.track_labels = track_labels
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1]-region[0]+1
            self.region = region
        super().__init__(num_samples, **kwargs)


    def plot_data(
            self, sequence, structure=None, structure2=None, interactions=None,
            interactions2=None, profile=None, annotations=None, domains=None,
            label='', ax=None, seqbar=True, title=True, panels=None,
            annotation_mode="track", track_height=None, profile_scale_factor=1,
            plot_error=False,
            ):
        ax = self.get_ax(ax)
        if panels is None:
            panels = {}
        panels = {
            'interactions': 'bottom',
            'interactions2': 'bottom',
            'structure': 'top',
            'profile': 'top'
        } | panels
        if annotations is None:
            annotations = []

        annotation_gap = 4 * len(annotations) * (annotation_mode == 'track')
        seqbar_height = 2 * seqbar
        domains_height = 6 * (domains is not None)
        if track_height is None:
            track_height = annotation_gap + seqbar_height + domains_height

        yvalues = {'bottom': 0, 'top': track_height}
        if structure is not None:
            structure = structure.as_interactions(structure2)
            plots.plot_interactions_arcs(
                ax=ax, interactions=structure, panel=panels['structure'],
                yvalue=yvalues[panels['structure']], region=self.region,
                )
            self.add_colorbar_args(structure.cmap)
        if interactions is not None:
            plots.plot_interactions_arcs(
                ax=ax, interactions=interactions,
                panel=panels['interactions'],
                yvalue=yvalues[panels['interactions']], region=self.region,
                )
            self.add_colorbar_args(interactions.cmap)
        if interactions2 is not None:
            plots.plot_interactions_arcs(
                ax=ax, interactions=interactions2,
                panel=panels['interactions2'],
                yvalue=yvalues[panels['interactions2']], region=self.region,
                )
            self.add_colorbar_args(interactions2.cmap)
        if profile is not None:
            if panels['profile'] == 'bottom':
                scale_factor = profile_scale_factor * -1
                bottom = 0
            else:
                scale_factor = profile_scale_factor
                bottom = track_height
            plots.plot_profile_bars(
                ax=ax, profile=profile, bottom=bottom,
                scale_factor=scale_factor, plot_error=plot_error,
                region=self.region,
                )
            self.add_colorbar_args(profile.cmap)

        yticks, ylabels = [], []
        if seqbar:
            plots.plot_sequence_track(
                ax, sequence.sequence, yvalue=0, height=2, ytrans="data"
                )
            yticks.append(1)
            ylabels.append('sequence')
            self.add_colorbar_args(styles.get_nt_cmap())
        if domains is not None:
            for domain in domains:
                plots.plot_domain_track(
                    ax=ax, spans=domain, yvalue=seqbar_height, height=6,
                    region=self.region,
                    )
            yticks.append(seqbar_height + 3)
            ylabels.append('domains')
        for i, annotation in enumerate(annotations):
            yvalue = seqbar_height + domains_height + 4 * (i+1) - 2
            plots.plot_annotation_track(
                ax, annotation=annotation, yvalue=yvalue, height=4,
                mode=annotation_mode, region=self.region,
                )
            yticks.append(yvalue)
            ylabels.append(annotation.name)
        if not self.track_labels:
            yticks, ylabels = [], []
        self.set_axis(
            ax=ax, track_height=track_height, yticks=yticks, ylabels=ylabels,
            )
        if title:
            ax.set_title(label)
        self.i += 1

    def set_axis(self, ax, track_height=0, xticks=20, xticks_minor=10,
                 max_height=300, yticks=None, ylabels=None):
        def get_ticks(x, mn, mx):
            return [tick for tick in range(x, mx+1, x) if mn <= tick <= mx]

        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set(position=('data', 0), visible=False)
        ax.spines['top'].set_color('none')
        height = min(max_height, self.nt_length/2)
        mn, mx = self.region
        ax.set(xlim=(mn - 0.5, mx + 0.5),
               ylim=(-height-1, height+1+track_height),
               xticks=get_ticks(xticks, mn, mx),
               yticks=yticks,
               yticklabels=ylabels,
               axisbelow=False)
        ax.set_xticks(get_ticks(xticks_minor, mn, mx), minor=True)
        for label in ax.get_xticklabels():
            label.set_bbox({
                "facecolor": "white",
                "edgecolor": "None",
                "alpha": 0.5,
                "boxstyle": "round,pad=0.1,rounding_size=0.2"
                })

    def set_figure_size(
            self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=0.03,
            width_ax_rel=0.03, width_ax_in=None, height_ax_in=None,
            height_gap_in=0.5, width_gap_in=0.5, top_in=1, bottom_in=1,
            left_in=1, right_in=1
            ):
        super().set_figure_size(
            fig=fig, ax=ax, rows=rows, cols=cols, height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel, width_ax_in=width_ax_in,
            height_ax_in=height_ax_in, height_gap_in=height_gap_in,
            width_gap_in=width_gap_in+self.track_labels*0.5, top_in=top_in,
            bottom_in=bottom_in, left_in=left_in, right_in=right_in
            )
