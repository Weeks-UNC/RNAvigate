from rnavigate import plots, data
import numpy as np


class Circle(plots.Plot):
    def __init__(self, num_samples, **kwargs):
        try:
            kwargs["subplot_kw"].update({"projection": "polar"})
        except KeyError:
            kwargs["subplot_kw"] = {"projection": "polar"}
        super().__init__(num_samples, **kwargs)
        self.pass_through = [
            "colors",
            "apply_color_to",
            "sequence",
            "title",
            "positions",
        ]
        self.zorder = {
            "annotations": 0,
            "data": 5,
            "nucleotide": 10,
            "sequence": 15,
            "position": 20,
        }

    def set_figure_size(
        self,
        fig=None,
        ax=None,
        rows=None,
        cols=None,
        height_ax_rel=0.035,
        width_ax_rel=0.035,
        width_ax_in=None,
        height_ax_in=None,
        height_gap_in=1,
        width_gap_in=1,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        super().set_figure_size(
            fig=fig,
            ax=ax,
            rows=rows,
            cols=cols,
            height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel,
            width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in,
            width_gap_in=width_gap_in,
            top_in=top_in,
            bottom_in=bottom_in,
            left_in=left_in,
            right_in=right_in,
        )

    def plot_data(
        self,
        sequence,
        structure=None,
        structure2=None,
        interactions=None,
        interactions2=None,
        profile=None,
        annotations=None,
        label=None,
        colors=None,
        gap=30,
        nt_ticks=(20, 5),
    ):
        if annotations is None:
            annotations = []
        seq_circle = data.SequenceCircle(sequence, gap=gap)
        ax = self.get_ax()
        if colors is None:
            colors = {}
        colors = {
            "sequence": None,
            "nucleotides": "sequence",
        } | colors
        for key in ["nucleotides", "sequence"]:
            if isinstance(colors[key], str) and colors[key] == "contrast":
                continue
            elif colors[key] is None:
                continue
            colors[key], colormap = seq_circle.get_colors(
                colors[key],
                profile=profile,
                structure=structure,
                annotations=annotations,
            )
            self.add_colorbar_args(colormap)
        if isinstance(colors["sequence"], np.ndarray):
            colors["nucleotides"], _ = sequence.get_colors("white")
        elif colors["sequence"] == "contrast":
            colors["sequence"] = plots.get_contrasting_colors(colors["nucleotides"])
        if structure is not None:
            structure = structure.as_interactions(structure2)
        if colors["sequence"] is not None:
            plots.plot_sequence_ss(ax, seq_circle, colors["sequence"])
        if colors["nucleotides"] is not None:
            plots.plot_nucleotides_ss(ax, seq_circle, colors["nucleotides"])
        if structure is not None:
            plots.plot_interactions_circle(ax, seq_circle, structure)
            self.add_colorbar_args(structure.cmap)
        if interactions is not None:
            plots.plot_interactions_circle(ax, seq_circle, interactions)
            self.add_colorbar_args(interactions.cmap)
        if interactions2 is not None:
            plots.plot_interactions_circle(ax, seq_circle, interactions2)
            self.add_colorbar_args(interactions2.cmap)
        for annotation in annotations:
            plots.plot_annotation_circle(ax, seq_circle, annotation)
        self.set_axis(
            ax=ax, label=label, seq_circle=seq_circle, gap=gap, nt_ticks=nt_ticks
        )

    def set_axis(self, ax, label, seq_circle, gap, nt_ticks):
        ax.set_title(label)
        ticks, labels = plots.get_nt_ticks(
            sequence=seq_circle.sequence, region=(1, seq_circle.length), gap=nt_ticks[0]
        )
        tick_idx = seq_circle.data["Nucleotide"].isin(ticks)
        ticks = seq_circle.data.loc[tick_idx, "Theta"]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        minor_ticks, _ = plots.get_nt_ticks(
            sequence=seq_circle.sequence,
            region=(1, seq_circle.length),
            gap=nt_ticks[1],
        )
        tick_idx = seq_circle.data["Nucleotide"].isin(minor_ticks)
        minor_ticks = seq_circle.data.loc[tick_idx, "Theta"]
        ax.set_xticks(minor=True, ticks=minor_ticks)
        ax.set(yticks=[], theta_zero_location="N", theta_direction=-1)
        ax.spines["polar"].set_bounds(np.pi * gap / 360, np.pi * (2 - gap / 360))
        # ax.spines["polar"].set_visible(False)
        ax.grid(False)
        self.i += 1
