
from rnavigate import plots, data
import numpy as np


class Circle(plots.Plot):

    def __init__(self, num_samples, **kwargs):
        try:
            kwargs['subplot_kw'].update({'projection': 'polar'})
        except KeyError:
            kwargs['subplot_kw'] = {'projection': 'polar'}
        super().__init__(num_samples, **kwargs)
        self.pass_through = ["colors", "apply_color_to", "sequence",
                             "title", "positions"]
        self.zorder = {"annotations": 0,
                       "data": 5,
                       "nucleotide": 10,
                       "sequence": 15,
                       "position": 20}

    def set_figure_size(
            self, fig=None, ax=None, rows=None, cols=None,
            height_ax_rel=1, width_ax_rel=1, width_ax_in=None,
            height_ax_in=None, height_gap_in=1, width_gap_in=0.5, top_in=1,
            bottom_in=0.5, left_in=0.5, right_in=0.5
            ):
        super().set_figure_size(
            fig=fig, ax=ax, rows=rows, cols=cols, height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel, width_ax_in=width_ax_in,
            height_ax_in=height_ax_in, height_gap_in=height_gap_in,
            width_gap_in=width_gap_in, top_in=top_in, bottom_in=bottom_in,
            left_in=left_in, right_in=right_in
            )

    def plot_data(
            self, sequence, structure=None, structure2=None, interactions=None,
            interactions2=None, profile=None, annotations=None, label=None,
            colors=None, title=True, positions=20, gap=30):
        if annotations is None:
            annotations = []
        seq_circle = data.SequenceCircle(sequence, gap=gap)
        ax = self.get_ax()
        if colors is None:
            colors = {}
        colors = {
            'sequence': None,
            'nucleotides': 'sequence',
        } | colors
        for key in ['nucleotides', 'sequence']:
            if isinstance(colors[key], str) and colors[key] == 'contrast':
                continue
            elif colors[key] is None:
                continue
            colors[key], colormap = seq_circle.get_colors(
                colors[key], profile=profile,structure=structure,
                annotations=annotations
                )
            self.add_colorbar_args(colormap)
        if isinstance(colors['sequence'], np.ndarray):
            colors['nucleotides'] = structure.get_colors('white')
        elif colors['sequence'] == 'contrast':
            colors['sequence'] = plots.get_contrasting_colors(
                colors['nucleotides']
                )
        if structure is not None:
            structure = structure.as_interactions(structure2)
        if colors['sequence'] is not None:
            plots.plot_sequence_ss(ax, seq_circle, colors['sequence'])
        if colors['nucleotides'] is not None:
            plots.plot_nucleotides_ss(ax, seq_circle, colors['nucleotides'])
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
        if title:
            ax.set_title(label)
        tick_set = seq_circle.data.query(
            f'Nucleotide % {positions} == 0 '
            f'| Nucleotide in [1, {seq_circle.length}]')
        ax.set(
            xticks=tick_set['Theta'],
            xticklabels=tick_set['Nucleotide'],
            yticks=[],
            theta_zero_location='N',
            theta_direction=-1)
        ax.spines['polar'].set_bounds(np.pi*gap/360, np.pi*(2-gap/360))
        ax.grid(False)
        self.i += 1
