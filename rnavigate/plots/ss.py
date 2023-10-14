from rnavigate import plots
import numpy as np


class SS(plots.Plot):
    def __init__(self, num_samples, **kwargs):
        super().__init__(num_samples, sharey=True, sharex=True, **kwargs)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.axis("off")
        self.xlims = [-2, 2]
        self.ylims = [-2, 2]

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=0.2, width_ax_rel=0.2,
                        width_ax_in=None, height_ax_in=None,
                        height_gap_in=0.5, width_gap_in=0.2,
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

    def plot_data(self, structure, interactions=None, interactions2=None,
                  profile=None, annotations=None, label='', colors=None,
                  positions=None, bp_style="dotted"):
        if annotations is None:
            annotations = []
        ax = self.get_ax()
        if colors is None:
            colors = {}
        colors = {
            'sequence': None,
            'nucleotides': 'sequence',
            'structure': 'grey',
            'basepairs': 'grey',
        } | colors
        for key in ['nucleotides', 'structure', 'basepairs', 'sequence']:
            if isinstance(colors[key], str) and colors[key] == 'contrast':
                continue
            elif colors[key] is None:
                continue
            colors[key] = structure.get_colors(
                colors[key], profile=profile,structure=structure,
                annotations=annotations
                )
        if isinstance(colors['sequence'], np.ndarray):
            colors['nucleotides'] = structure.get_colors('white')
        elif colors['sequence'] == 'contrast':
            colors['sequence'] = plots.get_contrasting_colors(
                colors['nucleotides']
                )
        if colors['sequence'] is not None:
            plots.plot_sequence_ss(ax, structure, colors['sequence'])
        if colors['nucleotides'] is not None:
            plots.plot_nucleotides_ss(ax, structure, colors['nucleotides'])
        if colors['structure'] is not None:
            plots.plot_structure_ss(ax, structure, colors['structure'])
        if colors['basepairs'] is not None:
            plots.plot_basepairs_ss(ax, structure, bp_style)
        if positions is not None:
            plots.plot_positions_ss(ax, structure, positions)
        if interactions is not None:
            plots.plot_interactions_ss(ax, structure, interactions)
        if interactions2 is not None:
            plots.plot_interactions_ss(ax, structure, interactions2)
        for annotation in annotations:
            plots.plot_annotation_ss(ax, structure, annotation)
        ax.set_title(label)
        self.set_axis(ax, structure)
        self.i += 1

    def set_axis(self, ax, structure):
        x = structure.xcoordinates
        y = structure.ycoordinates
        self.xlims = [
            min([self.xlims[0], min(x)-2]),
            max([self.xlims[1], max(x)+2])]
        self.ylims = [
            min([self.ylims[0], min(y)-2]),
            max([self.ylims[1], max(y)+2])]
        ax.set(ylim=self.ylims, xlim=self.xlims)
