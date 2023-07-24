from rnavigate import plots
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection


class AP(plots.Plot):
    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1]-region[0]+1
            self.region = region
        super().__init__(num_samples, **kwargs)
        self.pass_through = ["ax", "seqbar", "title", "profile_panel",
                             "interactions_panel", "interactions2_panel",
                             "ct_panel", "annotation_mode", "plot_error",
                             "profile_scale_factor", "annotation_gap"]

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=0.03, width_ax_rel=0.03,
                        width_ax_in=None, height_ax_in=None,
                        height_gap_in=0.1, width_gap_in=0.1,
                        top_in=1, bottom_in=1,
                        left_in=1, right_in=1):
        super().set_figure_size(fig=fig, ax=ax, rows=rows, cols=cols,
                                height_ax_rel=height_ax_rel,
                                width_ax_rel=width_ax_rel,
                                width_ax_in=width_ax_in,
                                height_ax_in=height_ax_in,
                                height_gap_in=height_gap_in,
                                width_gap_in=width_gap_in, top_in=top_in,
                                bottom_in=bottom_in, left_in=left_in,
                                right_in=right_in)

    def plot_data(self, seq, structure, structure2, interactions,
                  interactions2, profile, label, ax=None, seqbar=True,
                  title=True,
                  interactions_panel="bottom", interactions2_panel="bottom",
                  ct_panel="top", annotations=[], annotation_mode="track",
                  profile_panel="top", annotation_gap=None,
                  profile_scale_factor=1, plot_error=True):
        ax = self.get_ax(ax)
        if annotation_mode == "track" and annotation_gap is None:
            annotation_gap = 2*(2*len(annotations) + seqbar)
        elif annotation_gap is None:
            annotation_gap = 2*seqbar
        if structure is not None:
            structure = structure.as_interactions(structure2)
        self.set_axis(ax=ax, annotation_gap=annotation_gap)
        self.plot_arcs(
            ax=ax, data=structure, panel=ct_panel,
            annotation_gap=annotation_gap)
        self.plot_arcs(
            ax=ax, data=interactions, panel=interactions_panel,
            annotation_gap=annotation_gap)
        self.plot_arcs(
            ax=ax, data=interactions2, panel=interactions2_panel,
            annotation_gap=annotation_gap)
        self.plot_profile(
            ax=ax, profile=profile, annotation_gap=annotation_gap,
            panel=profile_panel, plot_error=plot_error,
            profile_scale_factor=profile_scale_factor)
        for i, annotation in enumerate(annotations):
            self.plot_annotation(ax, annotation=annotation, yvalue=-2-(4*i),
                                 mode=annotation_mode)
        if seqbar:
            self.add_sequence(ax, seq.sequence, yvalue=-annotation_gap,
                              ytrans="data")
        if title:
            self.add_title(ax, label)
        self.i += 1

    def set_axis(self, ax, annotation_gap=0, xticks=20, xticks_minor=10):
        def get_ticks(x, mn, mx):
            return [tick for tick in range(x, mx+1, x) if mn <= tick <= mx]

        ax.set_aspect('equal')
        ax.yaxis.set_visible(False)
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set(position=('data', -annotation_gap),
                                visible=False)
        ax.spines['top'].set_color('none')
        height = min(300, self.nt_length/2)
        mn, mx = self.region
        ax.set(xlim=(mn - 0.5, mx + 0.5),
               ylim=(-height-1-annotation_gap, height+1),
               xticks=get_ticks(xticks, mn, mx),
               axisbelow=False)
        ax.set_xticks(get_ticks(xticks_minor, mn, mx), minor=True)
        for label in ax.get_xticklabels():
            label.set_bbox({"facecolor": "white",
                            "edgecolor": "None",
                            "alpha": 0.5,
                            "boxstyle": "round,pad=0.1,rounding_size=0.2"})

    def plot_arcs(self, ax, data, panel, annotation_gap):
        if data is None:
            return
        ij_colors = data.get_ij_colors()
        self.add_colorbar_args(interactions=data)
        patches = []
        mn, mx = self.region
        for i, j, color in zip(*ij_colors):
            if j < i:  # flip the order
                i, j = j, i
            if ((i < mn) and (j < mn)) or ((i > mx) and (j > mx)):
                continue
            if panel == "top":
                center = ((i+j)/2., 0)
                theta1 = 0
                theta2 = 180
            elif panel == "bottom":
                center = ((i+j)/2., -annotation_gap)
                theta1 = 180
                theta2 = 360
            radius = 0.5+(j-i)/2.
            patches.append(Wedge(center, radius, theta1, theta2, color=color,
                                 width=1, ec='none'))
        ax.add_collection(PatchCollection(patches, match_original=True))

    def add_title(self, ax, label):
        ax.annotate(label, xy=(0.1, 0.9), xycoords="axes fraction",
                    fontsize=12, ha='left')

    def get_figsize(self):
        width = self.nt_length * 0.1 + 1
        height = min(self.nt_length, 602) * 0.1 + 1
        return (width*self.columns, height*self.rows)

    def plot_profile(self, ax, profile, annotation_gap, profile_scale_factor=1,
                     plot_error=True, panel="top"):
        if profile is None:
            return
        data = profile.get_plotting_dataframe()
        factor = profile_scale_factor
        values = data["Values"]
        colormap = data["Colors"]
        nts = data["Nucleotide"]
        mn, mx = self.region
        if panel == "top":
            bottom = 0
        elif panel == "bottom":
            values *= -1
            bottom = -annotation_gap
        if plot_error and ("Errors" in data.columns):
            yerr = data["Errors"]
            ax.bar(nts[mn-1:mx], values[mn-1:mx]*factor, align="center",
                   bottom=bottom, width=1, color=colormap[mn-1:mx],
                   edgecolor=colormap[mn-1:mx], linewidth=0.0,
                   yerr=yerr[mn-1:mx], ecolor=(0, 0, 1 / 255.0), capsize=0)
        else:
            ax.bar(nts[mn-1:mx], values[mn-1:mx]*factor, align="center",
                   width=1, color=colormap[mn-1:mx], bottom=bottom,
                   edgecolor=colormap[mn-1:mx], linewidth=0.0)

    def plot_annotation(self, ax, annotation, yvalue, mode):
        color = annotation.color
        modes = ["track", "vbar"]
        assert mode in modes, f"annotation mode must be one of: {modes}"
        a_type = annotation.annotation_type
        if a_type == "spans" or (a_type == "primers" and mode == "vbar"):
            for start, end in annotation:
                if start > end:
                    start, end = end, start
                if (self.region[0] > end) or (self.region[1] < start):
                    continue
                if (self.region[0] > start):
                    start = self.region[0]
                if (self.region[1] < end):
                    end = self.region[1]
                if mode == "track":
                    ax.plot([start, end], [yvalue]*2,
                            color=color, alpha=0.7, lw=11)
                elif mode == "vbar":
                    ax.axvspan(start-0.5, end+0.5,
                               fc=color, ec="none", alpha=0.1)
        elif a_type == "sites":
            sites = annotation[:]
            if mode == "track":
                ax.scatter(sites, [yvalue]*len(sites),
                           color=color, marker='*',
                           ec="none", alpha=0.7, s=20**2)
            elif mode == "vbar":
                for site in sites:
                    ax.axvline(site, color=color, ls=":")
        elif a_type == "groups":
            sites = annotation[:]
            span = [min(sites), max(sites)]
            if mode == "track":
                ax.plot(span, [yvalue, yvalue],
                        color=color, alpha=0.4, lw=11)
                ax.scatter(sites, [yvalue]*len(sites),
                            color=color, marker='o', ec="none", alpha=0.7,
                            s=11**2)
            elif mode == "vbar":
                ax.axvspan(span[0]-0.5, span[1]+0.5,
                            fc=color, ec="none", alpha=0.1)
                for site in sites:
                    ax.axvspan(site-0.5, site+0.5, fc=color, ec="none",
                               alpha=0.1)
        elif a_type == "primers":
            for start, end in annotation:
                if mode == "track":
                    if start < end:
                        xs = [start, end, end-1]
                        ys = [yvalue, yvalue, yvalue+1]
                    elif start > end:
                        xs = [start, end, end+1]
                        ys = [yvalue, yvalue, yvalue-1]
                    ax.plot(xs, ys, color=color, alpha=0.7)
