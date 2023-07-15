from scipy import stats
import seaborn as sns
import numpy as np
from rnavigate import plots


class LinReg(plots.Plot):
    def __init__(self, num_samples):
        super().__init__(num_samples)
        self.paired = None
        self.profiles = []
        self.sequences = []
        self.colors = []
        self.labels = []
        self.pass_through = ["colorby"]

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=None,
                        width_ax_in=7, height_ax_in=7,
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
        return (self.length, self.length)

    def get_figsize(self):
        return (7*self.columns, 7*self.rows)

    def plot_data(self, ct, profile, label, colorby="sequence"):
        self.labels.append(label)
        seq = profile.data["Sequence"].values
        seq_mask = [i for i, x in enumerate(seq) if x.isupper()]
        self.sequences.append(seq[seq_mask])
        colors = profile.get_colors(colorby, profile=profile, ct=ct)
        self.colors.append(colors)
        if ct is not None:
            self.paired = ct.ct != 0
        self.profiles.append(profile.data["Reactivity_profile"].values)
        if len(self.profiles) == self.length:
            for i in range(self.length):
                self.plot_kde(i)
                for j in range(self.length):
                    if i < j:
                        self.plot_regression(i, j)
                    elif i > j:
                        self.axes[i, j].set_axis_off()
            # handles, labels = self.axes[0, 1].get_legend_handles_labels()
            # self.axes[-1, 0].legend(handles, labels, loc=8,
            #                        title = colorby.capitalize())
            handles, labels = self.axes[0, 0].get_legend_handles_labels()
            self.axes[1, 0].legend(handles, labels, loc=9,
                                   title="Structure")

    def plot_regression(self, i, j):
        ax = self.axes[i, j]
        p1 = self.profiles[i]
        p2 = self.profiles[j]
        s1 = self.labels[i]
        s2 = self.labels[j]
        column = "Reactivity"
        colors = self.colors[i]

        ax.plot([0, 1], [0, 1], color='black')
        notNans = ~np.isnan(p1) & ~np.isnan(p2)
        p1 = p1[notNans]
        p2 = p2[notNans]
        colors = colors[notNans]
        gradient, _, r_value, _, _ = stats.linregress(p1, p2)
        ax.text(0.1, 0.8, f'R^2 = {r_value**2:.2f}\nslope = {gradient:.2f}',
                transform=ax.transAxes)
        ax.scatter(p1, p2, c=colors)
        ax.set(xscale='log',
               xlim=[0.00001, 0.3],
               yscale='log',
               ylim=[0.00001, 0.3],
               title=f'{s1} vs. {s2}: {column}')

    def plot_kde(self, i):
        if self.paired is None:
            return
        ax = self.axes[i, i]
        paired = self.paired
        profile = self.profiles[i]
        label = self.labels[i]
        valid = profile > 0
        sns.kdeplot(profile[paired & valid], bw_adjust=0.6, shade=True,
                    label='SS', ax=ax, log_scale=True)
        sns.kdeplot(profile[~paired & valid], bw_adjust=0.6, shade=True,
                    label='DS', ax=ax, log_scale=True)
        ax.annotate(label, (0.1, 0.9), xycoords="axes fraction")
        ax.set(xlim=(0.00001, 0.5))
