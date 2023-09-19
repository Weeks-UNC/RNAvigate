from scipy import stats
import seaborn as sns
import matplotlib as mpl
import numpy as np
from .plots import Plot
from  rnavigate import styles


class LinReg(Plot):
    def __init__(self, num_samples, scale='log', regression='pearson'):
        super().__init__(num_samples)
        linreg_axes = []
        kde_axes = []
        for i in range(num_samples):
            kde_axes.append(self.axes[i,i])
            self.axes[i, i].set_yticks([])
            for j in range(i+1, num_samples):
                linreg_axes.append(self.axes[i,j])
        self.axes[0, 1].get_shared_x_axes().join(*linreg_axes+kde_axes)
        self.axes[0, 1].get_shared_y_axes().join(*linreg_axes)
        self.lims = [0, 0]
        self.scale = scale
        self.regression = {
            'pearson': stats.pearsonr,
            'spearman': stats.spearmanr}[regression]
        self.profiles = []
        self.colors = []
        self.labels = []
        self.pass_through = ["colorby", "column"]

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=None,
                        width_ax_in=7, height_ax_in=7,
                        height_gap_in=1, width_gap_in=1,
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

    def get_rows_columns(self, rows=None, cols=None):
        return (self.length, self.length)

    def get_figsize(self):
        return (7*self.columns, 7*self.rows)

    def plot_data(self, sequence, ct, profile, annotations, label, column,
                  colorby="structure"):
        profile = profile.get_plotting_dataframe(column=column)
        self.labels.append(label)
        self.profiles.append(profile['Values'].to_numpy(copy=True))
        self.colors.append(
            sequence.get_colors(colorby, profile=profile, ct=ct,
                                annotations=annotations))
        self.add_legend(colorby, annotations)
        if len(self.profiles) == self.length:
            for i in range(self.length):
                self.plot_kde(i)
                for j in range(self.length):
                    if i < j:
                        self.plot_regression(i, j)
                    elif i > j:
                        self.axes[i, j].set_axis_off()
            if self.scale == 'linear':
                buffer = 0.05 * (self.lims[1] - self.lims[0])
                self.lims[0] -= buffer
                self.lims[1] += buffer
                self.axes[0, 1].set(
                    ylim=self.lims,
                    xlim=self.lims)

    def add_legend(self, colorby, annotations):
        handles = []
        if colorby == 'sequence':
            for label in ['A', 'U', 'C', 'G']:
                color = styles.get_nt_color(label)
                handles.append(mpl.lines.Line2D(
                    [], [], color=color, marker='.', linestyle='None',
                    markersize=10, label=label))
        elif colorby == 'structure':
            for label, color in zip(['paired', 'unpaired'], ['C0', 'C1']):
                handles.append(mpl.lines.Line2D(
                    [], [], color=color, marker='.', linestyle='None',
                    markersize=10, label=label))
        elif colorby == 'annotations':
            for annotation in annotations:
                label = annotation.name
                color = annotation.color
                handles.append(mpl.lines.Line2D(
                    [], [], color=color, marker='.', linestyle='None',
                    markersize=10, label=label))
        self.axes[1, 0].legend(handles=handles, loc=9, title=colorby)

    def plot_regression(self, i, j):
        ax = self.axes[i, j]
        p1 = self.profiles[i]
        p2 = self.profiles[j]
        s1 = self.labels[i]
        s2 = self.labels[j]
        colors = self.colors[i]

        # ax.plot([0, 1], [0, 1], color='black')
        if self.scale == 'log':
            notNans = (p1 > 0) & (p2 > 0)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            gradient, _, _, _, _ = stats.linregress(np.log10(p1), np.log10(p2))
            r_value, p_value = self.regression(np.log10(p1), np.log10(p2))
        if self.scale == 'linear':
            notNans = ~np.isnan(p1) & ~np.isnan(p2)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            gradient, _, _, _, _ = stats.linregress(p1, p2)
            r_value, p_value = self.regression(p1, p2)
        minimum, maximum = self.lims
        self.lims[0] = min([minimum, min(p1), min(p2)])
        self.lims[1] = max([maximum, max(p1), max(p2)])
        ax.text(
            0.1, 0.8,
            f'R^2: {r_value**2:.2f}\nslope: {gradient:.2f}\np: {p_value:.4f}',
            transform=ax.transAxes)
        ax.scatter(p1, p2, c=colors)
        ax.set(xscale=self.scale,
               yscale=self.scale,
               title=f'{s1} vs. {s2}')

    def plot_kde(self, i):
        if len(np.unique(self.colors[i])) > 5:
            return
        ax = self.axes[i, i]
        profile = self.profiles[i]
        label = self.labels[i]
        valid = profile > 0
        for color in np.unique(self.colors[i]):
            group = self.colors[i] == color
            sns.kdeplot(profile[group & valid], bw_adjust=0.6, shade=True,
                        ax=ax, log_scale=(self.scale=='log'), color=color)
        ax.annotate(label, (0.1, 0.9), xycoords="axes fraction")
