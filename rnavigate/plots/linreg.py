from scipy import stats
import seaborn as sns
import matplotlib as mpl
import numpy as np
from rnavigate import plots, styles


class LinReg(plots.Plot):
    def __init__(self, num_samples, scale='linear', regression='pearson'):
        super().__init__(num_samples)
        linreg_axes = []
        for row in range(num_samples-1):
            for col in range(row, num_samples-1):
                linreg_axes.append(self.axes[row, col])
        self.axes[0, 0].get_shared_x_axes().join(*linreg_axes)
        self.axes[0, 0].get_shared_y_axes().join(*linreg_axes)
        self.lims = [0, 0]
        self.scale = scale
        self.regression = {
            'pearson': stats.pearsonr,
            'spearman': stats.spearmanr}[regression]
        self.profiles = []
        self.colors = []
        self.labels = []

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=None,
                        width_ax_in=2, height_ax_in=2,
                        height_gap_in=0.3, width_gap_in=0.3,
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
        return (self.length-1, self.length-1)

    def plot_data(self, structure, profile, annotations, label,
                  column=None, colors="sequence"):
        if column is None:
            column = profile.metric
        self.labels.append(label)
        self.colors.append(
            profile.get_colors(colors, profile=profile, structure=structure))
        self.profiles.append(profile.data[column].to_numpy(copy=True))
        self.add_legend(colors, annotations)
        if len(self.profiles) == self.length:
            self.finalize()

    def finalize(self):
        # plot the regressions
        for row in range(self.length-1):
            for col in range(self.length-1):
                if row <= col:
                    self.plot_regression(i=row, j=col)
        # # change the plotting buffer for linear scale
        if self.scale == 'linear':
            buffer = 0.05 * (self.lims[1] - self.lims[0])
            self.lims[0] -= buffer
            self.lims[1] += buffer
            self.axes[0, 0].set(
                ylim=self.lims,
                xlim=self.lims)
        # format axis spines and labels
        ticks = self.axes[0, 0].get_yticks()[1:-1]
        for row in range(self.length-1):
            for col in range(self.length-1):
                ax = self.axes[row, col]
                if row <= col:
                    sns.despine(ax=ax)
                    ax.set_xticks(ticks)
                    ax.set_yticks(ticks)
                else:
                    self.axes[row, col].set_axis_off()
                if row < col:
                    ax.set_xticklabels([' ']*len(ticks))
                    ax.set_yticklabels([' ']*len(ticks))
                if row == 0:
                    ax.annotate(
                        self.labels[col+1], xy=(0.5, 1), xytext=(0, 5),
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='center', va='baseline')
                if col == (self.length - 2):
                    ax.annotate(
                        self.labels[row], xy=(1, 0.5), xytext=(5, 0),
                        rotation=-90,
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='left', va='center')

    def add_legend(self, colorby, annotations):
        handles = []
        if colorby == 'sequence':
            for label in ['A', 'U', 'C', 'G']:
                color = styles.get_nt_color(label)
                handles.append(mpl.lines.Line2D(
                    [], [], color=color, marker='o', linestyle='None',
                    markersize=10, label=label))
        elif colorby == 'structure':
            for label, color in zip(['paired', 'unpaired'], ['C0', 'C1']):
                handles.append(mpl.lines.Line2D(
                    [], [], color=color, marker='o', linestyle='None',
                    markersize=10, label=label))
        elif colorby == 'annotations':
            for annotation in annotations:
                label = annotation.name
                color = annotation.color
                handles.append(mpl.lines.Line2D(
                    [], [], color=color, marker='.', linestyle='None',
                    markersize=10, label=label))
        self.axes[1, 0].legend(handles=handles, loc=10, title=colorby)

    def plot_regression(self, i, j):
        ax = self.axes[i, j]
        p1 = self.profiles[i]
        p2 = self.profiles[j+1]
        colors = self.colors[i]

        # ax.plot([0, 1], [0, 1], color='black')
        if self.scale == 'log':
            notNans = (p1 > 0) & (p2 > 0)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            gradient, _, _, _, _ = stats.linregress(np.log10(p2), np.log10(p1))
            r_value, p_value = self.regression(np.log10(p2), np.log10(p1))
        if self.scale == 'linear':
            notNans = ~np.isnan(p1) & ~np.isnan(p2)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            gradient, _, _, _, _ = stats.linregress(p2, p1)
            r_value, p_value = self.regression(p2, p1)
        self.lims[0] = min([self.lims[0], min(p1), min(p2)])
        self.lims[1] = max([self.lims[1], max(p1), max(p2)])
        ax.text(
            0.1, 0.75,
            f'R^2: {r_value**2:.2f}\nslope: {gradient:.2f}\np: {p_value:.4f}',
            transform=ax.transAxes)
        ax.scatter(p2, p1, c=colors, marker='.')
        ax.set(xscale=self.scale,
               yscale=self.scale)

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
