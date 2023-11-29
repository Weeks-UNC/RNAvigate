from scipy import stats
import seaborn as sns
import matplotlib as mpl
import numpy as np
from rnavigate import plots, styles


class LinReg(plots.Plot):
    def __init__(self, num_samples, scale='linear', regression='pearson',
                 kde=False, region='all'):
        super().__init__(num_samples)
        self.region = region
        linreg_axes = []
        for row in range(num_samples-1):
            for col in range(row, num_samples-1):
                linreg_axes.append(self.axes[row, col])
        self.axes[0, 0].get_shared_x_axes().join(*linreg_axes)
        self.axes[0, 0].get_shared_y_axes().join(*linreg_axes)
        self.lims = [0, 0]
        self.scale = scale
        self.kde = kde
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
        if self.region == 'all':
            start, end = 1, profile.length
        else:
            start, end = self.region
        if column is None:
            column = profile.metric
        self.labels.append(label)
        colors, colormap = profile.get_colors(
            colors, profile=profile, structure=structure,
            annotations=annotations,
            )
        self.colors.append(colors[start-1:end])
        self.add_colorbar_args(colormap)
        self.profiles.append(
            profile.data[column].to_numpy(copy=True)[start-1:end]
            )
        if len(self.profiles) == self.length:
            for row in range(self.length-1):
                for col in range(self.length-1):
                    if row <= col:
                        self.plot_regression(i=row, j=col)
            self.finalize()

    def finalize(self):
        # change the plotting buffer for linear scale
        if self.scale == 'log' and self.kde:
            # change the tick labels to 'fake' a log scale plot
            ticks = []
            for power in range(-5, 5):
                power = 10**power
                if self.lims[0] > power or power > self.lims[1]:
                    continue
                ticks.append(power)
                # 10^n formatted tick labels
            for i in range(len(self.axes)):
                self.axes[i, i].set(xticks=ticks, yticks=ticks)
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

    def plot_regression(self, i, j):
        ax = self.axes[i, j]
        p1 = self.profiles[j+1]
        p2 = self.profiles[i]
        colors = self.colors[i]

        # ax.plot([0, 1], [0, 1], color='black')
        if self.scale == 'log':
            notNans = (p1 > 0) & (p2 > 0)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            r_value, _ = self.regression(np.log10(p2), np.log10(p1))
        if self.scale == 'linear':
            notNans = ~np.isnan(p1) & ~np.isnan(p2)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            r_value, _ = self.regression(p2, p1)
        self.lims[0] = min([self.lims[0], min(p1), min(p2)])
        self.lims[1] = max([self.lims[1], max(p1), max(p2)])
        ax.text(0.1, 0.95, f'r = {r_value:.2f}', transform=ax.transAxes,
                ha='left', va='top',
                bbox=dict(fc='white', alpha=0.5, ec='black'))
        if self.kde:
            if self.scale == 'log':
                sns.kdeplot(
                    ax=ax, x=p1, y=p2, fill=True, log_scale=True,
                    levels=np.arange(1, 11)/10
                    )
            elif self.scale == 'linear':
                sns.kdeplot(
                    ax=ax, x=p1, y=p2, fill=True,
                    levels=np.arange(1, 11)/10
                    )
        else:
            ax.scatter(p2, p1, c=colors, marker='.')
            ax.set(xscale=self.scale, yscale=self.scale)
