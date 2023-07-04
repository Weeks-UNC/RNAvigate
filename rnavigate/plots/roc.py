import numpy as np
from rnavigate import plots
from sklearn.metrics import roc_curve, auc


class ROC(plots.Plot):
    def __init__(self, num_samples):
        super().__init__(num_samples)
        self.a_ax = self.axes[0, 2]
        self.u_ax = self.axes[0, 3]
        self.g_ax = self.axes[1, 2]
        self.c_ax = self.axes[1, 3]
        gs = self.axes[0, 0].get_gridspec()
        for ax in self.axes[:, :2].flatten():
            ax.remove()
        self.main_ax = self.fig.add_subplot(gs[:, :2])
        self.pass_through = []

    def set_figure_size(self, fig=None, ax=None,
                        rows=None, cols=None,
                        height_ax_rel=None, width_ax_rel=None,
                        width_ax_in=6, height_ax_in=6,
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
        return (2, 4)

    def get_figsize(self):
        return (28, 14)

    def plot_data(self, ct, profile, label):
        self.i += 1

        metric = profile.metric
        valid = ~np.isnan(profile.data[metric])
        y = ct.ct[valid] == 0
        scores = profile.data.loc[valid, metric]
        tpr, fpr, _ = roc_curve(y, scores)
        auc_score = auc(tpr, fpr)
        self.main_ax.plot(tpr, fpr, label=f"{label}: AUC={auc_score:.2f}")
        self.main_ax.plot([0, 1], [0, 1], "k--")

        axes = [self.a_ax, self.u_ax, self.c_ax, self.g_ax]
        for ax, nt in zip(axes, 'AUCG'):
            ax.plot([0, 1], [0, 1], "k--")
            ax.set(title=nt,
                   aspect='equal')
            valid = ~np.isnan(profile.data[metric])
            valid &= profile.data["Sequence"] == nt
            y = ct.ct[valid] == 0
            scores = profile.data.loc[valid, metric]
            tpr, fpr, _ = roc_curve(y, scores)
            auc_score = auc(tpr, fpr)
            ax.plot(tpr, fpr, label=f"AUC={auc_score:.2f}")

        if self.i == self.length:
            self.main_ax.legend(loc=4)
            self.main_ax.set(title="Receiver Operator Characteristic",
                             ylabel="True Positive Rate",
                             xlabel="False Positive Rate",
                             aspect='equal')
            for ax in axes:
                ax.set(yticks=[], xticks=[])
                ax.legend()
