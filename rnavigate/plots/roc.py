import numpy as np
from .plots import Plot
from sklearn.metrics import roc_curve, auc


class ROC(Plot):
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

    def get_rows_columns(self, rows=None, cols=None):
        return (2, 4)

    def get_figsize(self):
        return (28, 14)

    def plot_data(self, ct, profile, label):
        self.i += 1

        valid = ~np.isnan(profile.data["Norm_profile"])
        y = ct.ct[valid] == 0
        scores = profile.data.loc[valid, "Norm_profile"]
        tpr, fpr, _ = roc_curve(y, scores)
        auc_score = auc(tpr, fpr)
        self.main_ax.plot(tpr, fpr, label=f"{label}: AUC={auc_score:.2f}")
        self.main_ax.plot([0, 1], [0, 1], "k--")

        axes = [self.a_ax, self.u_ax, self.c_ax, self.g_ax]
        for ax, nt in zip(axes, 'AUCG'):
            ax.plot([0, 1], [0, 1], "k--")
            ax.set(title=nt,
                   aspect='equal')
            valid = ~np.isnan(profile.data["Norm_profile"])
            valid &= profile.data["Sequence"] == nt
            y = ct.ct[valid] == 0
            scores = profile.data.loc[valid, "Norm_profile"]
            tpr, fpr, _ = roc_curve(y, scores)
            auc_score = auc(tpr, fpr)
            ax.plot(tpr, fpr, label=f"AUC={auc_score:.2f}")

        if self.i == self.length:
            self.main_ax.legend()
            self.main_ax.set(title="Receiver Operator Characteristic",
                             ylabel="True Positive Rate",
                             xlabel="False Positive Rate",
                             aspect='equal')
            for ax in axes:
                ax.set(yticks=[], xticks=[])
                ax.legend()
