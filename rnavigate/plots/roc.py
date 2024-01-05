import numpy as np
from rnavigate import plots
from sklearn.metrics import roc_curve, auc


class ROC(plots.Plot):
    """Plot ROC curves.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    **kwargs
        Keyword arguments passed to `rnavigate.plots.Plot`.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    a_ax : matplotlib.axes.Axes
        Axes object for the AUC plot of A nucleotides.
    u_ax : matplotlib.axes.Axes
        Axes object for the AUC plot of U nucleotides.
    g_ax : matplotlib.axes.Axes
        Axes object for the AUC plot of G nucleotides.
    c_ax : matplotlib.axes.Axes
        Axes object for the AUC plot of C nucleotides.
    main_ax : matplotlib.axes.Axes
        Axes object for the main ROC plot.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    i : int
        Index of the current plot.
    """

    def __init__(self, num_samples, **kwargs):
        """Initialize the plot."""
        super().__init__(num_samples, **kwargs)
        self.a_ax = self.axes[0, 2]
        self.u_ax = self.axes[0, 3]
        self.g_ax = self.axes[1, 2]
        self.c_ax = self.axes[1, 3]
        gs = self.axes[0, 0].get_gridspec()
        for ax in self.axes[:, :2].flatten():
            ax.remove()
        self.main_ax = self.fig.add_subplot(gs[:, :2])

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=None,
        width_ax_in=1.5,
        height_ax_in=1.5,
        height_gap_in=0.3,
        width_gap_in=0.2,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        """Set the figure size.

        Parameters
        ----------
        height_ax_rel : float
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float
            Width of the axes relative to the x-axis limits.
        width_ax_in : float
            Width of the axes in inches.
        height_ax_in : float
            Height of the axes in inches.
        height_gap_in : float
            Height of the gap between axes in inches.
        width_gap_in : float
            Width of the gap between axes in inches.
        top_in : float
            Height of the top margin in inches.
        bottom_in : float
            Height of the bottom margin in inches.
        left_in : float
            Width of the left margin in inches.
        right_in : float
            Width of the right margin in inches.
        """
        super().set_figure_size(
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

    def get_rows_columns(self, rows=None, cols=None):
        """Get the number of rows and columns.

        Parameters
        ----------
        rows : int
            Number of rows. This is ignored.
        cols : int
            Number of columns. This is ignored.

        Returns
        -------
        rows : int
            Number of rows. This is always 2.
        cols : int
            Number of columns. This is always 4.
        """
        return (2, 4)

    def plot_data(self, structure, profile, label, nts="AUCG"):
        """Plot the data.

        Parameters
        ----------
        structure : rnavigate.structure.Structure
            Structure object.
        profile : rnavigate.profile.Profile
            Profile object.
        label : str
            Sample name.
        nts : str, defaults to "AUCG"
            Which nucleotides to plot.
        """
        self.i += 1

        metric = profile.metric
        valid = ~profile.data[metric].isna()
        y = structure.pair_nts[valid] == 0
        scores = profile.data.loc[valid, metric]
        y = y.astype(int)
        scores = scores.astype(float)
        tpr, fpr, _ = roc_curve(y, scores)
        auc_score = auc(tpr, fpr)
        self.main_ax.plot(tpr, fpr, label=f"{label}: {auc_score:.2f}")
        self.main_ax.plot([0, 1], [0, 1], "k--")

        axes = {"A": self.a_ax, "U": self.u_ax, "C": self.c_ax, "G": self.g_ax}
        other_nts = {
            "A": ["a", "A"],
            "U": ["u", "U", "t", "T"],
            "C": ["c", "C"],
            "G": ["g", "G"],
        }
        for nt in nts:
            ax = axes[nt]
            ax.plot([0, 1], [0, 1], "k--")
            ax.set(title=nt, aspect="equal")
            valid = ~np.isnan(profile.data[metric])
            valid &= profile.data["Sequence"].isin(other_nts[nt])
            if sum(valid) == 0:
                print(f"RNAvigate warning: {profile} data is missing for {nt}")
                ax.plot([], [], label="N/A")
                continue
            y = structure.pair_nts[valid] == 0
            scores = profile.data.loc[valid, metric].to_numpy()
            y = y.astype(int)
            scores = scores.astype(float)
            tpr, fpr, _ = roc_curve(y, scores)
            auc_score = auc(tpr, fpr)
            ax.plot(tpr, fpr, label=f"{auc_score:.2f}")

        if self.i == self.length:
            self.main_ax.legend(title="Sample name: AUC", loc=4)
            self.main_ax.set(
                title="Receiver Operator Characteristic",
                ylabel="True Positive Rate",
                xlabel="False Positive Rate",
                aspect="equal",
            )
            for ax in axes.values():
                ax.set(yticks=[], xticks=[])
                ax.legend(title="AUC")
