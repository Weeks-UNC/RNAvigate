from rnavigate.plots import AP
from rnavigate.plots.plots import adjust_spines, clip_spines
from sklearn.metrics import roc_curve, auc
import numpy as np


class WindowedAUROC():
    """Class for computing and displaying windowed AUROC analysis. This
    analysis computes the ROC curve over a sliding window for the performance
    of per-nucleotide data (usually SHAPE-MaP or DMS-MaP Normalized reactivity)
    in predicting the base-pairing status of each nucleotide. The area under
    this curve (AUROC) is displayed compared to the median across the RNA.
    Below, an arc plot displays the secondary structure and per-nucleotide
    profile.

    Methods:
        __init__: Computes the AUROC array and AUROC median.
            By default, calls plot_auroc.
        plot_auroc: Displays the AUROC analysis over the given region.
            Returns Plot object

    Attributes:
        sample: an rnavigate.Sample to retrieve profile and secondary structure
        ct: sample.data[ct]
        profile: sample.data[profile]
        sequence: the sequence string of sample.data[ct]
        window: the size of the windows
        nt_length: the length of sequence string
        auroc: the auroc numpy array, length = nt_length, padded with np.nan
        median_auroc: the median of the auroc array
    """

    def __init__(self, sample, pad=40, region=None, profile="profile", ct="ct",
                 show=True):
        """Compute the AUROC for all windows. AUROC is a measure of how well a
        reactivity profile predicts paired vs. unpaired nucleotide status.

        Args:
            sample (rnav.Sample): Your rnavigate sample
            pad (int, optional): number of nucleotides on either side of a
                position to include in window
                Defaults to 40 (window=81).
            region (list of int: length 2, optional): Passed to self.plot_auroc
                Defaults to [1, RNA length].
            show (bool, optional): Creates a plot by calling self.plot_auroc.
                Defaults to True.
        """
        # ensure sample contains profile and ct data
        for data in [profile, ct]:
            assert data in sample.data.keys(), f"Sample missing {data} data"

        # store basic information
        self.sample = sample
        self.ct = sample.data[ct]
        self.profile = sample.data[profile]
        self.sequence = self.ct.sequence
        self.window = pad * 2 + 1
        self.nt_length = self.ct.length

        # get Norm_profile array and ct array
        profile = self.profile.data["Norm_profile"].values
        paired_nts = self.ct.ct

        # for each possible window: compute auroc and populate array
        self.auroc = np.full(len(profile), np.nan)
        for i in range(pad, len(profile)-pad):
            # get profile and ct values within window
            win_profile = profile[i-pad:i+pad+1]
            win_ct = paired_nts[i-pad:i+pad+1]
            # ignore positions where profile is nan
            valid = ~np.isnan(win_profile)
            # y: classification (paired or unpaired)
            y = win_ct[valid] == 0
            scores = win_profile[valid]
            # skip this window if there are less than 10 paired or unpaired nts
            if (sum(y) < 10) or (sum(~y) < 10):
                continue
            # add window auroc to array
            tpr, fpr, _ = roc_curve(y, scores)
            self.auroc[i] = auc(tpr, fpr)

        self.auroc_median = np.nanmedian(self.auroc)

        if show:
            self.plot_auroc(region=region)

    def plot_auroc(self, region=None):
        """Plot the result of the windowed AUROC analysis, with arc plot of
        structure and reactivity profile.

        Args:
            region (list of int: length 2, optional): Start and end nucleotide
                positions to plot. Defaults to [1, RNA length].
        """
        if region is None:
            start = 1
            stop = self.nt_length
            region = [start, stop]
            region_length = self.nt_length
        else:
            start, stop = region
            region_length = stop - start + 1

        plot = AP(1, region_length, cols=1, rows=1, region=region)
        ax = plot.axes[0, 0]

        # fill between auroc values and median, using secondary ax
        x_values = np.arange(start, stop + 1)
        auc_ax = ax.twinx()
        auc_ax.set_ylim(0.5, 1.6)
        auc_ax.set_yticks([0.5, self.auroc_median, 1.0])
        auc_ax.fill_between(x_values, self.auroc[start-1:stop],
                            self.auroc_median, fc='0.3', lw=0)
        adjust_spines(auc_ax, ["left"])
        clip_spines(auc_ax, ["left"])

        # add ct and reactivity profile track
        plot.plot_data(
            ct=self.ct,
            comp=None,
            interactions=None,
            interactions2=None,
            profile=self.profile,
            label="label",
            seqbar=False,
            title=False,
            annotations=[],
            plot_error=False)

        # Place Track Labels
        ax.set_title(f"{self.sample.sample}\n{start} - {stop}", loc='left',
                     fontdict={"fontsize": 48})
        ax.text(1.002, 6/8, "Secondary\nStructure",
                transform=ax.transAxes, fontsize=36, va='center')
        ax.text(1.002, 2/8, f"{self.window}-nt window\nAUROC",
                transform=ax.transAxes, va='center', fontsize=36)

        # limits, ticks, spines, and grid
        ax.set_ylim([-305, 315])
        ax.set_xticks(ticks=[x for x in range(500, stop, 500) if x > start])
        ax.set_xticks(ticks=[x for x in range(100, stop, 100) if x > start],
                      minor=True)
        ax.tick_params(ax='x', which='major', labelsize=36)
        adjust_spines(ax, ['bottom'])
        ax.grid(ax='x')

        # set figure size so that 100 ax units == 1 inch
        plot.set_figure_size(height_ax_rel=1/100, width_ax_rel=1/100)
        return plot
