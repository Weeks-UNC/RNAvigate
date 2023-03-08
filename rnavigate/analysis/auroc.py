from rnavigate.plots import AP
from rnavigate.plots.plots import adjust_spines, clip_spines
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import numpy as np


class WindowedAUROC():
    def __init__(self, sample, window=81, region=None, show=True):
        for data in ["profile", "ct"]:
            assert data in sample.data.keys(), f"Sample missing {data} data"
        assert window % 2 == 1, "Window must be an odd number."
        self.sample = sample
        self.sequence = self.sample.data["ct"].sequence
        self.window = window
        self.nt_length = sample.data["ct"].length

        profile = sample.data["profile"].data["Norm_profile"].values
        ct = sample.data["ct"].ct
        self.auroc = np.full(len(profile), np.nan)
        buffer = int((window - 1) / 2)

        for i in range(buffer, len(profile)-buffer):
            win_profile = profile[i-buffer:i+buffer+1]
            win_ct = ct[i-buffer:i+buffer+1]
            valid = ~np.isnan(win_profile)
            y = win_ct[valid] == 0
            scores = win_profile[valid]
            if (sum(y) < 10) or (sum(~y) < 10):
                continue
            tpr, fpr, _ = roc_curve(y, scores)
            self.auroc[i] = auc(tpr, fpr)
        self.auroc_median = np.nanmedian(self.auroc)

        if show:
            self.plot_LowSS(region=region)

    def plot_LowSS(self, region=None):
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

        x_values = np.arange(start, stop + 1)
        auc_ax = ax.twinx()
        auc_ax.set_ylim(0.5, 1.6)
        auc_ax.set_yticks([0.5, self.auroc_median, 1.0])
        auc_ax.fill_between(x_values, self.auroc[start-1:stop],
                            self.auroc_median, fc='0.3', lw=0)
        adjust_spines(auc_ax, ["left"])
        clip_spines(auc_ax, ["left"])

        # add ct and pairing probabilities track
        plot.plot_data(ct=self.sample.data["ct"], comp=None,
                       interactions=None, interactions2=None,
                       profile=self.sample.data["profile"], label="label",
                       colorbar=False, seqbar=False, title=False,
                       annotations=[], plot_error=False)

        # Place Track Labels
        ax.set_title(f"{self.sample.sample}\n{start} - {stop}", loc='left',
                     fontdict={"fontsize": 48})
        ax.text(1.002, 6/8, "Secondary\nStructure",
                transform=ax.transAxes, fontsize=36, va='center')
        ax.text(1.002, 2/8, f"{self.window}-nt window\nAUROC",
                transform=ax.transAxes, va='center', fontsize=36)

        ax.set_ylim([-305, 315])
        ax.set_xticks(ticks=[x for x in range(500, stop, 500) if x > start])
        ax.set_xticks(ticks=[x for x in range(100, stop, 100) if x > start],
                      minor=True)
        ax.tick_params(axis='x', which='major', labelsize=36)
        adjust_spines(ax, ['bottom'])
        ax.grid(axis='x')

        # set figure size so that 100 axis units == 1 inch
        plot.set_figure_size(yscale=1/100, xscale=1/100)
        return plot
