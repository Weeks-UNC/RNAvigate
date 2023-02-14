from rnavigate.plots import AP
from rnavigate.plots.plots import adjust_spines, clip_spines
from rnavigate.data import Annotation
import matplotlib.pyplot as plt
import numpy as np


class LowSS():
    def __init__(self, sample, window=55, region=None, show=True):
        for data in ["profile", "pairprob", "ct"]:
            assert data in sample.data.keys(), f"Sample missing {data} data"
        assert window % 2 == 1, "Window must be an odd number."
        self.sample = sample
        self.sequence = self.sample.data["ct"].sequence
        self.window = window
        self.nt_length = sample.data["ct"].length

        profile = sample.data["profile"].data["Norm_profile"].values
        self.median_profile = np.median(profile[~np.isnan(profile)])
        self.windowed_profile = self.windowed_median(profile)

        sample.data["pairprob"].set_entropy()
        entropies = sample.data["pairprob"].entropy
        self.median_entropy = np.median(entropies)
        self.windowed_entropy = self.windowed_median(entropies)

        self.in_lowss_region = np.zeros(self.nt_length, dtype=int)
        lowss_regions = []
        lowss_region = None
        for i, (entropy, profile) in enumerate(zip(self.windowed_entropy,
                                                   self.windowed_profile)):
            if (entropy < 0.08) & (profile < self.median_profile):
                start = max(0, i-(self.window//2))
                stop = min(self.nt_length, i+(self.window//2)+1)
                if lowss_region is None:
                    lowss_region = [start, stop]
                elif lowss_region[0] < start < lowss_region[1]:
                    lowss_region[1] = stop
                elif start > lowss_region[1]:
                    lowss_regions.append(lowss_region)
                    self.in_lowss_region[start:stop] = 1
                    lowss_region = [start, stop]
        self.lowss_regions = Annotation(span_list=lowss_regions, color="grey",
                                        sequence=self.sequence)
        sample.data["lowss"] = self.lowss_regions
        if show:
            self.plot_LowSS(region=region)

    def plot_LowSS(self, region=None):
        if region is None:
            start = 1
            stop = self.nt_length
            region = [start, stop]
            region_length = self.nt_length
        elif isinstance(region, int):
            region = self.lowss_regions[region-1]
            start, stop = region
            start = max(start - 150, 1)
            stop = min(stop + 150, self.nt_length)
            region_length = stop - start + 1
        else:
            start, stop = region
            region_length = stop - start + 1

        plot = AP(1, region_length, cols=1, rows=1, region=region)
        ax = plot.axes[0, 0]

        x_values = np.arange(start, stop + 1)
        prof_ax = ax.twinx()
        prof_ax.set_ylim(-3, 1)
        prof_ax.set_yticks([0.0, 0.4, 0.85])
        prof_ax.fill_between(x_values, [self.median_profile]*region_length,
                             self.windowed_profile[start-1:stop],
                             fc='0.3', lw=0)
        adjust_spines(prof_ax, ["left"])
        clip_spines(prof_ax, ["left"])

        ent_ax = ax.twinx()
        ent_ax.set_ylim(-1.5, 1.5)
        ent_ax.set_yticks(ticks=[0, 0.08, 0.5])
        ent_ax.fill_between(x_values, [0.08]*region_length,
                            self.windowed_entropy[start-1:stop],
                            fc='C1', lw=0)
        adjust_spines(ent_ax, ["left"])
        clip_spines(ent_ax, ["left"])

        # add ct and pairing probabilities track
        self.sample.filter_interactions("pairprob", "pairprob")
        plot.plot_data(ct=self.sample.data["ct"], comp=None,
                       interactions=self.sample.data["pairprob"],
                       interactions2=None, profile=None, label="label",
                       colorbar=False, seqbar=False, title=False,
                       annotations=[self.lowss_regions],
                       annotation_mode="vbar")

        # Place Track Labels
        ax.set_title(f"{self.sample.sample}\n{start} - {stop}", loc='left',
                     fontdict={"fontsize": 48})
        ax.text(1.002, 7/8, "SHAPE\nReactivity",
                fontsize=36, transform=ax.transAxes, va='center')
        ax.text(1.002, 5/8, "Shannon\nEntropy",
                fontsize=36, transform=ax.transAxes, va='center')
        ax.text(1.002, 3/8, "Secondary\nStructure",
                transform=ax.transAxes, fontsize=36, va='center')
        ax.text(1.002, 1/8, "Pairing\nProbability",
                transform=ax.transAxes, va='center', fontsize=36)

        # Place region labels
        for i, (lssr_start, lssr_stop) in enumerate(self.lowss_regions):
            middle = (lssr_stop + lssr_start) / 2
            if start < middle < stop:
                ax.text(middle, 550, str(i+1),
                        ha='center', va='center', fontsize=36)

        ax.set_ylim([-305, 915])
        ax.set_xticks(ticks=[x for x in range(500, stop+1, 500) if x > start])
        ax.set_xticks(ticks=[x for x in range(100, stop+1, 100) if x > start],
                      minor=True)
        ax.tick_params(axis='x', which='major', labelsize=36)
        adjust_spines(ax, ['bottom'])

        # set figure size by axis size + current figure margins
        l_ax, r_ax = ax.get_xlim()
        b_ax, t_ax = ax.get_ylim()
        w_ax = (r_ax - l_ax) / 100
        h_ax = (t_ax - b_ax) / 100
        l_fig = ax.figure.subplotpars.left
        r_fig = ax.figure.subplotpars.right
        t_fig = ax.figure.subplotpars.top
        b_fig = ax.figure.subplotpars.bottom
        figw = w_ax/(r_fig-l_fig)
        figh = h_ax/(t_fig-b_fig)
        ax.figure.set_size_inches(figw, figh)
        return plot

    def windowed_median(self, data):
        pads = [np.nan]*(self.window//2)
        win_median = []
        for i in range(len(data) - self.window+1):
            win_data = data[i:i+self.window]
            win_data = win_data[~np.isnan(win_data)]
            win_median.append(np.median(win_data))
        return np.array(pads + win_median + pads)
