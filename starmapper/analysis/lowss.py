from starmapper.plots import AP
import matplotlib.pyplot as plt
import numpy as np


class LowSS():
    def __init__(self, sample, window=55, region="all"):
        for data in ["profile", "probs", "ct"]:
            assert data in sample.data.keys(), f"Sample missing {data} data"
        assert window % 2 == 1, "Window must be an odd number."
        self.sample = sample
        self.window = window
        self.region = region
        self.nt_length = sample.data["ct"].length

        profile = sample.data["profile"].data["Norm_profile"].values
        self.median_profile = np.median(profile[~np.isnan(profile)])
        self.windowed_profile = self.windowed_median(profile)

        sample.data["probs"].set_entropy()
        entropies = sample.data["probs"].entropy
        self.windowed_entropy = self.windowed_median(entropies)

        self.in_lssr = np.zeros(self.nt_length, dtype=int)
        self.lowss_regions = []
        lssr = None
        for i, (ent, prof) in enumerate(zip(self.windowed_entropy,
                                            self.windowed_profile)):
            if (ent < 0.08) & (prof < self.median_profile):
                start = max(0, i-(self.window//2))
                stop = min(self.nt_length, i+(self.window//2)+1)
                if lssr is None:
                    lssr = [start, stop]
                elif lssr[0] < start < lssr[1]:
                    lssr[1] = stop
                elif start > lssr[1]:
                    self.lowss_regions.append(lssr)
                    lssr = [start, stop]
                self.in_lssr[start:stop] = 1

        self.plot_LowSS()

    def plot_LowSS(self):
        figsize = (self.nt_length/100, 12.06)
        plot = AP(1, self.nt_length, cols=1, rows=1, figsize=figsize,
                  region=self.region)
        region = np.s_[self.region[0]-1: self.region[1]]
        ax = plot.axes[0, 0]
        ax.set_ylim([-305, 901])
        ax.set_xticks(list(range(500, self.region[1], 500)))
        ax.tick_params(axis='x', which='major', labelsize=36)
        ax.spines['bottom'].set_position(('outward', 2))

        # Plot windowed profile and windowed entropy
        x_values = np.arange(self.region[0], self.region[1]+1)
        ax.fill_between(x_values, [self.median_profile*350+600]*self.nt_length,
                        self.windowed_profile[region]*350+600, fc='0.3')
        ax.fill_between(x_values, [300]*self.nt_length,
                        self.windowed_entropy[region]*600 + 300, fc='C1')

        # add shaded vertical bars over LowSS regions

        xvals = (np.arange(self.nt_length)+1)/self.nt_length
        ax.fill_between(xvals, [0]*self.nt_length, self.in_lssr, alpha=0.2,
                        fc='grey', transform=ax.transAxes)

        # add ct and pairing probabilities track
        self.sample.filter_ij("probs", "probs")
        plot.add_sample(self.sample, ct="ct", comp=None, ij="probs",
                        ij2=None, profile=None, label="label",
                        colorbar=False, seqbar=False, title=False)

        # Place Track Labels
        ax.text(0.002, 0.95, f"{self.sample.sample} {self.region}",
                transform=ax.transAxes, fontsize=48)
        ax.text(1.002, 0.85, f"{self.window} nt median SHAPE Reactivity",
                fontsize=36, transform=ax.transAxes,
                verticalalignment='center')
        ax.text(1.002, 0.55, f"{self.window} nt median Shannon entropy",
                fontsize=36, transform=ax.transAxes,
                verticalalignment='center')
        ax.text(1.002, 0.35, "Secondary Structure", transform=ax.transAxes,
                verticalalignment='center', fontsize=36)
        ax.text(1.002, 0.15, "Pairing Probabilities", transform=ax.transAxes,
                verticalalignment='center', fontsize=36)

        plt.tight_layout()

    def windowed_median(self, data):
        pads = [np.nan]*(self.window//2)
        win_median = []
        for i in range(len(data) - self.window+1):
            win_data = data[i:i+self.window]
            win_data = win_data[~np.isnan(win_data)]
            win_median.append(np.median(win_data))
        return np.array(pads + win_median + pads)
