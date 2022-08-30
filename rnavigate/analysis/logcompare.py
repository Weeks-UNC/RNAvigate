import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from ..plots import Plot


class LogCompare():
    def __init__(self, samples1, samples2, name1, name2, data="profile"):
        """Takes replicates of two samples for comparison. Replicates are
        required. Calculates the log division profile
        (log10(modified/untreated)) and minimizes the median of the absolute
        difference between these. Finally, creates a plot.

        Args:
            samples1 (list of Sample objects): Replicates of the first sample
            samples2 (list of Sample objects): Replicates of the second sample
            name1 (string): name of first sample
            name2 (string): name of second sample
            data (str, optional): Datatype to compare. Defaults to "profile".
        """
        self.data = data
        self.groups = {}
        self.load_replicates(*samples1, group=1)
        self.load_replicates(*samples2, group=2)
        sequences_match = self.groups[1]["seq"] == self.groups[2]["seq"]
        assert sequences_match, "Sample sequences do not match."
        self.groups[2][self.data] = self.rescale(self.groups[2][self.data],
                                                 self.groups[1][self.data])
        self.groups[1]["name"] = name1
        self.groups[2]["name"] = name2
        self.make_plots()

    def load_profile(self, sample):
        df = sample.data[self.data].data
        plus = df.Modified_rate.values.copy()
        minus = df.Untreated_rate.values.copy()
        profile = plus/minus
        profile = np.log(profile)
        profile[minus > 0.05] = np.nan
        return profile, sample.data[self.data].sequence

    def calc_scale_factor(self, profile, target_profile):
        def f(offset):
            return np.nanmedian(np.abs(profile + offset - target_profile))
        result = minimize_scalar(f, bounds=[-20, 20])
        offset = result.x
        return offset

    def rescale(self, profile, target_profile):
        offset = self.calc_scale_factor(profile, target_profile)
        return profile+offset

    def load_replicates(self, *samples, group):
        profile1, seq1 = self.load_profile(samples[0])
        primermask = np.array([c.islower() for c in seq1])
        profile1[primermask] = np.nan

        profiles = [profile1]
        for sample in samples[1:]:
            profile, seq = self.load_profile(sample)
            assert seq == seq1, "Replicate sequences do not match."
            profile[primermask] = np.nan
            profile = self.rescale(profile, profile1)
            profiles.append(profile)
        stacked = np.vstack(profiles)
        avgprofile = np.nanmean(stacked, axis=0)
        stderr = np.std(stacked, axis=0)
        self.groups[group] = {self.data: avgprofile,
                              "stderr": stderr,
                              "stacked": stacked,
                              "seq": seq}

    def plotseq(self, sequence, ax=None):
        if ax is None:
            ax = plt.gca()
        x = list(range(1, len(sequence)+1))
        trans = ax.get_xaxis_transform()
        for i in range(len(sequence)):
            c = sequence[i]
            ax.text(x[i]-0.5, 0.003, c, fontsize=8, transform=trans,
                    va='bottom', ha='center', clip_on=True)

    def make_plots(self):
        prof1, stderr1, stack1, seq, name1 = self.groups[1].values()
        prof2, stderr2, stack2, seq, name2 = self.groups[2].values()

        x = np.array(list(range(1, len(prof1)+1)))

        _, axes = plt.subplots(2, 1, figsize=(30, 14), sharex=True)
        ax = axes[0]
        ax.step(x, prof1, label=name1, color="C0")
        ax.step(x, prof2, label=name2, color="C1")
        ax.fill_between(x, prof1-stderr1, prof1+stderr1,
                        step='pre', color='C0', alpha=0.25, lw=0)
        ax.fill_between(x, prof2-stderr2, prof2+stderr2,
                        step='pre', color='C1', alpha=0.25, lw=0)
        ax.legend(loc='upper right')
        ax.axhline(0, color='black', lw=1, zorder=0)
        ax.set_ylabel('ln(Mod/BG)')
        Plot.add_sequence(ax=ax, sequence=seq)

        stack1 = stack1+10
        rescale_dms = self.rescale(stack1, stack2)
        stack1 = rescale_dms
        diff = stack2-stack1
        meandiff = np.nanmean(diff, axis=0)
        std_err = np.std(diff, axis=0)

        # get colormapping for differences
        # blue, white, red gradient from -5 to +5 stderrs
        z_scores = meandiff/std_err
        z_scores = (z_scores + 5) / 10
        z_scores[z_scores > 1] = 1
        z_scores[z_scores < 0] = 0
        colormap = plt.get_cmap("bwr")

        ax = axes[1]
        ax.bar(x-0.5, meandiff, color=colormap(z_scores), lw=0)
        ax.axhline(0, color='black', lw=1)
        ax.set_xlabel('nucleotide')
        ax.set_ylabel(f'{name2}-{name1}')
        ax.errorbar(x-0.5, meandiff, yerr=std_err,
                    fmt='none', color='black', lw=1)
        Plot.add_sequence(ax=ax, sequence=seq)

        axin1 = ax.inset_axes([0.8, 0.1, 0.15, 0.15])
        Plot.view_colormap(ax=axin1, ticks=[0, 5, 10], values=[-5, 0, 5],
                           title="Z-score", cmap="bwr")
        plt.tight_layout()
