from scipy import stats
import seaborn as sns
import numpy as np
from .plots import Plot


class LinReg(Plot):
    def __init__(self, num_samples):
        super().__init__(num_samples)
        self.paired = None
        self.profiles = []
        self.labels = []

    def get_rows_columns(self, number_of_samples):
        return (number_of_samples, number_of_samples)

    def get_figsize(self):
        return (7*self.columns, 7*self.rows)

    def plot_data(self, ct, profile, label, colorby="sequence"):
        self.labels.append(label)
        if self.paired is None:
            self.paired = np.array(ct.ct) != 0
        am = profile.get_alignment_map(ct)
        profile = profile.data["Reactivity_profile"].copy()
        prof = np.full(ct.length, np.nan)
        for i1, i2 in enumerate(am):
            if i2 != -1:
                prof[i1] = profile[i2]
        self.profiles.append(prof)
        if len(self.profiles) == self.length:
            for i in range(self.length):
                self.plot_kde(i)
                for j in range(self.length):
                    if i == j:
                        continue
                    self.plot_regression(i, j, colorby)

    def plot_regression(self, i, j, colorby):
        ax = self.axes[i, j]
        p1 = self.profiles[i]
        p2 = self.profiles[j]
        s1 = self.labels[i]
        s2 = self.labels[j]
        column = "Reactivity"

        ax.plot([0, 1], [0, 1], color='black')
        notNans = ~np.isnan(p1) & ~np.isnan(p2)
        p1 = p1[notNans]
        p2 = p2[notNans]
        gradient, _, r_value, _, _ = stats.linregress(p1, p2)
        ax.text(0.1, 0.8, f'R^2 = {r_value**2:.2f}\nslope = {gradient:.2f}',
                transform=ax.transAxes)
        if self.paired is not None:
            paired_mask = self.paired
            paired_mask = paired_mask[notNans]
            unpaired_mask = ~paired_mask
            ax.scatter(p1[paired_mask], p2[paired_mask], label="Paired")
            ax.scatter(p1[unpaired_mask], p2[unpaired_mask], label="Unpaired")
        elif colorby == "sequence":
            for nuc in "GUAC":
                sequence = self.profile["Sequence"][notNans]
                nuc_mask = [nt == nuc for nt in sequence]
                color = get_nt_color(nuc)
                ax.scatter(p1[nuc_mask], p2[nuc_mask], label=nuc, color=color)
        else:
            ax.scatter(p1, p2)
        ax.set(xscale='log',
               xlim=[0.00001, 0.3],
               xlabel=s1,
               yscale='log',
               ylim=[0.00001, 0.3],
               ylabel=s2,
               title=f'{s1} vs. {s2}: {column}')
        ax.legend(title=colorby, loc="lower right")

    def plot_kde(self, i):
        ax = self.axes[i, i]
        paired = self.paired
        profile = self.profiles[i]
        label = self.labels[i]
        valid = profile > 0
        sns.kdeplot(profile[paired & valid], bw_adjust=0.6, shade=True,
                    label='SS', ax=ax, log_scale=True)
        sns.kdeplot(profile[~paired & valid], bw_adjust=0.6, shade=True,
                    label='DS', ax=ax, log_scale=True)
        ax.annotate(label, (0.1, 0.9), xycoords="axes fraction")
        ax.set(xlim=(0.00001, 0.5))
