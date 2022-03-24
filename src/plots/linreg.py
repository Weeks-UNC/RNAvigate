from scipy import stats
import seaborn as sns
import numpy as np
from .plots import Plot


class LinReg(Plot):
    def __init__(self, num_samples):
        super().__init__(num_samples)
        self.paired = None
        self.profiles = []
        self.sequences = []
        self.labels = []
        self.pass_through = ["colorby"]

    def get_rows_columns(self):
        return (self.length, self.length)

    def get_figsize(self):
        return (7*self.columns, 7*self.rows)

    def plot_data(self, ct, profile, label, colorby="structure"):
        self.labels.append(label)
        seq = profile.data["Sequence"].values
        seq = [x for x in seq if x.isupper()]
        profile = profile.data["Reactivity_profile"].copy()
        self.paired = np.array(ct.ct) != 0
        if self.paired is None and colorby == "structure":
            am = profile.get_alignment_map(ct)
            prof = np.full(ct.length, np.nan)
            seq = ct.sequence
            for i1, i2 in enumerate(am):
                if i2 != -1:
                    prof[i1] = profile[i2]
        else:
            prof = profile.values
        self.profiles.append(prof)
        self.sequences.append(seq)
        if len(self.profiles) == self.length:
            for i in range(self.length):
                self.plot_kde(i)
                for j in range(self.length):
                    if i < j:
                        self.plot_regression(i, j, colorby)
                    elif i > j:
                        self.axes[i, j].set_axis_off()
            handles, labels = self.axes[0, 1].get_legend_handles_labels()
            self.axes[2, 0].legend(handles, labels, loc=10,
                                   title=colorby.capitalize())
            handles, labels = self.axes[0, 0].get_legend_handles_labels()
            self.axes[1, 0].legend(handles, labels, loc=10,
                                   title="Structure")

    def plot_regression(self, i, j, colorby):
        ax = self.axes[i, j]
        p1 = self.profiles[i]
        seq1 = self.sequences[i]
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
        if self.paired is not None and colorby == "structure":
            paired_mask = self.paired
            paired_mask = paired_mask[notNans]
            unpaired_mask = ~paired_mask
            ax.scatter(p1[paired_mask], p2[paired_mask], label="Paired")
            ax.scatter(p1[unpaired_mask], p2[unpaired_mask], label="Unpaired")
        elif colorby == "sequence":
            color_dict = {"A": "#f20000", "U": "#f28f00",
                          "G": "#00509d", "C": "#00c200"}
            for nuc in "GUAC":
                nuc_mask = [nt == nuc for nt in seq1]
                color = color_dict[nuc]
                ax.scatter(p1[nuc_mask], p2[nuc_mask], label=nuc, color=color)
        else:
            ax.scatter(p1, p2)
        ax.set(xscale='log',
               xlim=[0.00001, 0.3],
               yscale='log',
               ylim=[0.00001, 0.3],
               title=f'{s1} vs. {s2}: {column}')

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
