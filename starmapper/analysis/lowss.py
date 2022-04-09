from starmapper.plots import AP
import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics


def LowSS(sample):
    for data in ["profile", "probs", "ct"]:
        assert data in sample.data.keys(), f"Sample missing {data} data"
    plot = AP(2, sample.data["profile"].length, cols=1, rows=2)
    sample.data["probs"].filter(sample.data["probs"])
    plot.add_sample(sample, ct="ct", comp=None, ij="probs",
                    ij2=None, profile="profile", label="label", ax=1)
    plot.axes[0, 0].remove()
    entropy_ax = plot.fig.add_axes([0.02, 0.5, 0.96, 0.1666])
    entropy_ax.set_xlim(0.01, plot.nt_length+1)
    prof_ax = plot.fig.add_axes(
        [0.02, 0.6666, 0.96, 0.1666], sharex=entropy_ax)
    auroc_ax = plot.fig.add_axes(
        [0.02, 0.8333, 0.96, 0.1666], sharex=entropy_ax)

    w = 51

    def windowed_median(l, w):
        assert w % 2 == 1, "Window must be odd number"
        pads = [np.nan]*(w//2)
        win_median = []
        for i in range(len(l) - w+1):
            win_l = l[i:i+w]
            win_l = win_l[~np.isnan(win_l)]
            win_median.append(np.median(win_l))
        return np.array(pads + win_median + pads)

    x_values = range(sample.data["profile"].length)

    profile = sample.data["profile"].data["Norm_profile"].values
    values = profile[~np.isnan(profile)]
    prof_ax.set_ylim(0, 1)
    prof_ax.set_yticks([0, 0.4, 0.8])
    win_prof = windowed_median(profile, w)
    prof_ax.fill_between(x_values, [np.median(values)]
                         * len(profile), win_prof, fc='0.3')
    prof_ax.text("51 nt median SHAPE")

    unpaired = np.array([nt == 0 for nt in sample.data["ct"].ct])
    pads = [np.nan] * (w // 2)
    win_roc_auc = []
    for i in range(len(profile) - w+1):
        pred = profile[i:i+w]
        y = unpaired[i:i+w][~np.isnan(pred)]
        pred = pred[~np.isnan(pred)]
        win_roc_auc.append(metrics.roc_auc_score(y, pred))
    win_roc_auc = np.array(pads + win_roc_auc + pads)
    auroc_ax.fill_between(x_values, [0.5]*len(profile), win_roc_auc, fc='C0')
    auroc_ax.set_ylim(0.5, 1)

    sample.data["probs"].set_entropy()
    entropies = sample.data["probs"].entropy
    entropy_ax.set_ylim(0, 0.4)
    entropy_ax.set_yticks([0.0, 0.1, 0.2, 0.3])
    win_ents = windowed_median(entropies, w)
    entropy_ax.fill_between(x_values, [0]*len(entropies), win_ents, fc='C1')

    yvals = np.zeros(plot.nt_length, dtype=int)
    for i, (ent, prof) in enumerate(zip(win_ents, win_prof)):
        if (ent < 0.08) & (prof < 0.3):
            start = max(0, i-(w//2))
            stop = min(plot.nt_length, i+(w//2)+1)
            yvals[start:stop] = 1
    xvals = sample.data["profile"].data["Nucleotide"]/plot.nt_length
    for ax in [plot.axes[1, 0], prof_ax, entropy_ax, auroc_ax]:
        ax.bar(xvals, yvals, 1/plot.nt_length, linewidth=0,
               alpha=0.2, color='grey', transform=ax.transAxes)

    plt.tight_layout()
