from starmapper.plots import AP
import matplotlib.pyplot as plt
import numpy as np


def LowSS(sample):
    for data in ["profile", "probs", "ct"]:
        assert data in sample.data.keys(), f"Sample missing {data} data"
    nt_length = sample.data["profile"].length
    figsize = (nt_length/100, 12.06)
    plot = AP(1, nt_length, cols=1, rows=1, figsize=figsize)
    ax = plot.axes[0, 0]
    ax.set(ylim=(-305, 901))
    ax.spines['bottom'].set_position('center')

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
    win_prof = windowed_median(profile, w)
    ax.fill_between(x_values, [np.median(values)*350 + 600]
                    * len(profile), win_prof*350 + 600, fc='0.3')
    ax.text(5, 850, "51 nt median SHAPE")

    sample.data["probs"].set_entropy()
    entropies = sample.data["probs"].entropy
    win_ents = windowed_median(entropies, w)
    ax.fill_between(x_values, [300]*len(entropies),
                    win_ents*600 + 300, fc='C1')
    ax.text(5, 550, "51 nt median Shannon entropy")
    ax.text(5, 250, "DENV2 Genome - MFE & Pairing Probabilities")

    yvals = np.zeros(plot.nt_length, dtype=int)
    for i, (ent, prof) in enumerate(zip(win_ents, win_prof)):
        if (ent < 0.08) & (prof < 0.3):
            start = max(0, i-(w//2))
            stop = min(plot.nt_length, i+(w//2)+1)
            yvals[start:stop] = 1
    xvals = sample.data["profile"].data["Nucleotide"]/plot.nt_length
    ax.bar(xvals, yvals, 1/plot.nt_length, linewidth=0,
           alpha=0.2, color='grey', transform=ax.transAxes)

    sample.filter_ij("probs", "probs")
    plot.add_sample(sample, ct="ct", comp=None, ij="probs",
                    ij2=None, profile=None, label="label",
                    colormap=False, seqbar=False, title=False)

    plt.tight_layout()
