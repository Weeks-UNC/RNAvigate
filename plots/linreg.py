
def make_regression(self, comp_sample, ax=None, colorby="ct",
                    column="Reactivity_profile"):
    """Plots scatter plot of reactivity profile vs. reactivity profile and
    computes regression metrics R^2 and slope, which are annotated on the
    axis. Can color scatter plot by nucleotide or by paired status.

    Args:
        ax (pyplot axis): axis on which scatter plot appears
        comp_sample (plotmapper sample): sample to be compared
        colorby (str, optional): How to color the scatter plot. Options are
            "nucleotide" or "ct".
            Defaults to None
        column (str, optional): column from profile to be compared
            Defaults to "Reactivity_profile".
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(7, 7))
    p1 = self.profile[column].copy()
    p2 = comp_sample.profile[column].copy()

    ax.plot([0, 1], [0, 1], color='black')
    notNans = ~np.isnan(p1) & ~np.isnan(p2)
    p1 = p1[notNans]
    p2 = p2[notNans]
    gradient, _, r_value, _, _ = stats.linregress(p1, p2)
    ax.text(0.1, 0.8, f'R^2 = {r_value**2:.2f}\nslope = {gradient:.2f}',
            transform=ax.transAxes)
    if colorby == "ct":
        paired_list = self.ct.pairedResidueList()
        paired_mask = np.zeros(self.length["profile"], dtype=bool)
        for i in paired_list:
            paired_mask[i] = True
        paired_mask = paired_mask[notNans]
        unpaired_mask = ~paired_mask
        ax.scatter(p1[paired_mask], p2[paired_mask], label="Paired")
        ax.scatter(p1[unpaired_mask], p2[unpaired_mask], label="Unpaired")
    elif colorby == "nucleotide":
        for nuc in "GUAC":
            sequence = self.profile["Sequence"][notNans]
            nuc_mask = [nt == nuc for nt in sequence]
            color = get_nt_color(nuc)
            ax.scatter(p1[nuc_mask], p2[nuc_mask], label=nuc, color=color)
    else:
        ax.scatter(p1, p2)
    s1, s2 = self.sample, comp_sample.sample
    ax.set(xscale='log',
           xlim=[0.00001, 0.3],
           xlabel=s1,
           yscale='log',
           ylim=[0.00001, 0.3],
           ylabel=s2,
           title=f'{s1} vs. {s2}: {column}')
    ax.legend(title=colorby, loc="lower right")
