from rnavigate.plots import AP
from rnavigate.plots.plots import adjust_spines, clip_spines
from rnavigate.data import Annotation
import matplotlib.pyplot as plt
import numpy as np


class LowSS():
    """Class for computing and displaying Low SHAPE, low Shannon entropy
    regions (LowSS) given a sample containing SHAPE reactivities, pairing
    probabilities, and MFE structure.

    Methods:
        __init__: performs the analysis
        plot_lowss: displays the result and returns plot object

    Attributes:
        sample (rnavigate.Sample): the sample to retreive data from.
        sequence (str): sequence retreived from sample.data["ct"].sequence
        window (int): size of the windows, must be odd
        nt_length (int): length of sequence string
        windowed_profile (numpy.array): windowed median SHAPE reactivities
        median_profile (float): global median SHAPE reactivity
        windowed_entropy (numpy.array): windowed median Shannon entropy
        median_entropy (float): global median Shannon entropy
        in_lowss_region (numpy.array: bool): if nucleotide is in LowSS region
        lowss_regions (rnavigate.Annotations): defines LowSS regions
    """

    def __init__(self, sample, window=55, region=None, show=True):
        """Perform Low SHAPE, low Shannon entropy analysis given a sample with
        1) reactivities 2) MFE structure 3) pairing probabilities. Optionally,
        plot the result over a given region.

        Args:
            sample (rnavigate.Sample): sample with profile, pairprob and ct
            window (int, optional): Window size for calculating median SHAPE
                and median Shannon entropy. Defaults to 55.
            region (list of str: len 2, optional): start and end positions to
                plot. Defaults to None.
            show (bool, optional): Whether to create a plot. Defaults to True.
        """
        # Make sure this sample contains the necessary data
        for data in ["profile", "pairprob", "ct"]:
            assert data in sample.data.keys(), f"Sample missing {data} data"
        assert window % 2 == 1, "Window must be an odd number."
        # TODO: these values could change between calculations and plotting
        # store values
        self.sample = sample
        self.sequence = self.sample.data["ct"].sequence
        self.window = window
        self.nt_length = sample.data["ct"].length

        # Calculate overall median SHAPE and windowed median SHAPEs
        profile = sample.data["profile"].data["Norm_profile"].values
        self.median_profile = np.median(profile[~np.isnan(profile)])
        self.windowed_profile = self.windowed_median(profile)

        # Calculate median and windowed shannon entropies
        sample.data["pairprob"].set_entropy()
        entropies = sample.data["pairprob"].entropy
        self.median_entropy = np.median(entropies)
        self.windowed_entropy = self.windowed_median(entropies)

        # Find LowSS regions
        self.in_lowss_region = np.zeros(self.nt_length, dtype=int)
        lowss_regions = []
        lowss_region = None
        for i, (entropy, profile) in enumerate(zip(self.windowed_entropy,
                                                   self.windowed_profile)):
            if (entropy < 0.08) & (profile < self.median_profile):
                start = max(1, i+1 - (self.window//2))
                stop = min(self.nt_length, i+(self.window//2)+1)
                if lowss_region is None:
                    lowss_region = [start, stop]
                elif start <= lowss_region[1]+1:
                    lowss_region[1] = stop
                elif start > lowss_region[1]:
                    lowss_regions.append(lowss_region)
                    self.in_lowss_region[start:stop] = 1
                    lowss_region = [start, stop]
        self.lowss_regions = Annotation(
            input_data=lowss_regions, annotation_type="spans", color="grey",
            sequence=self.sequence)
        # TODO: this could overwrite data already stored in sample
        sample.data["lowss"] = self.lowss_regions
        if show:
            self.plot_LowSS(region=region)

    def plot_LowSS(self, region=None, colorbar=True):
        """Visualize LowSS analysis over the given region.

        Args:
            region (list of int: len 2, optional): start and end positions.
                Defaults to None. Can also be an integer, in which case that
                lowSS region +/- 150 nts are shown.

        Returns:
            rnavigate.plots.AP: LowSS visualization
        """
        # show entire RNA if region is not provided
        if region is None:
            start = 1
            stop = self.nt_length
            region = [start, stop]
        # if region is an integer, get that lowSS region +/- 150 nts
        elif isinstance(region, int):
            region = self.lowss_regions[region-1]
            start, stop = region
            start = max(start - 150, 1)
            stop = min(stop + 150, self.nt_length)
        # else region should be a list of 2 integers, start and end.
        else:
            start, stop = region
        region_length = stop - start + 1

        # create arc plot instance
        plot = AP(1, region_length, cols=1, rows=1, region=region)
        axis = plot.axes[0, 0]
        x_values = np.arange(start, stop + 1)
        # plot median SHAPE on secondary axis
        prof_ax = axis.twinx()
        prof_ax.set_ylim(-3, 1)
        prof_ax.set_yticks([0.0, 0.4, 0.85])
        prof_ax.fill_between(x_values, [self.median_profile]*region_length,
                             self.windowed_profile[start-1:stop],
                             fc='0.3', lw=0)
        adjust_spines(prof_ax, ["left"])
        clip_spines(prof_ax, ["left"])

        # plot median entropy on second secondary axis
        ent_ax = axis.twinx()
        ent_ax.set_ylim(-1.5, 1.5)
        ent_ax.set_yticks(ticks=[0, 0.08, 0.5])
        ent_ax.fill_between(x_values, [0.08]*region_length,
                            self.windowed_entropy[start-1:stop],
                            fc='C1', lw=0)
        adjust_spines(ent_ax, ["left"])
        clip_spines(ent_ax, ["left"])

        # add ct and pairing probabilities track
        self.sample.filter_interactions("pairprob", "pairprob")
        plot.plot_data(
            seq=self.sample.data["ct"],
            ct=self.sample.data["ct"],
            comp=None,
            interactions=self.sample.data["pairprob"],
            interactions2=None,
            profile=None,
            label="label",
            seqbar=False,
            title=False,
            annotations=[self.lowss_regions],
            annotation_mode="vbar")

        # Place Track Labels
        axis.set_title(f"{self.sample.sample}\n{start} - {stop}", loc='left',
                     fontdict={"fontsize": 48})
        axis.text(1.002, 7/8, "SHAPE\nReactivity",
                fontsize=36, transform=axis.transAxes, va='center')
        axis.text(1.002, 5/8, "Shannon\nEntropy",
                fontsize=36, transform=axis.transAxes, va='center')
        axis.text(1.002, 3/8, "Secondary\nStructure",
                transform=axis.transAxes, fontsize=36, va='center')
        axis.text(1.002, 1/8, "Pairing\nProbability",
                transform=axis.transAxes, va='center', fontsize=36)

        # Place region labels
        for i, (lssr_start, lssr_stop) in enumerate(self.lowss_regions):
            middle = (lssr_stop + lssr_start) / 2
            if start < middle < stop:
                axis.text(middle, 550, str(i+1),
                        ha='center', va='center', fontsize=36)

        axis.set_ylim([-305, 915])
        axis.set_xticks(ticks=[x for x in range(500, stop+1, 500) if x > start])
        axis.set_xticks(ticks=[x for x in range(100, stop+1, 100) if x > start],
                      minor=True)
        axis.tick_params(axis='x', which='major', labelsize=36)
        adjust_spines(axis, ['bottom'])

        # set figure size so that 100 axis units == 1 inch
        plot.set_figure_size(height_ax_rel=1/100, width_ax_rel=1/100)
        if colorbar:
            plot.plot_colorbars()
        return plot

    def windowed_median(self, data):
        """Creates an array of same length as data, with windowed median values
        and NaN pads at beginning and end

        Args:
            data (np.array): data for computing windowed medians

        Returns:
            np.array: windowed median array with NaN pads at either end
        """
        pads = [np.nan]*(self.window//2)
        win_median = []
        for i in range(len(data) - self.window+1):
            win_data = data[i:i+self.window]
            win_data = win_data[~np.isnan(win_data)]
            win_median.append(np.median(win_data))
        return np.array(pads + win_median + pads)
