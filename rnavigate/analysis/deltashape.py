#!/usr/bin/env python
#  deltaSHAPE software for detecting meaningful changes in SHAPE reactivity
#
#  - Requires two .map files as input (see README for details)
#  - See the README for required modules, installation, and execution help.
#  - Version 1.0
#  - Copyright Matthew J. Smola 2015
#    Refactored for RNAvigate by Patrick Irving 2023

###########################################################################
# GPL statement:                                                          #
#                                                                         #
# This program is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
###########################################################################

import warnings
import matplotlib.pyplot as plt
import numpy as np


class DeltaSHAPE():
    def __init__(self, sample1, sample2, pad=1, mask_5=0, mask_3=0, z_coeff=1.96,
                 z_thresh=0, ss_thresh=1, find_site=(2, 3)):
        """Detects meaningful differences in chemical probing reactivity according
        to the deltaSHAPE algorithm (doi:10.1021/acs.biochem.5b00977)

        Args:
            sample1 (rnavigate.Sample): First sample to compare
            sample2 (rnavigate.Sample): Second sample to compare
            pad (int, optional): Windows range = (nt-pad:nt+pad). Defaults to 1.
            mask_5 (int, optional): size of 5` primer site. Defaults to 0.
            mask_3 (int, optional): size of 3` primer site. Defaults to 0.
            z_coeff (float, optional): Ajust the Z-factor stringency by changing the equation coefficient. See the README for details. Defaults to 1.96.
            z_thresh (int, optional): Adjust the Z-factor stringency by changing the cutoff threshold. See the README for details. Defaults to 0.
            ss_thresh (int, optional): Set the cutoff threshold of standard score filtering. Defaults to 1.
            find_site (tuple, optional): Comma-separated pair of numbers indicating the window pad size and number of required hits when finding binding sites. Defaults to (2, 3).
        """

        # STEP ONE
        # extract relevant information from sample 1 and 2 profiles
        # sequence, normalized profile, normalized standard error
        columns = ["Nucleotide", "Sequence", "Norm_profile", "Norm_stderr"]
        sample1 = sample1.get_data_list("profile").data[columns].copy()
        sample2 = sample2.get_data_list("profile").data[columns].copy()
        assert all(sample1["Sequence"] == sample2["Sequence"]
                   ), "Sequences don't match"
        self.data = sample1.merge(sample2, how="left",
                                  on=["Nucleotide", "Sequence"],
                                  suffixes=("_1", "_2"))

        # STEP TWO
        # smooth data and errors
        s_data1, s_err1 = smooth(self.data["Norm_profile_1"],
                                 self.data["Norm_stderr_1"], pad)
        self.data["Smooth_profile_1"] = s_data1
        self.data["Smooth_stderr_1"] = s_err1

        s_data2, s_err2 = smooth(self.data["Norm_profile_2"],
                                 self.data["Norm_stderr_2"], pad)
        self.data["Smooth_profile_2"] = s_data2
        self.data["Smooth_stderr_2"] = s_err2

        # STEP THREE
        # subtract raw and smoothed data
        self.data.eval("diff = Norm_profile_1 - Norm_profile_2",
                       inplace=True)
        self.data.eval("Smooth_diff = Smooth_profile_1 - Smooth_profile_2",
                       inplace=True)
        self.data.eval("Positive = (Smooth_diff > 0)", inplace=True)

        # STEP FOUR
        # calculate Z-factors from smoothed data and smoothed errs
        self.data["z_factors"] = z_factor(
            s_data1, s_data2, s_err1, s_err2, z_coeff)

        # STEP FIVE
        # calculate Z-scores from difference of smoothed data1 and smoothed data2
        self.data["z_scores"] = calc_zScores(self.data["Smooth_diff"])

        # STEP SIX
        # identify (2x+1)-nt windows where y+ nts are sig. diff.
        # find_site=(x, y)
        site_pad, site_min = find_site
        self.data["Significant"] = np.full(len(self.data), False)
        for i in range(site_pad, len(self.data["diff"])-site_pad):
            win = range(i-site_pad, i+site_pad+1)
            count = 0
            maybes = []
            for j in win:
                if ((self.data.loc[j, "z_factors"] > z_thresh) and
                        (np.abs(self.data.loc[j, "z_scores"]) >= ss_thresh)):
                    count += 1
                    maybes.append(j)
            if count >= site_min:
                for k in maybes:
                    self.data.loc[k, "Significant"] = True

    def get_figsize(self):
        left_inches = 0.9
        right_inches = 0.4
        ax_width = len(self.data) * 0.1
        fig_height = 6
        fig_width = max(7, ax_width + left_inches + right_inches)
        return (fig_width, fig_height)

    def plot(self, xlims=None, ylims=None, front=0, back=0):
        """Creates and displays a plot of deltaSHAPE profile and highlights
        sites called.

        Args:
            xlims (tuple of 2 floats, optional): sets x-ax bounds.
                Defaults to None which uses Nucleotide range.
            ylims (tuple of 2 floats, optional): sets y-ax bounds.
                Defaults to None which uses smoothed difference range +/- 0.25

        Returns:
            matplotlib figure and ax objects
        """
        fig, ax = plt.subplots(1, figsize=self.get_figsize())

        ax.plot(self.data["Nucleotide"], self.data["Smooth_diff"],
                drawstyle='steps-mid', color='black')
        plt.axhline(0, color='black')

        # mask primer-binding regions
        plt.axvspan(0, front+0.5, color="grey", alpha=0.25)
        plt.axvspan(len(self.data)-back+0.5, len(self.data) + 0.5,
                    color="grey", alpha=0.25)

        # color deltaSHAPE sites
        ax.bar(self.data["Nucleotide"],
               self.data.eval("Smooth_diff * Positive * Significant"),
               width=1, ec=None, fc='#3EB452', lw=0)
        ax.bar(self.data["Nucleotide"],
               self.data.eval("Smooth_diff * (~Positive) * Significant"),
               width=1, ec=None, fc='#7F3B95', lw=0)

        # set axes limits
        # default xlim is min=1, max=length of RNA
        if xlims is None:
            xlims = (1, len(self.data))

        # default ylim is Smooth_diff +/- 0.25
        if ylims is None:
            ylims = (np.nanmin(self.data["Smooth_diff"])-0.25,
                     np.nanmax(self.data["Smooth_diff"])+0.25)
        ax.set(
            xlim=xlims,
            ylim=ylims,
            xlabel="Nucleotide",
            ylabel=r'$\Delta$SHAPE')
        ax.tick_params(
            which='both',
            direction='out',
            top='off',
            right='off')

        # turn off UserWarnings temporarily so that plt.tight_layout() doesn't print a warning to the screen.
        warnings.simplefilter("ignore", UserWarning)
        # set the plot layout
        plt.tight_layout()
        # turn warnings back on in case something terrible happens.
        warnings.resetwarnings()

        return fig, ax

# "mask_5 int: Specify the number of nucleotides at the 5' end to ignore. Default: 0"
# "mask_3 int: Specify the number of nucleotides at the 3' end to ignore. Default: 0"
# 'pad int: Indicate the smoothing window size. Window = 2*pad+1. To turn off smoothing, set PAD = 0. Default: 1'
# 'z_coeff float: Ajust the Z-factor stringency by changing the equation coefficient. See the README for details. Default: 1.96'
# 'z_thresh float: Adjust the Z-factor stringency by changing the cutoff threshold. See the README for details. Default: 0'
# 'ss_thres float: Set the cutoff threshold of standard score filtering. Default: 1.0'
# 'find_site tuple: Comma-separated pair of numbers indicating the window pad size and number of required hits when finding binding sites. Default settings look for 3+ nucleotides within a 5-nucleotide window. See the README for details. Default: 2,3'

# ('-o', '--out', type=str, default="differences.txt", help='Name and location of output file to be written. Default: ./differences.txt')
# ('--magrank', action='store_true', help='Sort output file by decreasing deltaSHAPE magnitude. Default: OFF')
# ('--all', action='store_true', help='Output data for all nucleotides. Insignificant changes are listed as zero. Default: OFF')
# ('--pdf', action='store_true', help='Save plot as PDF. If output file is given, PDF will have same prefix. Default: OFF')
# ('--noshow', action='store_true', help='Generate the plot but do not show it. Typically used with --pdf. Default: display plot')
# ('--noplot', action='store_true', help='Skip plotting completely. Default: OFF')
# ('--dots', action='store_true', help='Plot markers indicating nucleotides that pass Z-factor and standard score filtering. This can get unweildy for large RNAs (>1000). Standard score (open) dots are plotted above Z-factor (filled) dots. Default: OFF')
# ('--Zdots', action='store_true', help='Plot markers indicating only nucleotides that pass Z-factor filtering. Default: OFF')
# ('--SSdots', action='store_true', help='Plot markers indicating only nucleotides that pass standard score filtering. Default: OFF')
# ('--colorfill', action='store_true', help='Highlight deltaSHAPE sites with coloration beneath the plot line for "prettier" figures. Default: OFF')
# ('--ymin', type=float, default=-999, help='Set plot y-ax minimum. Default: Determined automatically')
# ('--ymax', type=float, default=-999, help='Set plot y-ax maximum. Default: Determined automatically')
# ('--xmin', type=float, default=-999, help='Set plot x-ax minimum. Default: Determined automatically')
# ('--xmax', type=float, default=-999, help='Set plot x-ax maximum. Default: Determined automatically')


def smooth(data, err, pad):
    new_data, new_err = [], []
    # eventually we want to exclude no-data nucleotides.
    # create a list ("mask") to store which positions to ignore later.
    mask = []
    for i in range(len(data)):
        if data[i] == -999 or np.isnan(data[i]) == True:
            mask.append(i)
    # you can't center a window at the first nucleotide so mask until a full centered window can be placed
    for i in range(pad):
        new_data.append(np.nan)
        new_err.append(np.nan)
    # proceed by windows to smooth the data
    for i in range(pad, len(data)-pad):
        # use numpy masked array to calculate average without including no-data (nan) nucleotides.
        new_data.append(np.mean(np.ma.MaskedArray(
            [j for j in data[i-pad:i+pad+1]], np.isnan([j for j in data[i-pad:i+pad+1]]))))

        # use stats.nanmean to calculate average without including no-data (nan) nucleotides. This causes long_scalars runtime warnings.
        #new_data.append(stats.nanmean([j for j in data[i-pad:i+pad+1] if np.isnan(j) != True]))
        errs = np.array(err[i-pad:i+pad+1])
        squerrs = np.power([j for j in errs if np.isnan(j) != True], 2)
        total = np.sum(squerrs)
        sqrt = np.sqrt(total)
        new_err.append(sqrt/len(data[i-pad:i+pad+1]))
    for i in range(pad):
        new_data.append(np.nan)
        new_err.append(np.nan)
    for i in mask:
        new_data[i] = np.nan
        new_err[i] = np.nan
    return np.array(new_data), np.array(new_err)


def z_factor(data1, data2, err1, err2, factor=1.96):
    z_factors = []
    for i in range(len(data1)):
        if data1[i] == 'nan' or data2[i] == 'nan':
            z_factors.append(float('nan'))
        else:
            # 1.645 = 90% confidence interval
            top = factor * (err2[i] + err1[i])
            bot = abs(data2[i] - data1[i])
            if bot == 0:
                z_factors.append(float('nan'))
            else:
                z = (1 - (top / bot))
                z_factors.append(z)
    return z_factors


def calc_zScores(diffs):
    mean = np.nanmean(diffs)
    sigma = np.nanstd(diffs)
    # calc Z-score
    z_scores = (diffs - mean) / sigma
    return np.array(z_scores)
