import pandas as pd
import plotmapper as MaP
from data.data import Data
from data.ct import CT
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np


def get_default_min_max(metric):
    min_max = {'Percentile': [0.98, 1.0],
               "Statistic": [-100, 100],
               'Zij': [-50, 50],
               "Class": [0, 2],
               "Metric": [0, 0.001],
               'Distance': [10, 80],
               'Probability': [0, 1]
               }[metric]
    return min_max


def get_default_fill(metric):
    fill = {'Class': -1,
            'Statistic': 0.0,
            'Zij': 0.0,
            'Metric': 0.0,
            'Distance': 1000,
            'Percentile': 0.0}[metric]
    return fill


def get_default_cmap(metric):
    cmap = {'Class': mp.colors.ListedColormap([[0.3, 0.3, 0.3, 0.2],
                                               [0.0, 0.0, 0.95, 0.6],
                                               [0.12, 0.76, 1.0, 0.6]]),
            'Statistic': 'bwr',
            'Zij': 'bwr',
            'Metric': 'YlGnBu',
            'Distance': 'jet',
            'Percentile': 'YlGnBu',
            'Probability': 'rainbow_r'
            }[metric]
    cmap = plt.get_cmap(cmap)
    cmap = cmap(np.arange(cmap.N))
    cmap[:, -1] = np.full((len(cmap)), 0.6)  # set default alpha to 0.6
    if metric == 'Class':
        cmap[0, -1] = 0.2  # alpha of non primary and secondary pairs to 0.2
    if metric == 'Distance':
        # set color of max distance and no data distances to gray
        cmap[-1, :] = np.array([80/255., 80/255., 80/255., 0.2])
    cmap = mp.colors.ListedColormap(cmap)
    return cmap


def get_default_metric(ij_data):
    metric = {'rings': 'Statistic',
              'pairs': 'Class',
              'deletions': 'Percentile',
              'probs': 'Probability'
              }[ij_data]
    return metric


class IJ(Data):

    def __init__(self, datatype, filepath, sequence=None, fasta=None):
        super().__init__(sequence, fasta)
        self.datatype = datatype
        self.default_metric = {'rings': 'Statistic',
                               'pairs': 'Class',
                               'deletions': 'Percentile',
                               'probs': 'Probability'
                               }[self.datatype]
        self.path = filepath
        read_file = {"probs": self.read_probs,
                     "rings": self.read_rings,
                     "deletions": self.read_deletions,
                     "pairs": self.read_pairs
                     }[datatype]
        read_file(filepath)

    def read_probs(self, probs):
        with open(probs, 'r') as file:
            self.header = file.readline()
        self.window = 1
        data = pd.read_csv(probs, sep='\t', header=1)
        data["Probability"] = 10 ** (-data["-log10(Probability)"])
        self.data = data

    def read_rings(self, rings):
        with open(rings, 'r') as file:
            self.header = file.readline()
        split_header = self.header.split('\t')
        window = split_header[1].split('=')[1]
        self.window = int(window)
        self.data = pd.read_csv(rings, sep='\t', header=1)
        self.data.rename(columns={"+/-": "Sign"}, inplace=True)

    def read_deletions(self, deletions):
        column_names = ['Gene', 'i', 'j', 'Metric']
        data = pd.read_csv(deletions, sep='\t', names=column_names, header=0)
        data["Percentile"] = data['Metric'].rank(method='max', pct=True)
        self.data = data
        self.window = 1

    def read_pairs(self, pairs):
        with open(pairs, 'r') as file:
            self.header = file.readline()
        split_header = self.header.split('\t')
        window = split_header[1].split('=')[1]
        self.window = int(window)
        self.data = pd.read_csv(pairs, sep='\t', header=1)

    def mask_on_ct(self, ct, cdAbove=None, cdBelow=None,
                   ss_only=False, ds_only=False, paired_only=False):
        mask = []
        i_j_keep = ["i_offset", "j_offset", "mask"]
        for _, i, j, keep in self.data[i_j_keep].itertuples():
            true_so_far = keep
            if paired_only and true_so_far:
                true_so_far = ct.ct[i-1] == j
            if ss_only and true_so_far:
                true_so_far = (ct.ct[i-1] == 0) and (ct.ct[j-1] == 0)
            if ds_only and true_so_far:
                true_so_far = (ct.ct[i-1] != 0) and (ct.ct[j-1] != 0)
            if cdAbove is not None and true_so_far:
                true_so_far = ct.contactDistance(i, j) > cdAbove
            if cdBelow is not None and true_so_far:
                true_so_far = ct.contactDistance(i, j) < cdBelow
            mask.append(true_so_far)
        self.update_mask(mask)

    def mask_on_profile(self, profile, profAbove=None, profBelow=None):
        clip, pad = MaP.get_clip_pad(self.sequence, profile.sequence)
        norm_prof = self.profile["Norm_profile"]
        mask = []
        for _, i, j in self.data[["i", "j"]].itertuples():
            keep_ij = (clip[0] < i < clip[1]) and (clip[0] < j < clip[1])
            prof_i = norm_prof[i-1+pad[0]]
            prof_j = norm_prof[j-1+pad[0]]
            if profAbove is not None and keep_ij:
                keep_ij = (prof_i >= profAbove) and (prof_j >= profAbove)
            if profBelow is not None and keep_ij:
                keep_ij = (prof_i <= profBelow) and (prof_j <= profBelow)
            mask.append(keep_ij)
        self.update_mask(mask)

    def set_mask_offset(self, fit_to):
        alignment_map = self.get_alignment_map(fit_to)
        i = [alignment_map[i-1] for i in self.data["i"].values]
        j = [alignment_map[i-1] for i in self.data["j"].values]
        self.data['mask'] = (i != 0) & (j != 0)
        self.data['i_offset'] = i
        self.data['j_offset'] = j

    def update_mask(self, mask):
        self.data["mask"] = self.data["mask"] & mask

    def filter(self, fit_to, profile=None, ct=None, cdAbove=None,
               cdBelow=None, ss_only=False, ds_only=False,
               profAbove=None, profBelow=None, all_pairs=False,
               **kwargs):
        self.set_mask_offset(fit_to)
        if fit_to.datatype == 'pdb':
            mask = []
            for i, j in zip(self.data["i_offset"], self.data["j_offset"]):
                mask.append(i in fit_to.validres and j in fit_to.validres)
            mask = np.array(mask, dtype=bool)
            self.data["mask"] = self.data["mask"] & mask
        if profAbove is not None or profBelow is not None:
            message = "Profile filters require a profile object."
            assert isinstance(profile, "profile"), message
            self.mask_on_profile(profile, profAbove, profBelow)
        if cdAbove is not None or cdBelow is not None or ss_only or ds_only:
            message = "CT filtering requires a ct object."
            assert isinstance(ct, CT), message
            self.mask_on_ct(ct, cdAbove, cdBelow, ss_only, ds_only)
        if not all_pairs and self.datatype == 'pairs':
            self.update_mask(self.data["Class"] != 0)
        if self.datatype == 'probs':
            self.update_mask(self.data["Probability"] >= 0.03)
        for key in kwargs.keys():
            try:
                self.update_mask(self.data[key] > kwargs[key])
            except KeyError:
                print(f"{key} is not a valid column of {self.datatype} dataFrame")

    def get_ij_colors(self, metric=None, min_max=None, cmap=None):
        if metric is None:
            metric = self.default_metric
        if min_max is None:
            minimum, maximum = get_default_min_max(metric)
        else:
            minimum, maximum = min_max
        columns = ["i_offset", "j_offset", metric]
        if self.datatype == 'rings':
            columns.append("Sign")
            data = self.data.loc[self.data["mask"], columns].copy()
            data[metric] = data[metric]*data["Sign"]
        else:
            data = self.data.loc[self.data["mask"], columns].copy()
        if metric != 'Class':  # Normalize to between 0 and 1
            data.loc[data[metric] < minimum, metric] = minimum
            data.loc[data[metric] > maximum, metric] = maximum
            ascending = metric not in ['Distance']  # high distance is bad
            data = data.sort_values(by=metric, ascending=ascending)
            data[metric] = (data[metric]-minimum) / (maximum-minimum)
        else:  # Class data have a weird order
            data = pd.concat([data[data[metric] == 0],
                              data[data[metric] == 2],
                              data[data[metric] == 1]])
        i = data["i_offset"].values
        j = data["j_offset"].values
        if cmap is None:
            cmap = get_default_cmap(metric)
        else:
            cmap = plt.get_cmap(cmap)
        colors = cmap(data[metric].values)
        return i, j, colors

    def print_new_file(self, outfile=None, **kwargs):
        self.filter_ij_data(**kwargs)
        columns = [c if c != "Sign" else "+/-" for c in self.data.columns]
        exclude_columns = ["i_offset", "j_offset",
                           "mask", "Distance", "Percentile"]
        for col in exclude_columns:
            if col in columns:
                columns.remove(col)
        csv = self.data.to_csv(columns=columns, sep='\t', index=False,
                               line_terminator='\n')
        if outfile is not None:
            with open(outfile, 'w') as out:
                out.write(self.header)
                out.write(csv)
        else:
            print(self.header, csv)

    # def set_3d_distances(self, ij_data):
    #     self.filter_ij_data(ij_data, "pdb")
    #     data = self.ij_data[ij_data].copy()
    #     if "Distance" not in data:
    #         distances = []
    #         for _, i, j in data[["i_offset", "j_offset"]].itertuples():
    #             distances.append(self.get_3d_distance(i, j))
    #         data["Distance"] = distances
    #     self.ij_data[ij_data] = data
