import pandas as pd
from .data import Data
from .ct import CT
from .profile import Profile
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from operator import ge, le, gt, lt


class IJ(Data):

    def __init__(self, datatype, filepath, sequence=None, fasta=None):
        super().__init__(sequence, fasta)
        self.datatype = datatype
        self.path = filepath
        self.default_metric = {'rings': 'Statistic',
                               'pairs': 'Class',
                               'deletions': 'Percentile',
                               'probs': 'Probability'
                               }[self.datatype]
        self._metric = self.default_metric
        read_file = {"probs": self.read_probs,
                     "rings": self.read_rings,
                     "deletions": self.read_deletions,
                     "pairs": self.read_pairs
                     }[datatype]
        read_file(filepath)

    @property
    def min_max(self):
        min_max = {'Percentile': [0.98, 1.0],
                   "Statistic": [-100, 100],
                   'Zij': [-20, 20],
                   "Class": [0, 2],
                   "Metric": [0, 0.001],
                   'Distance': [10, 80],
                   'Probability': [0, 1]
                   }[self.metric]
        return min_max

    @property
    def fill(self):
        fill = {'Class': -1,
                'Statistic': 0.0,
                'Zij': 0.0,
                'Metric': 0.0,
                'Distance': 1000,
                'Percentile': 0.0}[self.metric]
        return fill

    @property
    def cmap(self):
        cmap = {'Class': mp.colors.ListedColormap([[0.3, 0.3, 0.3, 0.2],
                                                   [0.0, 0.0, 0.95, 0.6],
                                                   [0.12, 0.76, 1.0, 0.6]]),
                'Statistic': 'bwr',
                'Zij': 'bwr',
                'Metric': 'YlGnBu',
                'Distance': 'jet',
                'Percentile': 'YlGnBu',
                'Probability': 'inferno_r'
                }[self.metric]
        cmap = plt.get_cmap(cmap)
        cmap = cmap(np.arange(cmap.N))
        cmap[:, -1] = np.full((len(cmap)), 0.6)  # set default alpha to 0.6
        if self.metric == 'Class':
            cmap[0, -1] = 0.2  # alpha of non 1ary and 2ary pairs to 0.2
        if self.metric == 'Distance':
            # set color of max distance and no data distances to gray
            cmap[-1, :] = np.array([80/255., 80/255., 80/255., 0.2])
        cmap = mp.colors.ListedColormap(cmap)
        return cmap

    @property
    def metric(self):
        return self._metric

    @metric.setter
    def metric(self, value):
        if value in self.data.keys():
            self._metric = value
        elif value == "Distance":
            print("Please set 3D distances using IJ.set_3d_distances(PDB)")
        elif value is None:
            self._metric = self.default_metric
        else:
            print(f"{value} is not a valid metric of {self.datatype}")
            self._metric = self.default_metric

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

    def mask_on_sequence(self, compliment_only, nts):
        mask = []
        comp = {'A': 'U', 'U': 'A', 'G': 'CT', 'C': 'G'}
        for _, i, j in self.data[["i", "j"]].itertuples():
            t = True
            for w in range(self.window):
                seq_iw = self.sequence[i+w-1].upper()
                seq_jw = self.sequence[j-w+2].upper()
                if compliment_only:
                    t = t and seq_iw in comp[seq_jw]
                if nts is not None:
                    t = t and seq_iw in nts and seq_jw in nts
            mask.append(t)
        self.update_mask(mask)

    def mask_on_ct(self, ct, cdAbove=None, cdBelow=None,
                   ss_only=False, ds_only=False, paired_only=False):
        if isinstance(ct, list):
            for each in ct:
                self.mask_on_ct(each, cdAbove, cdBelow, ss_only, ds_only,
                                paired_only)
                return
        message = "CT filtering requires a ct object."
        assert isinstance(ct, CT), message
        am = self.get_alignment_map(ct)
        mask = []
        i_j_keep = ["i", "j", "mask"]
        for _, i, j, keep in self.data[i_j_keep].itertuples():
            i = am[i-1]+1
            j = am[j-1]+1
            true_so_far = keep
            if paired_only and true_so_far:
                for w in range(self.window):
                    true_so_far = ct.ct[i+w-1] == j+self.window-w-1
                    if not true_so_far:
                        break
            if ss_only and true_so_far:
                # TODO: which windows are ss vs. ds? At least 2 per window.
                true_so_far = (ct.ct[i-1] == 0) and (ct.ct[j-1] == 0)
            if ds_only and true_so_far:
                true_so_far = (ct.ct[i-1] != 0) and (ct.ct[j-1] != 0)
            if (cdAbove is not None or cdBelow is not None) and true_so_far:
                cd = []
                for iw in range(self.window):
                    for jw in range(self.window):
                        cd.append(ct.contactDistance(i+iw, j+jw))
                cd = min(cd)
                if cdAbove is not None:
                    true_so_far = cd > cdAbove
                if cdBelow is not None:
                    true_so_far = cd < cdBelow
            mask.append(true_so_far)
        self.update_mask(mask)

    def mask_on_profile(self, profile, profAbove=None, profBelow=None):
        alignment_map = self.get_alignment_map(profile)
        norm_prof = profile.data["Norm_profile"]
        mask = []
        for _, i, j in self.data[["i", "j"]].itertuples():
            index_i = alignment_map[i-1]
            index_j = alignment_map[j-1]
            if (index_i != -1) and (index_j != -1):
                prof_i = norm_prof[index_i]
                prof_j = norm_prof[index_j]
                if profAbove is not None:
                    keep_ij = (prof_i >= profAbove) and (prof_j >= profAbove)
                if profBelow is not None:
                    keep_ij = (prof_i <= profBelow) and (prof_j <= profBelow)
            mask.append(keep_ij)
        self.update_mask(mask)

    def mask_nts(self, exclude, isolate):
        mask = np.full(len(self.data), True)
        for index, i, j in self.data[["i", "j"]].itertuples():
            if exclude is not None:
                if i in exclude or j in exclude:
                    mask[index] = False
            if isolate is not None:
                if i not in isolate and j not in isolate:
                    mask[index] = False
        self.update_mask(mask)

    def set_mask_offset(self, fit_to):
        am = self.get_alignment_map(fit_to)
        i = np.array([am[i-1]+1 for i in self.data["i"].values])
        j = np.array([am[j-1]+1 for j in self.data["j"].values])
        mask = np.ones(len(i), dtype=bool)
        for w in range(self.window):
            mask = mask & ((i+w) != 0) & ((j+w) != 0)
        self.data["mask"] = mask
        self.data['i_offset'] = i
        self.data['j_offset'] = j

    def update_mask(self, mask):
        self.data["mask"] = self.data["mask"] & mask

    def filter(self, fit_to,
               ct=None,  # required if any of next are passed
               cdAbove=None, cdBelow=None,
               paired_only=False, ss_only=False, ds_only=False,
               profile=None, profAbove=None, profBelow=None,
               compliments_only=False, nts=None,
               all_pairs=False,
               exclude_nts=None, isolate_nts=None,
               positive_only=False, negative_only=False,
               **kwargs):
        self.set_mask_offset(fit_to)
        if fit_to.datatype == 'pdb':
            mask = []
            for i, j in zip(self.data["i_offset"], self.data["j_offset"]):
                mask.append(i in fit_to.validres and j in fit_to.validres)
            mask = np.array(mask, dtype=bool)
            self.data["mask"] = self.data["mask"] & mask
        if exclude_nts is not None or isolate_nts is not None:
            self.mask_nts(exclude_nts, isolate_nts)
        if profAbove is not None or profBelow is not None:
            message = "Profile filters require a profile object."
            assert isinstance(profile, Profile), message
            self.mask_on_profile(profile, profAbove, profBelow)
        if cdAbove is not None or cdBelow is not None or ss_only or ds_only or paired_only:
            self.mask_on_ct(ct, cdAbove, cdBelow,
                            ss_only, ds_only, paired_only)
        if compliments_only or nts is not None:
            self.mask_on_sequence(compliments_only, nts)
        if not all_pairs and self.datatype == 'pairs':
            self.update_mask(self.data["Class"] != 0)
        if self.datatype == 'probs':
            self.update_mask(self.data["Probability"] >= 0.03)
        if positive_only:
            self.update_mask(self.data["Sign"] == 1)
        if negative_only:
            self.update_mask(self.data["Sign"] == -1)
        for key in kwargs.keys():
            try:
                self.update_mask(self.data[key] > kwargs[key])
            except KeyError:
                print(f"{key} is invalid column of {self.datatype} dataFrame")

    def get_ij_colors(self, min_max=None, cmap=None):
        metric = self.metric
        if min_max is None:
            minimum, maximum = self.min_max
        else:
            minimum, maximum = min_max
        columns = ["i_offset", "j_offset", metric]
        if self.metric in ['Zij', 'Statistic']:
            data = self.data.loc[self.data["mask"], columns+['Sign']].copy()
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
        if cmap is None:
            cmap = self.cmap
        else:
            cmap = plt.get_cmap(cmap)
        i_list, j_list, colors = [], [], []
        if len(data[metric].values) == 0:
            return i_list, j_list, colors
        for _, i, j, datum in data[columns].itertuples():
            for w in range(self.window):
                i_list.append(i + w)
                j_list.append(j + self.window - 1 - w)
                colors.append(cmap(datum))
        return i_list, j_list, colors

    def print_new_file(self, outfile=None, **kwargs):
        data = self.data.copy()
        data = data[data["mask"]]
        data["i"] = data["i_offset"]
        data["j"] = data["j_offset"]
        if "Sign" in data.columns:
            data.rename({"Sign": "+/-"}, inplace=True)
        exclude_columns = ["i_offset", "j_offset",
                           "mask", "Distance", "Percentile"]
        columns = [col for col in data.columns if col not in exclude_columns]
        csv = data.to_csv(columns=columns, sep='\t', index=False,
                          line_terminator='\n')
        if outfile is not None:
            with open(outfile, 'w') as out:
                out.write(self.header)
                out.write(csv)
        else:
            print(self.header, csv)

    def set_3d_distances(self, pdb):
        alignment_map = self.get_alignment_map(pdb)
        distances = []
        for _, i, j in self.data[["i", "j"]].itertuples():
            # for windows calculate average of window*window distances
            # e.g. window==3 results in average of nine distances
            distance = 0
            for iw in range(self.window):
                for jw in range(self.window):
                    io = alignment_map[i+iw-1]+1
                    jo = alignment_map[j+jw-1]+1
                    distance += pdb.get_distance(io, jo)/(self.window**2)
            distances.append(distance)
        self.data["Distance"] = distances
