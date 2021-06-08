import Bio.SeqIO
import pandas as pd
import plotmapper as MaP


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


class IJ():

    def __init__(self, datatype, filepath, sequence=None, fasta=None):
        self.datatype = datatype
        self.default_metric = {'rings': 'Statistic',
                               'pairs': 'Class',
                               'deletions': 'Percentile',
                               'probs': 'Probability'
                               }[self.datatype]
        self.path = filepath
        sequence_passed = sequence is not None or fasta is not None
        assert sequence_passed, "ij_data object requires a sequence"
        if sequence is not None:
            self.sequence = sequence
            self.length = len(sequence)
        if fasta is not None:
            fasta = list(Bio.SeqIO.parse(open(fasta), 'fasta'))
            self.sequence = str(fasta[0].seq).upper().replace("T", "U")
            self.length = len(self.sequence)
            self.gene = fasta[0].id
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
        self.data.rename(columns={"+/-": "Sign"})

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
        clip, pad = MaP.get_clip_pad(self.sequence, fit_to.sequence)
        start = clip[0]
        end = clip[1]
        i = self.data["i"].values
        j = self.data["j"].values
        iinrange = (start < i) & (i < end)
        jinrange = (start < j) & (j < end)
        self.data['mask'] = iinrange & jinrange
        offset = pad[0]-start
        self.data['i_offset'] = i+offset
        self.data['j_offset'] = j+offset

    def update_mask(self, mask):
        self.data["mask"] = self.data["mask"] & mask

    def filter_ij_data(self, fit_to, profile=None, ct=None, cdAbove=None,
                       cdBelow=None, ss_only=False, ds_only=False,
                       profAbove=None, profBelow=None, all_pairs=False,
                       **kwargs):
        self.set_mask_offset(fit_to)
        if fit_to.type == 'pdb':
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
            assert isinstance(ct, "CT"), message
            self.mask_on_ct(ct, cdAbove, cdBelow, ss_only, ds_only)
        if not all_pairs and ij_data == 'pairs':
            self.update_mask(self.data["Class"] != 0)
        if ij_data == 'probs':
            self.update_mask(self.data["Probability"] >= 0.03)
        for key in kwargs.keys():
            try:
                self.update_mask(self.data[key] > kwargs[key])
            except KeyError:
                print(f"{key} is not a valid column of {self.type} dataFrame")

    def get_ij_colors(self, metric=None, min_max=None, cmap=None):
        if metric is None:
            metric = get_default_metric(self.type)
        if min_max is None:
            minimum, maximum = get_default_min_max(metric)
        else:
            minimum, maximum = min_max
        columns = ["i_offset", "j_offset", metric]
        if self.type == 'rings':
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

    def print_new_ij_file(self, outfile=None, **kwargs):
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
