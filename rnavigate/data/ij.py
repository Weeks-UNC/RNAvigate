import pandas as pd
from .data import Data
from .ct import CT
from .profile import Profile
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from operator import ge, le, gt, lt, eq, ne


class IJ(Data):

    def __init__(self, datatype="ij", dataframe=None, default_metric=None,
                 filepath=None, sep='\t', read_csv_kw={}, window=1,
                 sequence=None, fasta=None,
                 fill={}, cmaps={}, mins_maxes={}):
        super().__init__(sequence=sequence, filepath=fasta)
        self.window = window
        self.datatype = datatype
        self.filepath = filepath
        if dataframe is not None:
            self.data = dataframe
        elif filepath is not None:
            self.read_file(filepath=filepath, sep=sep, read_csv_kw=read_csv_kw)
        self.columns = self.data.columns
        self._fill_values = {'Distance': 1000}
        self._fill_values.update(fill)
        self._cmaps = {'Distance': 'jet'}
        self._cmaps.update(cmaps)
        self._mins_maxes = {'Distance': [10, 80]}
        self._mins_maxes.update(mins_maxes)
        self.default_metric = default_metric
        self.metric = self.default_metric

    def read_file(self, filepath, sep, read_csv_kw):
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)

    @property
    def fill(self):
        return self._fill_values[self.metric]

    @property
    def cmap(self):
        return self._cmap

    @cmap.setter
    def cmap(self, cmap):
        if mp.colors.is_color_like(cmap):
            cmap = mp.colors.ListedColormap([cmap])
        elif isinstance(cmap, list) and all(mp.colors.is_color_like(c) for c in cmap):
            cmap = mp.colors.ListedColormap(cmap)
        cmap = plt.get_cmap(cmap)
        cmap = cmap(np.arange(cmap.N))
        cmap[:, -1] = np.full((len(cmap)), 0.6)  # set default alpha to 0.6
        if self.metric == 'Distance':
            # set color of max distance and no data distances to gray
            cmap[-1, :] = np.array([80/255., 80/255., 80/255., 0.2])
        cmap = self.modify_cmap(cmap)
        cmap = mp.colors.ListedColormap(cmap)
        self._cmap = cmap

    def modify_cmap(self, cmap):
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
        self.cmap = self._cmaps[self._metric]
        self.min_max = self._mins_maxes[self._metric]

    def mask_on_sequence(self, compliment_only, nts, return_mask=False):
        mask = []
        comp = {'A': 'U', 'U': 'AG', 'G': 'CU', 'C': 'G'}
        for _, i, j in self.data[["i", "j"]].itertuples():
            keep = []
            for w in range(self.window):
                i_nt = self.sequence[i+w-1].upper()
                j_nt = self.sequence[j-w+1].upper()
                if compliment_only:
                    keep.append(i_nt in comp[j_nt])
                if nts is not None:
                    keep.append(i_nt in nts and j_nt in nts)
            mask.append(all(keep))
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_ct(self, ct, cdAbove=None, cdBelow=None, ss_only=False,
                   ds_only=False, paired_only=False, return_mask=False):
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
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_profile(self, profile, profAbove=None, profBelow=None):
        alignment_map = self.get_alignment_map(profile)
        norm_prof = profile.data["Norm_profile"]
        mask = []
        for _, i, j in self.data[["i", "j"]].itertuples():
            index_i = alignment_map[i-1]
            index_j = alignment_map[j-1]
            if (index_i != -1) and (index_j != -1):
                prof_i = np.median(norm_prof[index_i:index_i+self.window])
                prof_j = np.median(norm_prof[index_j:index_j+self.window])
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
        if fit_to.datatype == 'pdb':
            mask = []
            for i, j in zip(self.data["i_offset"], self.data["j_offset"]):
                mask.append(i in fit_to.validres and j in fit_to.validres)
            mask = np.array(mask, dtype=bool)
            self.update_mask(mask)

    def update_mask(self, mask):
        self.data["mask"] = self.data["mask"] & mask

    def filter(self, fit_to,
               ct=None,  # required if any of next are passed
               cdAbove=None, cdBelow=None,
               paired_only=False, ss_only=False, ds_only=False,
               profile=None, profAbove=None, profBelow=None,
               compliments_only=False, nts=None,
               exclude_nts=None, isolate_nts=None,
               resolve_conflicts=None,
               **kwargs):
        self.set_mask_offset(fit_to)
        if exclude_nts is not None or isolate_nts is not None:
            self.mask_nts(exclude_nts, isolate_nts)
        if profAbove is not None or profBelow is not None:
            message = "Profile filters require a profile object."
            assert isinstance(profile, Profile), message
            self.mask_on_profile(profile, profAbove, profBelow)
        mask_on_ct = any([cdAbove is not None, cdBelow is not None,
                          ss_only, ds_only, paired_only])
        if mask_on_ct:
            self.mask_on_ct(ct, cdAbove, cdBelow,
                            ss_only, ds_only, paired_only)
        if compliments_only or nts is not None:
            self.mask_on_sequence(compliments_only, nts)
        kwargs = self.data_specific_filter(**kwargs)
        for key in kwargs.keys():
            if key in self.data.keys():
                self.update_mask(self.data[key] > kwargs[key])
            elif "_" in key:
                key2, comparison = key.rsplit("_", 1)
                operators = {"ge": ge, "le": le,
                             "gt": gt, "lt": lt,
                             "eq": eq, "ne": ne}
                if key2 in self.data.keys() and comparison in operators.keys():
                    operator = operators[comparison]
                    self.update_mask(operator(self.data[key2], kwargs[key]))
                else:
                    print(f"{key}={kwargs[key]} is not a valid filter.")
            else:
                print(f"{key}={kwargs[key]} is not a valid filter.")
        if resolve_conflicts is not None:
            self.resolve_conflicts(resolve_conflicts)

    def data_specific_filter(self, **kwargs):
        return kwargs

    def get_ij_colors(self, min_max=None, cmap=None):
        if cmap is None:
            cmap = self.cmap
        else:
            cmap = plt.get_cmap(cmap)
        if min_max is None:
            min_max = self.min_max
        data = self.get_normalized_ij_data(min_max=min_max)
        i_list, j_list, colors = [], [], []
        if len(data) == 0:
            return i_list, j_list, colors
        for _, i, j, datum in data.itertuples():
            for w in range(self.window):
                i_list.append(i + w)
                j_list.append(j + self.window - 1 - w)
                colors.append(cmap(datum))
        return i_list, j_list, colors

    def get_normalized_ij_data(self, min_max):
        metric = self.metric
        columns = ["i_offset", "j_offset", metric]
        data = self.data.loc[self.data["mask"], columns].copy()
        ascending = metric not in ['Distance']  # high distance is bad
        data.sort_values(by=metric, ascending=ascending, inplace=True)
        norm = plt.Normalize(min_max[0], min_max[1], clip=True)
        data[metric] = norm(data[metric])
        return data

    def print_new_file(self, outfile=None):
        data = self.data.copy()
        data = data[data["mask"]]
        data["i"] = data["i_offset"]
        data["j"] = data["j_offset"]
        if "Sign" in data.columns:
            data.rename({"Sign": "+/-"}, inplace=True)
        csv = data.to_csv(columns=self.columns, sep='\t', index=False,
                          line_terminator='\n')
        if outfile is not None:
            with open(outfile, 'w') as out:
                out.write(self.header)
                out.write(csv)
        else:
            print(self.header, csv)

    def set_3d_distances(self, pdb, atom):
        alignment_map = self.get_alignment_map(pdb)
        distance_matrix = pdb.get_distance_matrix(atom=atom)
        i = self.data["i"].values
        j = self.data["j"].values
        distances = np.zeros(len(i))
        pairs = self.window**2
        for iw in range(self.window):
            for jw in range(self.window):
                io = alignment_map[i+iw-1]
                jo = alignment_map[j+jw-1]
                distances += distance_matrix[io, jo]/pairs
        self.data["Distance"] = distances

    def resolve_conflicts(self, metric=None):
        """Resolves conflicting windows using the Maximal Weighted Independent
        Set. The weights are taken from the metric value. The graph is first
        broken into components to speed up the identification of the MWIS. Then
        the mask is updated to only include the MWIS.
        """
        def get_components(graph):
            """Using DFS, returns a list of component subgraphs of graph."""
            def dfs(key, value, component):
                """Depth first seach to identify a component."""
                visited[key] = True
                component[key] = value
                for v in value:
                    if not visited[v]:
                        component = dfs(v, graph[v], component)
                return dict(component)

            components = []
            visited = {key: False for key in graph.keys()}
            for key, value in graph.items():
                if not visited[key]:
                    components.append(dfs(key, value, {}))
            return components

        def graphSets(graph, node_weights):
            """Finds the Maximal Weighted Independent Set of a graph."""
            # Base Case - Given Graph has no nodes
            if(len(graph) == 0):
                return []
            # Select a vertex from the graph
            current_vertex = list(graph.keys())[0]
            # Create a copy of the graph and remove the current vertex
            new_graph = dict(graph)
            del new_graph[current_vertex]
            # Case 1 - excluding current Vertex
            case_1 = graphSets(new_graph, node_weights)
            # Case 2 - including current Vertex, and excluding neighbors
            for neighbor in graph[current_vertex]:
                if(neighbor in new_graph):
                    del new_graph[neighbor]
            case_2 = [current_vertex] + graphSets(new_graph, node_weights)
            # Our final result is the one which has highest weight, return it
            weight_1 = sum(node_weights[node] for node in case_1)
            weight_2 = sum(node_weights[node] for node in case_2)
            if(weight_1 > weight_2):
                return case_1
            return case_2

        # Building the graph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if metric is None:
            metric = self.metric
        # a list of indices relating back to the original mask
        indices = np.where(self.data["mask"])[0]
        # the correlation data to be considered
        data = self.data.loc[indices].copy()
        window = self.window
        # node_weight[i from indices] = metric value for that correlation
        node_weights = {}
        # graph[i from indices] = list of indices of conflicting correlations
        graph = {}
        # for each considered correlation, record conflicts and node weights
        for v1, i, j, g in data[["i", "j", metric]].itertuples():
            if v1 not in graph:
                graph[v1] = []
            node_weights[v1] = g  # recording the node weight
            # conflicts = indices of all overlapping, but not parallel corrs
            conflicts = np.where(((abs(data["i"] - i) < window) |
                                  (abs(data["j"] - i) < window) |
                                  (abs(data["i"] - j) < window) |
                                  (abs(data["j"] - j) < window)) &
                                 ((data["j"] - j) != (i - data["i"])))[0]
            # adding conflicts to graph dictionary
            for conflict in conflicts:
                v2 = indices[conflict]
                if v2 > v1:  # prevents self-loop and parallel edges
                    if v2 not in graph:
                        graph[v2] = []
                    graph[v1].append(v2)
                    graph[v2].append(v1)
        # Finding the maximum set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Split the graph into components. Get the max set for each component.
        # The graph's max set is the same as the union of these max sets.
        max_set = []
        for component in get_components(graph):
            max_set += graphSets(component, node_weights)
        # set these indices to true and update the original mask
        new_mask = np.zeros(len(self.data), dtype=bool)
        new_mask[max_set] = True
        self.update_mask(new_mask)


class SHAPEJuMP(IJ):
    def __init__(self, filepath, datatype="shapejump", sequence=None):
        default_metric = 'Percentile'
        fill = {'Metric': 0.0, 'Percentile': 0.0}
        cmaps = {'Metric': 'YlGnBu', 'Percentile': 'YlGnBu'}
        mins_maxes = {'Percentile': [0.98, 1.0], "Metric": [0, 0.001]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, sequence=sequence,
                         fill=fill, cmaps=cmaps, mins_maxes=mins_maxes)

    def read_file(self, filepath, **kwargs):
        column_names = ['Gene', 'i', 'j', 'Metric']
        data = pd.read_csv(filepath, sep='\t', names=column_names, header=0)
        data["Percentile"] = data['Metric'].rank(method='max', pct=True)
        self.data = data
        self.window = 1


class RINGMaP(IJ):
    def __init__(self, filepath, datatype="ringmap", sequence=None):
        default_metric = 'Statistic'
        fill = {'Statistic': 0.0, 'Zij': 0.0}
        cmaps = {'Statistic': 'bwr', 'Zij': 'bwr'}
        mins_maxes = {"Statistic": [-100, 100], 'Zij': [-8, 8]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, sequence=sequence,
                         fill=fill, cmaps=cmaps, mins_maxes=mins_maxes)

    def read_file(self, filepath, **kwargs):
        with open(filepath, 'r') as file:
            self.header = file.readline()
        split_header = self.header.split('\t')
        window = split_header[1].split('=')[1]
        self.window = int(window)
        self.data = pd.read_csv(filepath, sep='\t', header=1)
        self.data.rename(columns={"+/-": "Sign"}, inplace=True)

    def data_specific_filter(self, positive_only, negative_only, **kwargs):
        if positive_only:
            self.update_mask(self.data["Sign"] == 1)
        if negative_only:
            self.update_mask(self.data["Sign"] == -1)
        return kwargs

    def get_normalized_ij_data(self, min_max):
        if self.metric == 'Distance':
            return super().get_normalized_ij_data(min_max)
        if min_max is None:
            min_max = self.min_max
        metric = self.metric
        columns = ["i_offset", "j_offset", metric]
        data = self.data.loc[self.data["mask"], columns+['Sign']].copy()
        data[metric] = data[metric]*data["Sign"]
        data.sort_values(by=metric, ascending=True, inplace=True)
        norm = plt.Normalize(min_max[0], min_max[1], clip=True)
        data[metric] = norm(data[metric])
        return data[columns]


class PAIRMaP(IJ):
    def __init__(self, filepath, datatype="pairmap", sequence=None):
        default_metric = 'Class'
        fill = {'Class': -1, 'Statistic': 0.0, 'Zij': 0.0}
        cmaps = {'Class': mp.colors.ListedColormap([[0.3, 0.3, 0.3, 0.2],
                                                    [0.0, 0.0, 0.95, 0.6],
                                                    [0.12, 0.76, 1.0, 0.6]]),
                 'Statistic': 'bwr', 'Zij': 'bwr'}
        mins_maxes = {"Statistic": [-100, 100],
                      'Zij': [-8, 8], "Class": [0, 2]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, sequence=sequence,
                         fill=fill, cmaps=cmaps, mins_maxes=mins_maxes)

    def read_file(self, filepath, **kwargs):
        with open(filepath, 'r') as file:
            self.header = file.readline()
        self.window = int(self.header.split('\t')[1].split('=')[1])
        self.data = pd.read_csv(filepath, sep='\t', header=1)
        self.data.rename(columns={"Sig.": "Statistic"}, inplace=True)

    def modify_cmap(self, cmap):
        if self.metric == 'Class':
            cmap[0, -1] = 0.2  # alpha of non 1ary and 2ary pairs to 0.2
        return cmap

    def data_specific_filter(self, all_pairs=False, **kwargs):
        if not all_pairs and self.datatype == 'pairs':
            self.update_mask(self.data["Class"] != 0)
        return kwargs

    def get_normalized_ij_data(self, min_max):
        if self.metric != 'Class':
            return super().get_normalized_ij_data(min_max)
        if min_max is None:
            min_max = self.min_max
        metric = self.metric
        columns = ["i_offset", "j_offset", metric]
        data = self.data.loc[self.data["mask"], columns].copy()
        # Class data have a weird order
        data = pd.concat([data[data[metric] == 0],
                          data[data[metric] == 2],
                          data[data[metric] == 1]])
        return data


class PairProb(IJ):
    def __init__(self, filepath, datatype="pairprob", sequence=None):
        default_metric = 'Probability'
        fill = {'Probability': 0}
        cmaps = {'Probability': 'inferno_r'}
        mins_maxes = {'Probability': [0, 1]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, sequence=sequence,
                         fill=fill, cmaps=cmaps, mins_maxes=mins_maxes)

    def read_file(self, filepath, **kwargs):
        with open(filepath, 'r') as file:
            self.header = file.readline()
        lengths_match = int(self.header.strip()) == self.length
        assert lengths_match, "DP and sequence are different lengths"
        self.window = 1
        data = pd.read_table(filepath, header=1, names=['i', 'j', 'log10p'])
        data["Probability"] = 10 ** (-data["log10p"])
        self.data = data

    def set_entropy(self, printOut=False, toFile=None):
        assert self.datatype == "probs", "IJ must be a pairing probability .dp"
        self.data.eval('nlogn = log10p * 10 ** ( - log10p )', inplace=True)
        entropy = np.zeros(self.length)
        for i in range(self.length):
            mask = (self.data["i"] == i+1) | (self.data["j"] == i+1)
            entropy[i] = self.data.loc[mask, "nlogn"].sum()
        # catch rounding errors:
        entropy[np.where(entropy < 0)] = 0
        if printOut:
            print(*[f"{i+1} {s}" for i, s in enumerate(entropy)], sep="\n")
        if toFile:
            with open(toFile) as outf:
                for i, s in enumerate(entropy):
                    outf.write(f"{i+1}\t{s}\n")
        self.entropy = entropy

    def data_specific_filter(self, **kwargs):
        self.update_mask(self.data["Probability"] >= 0.03)
        return kwargs
