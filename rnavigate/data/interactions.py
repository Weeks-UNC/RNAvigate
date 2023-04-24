import pandas as pd
from .data import Data
from .ct import CT
from .profile import Profile
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from operator import ge, le, gt, lt, eq, ne
import seaborn as sns


class Interactions(Data):
    def __init__(self, datatype="interactions", dataframe=None,
                 default_metric=None,
                 filepath=None, sep='\t', read_csv_kw={}, window=1,
                 sequence=None, fasta=None,
                 fill={}, cmaps={}, mins_maxes={}):
        """Given a dataframe or a data file, construct the interactions object

        Args:
            datatype (str, optional): becomes self.datatype, indicates the
                datatype. Defaults to "interactions".
            dataframe (pandas DataFrame, optional): a dataframe containing
                interactions data. Must have at least "i" and "j" columns
                indicating the 5' and 3' ends of the interactions. Defaults to
                None.
            default_metric (str, optional): column name to use as the default
                metric. Defaults to None.
            filepath (str, optional): path to a file containing interactions
                data. Defaults to None.
            sep (str, optional): passed to pandas read_csv. Defaults to '\t'.
            read_csv_kw (dict, optional): other options for read_csv. Defaults to {}.
            window (int, optional): 5' and 3' interactions windows. Defaults to 1.
            sequence (str, optional): sequence string. Defaults to None.
            fasta (str, optional): path to fasta file. Defaults to None.
            fill (dict, optional): dictionary specifying a fill value (values)
                to use with a metric (keys). Defaults to {}.
            cmaps (dict, optional): specifies cmaps (values) to use with
                metrics (keys). Defaults to {}.
            mins_maxes (dict, optional): specifies minimum and maximum values
                as a list of floats (values) to use with given metric (keys).
                Defaults to {}.
        """
        super().__init__(sequence=sequence, filepath=fasta)
        self.window = window
        self.datatype = datatype
        self.filepath = filepath
        if dataframe is not None:
            self.data = dataframe
        elif filepath is not None:
            self.read_file(filepath=filepath, sep=sep, read_csv_kw=read_csv_kw)
        self.columns = self.data.columns
        self._fill_values = {'Distance': np.nan}
        self._fill_values.update(fill)
        self._cmaps = {'Distance': 'jet'}
        self._cmaps.update(cmaps)
        self._mins_maxes = {'Distance': [10, 80]}
        self._mins_maxes.update(mins_maxes)
        self.default_metric = default_metric
        self.metric = self.default_metric

    def read_file(self, filepath, sep, read_csv_kw):
        """Convert data file to pandas dataframe and store as self.data

        Args:
            filepath (str): path to data file containing interactions
            sep (str): field separator character
            read_csv_kw (dict): kwargs dictionary passed to pd.read_csv
        """
        self.data = pd.read_csv(filepath, sep=sep, **read_csv_kw)

    @property
    def fill(self):
        """retreive fill value for currently set metric

        Returns:
            number-like: value used to fill in missing values
        """
        return self._fill_values[self.metric]

    @property
    def cmap(self):
        """Get the currect colormap

        Returns:
            matplotlib colormap: colormap for mapping data to colors
        """
        return self._cmap

    @cmap.setter
    def cmap(self, cmap):
        """Sets the colormap to be used for mapping data to colors

        Args:
            cmap (str | list): a valid matplotlib color-like, list of colors,
            colormap name or colormap object
        """
        if mp.colors.is_color_like(cmap):
            cmap = mp.colors.ListedColormap([cmap])
        elif (isinstance(cmap, list) and
              all(mp.colors.is_color_like(c) for c in cmap)):
            cmap = mp.colors.ListedColormap(cmap)
        cmap = plt.get_cmap(cmap)
        cmap = cmap(np.arange(cmap.N))
        cmap[:, -1] = np.full((len(cmap)), 0.6)  # set default alpha to 0.6
        cmap = self.modify_cmap(cmap)
        cmap = mp.colors.ListedColormap(cmap)
        self._cmap = cmap

    def modify_cmap(self, cmap):
        """Subclass specific function to modify cmaps

        Args:
            cmap (matplotlib colormap): colormap to be modified

        Returns:
            matplotlib colormap: modified colormap
        """
        return cmap

    @property
    def metric(self):
        """Retreive the currently set metric

        Returns:
            str: column name of self.data
        """
        return self._metric

    @metric.setter
    def metric(self, value):
        """Sets the metric, cmap, and min_max values. If metric is "Distance"
        or "Distance_atom", Distances are calculated.

        Args:
            value (str | tuple): valid column name of self.data. If "Distance",
                must be provided as a tuple along with the PDB object to
                compute distances
        """
        if value in self.data.keys():
            self._metric = value
        elif isinstance(value, tuple):
            value, pdb = value
            if value.startswith("Distance_"):
                value, atom = value.split("_")
                self.set_3d_distances(pdb, atom)
                self._metric = value
            elif value == "Distance":
                self.set_3d_distances(pdb, "O2'")
                self._metric = value
        elif value is None:
            self._metric = self.default_metric
        else:
            print(f"{value} is not a valid metric of {self.datatype}")
            self._metric = self.default_metric
        if self._metric in self._cmaps:
            self.cmap = self._cmaps[self._metric]
        else:
            self.cmap = "gray"
        if self._metric in self._mins_maxes:
            self.min_max = self._mins_maxes[self._metric]
        else:
            self.min_max = [min(self.data[self.metric]),
                            max(self.data[self.metric])]

    def mask_on_sequence(self, compliment_only, nts, return_mask=False):
        """Mask interactions based on sequence content

        Args:
            compliment_only (bool): require that i and j windows are reverse
                complimentary
            nts (str): require that all nucleotides in i and j windows are in
                nts
            return_mask (bool, optional): whether to return the mask instead of
                updating the current filter. Defaults to False.

        Returns:
            numpy array: the mask array (only if return_mask==True)
        """
        mask = []
        comp = {'A': 'U', 'U': 'AG', 'G': 'CU', 'C': 'G'}
        for _, i, j in self.data[["i", "j"]].itertuples():
            keep = []
            for w in range(self.window):
                i_nt = self.sequence[i+w-1].upper()
                j_nt = self.sequence[j+self.window-w-2].upper()
                if compliment_only:
                    keep.append(i_nt in comp[j_nt])
                if nts is not None:
                    keep.append(i_nt in nts and j_nt in nts)
            mask.append(all(keep))
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_ct(self, ct, min_cd=None, max_cd=None, ss_only=False,
                   ds_only=False, paired_only=False, return_mask=False):
        """Mask interactions based on secondary structure

        Args:
            ct (CT or subclass): a data object containing secondary structure
                information, or a list of these data objects, filters are applied
                based on all structures.
            min_cd (int, optional): minimum allowable contact distance. Defaults to None.
            max_cd (int, optional): maximum allowable contact distance. Defaults to None.
            ss_only (bool, optional): whether to require i and j to be single-
                stranded. Defaults to False.
            ds_only (bool, optional): whether to require i and j to be double-
                stranded. Defaults to False.
            paired_only (bool, optional): whether to require that i and j are
                base paired. Defaults to False.
            return_mask (bool, optional): whether to return the new mask array
                instead of updating the current filter. Defaults to False.

        Returns:
            numpy array: the mask array (only if return_mask==True)
        """
        if isinstance(ct, list):
            for each in ct:
                self.mask_on_ct(each, min_cd, max_cd, ss_only, ds_only,
                                paired_only)
                return
        assert ct.datatype == "ct", "CT filtering requires a ct object."
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
            if (min_cd is not None or max_cd is not None) and true_so_far:
                cd = []
                for iw in range(self.window):
                    for jw in range(self.window):
                        cd.append(ct.contactDistance(i+iw, j+jw))
                cd = min(cd)
                if min_cd is not None:
                    true_so_far = cd >= min_cd
                if max_cd is not None:
                    true_so_far = cd <= max_cd
            mask.append(true_so_far)
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_profile(self, profile, min_profile=None, max_profile=None,
                        return_mask=False):
        """Masks interactions based on per-nucleotide information. Positions
        that are not mapped to profile or are np.nan values in profile are not
        masked.

        Args:
            profile (Profile or subclass): a data object containing
                per-nucleotide information
            min_profile (float, optional): minimum allowable per-nucleotide
                value. Defaults to None.
            max_profile (float, optional): maximum allowable per-nucleotide
                value. Defaults to None.
            return_mask (bool, optional): whether to return mask instead of
                updating the filter. Defaults to False.

        Returns:
            numpy array: the mask array (only if return_mask==True)
        """
        alignment_map = self.get_alignment_map(profile)
        norm_prof = profile.data["Norm_profile"]
        mask = np.full(len(self.data), True)
        for idx, i, j in self.data[["i", "j"]].itertuples():
            index_i = alignment_map[i-1]
            index_j = alignment_map[j-1]
            keep_ij = True
            if (index_i != -1) and (index_j != -1):
                prof_i = np.nanmedian(norm_prof[index_i:index_i+self.window])
                prof_j = np.nanmedian(norm_prof[index_j:index_j+self.window])
                if min_profile is not None:
                    keep_ij &= (prof_i >= min_profile) | np.isnan(prof_i)
                    keep_ij &= (prof_j >= min_profile) | np.isnan(prof_j)
                if max_profile is not None:
                    keep_ij &= (prof_i <= max_profile) | np.isnan(prof_i)
                    keep_ij &= (prof_j <= max_profile) | np.isnan(prof_j)
            mask[idx] &= keep_ij
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_position(self, exclude, isolate, return_mask=False):
        """Masks interactions based on position in sequence

        Args:
            exclude (list of int, optional): a list of nucleotide positions to
                exclude if i or j is in list
            isolate (list of int, optional): a list of nucleotide positions to
                isolate if i and j are not in list
            return_mask (bool, optional): whether to return mask instead of
                updating the filter. Defaults to False.

        Returns:
            numpy array: the mask array (only if return_mask==True)
        """
        mask = np.full(len(self.data), True)
        for index, (i, j) in self.data[["i", "j"]].iterrows():
            if exclude is not None:
                if (i in exclude) or (j in exclude):
                    mask[index] = False
            if isolate is not None:
                if (i not in isolate) and (j not in isolate):
                    mask[index] = False
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_distance(self, max_dist, min_dist, return_mask=False):
        """Mask interactions based on their primary sequence distance (j-i).

        Args:
            max_dist (int): maximum allowable distance
            min_dist (int): minimum allowable distance
            return_mask (bool, optional): whether to return mask instead of
                updating the filter. Defaults to False.

        Returns:
            numpy array: the mask array (only if return_mask==True)
        """
        primary_distances = np.absolute(self.data.eval("i - j"))
        mask = self.data['mask']
        if min_dist is not None:
            mask &= primary_distances >= min_dist
        if max_dist is not None:
            mask &= primary_distances <= max_dist
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_values(self, **kwargs):
        """Mask interactions on values in self.data. Each keyword should have
        the format "column_operator" where column is a valid column name of
        the dataframe and operator is one of:
            "ge": greater than or equal to
            "le": less than or equal to
            "gt": greater than
            "lt": less than
            "eq": equal to
            "ne": not equal to
        The values given to these keywords are then used in the comparison and
        False comparisons are filtered out. e.g.:
            self.mask_on_values(Statistic_ge=23) evaluates to:
            self.update_mask(self.data["Statistic"] >= 23)
        """
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

    def set_mask_offset(self, fit_to, prefiltered=False):
        """Define new i and j positions that map to the provided fit_to
        sequence. Interactions in which i or j does not map are masked. If
        fit_to object is a PDB and i and j are not in the structure, they are
        masked.

        Args:
            fit_to (Data or subclass): a data object containing a sequence
            prefiltered (bool, optional): if True, original mask values are
                kept. Defaults to False.
        """
        am = self.get_alignment_map(fit_to)
        i = np.array([am[i-1]+1 for i in self.data["i"].values])
        j = np.array([am[j-1]+1 for j in self.data["j"].values])
        if prefiltered:
            mask = self.data["mask"].copy()
        else:
            mask = np.ones(len(i), dtype=bool)
        for w in range(self.window):
            mask = mask & ((i+w) != 0) & ((j+w) != 0)
        self.data["mask"] = mask
        self.data['i_offset'] = i
        self.data['j_offset'] = j
        if fit_to.datatype == 'pdb':
            mask = []
            for i, j in zip(self.data["i_offset"], self.data["j_offset"]):
                mask.append(fit_to.is_valid_idx(seq_idx=i)
                            and fit_to.is_valid_idx(seq_idx=j))
            mask = np.array(mask, dtype=bool)
            self.update_mask(mask)

    def update_mask(self, mask):
        """Given a new masking array, the mask is updated

        Args:
            mask (numpy array of bool): the new mask to update on
        """
        self.data["mask"] = self.data["mask"] & mask

    def filter(self, fit_to,
               prefiltered=False,
               ct=None,  # required if any of next are passed
               min_cd=None, max_cd=None,
               paired_only=False, ss_only=False, ds_only=False,
               profile=None, min_profile=None, max_profile=None,
               compliments_only=False, nts=None,
               max_distance=None, min_distance=None,
               exclude_nts=None, isolate_nts=None,
               resolve_conflicts=None,
               **kwargs):
        """Convenience function that applies the above filters simultaneously.

        Args:
            fit_to (Data or subclass): passed to self.set_mask_update()
            prefiltered (bool, optional): passed to self.set_mask_update().
                Defaults to False.
            ct (CT or subclass, optional): passed to self.mask_on_ct().
                Defaults to None.
            min_cd (int, optional): passed to self.mask_on_ct(). Defaults to
                None.
            max_cd (int, optional): passed to self.mask_on_ct(). Defaults to
                None.
            paired_only (bool, optional): passed to self.mask_on_ct(). Defaults
                to False.
            ss_only (bool, optional): passed to self.mask_on_ct(). Defaults to
                False.
            ds_only (bool, optional): passed to self.mask_on_ct(). Defaults to
                False.
            profile (Profile or subclass, optional): passed to
                self.mask_on_profile(). Defaults to None.
            min_profile (float, optional): passed to self.mask_on_profile().
                Defaults to None.
            max_profile (float, optional): passed to self.mask_on_profile().
                Defaults to None.
            compliments_only (bool, optional): passed to
                self.mask_on_sequence(). Defaults to False.
            nts (str, optional): passed to self.mask_on_sequence(). Defaults
                to None.
            max_distance (int, optional): passed to self.mask_on_distance().
                Defaults to None.
            min_distance (int, optional): passed to self.mask_on_distance().
                Defaults to None.
            exclude_nts (list of int, optional): passed to
                self.mask_on_position(). Defaults to None.
            isolate_nts (list of int, optional): passed to
                self.mask_on_position(). Defaults to None.
            resolve_conflicts (str, optional): passed to
                self.resolve_conflicts(). Defaults to None.
            **kwargs: additional arguments are first passed to
                self.data_specific_filter(), remaining kwargs are passed to
                self.mask_on_values()
        """
        def filters_are_on(*filters):
            return any(f not in [None, False] for f in filters)

        self.set_mask_offset(fit_to, prefiltered=prefiltered)
        # TODO: decide what prefiltered does exactly...
        if prefiltered:
            return
        if filters_are_on(exclude_nts, isolate_nts):
            self.mask_on_position(exclude_nts, isolate_nts)
        if filters_are_on(max_distance, min_distance):
            self.mask_distance(max_dist=max_distance, min_dist=min_distance)
        if filters_are_on(min_profile, max_profile):
            self.mask_on_profile(profile, min_profile, max_profile)
        if filters_are_on(min_cd, max_cd, ss_only, ds_only, paired_only):
            self.mask_on_ct(ct, min_cd, max_cd, ss_only, ds_only, paired_only)
        if filters_are_on(compliments_only, nts):
            self.mask_on_sequence(compliments_only, nts)
        kwargs = self.data_specific_filter(**kwargs)
        self.mask_on_values(**kwargs)
        if resolve_conflicts is not None:
            self.resolve_conflicts(resolve_conflicts)

    def data_specific_filter(self, **kwargs):
        """Does nothing for the base Interactions class, can be overwritten in
        subclasses.

        Returns:
            dict: dictionary of keyword argument pairs
        """
        return kwargs

    def get_ij_colors(self, min_max=None, cmap=None):
        """"Gets i, j, and colors lists for plotting interactions. i and j are
        the 5' and 3' ends of each interaction, and colors is the color to use
        for each interaction. Values of self.data[self.metric] are normalized
        to 0 to 1, which correspond to self.min_max values. These are then
        mapped to a color using self.cmap.

        Args:
            min_max (list of float, optional): overrides self.min_max.
                Defaults to None.
            cmap (str, optional): a colormap name which overrides self.cmap.
                Defaults to True.

        Returns:
            list, list, list: 5' and 3' ends of each pair, color for each pair
        """
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
                if np.isnan(datum):
                    colors.append('grey')
                else:
                    colors.append(cmap(datum))
        return i_list, j_list, colors

    def get_normalized_ij_data(self, min_max):
        """Retreives values in self.data[self.metric] that are not masked,
        normalized between 0 and 1, which correspond to min_max values.

        Args:
            min_max (list of float): min and max used to normalize values

        Returns:
            numpy array: filtered and normalized interactions values
        """
        metric = self.metric
        columns = ["i_offset", "j_offset", metric]
        data = self.data.loc[self.data["mask"], columns].copy()
        ascending = metric not in ['Distance']  # high distance is bad
        data.sort_values(by=metric, ascending=ascending, inplace=True,
                         na_position='first')
        norm = plt.Normalize(min_max[0], min_max[1], clip=True)
        data[metric] = norm(data[metric])
        return data

    def print_new_file(self, outfile=None):
        """Prints a new file containing repositioned and filtered interactions
        in the original format

        Args:
            outfile (str, optional): path to an output file. If None, file
                string is printed to console. Defaults to None.
        """
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
        """Creates or overwrites values in self.data["Distance"] by calculating
        the distance between atoms in i and j in the PDB structure.

        Args:
            pdb (PDB): a data object containing atomic coordinates
            atom (str): an atom id
        """
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


class SHAPEJuMP(Interactions):
    def __init__(self, filepath, datatype="shapejump", fasta=None,
                 sequence=None):
        """Constructs an Interactions object from SHAPEJuMP data

        Args:
            filepath (str): path to ShapeJumper deletions file
            datatype (str, optional): stored as self.datatype. Defaults to
                "shapejump".
            fasta (str, optional): path to fasta file. Defaults to None.
            sequence (str, optional): a sequence string. Defaults to None.
        """
        default_metric = 'Percentile'
        fill = {'Metric': 0.0, 'Percentile': 0.0}
        cmaps = {'Metric': 'YlGnBu', 'Percentile': 'YlGnBu'}
        mins_maxes = {'Percentile': [0.98, 1.0], "Metric": [0, 0.001]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, fasta=fasta,
                         sequence=sequence, fill=fill, cmaps=cmaps,
                         mins_maxes=mins_maxes)

    def read_file(self, filepath, sep=None, read_csv_kw=None):
        """Parses a deletions.txt file and stores data as a dataframe at
        self.data, sets self.window=1, and calculates a "Percentile" column.

        Args:
            filepath (str): path to deletions.txt file
            sep (str, optional): passed to pandas.read_csv(). Defaults to None.
            read_csv_kw (dict, optional): kwargs passed to pandas.read_csv().
                Defaults to None.
        """
        column_names = ['Gene', 'i', 'j', 'Metric']
        data = pd.read_csv(filepath, sep='\t', names=column_names, header=0)
        data["Percentile"] = data['Metric'].rank(method='max', pct=True)
        self.data = data
        self.window = 1


class RINGMaP(Interactions):
    def __init__(self, filepath, datatype="ringmap", sequence=None):
        """Constructs an Interactions object from RING-MaP data

        Args:
            filepath (str): path to RingMapper correlations file
            datatype (str, optional): stored as self.datatype. Defaults to
                "ringmap".
            fasta (str, optional): path to fasta file. Defaults to None.
            sequence (str, optional): a sequence string. Defaults to None.
        """
        default_metric = 'Statistic'
        fill = {'Statistic': 0.0, 'Zij': 0.0}
        cmaps = {'Statistic': 'bwr', 'Zij': 'bwr'}
        mins_maxes = {"Statistic": [-100, 100], 'Zij': [-8, 8]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, sequence=sequence,
                         fill=fill, cmaps=cmaps, mins_maxes=mins_maxes)

    def read_file(self, filepath, sep=None, read_csv_kw=None):
        """Parses a correlations file and stores data as a dataframe at
        self.data, sets self.window=1, and renames "+/-" column to "Sign".

        Args:
            filepath (str): path to a correlations file
            sep (str, optional): passed to pandas.read_csv(). Defaults to None.
            read_csv_kw (dict, optional): kwargs passed to pandas.read_csv().
                Defaults to None.
        """
        with open(filepath, 'r') as file:
            self.header = file.readline()
        split_header = self.header.split('\t')
        window = split_header[1].split('=')[1]
        self.window = int(window)
        self.data = pd.read_csv(filepath, sep='\t', header=1)
        self.data.rename(columns={"+/-": "Sign"}, inplace=True)

    def data_specific_filter(self, positive_only=False, negative_only=False,
                             **kwargs):
        """Adds filters for "Sign" column to parent filter() function

        Args:
            positive_only (bool, optional): whether to require that sign is 1.
                Defaults to False.
            negative_only (bool, optional): whether to require that sign is -1.
                Defaults to False.

        Returns:
            dict: any additional keyword-argument pairs are returned
        """
        if positive_only:
            self.update_mask(self.data["Sign"] == 1)
        if negative_only:
            self.update_mask(self.data["Sign"] == -1)
        return kwargs

    def get_normalized_ij_data(self, min_max):
        """Overwrites parent get_normalized_ij_data. Uses these values instead:
            self.data[self.metric] * self.data["Sign"]
        Except when self.metric is "Distance".

        Args:
            min_max (list of float): min and max values for normalization

        Returns:
            numpy array: normalized values for color mapping
        """
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


class PAIRMaP(Interactions):
    def __init__(self, filepath, datatype="pairmap", sequence=None):
        """Constructs an Interactions object from PAIR-MaP data

        Args:
            filepath (str): path to PAIR-MaP pairmap.txt file
            datatype (str, optional): stored as self.datatype. Defaults to
                "ringmap".
            fasta (str, optional): path to fasta file. Defaults to None.
            sequence (str, optional): a sequence string. Defaults to None.
        """
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

    def read_file(self, filepath, sep=None, read_csv_kw=None):
        """Parses a pairmap.txt file and stores data as a dataframe at
        self.data, sets self.window (usually 3, from header).

        Args:
            filepath (str): path to a PAIR-MaP pairmap.txt file
            sep (str, optional): passed to pandas.read_csv(). Defaults to None.
            read_csv_kw (dict, optional): kwargs passed to pandas.read_csv().
                Defaults to None.
        """
        with open(filepath, 'r') as file:
            self.header = file.readline()
        self.window = int(self.header.split('\t')[1].split('=')[1])
        self.data = pd.read_csv(filepath, sep='\t', header=1)
        self.data.rename(columns={"Sig.": "Statistic"}, inplace=True)

    def modify_cmap(self, cmap):
        """Changes the alpha value of cmap for non-primary and -secondary PAIRs

        Args:
            cmap (mpl colormap): colormap from Interactions.cmap setter

        Returns:
            mpl colormap: modified colormap
        """
        if self.metric == 'Class':
            cmap[0, -1] = 0.2  # alpha of non 1ary and 2ary pairs to 0.2
        return cmap

    def data_specific_filter(self, all_pairs=False, **kwargs):
        """Used by Interactions.filter(). By default, non-primary and
        -secondary pairs are removed. all_pairs=True changes this behavior.

        Args:
            all_pairs (bool, optional): whether to include all PAIRs.
                Defaults to False.

        Returns:
            dict: remaining kwargs are passed back to Interactions.filter()
        """
        if not all_pairs:
            self.update_mask(self.data["Class"] != 0)
        return kwargs

    def get_normalized_ij_data(self, min_max):
        """Same as parent function, unless metric is set to "Class", in which
        case ij pairs are returned in a different order.

        Args:
            min_max (list of int, length 2): minimum and maximum bounds for
                colormapping

        Returns:
            pandas DataFrame: Dataframe providing i, j, and normalized data
                values for plotting
        """
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


class PairProb(Interactions):
    def __init__(self, filepath, datatype="pairprob", sequence=None):
        """Constructs Interactions data from a pairing probability text file
        containing i, j, and -log10(P) values. Can be obtained using partition
        and ProbabilityPlot functions from RNAStructure (Matthews Lab).

        Args:
            filepath (str): path to pairing probability text file
            datatype (str, optional): "pairprob". Defaults to "pairprob".
            sequence (str, optional): Sequence string. Defaults to None.
        """
        default_metric = 'Probability'
        fill = {'Probability': 0}
        cmaps = {'Probability': sns.cubehelix_palette(10, 0.7, 0.9, 1.5, 2.5,
                                                      1, 0.4, False, True)}
        mins_maxes = {'Probability': [0, 1]}
        super().__init__(filepath=filepath, datatype=datatype,
                         default_metric=default_metric, sequence=sequence,
                         fill=fill, cmaps=cmaps, mins_maxes=mins_maxes)

    def read_file(self, filepath, sep=None, read_csv_kw=None):
        """Parses a pairing probability text file to create a DataFrame
        containing i, j, -log10(P) and Probability (0-1).

        Args:
            filepath (str): path to pairing probability text file
            sep (None, optional): ignored. Defaults to None.
            read_csv_kw (None, optional): ignored. Defaults to None.
        """
        with open(filepath, 'r') as file:
            self.header = file.readline()
        lengths_match = int(self.header.strip()) == self.length
        assert lengths_match, "DP and sequence are different lengths"
        self.window = 1
        data = pd.read_table(filepath, header=1, names=['i', 'j', 'log10p'])
        data["Probability"] = 10 ** (-data["log10p"])
        self.data = data

    def set_entropy(self, printOut=False, toFile=None):
        """Calculates per-nucleotide Shannon entropy and stores as self.entropy

        Args:
            printOut (bool, optional): whether to print the result.
                Defaults to False.
            toFile (str, optional): file to write result to. Defaults to None.
        """
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
        """Used by parent filter function. By default, filters out pairs with
        probability less that 3%

        Returns:
            dict: keyword arguments are passed back to Interactions.filter()
        """
        self.update_mask(self.data["Probability"] >= 0.03)
        return kwargs


class AllPossible(Interactions):
    def __init__(self, filepath, sequence=None, window=1):
        """Constructs Interactions data from a sequence. One interaction will
        be made for every possible pair of nucleotides.

        Args:
            filepath (None): Ignored.
            sequence (str, optional): Sequence string. Defaults to None.
            window (int, optional): Window size of each interaction.
                Defaults to 1.
        """
        data = {'i': [], 'j': [], 'data': []}
        for i in range(1, len(sequence)):
            for j in range(i+1, len(sequence)+1):
                data['i'].append(i)
                data['j'].append(j)
                data['data'].append(1)
        data = pd.DataFrame(data)
        super().__init__(datatype="interactions", dataframe=data,
                         filepath=filepath,
                         default_metric='data', window=window,
                         sequence=sequence, fill={'data': 1},
                         cmaps={'data': 'magenta'},
                         mins_maxes={'data': [0, 2]})
