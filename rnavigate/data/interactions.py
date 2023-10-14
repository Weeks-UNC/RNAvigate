from operator import ge, le, gt, lt, eq, ne
import pandas as pd
import matplotlib.colors as mpc
import numpy as np
import seaborn as sns
from rnavigate import data


class Interactions(data.Data):
    def __init__(self, input_data, sequence, metric, metric_defaults,
                 read_table_kw=None, window=1):
        """Given a dataframe or a data file, construct the interactions object

        Args:
            input_data (str | pandas.DataFrame):
                path to a file or dataframe containing interactions data.
                Must have at least "i" and "j" columns indicating the 5' and 3'
                ends of the interactions.
            sequence (str | pandas.DataFrame):
                sequence string, fasta file, or a pandas dataframe containing
                a "Sequence" column.
            default_metric (str, optional): column name to use as the default
                metric. Defaults to None.
            read_table_kw (dict, optional): other options for read_table.
                Defaults to {}.
            window (int, optional): 5' and 3' interactions windows.
                Defaults to 1.
            fasta (str, optional): path to fasta file. Defaults to None.
            fill (dict, optional): dictionary specifying a fill value (values)
                to use with a metric (keys). Defaults to {}.
            cmaps (dict, optional): specifies cmaps (values) to use with
                metrics (keys). Defaults to {}.
            mins_maxes (dict, optional): specifies minimum and maximum values
                as a list of floats (values) to use with given metric (keys).
                Defaults to {}.
        """
        self.window = window
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw)
        self.reset_mask()

    def mask_on_sequence(self, compliment_only=None, nts=None,
                         return_mask=False):
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
        compliment = {"A": "UT", "U": "AG", "T": "AG", "G": "CU", "C": "G"}
        for _, i, j in self.data[["i", "j"]].itertuples():
            keep = []
            for w in range(self.window):
                i_nt = self.sequence[i + w - 1].upper()
                j_nt = self.sequence[j + self.window - w - 2].upper()
                if compliment_only:
                    keep.append(i_nt in compliment[j_nt])
                if nts is not None:
                    keep.append(i_nt in nts and j_nt in nts)
            mask.append(all(keep))
        if return_mask:
            return mask
        else:
            self.update_mask(mask)

    def mask_on_structure(self, structure, min_cd=None, max_cd=None,
                          ss_only=False, ds_only=False, paired_only=False,
                          return_mask=False):
        """Mask interactions based on secondary structure

        Args:
            structure (SecondaryStructure):
                A SecondaryStructure object or list of these objects
                Filters are applied based on all structures.
            min_cd (int, optional): minimum allowable contact distance.
                Defaults to None.
            max_cd (int, optional): maximum allowable contact distance.
                Defaults to None.
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
        if isinstance(structure, list):
            for each in structure:
                self.mask_on_structure(
                    each, min_cd, max_cd, ss_only, ds_only, paired_only)
                return
        alignment = data.SequenceAlignment(self, structure)
        mask = []
        i_j_keep = ["i", "j", "mask"]
        for _, i, j, keep in self.data[i_j_keep].itertuples():
            i = alignment.map_positions(i)
            j = alignment.map_positions(j)
            true_so_far = keep and (-1 not in [i, j])
            if paired_only and true_so_far:
                for w in range(self.window):
                    true_so_far = structure.ct[i+w-1] == j+self.window-w-1
                    if not true_so_far:
                        break
            if (ss_only or ds_only) and true_so_far:
                percentage = 0
                for w in range(self.window):
                    percentage += int(structure.ct[i-1+w] == 0)
                    percentage += int(structure.ct[j-1+w] == 0)
                percentage /= (self.window * 2)
                if ss_only:
                    true_so_far = percentage > 0.501
                if ds_only:
                    true_so_far = percentage < 0.501
            if (min_cd is not None or max_cd is not None) and true_so_far:
                cd = []
                for iw in range(self.window):
                    for jw in range(self.window):
                        cd.append(structure.contact_distance(i + iw, j + jw))
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

    def mask_on_profile(
            self, profile, min_profile=None, max_profile=None,
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
        mapping = data.SequenceAlignment(self, profile).mapping
        norm_prof = profile.data["Norm_profile"]
        mask = np.full(len(self.data), True)
        for idx, i, j in self.data[["i", "j"]].itertuples():
            index_i = mapping[i - 1]
            index_j = mapping[j - 1]
            keep_ij = True
            if (index_i != -1) and (index_j != -1):
                prof_i = np.nanmedian(norm_prof[index_i:index_i + self.window])
                prof_j = np.nanmedian(norm_prof[index_j:index_j + self.window])
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

    def mask_on_position(
            self, exclude=None, isolate=None, return_mask=False):
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

    def mask_on_distance(
            self, max_dist=None, min_dist=None, return_mask=False):
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
        mask = self.data["mask"]
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
        for key, value in kwargs.items():
            if key in self.data.keys():
                self.update_mask(self.data[key] > value)
            elif "_" in key:
                key2, comparison = key.rsplit("_", 1)
                operators = {
                    "ge": ge,
                    "le": le,
                    "gt": gt,
                    "lt": lt,
                    "eq": eq,
                    "ne": ne}
                if key2 in self.data.columns and comparison in operators:
                    operator = operators[comparison]
                    self.update_mask(operator(self.data[key2], value))
                else:
                    print(f"{key}={value} is not a valid filter.")
            else:
                print(f"{key}={value} is not a valid filter.")

    def reset_mask(self):
        """Resets the mask to all True (removes previous filters)"""
        self.data['mask'] = np.ones(len(self.data), dtype=bool)

    def get_aligned_data(self, alignment):
        """Get a new copy of the data with i and j mapped to new positions
        using an alignment. Interactions in which i or j does not map are
        dropped.
        """
        new_data = alignment.map_dataframe(
            dataframe = self.data[self.data['mask']],
            position_columns=['i', 'j'])
        return self.__class__(
            input_data=new_data,
            sequence=alignment.target_sequence,
            metric=self._metric,
            metric_defaults=self.metric_defaults,
            window=self.window)

    def update_mask(self, mask):
        """Given a new masking array, the mask is updated

        Args:
            mask (numpy array of bool): the new mask to update on
        """
        self.data["mask"] = self.data["mask"] & mask

    def filter(
            self, prefiltered=False,
            # mask on structure
            structure=None, min_cd=None, max_cd=None, paired_only=False,
            ss_only=False, ds_only=False,
            # mask on profile
            profile=None, min_profile=None, max_profile=None,
            # mask on sequence
            compliments_only=False, nts=None,
            # mask on position
            max_distance=None, min_distance=None, exclude_nts=None,
            isolate_nts=None,
            # others
            resolve_conflicts=None, **kwargs):
        """Convenience function that applies the above filters simultaneously.

        Args:
            prefiltered (bool, optional): passed to self.set_mask_update().
                Defaults to False.
            structure (SecondaryStructure, optional):
                passed to self.mask_on_structure(). Defaults to None.
            min_cd (int, optional): passed to self.mask_on_structure(). Defaults to
                None.
            max_cd (int, optional): passed to self.mask_on_structure(). Defaults to
                None.
            paired_only (bool, optional): passed to self.mask_on_structure(). Defaults
                to False.
            ss_only (bool, optional): passed to self.mask_on_structure(). Defaults to
                False.
            ds_only (bool, optional): passed to self.mask_on_structure(). Defaults to
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

        if prefiltered:
            return
        self.reset_mask()
        if filters_are_on(exclude_nts, isolate_nts):
            self.mask_on_position(exclude_nts, isolate_nts)
        if filters_are_on(max_distance, min_distance):
            self.mask_on_on_distance(max_dist=max_distance, min_dist=min_distance)
        if filters_are_on(min_profile, max_profile):
            self.mask_on_profile(profile, min_profile, max_profile)
        if filters_are_on(min_cd, max_cd, ss_only, ds_only, paired_only):
            self.mask_on_structure(
                structure, min_cd, max_cd, ss_only, ds_only, paired_only)
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

    def get_ij_colors(self):
        """Gets i, j, and colors lists for plotting interactions. i and j are
        the 5' and 3' ends of each interaction, and colors is the color to use
        for each interaction. Values of self.data[self.metric] are normalized
        to 0 to 1, which correspond to self.min_max values. These are then
        mapped to a color using self.cmap.

        Returns:
            list, list, list: 5' and 3' ends of each pair, color for each pair
        """
        if len(self.data) == 0:
            return [], [], []

        dataframe = self.get_sorted_data()
        i_list, j_list, colors = [], [], []
        for _, i, j, datum in dataframe[['i','j', self.metric]].itertuples():
            for w in range(self.window):
                i_list.append(i + w)
                j_list.append(j + self.window - 1 - w)
                colors.append(datum)
        colors = self.cmap.values_to_hexcolors(colors, 0.6)
        return i_list, j_list, colors

    def get_sorted_data(self):
        """Returns a sorted copy of interactions data.

        Returns:
            pandas.DataFrame: i, j, and metric values, sorted
        """
        ascending = self.metric not in ["Distance"]  # high distance is bad
        dataframe = self.data[['i', 'j', self.metric]].copy()
        dataframe.sort_values(
            by=self.metric,
            ascending=ascending,
            inplace=True,
            na_position="first")
        return dataframe

    def print_new_file(self, outfile=None):
        """Prints a new file containing repositioned and filtered interactions
        in the original format

        Args:
            outfile (str, optional): path to an output file. If None, file
                string is printed to console. Defaults to None.
        """
        data = self.data.copy()
        data = data[data["mask"]]
        if "Sign" in data.columns:
            data.rename({"Sign": "+/-"}, inplace=True)
        csv = data.to_csv(
            columns=self.columns, sep="\t", index=False, line_terminator="\n"
        )
        if outfile is not None:
            with open(outfile, "w") as out:
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
        mapping = data.SequenceAlignment(self, pdb).mapping
        distance_matrix = pdb.get_distance_matrix(atom=atom)
        i = self.data["i"].values
        j = self.data["j"].values
        distances = np.zeros(len(i))
        pairs = self.window**2
        for iw in range(self.window):
            for jw in range(self.window):
                io = mapping[i + iw - 1]
                jo = mapping[j + jw - 1]
                distances += distance_matrix[io, jo] / pairs
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
            if len(graph) == 0:
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
                if neighbor in new_graph:
                    del new_graph[neighbor]
            case_2 = [current_vertex] + graphSets(new_graph, node_weights)
            # Our final result is the one which has highest weight, return it
            weight_1 = sum(node_weights[node] for node in case_1)
            weight_2 = sum(node_weights[node] for node in case_2)
            if weight_1 > weight_2:
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
            conflicts = np.where(
                (
                    (abs(data["i"] - i) < window)
                    | (abs(data["j"] - i) < window)
                    | (abs(data["i"] - j) < window)
                    | (abs(data["j"] - j) < window)
                )
                & ((data["j"] - j) != (i - data["i"]))
            )[0]
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
    def __init__(self, input_data, sequence=None, metric='Percentile',
                 metric_defaults=None, read_table_kw=None, window=1):
        """Constructs an Interactions object from SHAPEJuMP data"""
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            'Percentile': {
                'metric_column': 'Percentile',
                'cmap': 'YlGnBu',
                'normalization': 'min_max',
                'values': [0.98, 1.0],
                'title': 'SHAPE-JuMP: percentile',
                'extend': 'min'},
            'Metric': {
                'metric_column': 'Metric',
                'cmap': 'YlGnBu',
                'normalization': 'min_max',
                'values': [0, 0.001],
                'extend': 'max',
                'title': 'SHAPE-JuMP: rate'},
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=1)

    def read_file(self, input_data, read_table_kw=None):
        """Parses a deletions.txt file and stores data as a dataframe at
        self.data, sets self.window=1, and calculates a "Percentile" column.

        Args:
            input_data (str): path to deletions.txt file
            read_table_kw (dict, optional): kwargs passed to pandas.read_table().
                Defaults to None.
        """
        column_names = ["Gene", "i", "j", "Metric"]
        data = pd.read_table(input_data, names=column_names, header=0,
                           **read_table_kw)
        data["Percentile"] = data["Metric"].rank(method="max", pct=True)
        return data


class RINGMaP(Interactions):
    def __init__(self, input_data, sequence=None, metric='Statistic',
                 metric_defaults=None, read_table_kw=None, window=1):
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            'Statistic': {
                'metric_column': 'Statistic',
                'cmap': 'bwr',
                'normalization': 'min_max',
                'values': [-100, 100],
                'title': 'RING-MaP: Gapc',
                'extend': 'both'},
            'Zij': {
                'metric_column': 'Zij',
                'cmap': 'bwr',
                'normalization': 'min_max',
                'values': [-8, 8],
                'title': 'RING-MaP: Zij',
                'extend': 'both'}
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window)

    def read_file(self, filepath, read_table_kw=None):
        """Parses a correlations file and stores data as a dataframe at
        self.data, sets self.window=1, and renames "+/-" column to "Sign".

        Args:
            filepath (str):
                path to a RingMapper correlations file
            read_table_kw (dict, optional):
                kwargs passed to pandas.read_table().
                Defaults to None.
        """
        with open(filepath, "r") as file:
            self.header = file.readline()
        split_header = self.header.split("\t")
        window = split_header[1].split("=")[1]
        self.window = int(window)
        dataframe = pd.read_table(filepath, header=1, **read_table_kw)
        dataframe.rename(columns={"+/-": "Sign"}, inplace=True)
        return dataframe

    def data_specific_filter(
            self, positive_only=False, negative_only=False, **kwargs):
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

    def get_sorted_data(self):
        """Overwrites parent get_normalized_ij_data. Uses these values instead:
            self.data[self.metric] * self.data["Sign"]
        Except when self.metric is "Distance".

        Args:
            min_max (list of float): min and max values for normalization

        Returns:
            numpy array: normalized values for color mapping
        """
        if self.metric == "Distance":
            return super().get_sorted_data()
        columns = ['i', 'j', self.metric]
        dataframe = self.data[columns+["Sign"]].copy()
        dataframe.eval(f"{self.metric} = {self.metric} * Sign", inplace=True)
        dataframe.sort_values(by=self.metric, ascending=True, inplace=True)
        return dataframe[columns]


class PAIRMaP(RINGMaP):
    def __init__(self, input_data, sequence=None, metric='Class',
                 metric_defaults=None, read_table_kw=None, window=1):
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            'Class': {
                'metric_column': 'Class',
                'cmap': mpc.ListedColormap([
                    [0.7, 0.7, 0.7],
                    [0.0, 0.0, 0.95],
                    [0.12, 0.76, 1.0]]),
                'normalization': 'none',
                'values': None,
                'title': 'PAIR-MaP',
                'extend': 'both',
                'ticks': [0, 1, 2],
                'tick_labels': ["Complimentary", "Primary", "Secondary"]}
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window)

    def read_file(self, filepath, read_table_kw=None):
        """Parses a pairmap.txt file and stores data as a dataframe at
        self.data, sets self.window (usually 3, from header).

        Args:
            filepath (str): path to a PAIR-MaP pairmap.txt file
            read_table_kw (dict, optional): kwargs passed to pandas.read_table().
                Defaults to None.
        """
        with open(filepath, "r") as file:
            self.header = file.readline()
        self.window = int(self.header.split("\t")[1].split("=")[1])
        dataframe = pd.read_table(filepath, header=1)
        dataframe.rename(columns={"Sig.": "Statistic"}, inplace=True)
        dataframe['Sign'] = 1
        return dataframe

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

    def get_sorted_data(self):
        """Same as parent function, unless metric is set to "Class", in which
        case ij pairs are returned in a different order.

        Args:
            min_max (list of int, length 2): minimum and maximum bounds for
                colormapping

        Returns:
            pandas DataFrame: Dataframe providing i, j, and normalized data
                values for plotting
        """
        if self.metric != "Class":
            return super().get_sorted_data()
        metric = self.metric
        columns = ["i", "j", metric]
        dataframe = self.data[columns].copy()
        # Class data have a weird order
        dataframe = pd.concat(
            [dataframe[dataframe[metric] == 0],
             dataframe[dataframe[metric] == 2],
             dataframe[dataframe[metric] == 1]])
        return dataframe


class PairingProbability(Interactions):
    def __init__(self, input_data, sequence=None, metric='Probability',
                 metric_defaults=None, read_table_kw=None, window=1):
        """Constructs Interactions data from a pairing probability text file
        containing i, j, and -log10(P) values. Can be obtained using partition
        and ProbabilityPlot functions from RNAStructure (Matthews Lab).

        Args:
            input_data (str): path to pairing probability text file
            sequence (str, optional): Sequence string. Defaults to None.
        """
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            'Probability': {
                'metric_column': 'Probability',
                'cmap': sns.cubehelix_palette(
                    10, 0.7, 0.9, 1.5, 2.5, 1, 0.4, False, True),
                'normalization': 'min_max',
                'values': [0, 1],
                'title': 'Pairing probability',
                'extend': 'neither'},
            'Probability_old': {
                'metric_column': 'Probability',
                'cmap': mpc.ListedColormap([
                    (150,150,150),
                    (255,204,0),
                    (72,143,205),
                    (81, 184, 72)]),
                'normalization': 'bins',
                'values': [0.1, 0.3, 0.8],
                'extend': 'neither',
                'title': 'Pairing probability'},
            'Probability_continuous': {
                'metric_column': 'Probability',
                'cmap': 'plasma_r',
                'normalization': 'min_max',
                'values': [0.0, 1.0],
                'extend': 'neither',
                'title': 'Pairing probability'}
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window)

    def read_file(self, filepath, read_table_kw=None):
        """Parses a pairing probability text file to create a DataFrame
        containing i, j, -log10(P) and Probability (0-1).

        Args:
            filepath (str): path to pairing probability text file
            read_table_kw (None, optional): ignored. Defaults to None.
        """
        with open(filepath, "r") as file:
            self.header = file.readline()
        self.window = 1
        dataframe = pd.read_table(filepath, header=1, names=["i", "j", "log10p"])
        dataframe["Probability"] = 10 ** (-dataframe["log10p"])
        return dataframe

    def get_entropy_profile(self, print_out=False, save_file=None):
        """Calculates per-nucleotide Shannon entropy and stores as self.entropy

        Args:
            print_out (bool, optional): whether to print the result.
                Defaults to False.
            save_file (str, optional): file to write result to.
                Defaults to None.
        """
        self.data.eval("nlogn = log10p * 10 ** ( - log10p )", inplace=True)
        entropy = np.zeros(self.length)
        args = {'labels': np.arange(self.length)+1, 'fill_value': 0}
        i_sum = self.data[['i', 'nlogn']].groupby('i').sum().reindex(**args)
        j_sum = self.data[['j', 'nlogn']].groupby('j').sum().reindex(**args)
        ij_sum = i_sum + j_sum
        entropy_df = pd.DataFrame({
            'Sequence':list(self.sequence),
            'Nucleotide': np.arange(self.length)+1,
            'Entropy': ij_sum['nlogn']
            })
        entropy_df = entropy_df.set_index('Nucleotide')
        # catch rounding errors:
        entropy_df.loc[entropy_df['Entropy'] < 0, 'Entropy'] = 0
        if print_out:
            print(*[f"{i+1} {s}" for i, s in enumerate(entropy)], sep="\n")
        if save_file is not None:
            with open(save_file) as outf:
                for i, s in enumerate(entropy):
                    outf.write(f"{i+1}\t{s}\n")
        return data.Profile(
            input_data=entropy_df,
            metric='Entropy',
            metric_defaults={
                'Entropy': {
                    'metric_column': 'Entropy',
                    'error_column': None,
                    'color_column': None,
                    'cmap': 'rainbow',
                    'normalization': 'norm',
                    'values': [0.0, 0.2],
                    'extend': 'right',
                    'title': 'Shannon entropy',
                    'alpha': 1
            }})



    def data_specific_filter(self, **kwargs):
        """Used by parent filter function. By default, filters out pairs with
        probability less that 3%

        Returns:
            dict: keyword arguments are passed back to Interactions.filter()
        """
        self.update_mask(self.data["Probability"] >= 0.03)
        return kwargs


class AllPossible(Interactions):
    def __init__(self, sequence, metric='data', input_data=None,
                 metric_defaults=None, read_table_kw=None, window=1):
        if isinstance(sequence, data.Sequence):
            sequence = sequence.sequence
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            'data': {
                'metric_column': 'data',
                'cmap': 'magenta',
                'normalization': 'none',
                'values': None,
                'title': 'Hypthetical pairs',
                'extend': 'neither'}
        } | metric_defaults
        if input_data is not None:
            dataframe=input_data
        else:
            dataframe = {"i": [], "j": [], "data": []}
            for i in range(1, len(sequence)):
                for j in range(i + 1, len(sequence) + 1):
                    dataframe["i"].append(i)
                    dataframe["j"].append(j)
                    dataframe["data"].append(0)
            dataframe = pd.DataFrame(dataframe)
        super().__init__(
            input_data=dataframe,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window)


class StructureInteractions(Interactions):
    def __init__(self, input_data, sequence, structure2=None):
        metric = "Structure"
        metric_defaults = {
            'Structure': {
                'metric_column': 'Structure',
                'cmap': 'grey',
                'normalization': 'none',
                'ticks': [],
                'title': 'Base-pairs',
                'extend': 'neither'}}
        if structure2 is not None:
            input_data = input_data.merge(
                structure2,
                how="left",
                on=["i", "j"],
                indicator="Which_structure",
                suffixes=["_left", "_right"])
            categories = {'both': 0, 'left_only': 1, 'right_only': 2}
            input_data["Which_structure"] = [
                categories[c] for c in input_data['Which_structure']
                ]
            input_data['Which_structure'].astype(int)
            metric = "Which_structure"
            metric_defaults = {
                'Structure_left': {
                    'metric_column': 'Structure_left',
                    'cmap': 'grey',
                    'normalization': 'none',
                    'ticks': [],
                    'title': 'Base-pairs',
                    'extend': 'neither'},
                'Structure_right': {
                    'metric_column': 'Structure_right',
                    'cmap': 'grey',
                    'normalization': 'none',
                    'ticks': [],
                    'title': 'Base-pairs',
                    'extend': 'neither'},
                'Which_structure': {
                    'metric_column': 'Which_structure',
                    'cmap': [
                        (150/255., 150/255., 150/255.), # shared
                        (38/255., 202/255., 145/255.),  # left
                        (153/255., 0.0, 1.0),           # right
                    ],
                    'normalization': 'none',
                    'title': 'Base-pairs by structure',
                    'extend': 'neither'
                }
            }
        super().__init__(input_data, sequence, metric, metric_defaults)
