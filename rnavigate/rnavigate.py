import os.path
from rnavigate import data
from rnavigate.data_loading import create_data

class Sample:
    """
    The main RNAvigate object, representing an RNA structure with experimental
    and/or computational data.
    """

    def __init__(self,
                 sample=None,
                 inherit=None,
                 dance_prefix=None,
                 overwrite_inherited_defaults=False,
                 **data_keywords):
        """Creates a sample object which connects all chemical probing and
        structural data for a single experiment. Contains convenience methods
        to plot, filter, compare and retrieve this data. Every argument is
        optional and defaults to None.

        Args:
            sample (str, optional):
                Label to be used in plot legends.
            inherit (Sample, optional):
                inherit all data from this Sample. Does NOT make copies:
                operations on inherit change self and vice versa. Saves time on
                expensive operations and memory on large data structures.
            dance_prefix (str, optional): "path/to/file_prefix"
                path and prefix for DanceMapper data, automatically locates
                data files and stores DANCE components as a list of Sample
                objects at self.dance
            **data_keywords:
        """
        self.sample = sample
        self.parent = None
        self.data = {}
        self.inputs = {}

        # inherit all data from another Sample
        if hasattr(inherit, "data") and isinstance(inherit.data, dict):
            self.data |= inherit.data
            self.inputs |= inherit.inputs

        # getting ready for default data keyword setup
        default_data = {
            data.Profile: "default_profile",
            data.SecondaryStructure: "default_structure",
            data.PDB: "default_pdb",
        }
        # get data and store inputs for all data_keywords
        for each_keyword in default_data.values():
            if each_keyword not in self.data or overwrite_inherited_defaults:
                self.data[each_keyword] = None
        for data_keyword, kwargs in data_keywords.items():
            self.set_data(data_keyword, kwargs=kwargs)
            for each_type, each_keyword in default_data.items():
                if (isinstance(self.data[data_keyword], each_type)
                        and self.data[each_keyword] is None):
                    self.data[each_keyword] = self.data[data_keyword]

        # for all DANCE-MaP data
        if dance_prefix is not None:
            self.init_dance(dance_prefix)

    @property
    def annotations(self):
        annotations = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.Annotation):
                annotations.append(data_keyword)
        return annotations

    @property
    def profiles(self):
        profiles = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.Profile):
                profiles.append(data_keyword)
        return profiles

    @property
    def structure(self):
        structure = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.SecondaryStructure):
                structure.append(data_keyword)
        return structure

    @property
    def interactions(self):
        interactions = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.Interactions):
                interactions.append(data_keyword)
        return interactions

    @property
    def pdbs(self):
        pdbs = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.PDB):
                pdbs.append(data_keyword)
        return pdbs

    def set_data(self, data_keyword, kwargs, overwrite_keyword=False):
        """This method is used internally to parse data_keywords passed to
        rnavigate.Sample initialization.

        Args:
            data_keyword (str):
                a keyword for the data dictionary
            kwargs (dict):
                a dictionary used to create the rnavigate.data object
            overwrite_keyword (bool, optional):
                whether to overwrite a pre-existing data_keyword.
                Defaults to False.

        Raises:
            ValueError: if the data keyword already exists and
                overwrite_keyword is False
        """
        if (data_keyword in self.inputs) and not overwrite_keyword:
            raise ValueError(
                f"'{data_keyword}' is already a data keyword. "
                "Choose a different one.")
        try:
            self.data[data_keyword] = create_data(
                sample=self, **{data_keyword: kwargs})
        except BaseException as exception:
            raise ValueError(f"issue while loading {data_keyword}") from exception
        self.inputs[data_keyword] = kwargs

    def init_dance(self, filepath):
        """Initializes a list of Sample objects which each represent a
        component of the DANCE model, available as self.dance.

        Args:
            prefix (str): Path to DanceMapper output file prefixes. Finds
                profiles, rings, pairs, structures, and pairing probabilities
                if present.
        """
        reactivityfile = f"{filepath}-reactivities.txt"
        # read in 2 line header
        with open(reactivityfile) as inf:
            header1 = inf.readline().strip().split()
            header2 = inf.readline().strip().split()
        # number of components
        self.dance_components = int(header1[0])
        # population percentage of each component
        self.dance_percents = header2[1:]
        # dance is a list containing one sample for each component
        self.dance = []
        # build column names for reading in BM file
        for i in range(self.dance_components):
            kwargs = {
                "sample": f"{self.sample}: {i} - {self.dance_percents[i]}",
                "dancemap": {"filepath": reactivityfile,
                             "component": i},
                "ringmap": f"{filepath}-{i}-rings.txt",
                "pairmap": f"{filepath}-{i}-pairmap.txt",
                "ss": [f"{filepath}-{i}.f.ct",  # if using --pk
                       f"{filepath}-{i}.ct"],  # if regular fold used
                "pairprob": f"{filepath}-{i}.dp"
            }
            for key in ["ringmap", "pairmap", "pairprob"]:
                if not os.path.isfile(kwargs[key]):
                    kwargs.pop(key)
            for structure_file in kwargs["ss"]:
                if os.path.isfile(structure_file):
                    kwargs["ss"] = structure_file
            if isinstance(kwargs["ss"], list):
                kwargs.pop("ss")
            sample = Sample(**kwargs)
            sample.parent = self
            self.dance.append(sample)

###############################################################################
# filtering and data retrieval
#     get_data
#     filter_interactions
#     filter_dance
###############################################################################

    def get_data(self, data_keyword, data_class=None):
        """Replaces data keyword with data object, even if nested.

        Args:
            data_keyword (str | rnavigate.data | list | dict):
                If a string, returns `self.data[data_keyword]`
                If `None`, returns `None`.
                If an existing rnavigate.data object, returns `data_keyword`
                If a dictionary, returns:
                    `{key: self.get_data(value) for each item in data_keyword}`
                If a list, returns:
                    `[self.get_data(item) for each item in data_keyword]`

        Returns:
            same type as data_keyword: data keywords are replaced with
                matching rnavigate.data objects from self
        """
        # handle special cases
        not_in_sample = ValueError(f"{data_keyword} not in {self.sample}.")
        wrong_class = ValueError(f"{data_keyword} is not {data_class}")
        if data_class is None:
            data_class = data.Sequence
        if isinstance(data_keyword, dict):
            return {k: self.get_data(v, data_class) for k, v in data_keyword.items()}
        elif isinstance(data_keyword, list):
            return [self.get_data(v, data_class) for v in data_keyword]
        elif isinstance(data_keyword, data.Sequence):
            if not isinstance(data_keyword, data_class):
                raise wrong_class
            return data_keyword
        elif data_keyword is None:
            return None
        # if keyword is in sample.data, retreive Sequence object
        try:
            return self.data[data_keyword]
        except KeyError:
            try:
                return self.parent.data[data_keyword]
            except (KeyError, AttributeError) as exception:
                raise not_in_sample from exception

    def filter_interactions(self, interactions, metric=None, cmap=None,
                            normalization=None, values=None, **kwargs):
        """sets coloring properties and filters interactions data.

        Args:
            interactions (rnavigate.data.Interactions | str):
                Interactions object to be filtered. If a string, value is
                replaced with self.get_data(interactions)
            metric (str, optional):
                column of interactions data to be used as metric for coloring
                interactions.
                "Distance" will compute 3D distance in "pdb", defaulting to
                2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
                those atoms to compute distance.
            cmap (str | list, optional):
                sets the interactions colormap, used to color interactions
                according to metric values.
            normalization (str, optional): 
                `'norm'`: extreme values in colormap are given to the extreme
                    values of interactions metric data
                `'bins'`: data are colored according to which bin they fall into
                    `values` defines bins (list, length = 2 less than cmap)
                `'min_max'`: extreme values in cmap are given to values beyond
                    minimum and maximum, defined by `values`
            values:
                behavior depends on normalization
                `'norm'`: values are not needed
                `'bins'`: values should be 2 less than the number of categories
                `'min_max'`: list or tuple containing the minimum and maximum
            **kwargs: Other arguments are passed to interactions.filter()
        """
        # check for valid interactions data
        if interactions is None:
            return
        try:
            interactions = self.data[interactions]
        except KeyError as exception:
            raise KeyError(
                f'{interactions} is not in {self.sample}') from exception
        if not isinstance(interactions, data.Interactions):
            raise ValueError(f'{interactions} is not interactions data')

        if metric is None:
            metric = interactions.default_metric
        elif metric is not None and metric.startswith("Distance"):
            if len(metric.split('_')) == 2:
                metric, atom = metric.split('_')
            else:
                atom = "O2'"
            interactions.set_3d_distances(self.data["pdb"], atom)
        try:
            metric = interactions.metric_defaults[metric]
        except KeyError:
            metric = {'metric_column': metric}
        if cmap is not None:
            metric['cmap'] = cmap
        if normalization is not None:
            metric['normalization'] = normalization
        if values is not None:
            metric['values'] = values
        interactions.metric = metric
        for each_type in ["profile", "structure"]:
            if each_type in kwargs:
                kwargs[each_type] = self.data[kwargs[each_type]]
            elif f"default_{each_type}" in self.data.keys():
                kwargs[each_type] = self.data[f"default_{each_type}"]
        interactions.filter(**kwargs)

    def dance_filter(self, positive_only=True, min_cd=15, Statistic_ge=23,
                     ss_only=True, **kwargs):
        """Applies a standard filter to plot DANCE rings, pairs, and predicted
        structures together.

        Args:
            positive_only (bool, optional): Remove negative correlations.
                Defaults to True.
            min_cd (int, optional): Filters rings by contact distance based
                on predicted structures for *ALL* DANCE components.
                Defaults to 15.
            Statistic_ge (int, optional): Lower bound for MI.
                Defaults to 23.
            ss_only (bool, optional): Filters out rings with at least one leg
                in a double stranded region based on that DANCE component.
                Defaults to True.
            **kwargs: additional arguments are used to filter "ringmap" data.
        """
        kwargs = {}
        if positive_only:
            kwargs["positive_only"] = True
        if ss_only:
            kwargs["ss_only"] = True
        kwargs["Statistic_ge"] = Statistic_ge
        ctlist = [dance.data["default_structure"] for dance in self.dance]
        for dance in self.dance:
            dance_ct = dance.data["default_structure"]
            fit_to = get_sequence(
                sequence=fit_to, sample=dance, default='default_structure')
            dance.data["ringmap"].filter(structure=dance_ct, **kwargs)
            dance.data["ringmap"].mask_on_ct(ctlist, min_cd=min_cd)
            dance.data["pairmap"].filter(structure=dance_ct, paired_only=True)
