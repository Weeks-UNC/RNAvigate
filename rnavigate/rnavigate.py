import os.path
from rnavigate import data
from rnavigate.data_loading import create_data, get_sequence

class Sample:
    """Loads and organizes RNA structural data for use with plotting functions.

    The Sample class stores all of the relevant experimental and computational
    structural data for a single RNA experiment. Between samples, common data
    types should be given a common data keyword so they can be easily compared.
    """
    def __init__(
            self, sample, inherit=None, dance_prefix=None,
            overwrite_inherited_defaults=False, **data_keywords):
        """Creates a Sample.

        Required arguments:
            sample (string)
                An arbitrary name. This will be used as a label in plot legends
                and titles to differentiate it from other samples

        Optional arguments:
            inherit (Sample)
                Another Sample from which to inherit stored data. This does not
                make additional copies of the data: i.e. operations that make
                changes to inherited data change the original sample, and any
                other samples that inherited that data. This can be useful to
                save time and memory on operations and large data structures.
            dance_prefix (path, e.g. "path/to/file_prefix")
                path and prefix for DanceMapper data, automatically locates
                data files and stores DANCE components as a list of Sample
                objects at self.dance

        Data keywords:
            There are many built-in data keywords with different expectations
            and behaviors. For a full list with expected input formats and
            output behavior, type:

            rnavigate.standard_data_keywords

            Data keywords are used in the following contexts, I'll use an
            arbitrary data keyword ("keyword") to illustrate:
                rnav.Sample initialization:
                    my_sample = rnav.Sample(
                        sample="my sample name",
                        keyword="expected_input")
                rnav plotting functions:
                    rnav.plot_skyline(
                        samples=[my_sample],
                        profile="keyword")
                data retrieval:
                    my_data = my_sample.get_data("keyword")

            Data keywords can either be a standard keyword or an arbitrary one.
            During rnav.Sample initialization, if an arbitrary keyword is
            desired, a standard data keyword must be used to specify how to
            parse the inputs.

            For the example above, if we wanted to give our pretend standard
            data keyword ("keyword") a different name ("arbitrary"):
                my_sample = rnav.Sample(
                    sample="my sample name",
                    arbitrary={"keyword": "expected_input"})
                rnav.plot_skyline(samples=[my_sample], profile="arbitrary")
                my_data = my_sample.get_data("arbitrary")

            The arbitrary keywords must follow some simple rules:
                1. Cannot already be a data keyword
                2. Cannot consist only of valid nucleotides: AUCGTaucgt
                3. Cannot start with a number: 0123456789
                4. Cannot contain any symbols accept underscore:
                   Not allowed: !@#$%^&*()-+=`~|}{]['";:/?.><
                5. Should avoid "sequence"
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
        for data_keyword, inputs in data_keywords.items():
            self.set_data(data_keyword, inputs=inputs)
            for each_type, each_keyword in default_data.items():
                if (isinstance(self.data[data_keyword], each_type)
                        and self.data[each_keyword] is None):
                    self.data[each_keyword] = self.data[data_keyword]

        # for all DANCE-MaP data
        if dance_prefix is not None:
            self.init_dance(dance_prefix)

    @property
    def annotations(self):
        """Print a list of data keywords associated with annotations"""
        annotations = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.Annotation):
                annotations.append(data_keyword)
        return annotations

    @property
    def profiles(self):
        """list data keywords associated with per-nucleotide data"""
        profiles = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.Profile):
                profiles.append(data_keyword)
        return profiles

    @property
    def structures(self):
        """list data keywords associated with secondary structures"""
        structure = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.SecondaryStructure):
                structure.append(data_keyword)
        return structure

    @property
    def interactions(self):
        """list data keywords associated with inter-nucleotide data"""
        interactions = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.Interactions):
                interactions.append(data_keyword)
        return interactions

    @property
    def pdbs(self):
        """list data keywords associated with 3D molecular structures"""
        pdbs = []
        for data_keyword, data_object in self.data.items():
            if isinstance(data_object, data.PDB):
                pdbs.append(data_keyword)
        return pdbs

    def set_data(self, data_keyword, inputs, overwrite_keyword=False):
        """Add data to Sample using the given data keyword and inputs

        This methods works similarly to the data keywords arguments used
        during Sample initialization:

            my_sample = rnavigate.Sample(
                sample='name',
                data_keyword=inputs)

            is equivalent to:

            my_sample = rnavigate.Sample(
                sample='name')
            my_sample.add_data(
                'data_keyword', inputs)

        Required arguments:
            data_keyword (string)
                a data keyword (arbitrary or standard) used to store and/or
                parse the inputs
            inputs (dictionary or RNAvigate Data)
                a dictionary used to create the data object

        Optional arguments:
            overwrite_keyword (bool)
                whether to overwrite a pre-existing data_keyword
                Defaults to False.

        Raises:
            ValueError:
                the data keyword already exists and overwrite_keyword is False
            ValueError:
                there was an issue parsing the data
            
        """
        if (data_keyword in self.inputs) and not overwrite_keyword:
            raise ValueError(
                f"'{data_keyword}' is already a data keyword. "
                "Choose a different one.")
        try:
            self.data[data_keyword] = create_data(
                sample=self, **{data_keyword: inputs})
        except BaseException as exception:
            raise ValueError(f"issue while loading {data_keyword}") from exception
        self.inputs[data_keyword] = inputs

    def init_dance(self, filepath):
        """Create a list of Samples, one for each component of the DANCE model

        The new list is stored in the .dance attribute.

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

        Required arguments:
            data_keyword (Data or data keyword or list/dict of such types)
                If None, returns None.
                If a data keyword, returns associated data from sample
                    or sample.parent
                If Data, returns that data.
                If a list or dictionary, returns list or dictionary with
                    data keyword values replaced with associated Data
            data_class (RNAvigate Data class)
                If provided, ensures that returned data is of this type.

        Returns:
            Same type as data_keyword argument, but data keywords are replaced
                with associated data

        Raises:
            ValueError:
                if data is not found in sample or sample.parent sample
            ValueError:
                if the data retrieved is not of the specified data_class
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
            data_obj = self.data[data_keyword]
        except KeyError:
            try:
                data_obj = self.parent.data[data_keyword]
                print(f'"{data_keyword}" not in "{self.sample}". '
                      f'Using data from parent sample: "{self.parent.sample}"')
            except (KeyError, AttributeError) as exception:
                raise not_in_sample from exception
        if not isinstance(data_obj, data_class):
            raise wrong_class
        return data_obj

    def filter_interactions(
            self, interactions, metric=None, cmap=None, normalization=None,
            values=None, **kwargs):
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
                `'bins'`: list of floats containing the boundaries between bins
                    One fewer than the number of categories
                `'min_max'`: list of floats containing the minimum and maximum
            **kwargs: Other arguments are passed to interactions.filter()
        """
        # check for valid interactions data
        if interactions is None:
            return
        try:
            interactions = self.get_data(interactions)
        except KeyError as exception:
            raise KeyError(
                f'{interactions} is not in {self.sample}') from exception
        if not isinstance(interactions, data.Interactions):
            raise ValueError(f'{interactions} is not interactions data')

        if metric is None:
            metric = interactions.default_metric
        elif metric.startswith("Distance"):
            if len(metric.split('_')) == 2:
                metric, atom = metric.split('_')
            else:
                atom = "O2'"
            interactions.set_3d_distances(self.get_data("default_pdb"), atom)
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

    def dance_filter(
            self, fit_to='default_structure', positive_only=True, min_cd=15,
            Statistic_ge=23, ss_only=True, **kwargs
            ):
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
        ctlist = [dance.get_data("default_structure") for dance in self.dance]
        for dance in self.dance:
            dance_ct = dance.data["default_structure"]
            fit_to = get_sequence(
                sequence=fit_to, sample=dance, default='default_structure')
            dance.data["ringmap"].filter(structure=dance_ct, **kwargs)
            dance.data["ringmap"].mask_on_structure(ctlist, min_cd=min_cd)
            dance.data["pairmap"].filter(structure=dance_ct, paired_only=True)
