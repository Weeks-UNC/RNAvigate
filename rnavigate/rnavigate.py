from rnavigate import data
from rnavigate.data_loading import create_data

class Sample:
    """Loads and organizes RNA structural data for use with plotting functions.

    The Sample class stores all of the relevant experimental and computational
    structural data for a single RNA experiment. Between samples, common data
    types should be given a common data keyword so they can be easily compared.
    """
    def __init__(
            self, sample, inherit=None, keep_inherited_defaults=True,
            **data_keywords
            ):
        """Creates a Sample.

        Required arguments:
            sample (string)
                An arbitrary name. This will be used as a label in plot legends
                and titles to differentiate it from other samples

        Optional arguments:
            inherit (Sample or list of Samples)
                Data keywords and associated data from other samples become the
                data keywords and associated data from this sample. This does
                not make additional copies of the data: i.e. operations that
                make changes to inherited data change the original sample, and
                any other samples that inherited that data. This can be useful
                to save time and memory on operations and large data structures
                that are shared between samples.
            keep_inherited_defaults (True or False)
                whether to keep inherited default keywords
                defaults to True

        Data keywords:
            There are many built-in data keywords with different expectations
            and behaviors. For a full list with expected input formats and
            output behavior, visit:

            https://rnavigate.readthedocs.io/en/latest/loading-data/
        """
        self.sample = sample
        self.inputs = {}
        self.data = {}
        self.defaults = {
                'default_annotation': None,
                'default_profile': None,
                'default_structure': None,
                'default_interactions': None,
                'default_pdb': None
                }

        # inherit data from other Samples
        self.inherit_data(
            inherit=inherit,
            keep_inherited_defaults=keep_inherited_defaults,
            overwrite=False)

        # get data and store inputs for all data_keywords
        for data_keyword, inputs in data_keywords.items():
            self.set_data(data_keyword=data_keyword, inputs=inputs)
            self.set_as_default(data_keyword=data_keyword, overwrite=False)

    def inherit_data(self, inherit, keep_inherited_defaults, overwrite):
        if isinstance(inherit, (list, tuple)):
            for inherit_sample in inherit[::-1]:
                self.inherit_data(
                    inherit_sample, keep_inherited_defaults, overwrite
                    )
        elif isinstance(inherit, Sample):
            if overwrite:
                self.data |= inherit.data
                self.inputs |= inherit.inputs
            else:
                self.data = inherit.data | self.data
                self.inputs = inherit.inputs | self.inputs
            if keep_inherited_defaults and overwrite:
                self.defaults |= inherit.defaults
            elif keep_inherited_defaults:
                self.defaults = inherit.defaults | self.defaults
        elif inherit is not None:
            raise ValueError(
                "inherit only accepts rnav.Sample or list of rnav.Sample"
                )

    def set_data(self, data_keyword, inputs, overwrite=False):
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
            overwrite (bool)
                whether to overwrite a pre-existing data_keyword
                Defaults to False.

        Raises:
            ValueError:
                the data keyword already exists and overwrite is False
            ValueError:
                there was an issue parsing the data
            
        """
        if (data_keyword in self.data) and not overwrite:
            raise ValueError(
                f"'{data_keyword}' is already a data keyword. "
                "Choose a different one.")
        try:
            self.data[data_keyword] = create_data(
                sample=self, data_keyword=data_keyword, inputs=inputs)
        except BaseException as e:
            raise ValueError(f"issue while loading {data_keyword}") from e
        self.inputs[data_keyword] = inputs

    def get_data(self, data_keyword, data_class=None):
        """Replaces data keyword with data object, even if nested.

        Required arguments:
            data_keyword (Data or data keyword or list/dict of such types)
                If None, returns None.
                If a data keyword, returns associated data from sample
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
                if data is not found in sample
            ValueError:
                if the data retrieved is not of the specified data_class
        """
        # handle special cases
        not_in_sample = ValueError(f"{data_keyword} not in {self.sample}.")
        wrong_class = ValueError(f"{data_keyword} is not {data_class}")
        if data_class is None:
            data_class = data.Sequence
        if isinstance(data_keyword, dict):
            return {
                k: self.get_data(v, data_class)
                for k, v in data_keyword.items()
                }
        elif isinstance(data_keyword, list):
            return [self.get_data(v, data_class) for v in data_keyword]
        elif isinstance(data_keyword, data.Sequence):
            if not isinstance(data_keyword, data_class):
                raise wrong_class
            return data_keyword
        elif data_keyword is None:
            return None
        elif (isinstance(data_keyword, str)
                and data_keyword.startswith("default_")):
            data_keyword = self.defaults[data_keyword]
            if data_keyword is None:
                raise not_in_sample
        # if keyword is in sample.data, retreive Sequence object
        try:
            data_obj = self.data[data_keyword]
        except KeyError as exception:
            raise not_in_sample from exception
        if not isinstance(data_obj, data_class):
            raise wrong_class
        return data_obj

    def set_as_default(self, data_keyword, overwrite=True):
        data_object = self.data[data_keyword]
        if isinstance(data_object, data.Annotation):
            default_keyword = "default_annotation"
        elif isinstance(data_object, data.Profile):
            default_keyword = "default_profile"
        elif isinstance(data_object, data.SecondaryStructure):
            default_keyword = "default_structure"
        elif isinstance(data_object, data.Interactions):
            default_keyword = "default_interactions"
        elif isinstance(data_object, data.PDB):
            default_keyword = "default_pdb"
        else:
            return
        if overwrite or self.defaults[default_keyword] is None:
            self.defaults[default_keyword] = data_keyword

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
            interactions = self.get_data(interactions, data.Interactions)
        except KeyError as exception:
            raise KeyError(
                f'{interactions} is not in {self.sample}') from exception

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
        if 'profile' in kwargs:
            kwargs['profile'] = self.data[kwargs['profile']]
        else:
            kwargs['profile'] = self.get_data(f"default_{'profile'}")
        if 'structure' in kwargs:
            kwargs['structure'] = self.data[kwargs['structure']]
        else:
            kwargs['structure'] = self.get_data(f"default_{'structure'}")
        interactions.filter(**kwargs)

    def print_data_keywords(self):
        data_keywords={
            'annotations': [],
            'profiles': [],
            'structures': [],
            'interactions': [],
            'pdbs': []}
        for data_keyword, data_object in self.data.items():
            if data_keyword in self.defaults.values():
                data_keyword += ' (default)'
            if isinstance(data_object, data.Annotation):
                data_keywords['annotations'].append(data_keyword)
            if isinstance(data_object, data.Profile):
                data_keywords['profiles'].append(data_keyword)
            if isinstance(data_object, data.SecondaryStructure):
                data_keywords['structures'].append(data_keyword)
            if isinstance(data_object, data.Interactions):
                data_keywords['interactions'].append(data_keyword)
            if isinstance(data_object, data.PDB):
                data_keywords['pdbs'].append(data_keyword)
        print(f'{self.sample} data keywords:')
        for k, v in data_keywords.items():
            print(f'\n  {k}:')
            print('    '+"\n    ".join(v))
