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
            self, sample, inherit=None, keep_inherited_defaults=False,
            **data_keywords
            ):
        """Creates a Sample.

        Required arguments:
            sample (string)
                An arbitrary name. This will be used as a label in plot legends
                and titles to differentiate it from other samples

        Optional arguments:
            inherit (Sample or list of Samples)
                Other Samples from which to inherit stored data. This does not
                make additional copies of the data: i.e. operations that make
                changes to inherited data change the original sample, and any
                other samples that inherited that data. This can be useful to
                save time and memory on operations and large data structures.
            keep_inherited_defaults (True or False)
                whether to keep inherited default keywords

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
        self.inputs = {}
        self.data = {}
        self.defaults = {
                'default_annotations': None,
                'default_profile': None,
                'default_structure': None,
                'default_interactions': None,
                'default_pdb': None
                }

        # inherit data from other Samples
        if isinstance(inherit, Sample):
            inherit = [inherit]
        for inherit_sample in inherit[::-1]:
            if isinstance(inherit, Sample):
                self.data |= inherit_sample.data
                self.inputs |= inherit_sample.inputs
                if keep_inherited_defaults:
                    self.defaults |= inherit_sample.defaults
            else:
                raise ValueError(
                    "inherit only accepts rnav.Sample or list of rnav.Sample"
                    )

        # get data and store inputs for all data_keywords
        for data_keyword, inputs in data_keywords.items():
            self.set_data(data_keyword=data_keyword, inputs=inputs)
            self.set_as_default(data_keyword=data_keyword, overwrite=False)

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

    def print_data_keywords(self):
        data_keywords={
            'annotations': [],
            'profiles': [],
            'structures': [],
            'interactions': [],
            'pdbs': [],
            'other': [],
            'missing': []}
        for data_keyword, data_object in self.data.items():
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
            if data_object is None:
                data_keywords['missing'].append(data_keyword)
            else:
                data_keywords['other'].append(data_keyword)
        print(f"{self.sample} Data keywords:\n"
              "  Annotations:\n"
              f"    {'   '.join(data_keywords['annotations'])}\n"
              "  Profiles:\n"
              f"    {'   '.join(data_keywords['profiles'])}\n"
              "  Structures:\n"
              f"    {'   '.join(data_keywords['structures'])}\n"
              "  Interactions:\n"
              f"    {'   '.join(data_keywords['interactions'])}\n"
              "  PDBs:\n"
              f"    {'   '.join(data_keywords['pdbs'])}\n"
              "  Other:\n"
              f"    {'   '.join(data_keywords['other'])}\n"
              "  Missing:"
              f"    {'   '.join(data_keywords['missing'])}")
