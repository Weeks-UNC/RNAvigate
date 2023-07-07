"""Contains rnavigate.Sample object"""

# external python packages
import os.path

# modules in RNAvigate
from rnavigate import plots
from rnavigate import data
from rnavigate.plotting_functions import get_sequence

_required = object()
data_keyword_defaults = {
    "fasta": {
        "data_class": data.Sequence},
    "log": {
        "data_class": data.Log},
    "shapemap": {
        "data_class": data.SHAPEMaP},
    "dmsmap": {
        "data_class": data.SHAPEMaP,
        "dms": True},
    "dancemap": {
        "data_class": data.DanceMaP,
        "component": _required},
    "rnpmap": {
        "data_class": data.RNPMaP},
    "ringmap": {
        "data_class": data.RINGMaP,
        "sequence": "profile"},
    "pairmap": {
        "data_class": data.PAIRMaP,
        "sequence": "profile"},
    "allcorrs": {
        "data_class": data.RINGMaP,
        "sequence": "profile"},
    "shapejump": {
        "data_class": data.SHAPEJuMP,
        "sequence": _required},
    "pairprob": {
        "data_class": data.PairingProbability,
        "sequence": "profile"},
    "ss": {
        "data_class": data.SecondaryStructure},
    "pdb": {
        "data_class": data.PDB,
        "chain": _required},
    "allpossible": {
        "data_class": data.AllPossible,
        "sequence": _required},
    "motif": {
        "data_class": data.Motif,
        "sequence": _required,
        "color": _required},
    "orfs": {
        "data_class": data.ORFs,
        "sequence": _required,
        "color": _required},
    "spans": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "sites": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "groups": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "primers": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
}

def create_data(sample=None, **data_keyword):
    """Convenience function for creating rnavigate.data objects. This function
    is used to parse **data_keywords passed to rnavigate.Sample, but can also
    be used on it's own using the same syntax.

    Args:
        sample (rnavigate.Sample, optional):
            Requried only if 'sequence' is provided as a data keyword.
                {'sequence': 'data_keyword'}
                    is replaced with:
                {'sequence': sample.data['data_keyword'].sequence}
            Defaults to None.
        **data_keyword:
            There must only one additional argument provided.
            This argument is used to create the rnavigate.data object, it's
            syntax is flexible, see the example below.

    Usage:
        In this example we create a data object from a fasta file using
        three different syntaxes, starting with the most explicit:

            get_data(arbitrary_name={
                        'data_class':'fasta',
                        'input_data': 'my_sequence.fa'})

        This can be simplified by using the value of 'data_class' to replace
        the key 'input_data':

            get_data(arbitrary_name={
                        'fasta': 'my_sequence.fa'})

        Above, an arbitrary_name is used. This allows Data to be assigned to
        arbitrary data keywords of the given Sample. This name can be replaced
        by our data class:

            get_data(fasta='my_sequence.fa')

        In all three cases, the result returned is equivalent to:

            rnavigate.data.Sequence(input_data="my_sequence.fa")
    """
    # Parse data_keyword for name and inputs
    if len(data_keyword) > 1:
        raise ValueError('Only one data_keyword can be provided')

    name = list(data_keyword.keys())[0]
    inputs = data_keyword[name]

    # if given existing rnavigate.data object return
    if isinstance(inputs, data.Data):
        return inputs

    # convert case 1 or 2 into case 3:
    # 1) name='A', inputs='B'
    # 2) inputs={'A': 'B'}
    # 3) inputs={'data_keyword': 'A', 'input_data': 'B'}
    if not isinstance(inputs, dict):
        if name in data_keyword_defaults:
            inputs = {
                'data_keyword': name,
                'input_data': inputs}
    elif all([name in data_keyword_defaults,
              name not in inputs,
              'data_keyword' not in inputs]):
        inputs['data_keyword'] = name
    else:
        for key in list(inputs.keys()):
            if key in data_keyword_defaults:
                inputs['data_keyword'] = key
                inputs['input_data'] = inputs.pop(key)

    # handle special cases
    if inputs['data_keyword'] == 'allpossible':
        inputs['sequence'] = inputs.pop('input_data')
    elif inputs['data_keyword'] in ['sites', 'spans', 'primers', 'groups']:
        inputs['annotation_type'] = inputs['data_keyword']

    # 'filepath' or 'dataframe' may be used in place of 'input_data'
    for alt_input_data in ['filepath', 'dataframe']:
        if alt_input_data in inputs:
            inputs['input_data'] = inputs.pop(alt_input_data)
    if 'fasta' in inputs:
        inputs['sequence'] = inputs.pop('fasta')
    elif 'seq_source' in inputs:
        inputs['sequence'] = inputs.pop('seq_source')

    # retrieve defaults for given data_class
    data_class = inputs.pop('data_keyword')
    try:
        inputs = data_keyword_defaults[data_class] | inputs
        data_class = inputs.pop('data_class')
    except KeyError as exception:
        raise KeyError(
            f"{data_class} is not a valid data keyword.") from exception

    # get sequence from another object if appropriate
    if 'sequence' in inputs:
        inputs["sequence"] = get_sequence(inputs["sequence"], sample)

    # check if any required arguments are not provided
    required_arguments = []
    for key, value in inputs.items():
        if value is _required:
            required_arguments.append(key)
    if len(required_arguments) > 0:
        raise ValueError(
            f"Required arguments for {data_class} were not provided:\n"
            ", ".join(required_arguments))

    # instantiate and return the data class
    try:
        return data_class(**inputs)
    except BaseException as exception:
        raise ValueError(
            f"data_class={data_class}\ninputs={inputs}") from exception


class Sample:
    """
    The main RNAvigate object, representing an RNA structure with experimental
    and/or computational data.
    """

    def __init__(self,
                 sample=None,
                 inherit=None,
                 dance_prefix=None,
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

        # get data and store inputs for all data_keywords
        for data_keyword, kwargs in data_keywords.items():
            self.set_data(data_keyword, kwargs=kwargs)

        # for all DANCE-MaP data
        if dance_prefix is not None:
            self.init_dance(dance_prefix)

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

        self.data[data_keyword] = create_data(
            sample=self, **{data_keyword: kwargs})
        self.inputs[data_keyword] = kwargs
        default_profiles = ["shapemap", "dmsmap", "dancemap", "rnpmap"]
        if (data_keyword in default_profiles) and ("profile" not in self.data):
            self.data["profile"] = self.data[data_keyword]


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
#     get_data_list
#     filter_interactions
#     filter_dance
###############################################################################

    def get_data(self, data_keyword):
        """Returns data based on data_keyword.

        Args:
            data_keyword (str | rnavigate.data | list | dict):
                If a string, returns self.data[data_keyword]
                If an existing rnavigate.data object, returns data_keyword
                If a dictionary, returns:
                    {key: self.get_data(value) for each item in data_keyword}
                If a list, returns:
                    [self.get_data(item) for each item in data_keyword]

        Returns:
            same type as data_keyword: data keywords are replaced with
                matching rnavigate.data objects from self
        """
        # handle special cases
        if isinstance(data_keyword, dict):
            return {k: self.get_data(v) for k, v in data_keyword.items()}
        elif isinstance(data_keyword, list):
            return [self.get_data(v) for v in data_keyword]
        elif isinstance(data_keyword, data.Sequence):
            return data_keyword
        elif data_keyword == "label":
            return self.sample
        elif data_keyword is None:
            return None
        # if keyword is in sample.data, retreive Sequence object
        try:
            return self.data[data_keyword]
        except (KeyError, TypeError):
            try:
                return self.parent.data[data_keyword]
            except (KeyError, AttributeError):
                print(f"{data_keyword} data not found in {self.sample}")
                return None

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
                values between given minimum and maximum, values outside of
                these limits will be given the most extreme color values.
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
        metric = {'metric_column': metric}
        if cmap is not None:
            metric['cmap'] = cmap
        if normalization is not None:
            metric['normalization'] = normalization
        if values is not None:
            metric['values'] = values
        interactions.metric = metric
        for datatype in ["profile", "ct"]:
            if datatype in kwargs:
                kwargs[datatype] = self.data[kwargs[datatype]]
            elif datatype in self.data.keys():
                kwargs[datatype] = self.data[datatype]
        interactions.filter(**kwargs)

    def dance_filter(self, fit_to=None, filterneg=True, cdfilter=15,
                     sigfilter=23, ssfilter=True, **kwargs):
        """Applies a standard filter to plot DANCE rings, pairs, and predicted
        structures together.

        Args:
            fit_to (str or Data object, optional): Data will be fitted to this.
            filterneg (bool, optional): Remove negative correlations.
                Defaults to True.
            cdfilter (int, optional): Filters rings by contact distance based
                on predicted structures for *ALL* DANCE components.
                Defaults to 15.
            sigfilter (int, optional): Lower bound for MI.
                Defaults to 23.
            ssfilter (bool, optional): Filters out rings with at least one leg
                in a double stranded region based on that DANCE component.
                Defaults to True.
            **kwargs: additional arguments are used to filter "ringmap" data.
        """
        kwargs = {}
        if filterneg:
            kwargs["positive_only"] = True
        if ssfilter:
            kwargs["ss_only"] = True
        kwargs["Statistic_ge"] = sigfilter
        ctlist = [dance.data["ct"] for dance in self.dance]
        for dance in self.dance:
            dance_ct = dance.data["ct"]
            fit_to = get_sequence(
                sequence=fit_to, sample=dance, default='ct')
            dance.data["ringmap"].filter(fit_to=fit_to, ct=dance_ct, **kwargs)
            dance.data["ringmap"].mask_on_ct(ctlist, min_cd=cdfilter)
            dance.data["pairmap"].filter(fit_to=fit_to, ct=dance_ct,
                                         paired_only=True)

###############################################################################
# sample plotting functions
#     plot_shapemapper
###############################################################################

    def plot_shapemapper(self, panels=None):
        """Makes a standard ShapeMapper2 profile plot with 3 panels: Normalized
        Reactivities, modified and untreated mutation rates, and modified and
        untreated read depths.

        Args:
            panels (list, optional): Which of the three panels to include.
                Defaults to ["profile", "rates", "depth"].

        Returns:
            Plot object:
        """
        if panels is None:
            panels = ["profile", "rates", "depth"]
        plot = plots.SM(self.data["profile"].length, panels=panels)
        plot.add_sample(self, profile="profile", label="label")
        plot.set_figure_size()
        return plot
