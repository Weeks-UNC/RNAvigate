#!/usr/bin/env python

# general python packages
import os.path
import numpy as np
import json

# modules in RNAvigate
from .data import get_ss_class
from . import plots
from . import data
from . import styles


def create_code_button():
    """Not functioning. Meant to create a button for HTML exported notebooks
    which toggles hiding or showing code cells for prettier reports."""
    # TODO: Figure out why this doesn't work or find a new solution.
    from IPython.display import display, HTML
    display(HTML('''<script>
                 code_show=true;
                 function code_toggle() {
                 if (code_show) {$('div.input').hide();}
                 else {$('div.input').show();}
                 code_show = !code_show
                 }
                 $( document ).ready(code_toggle);
                 </script>
                 <form action="javascript:code_toggle()">
                 <input type="submit" value="Hide/show raw code.">
                 </form>'''))


def get_color_list(length, default, color_regions):
    """Creates a list of colors for use with colors arguments

    Args:
        length (int): length of sequence color list will be used on.
        default (matplotlib color): default color
        color_regions (dict): dictionary: keys are matplot lib colors and
            values are lists of 1-indexed pairs of start and end positions.

    Returns:
        list: a list of colors
            e.g. get_color_list(
                length=10,
                default="grey",
                color_regions={
                    "red": [[2, 4]],
                    "black": [[5, 6], [8, 10]]
                }
            )
            returns:
                ["grey", "red", "red", "red", "black",
                 "black", "grey", "black", "black", "black"]
    """
    # TODO: somehow implement this using annotations?
    c_list = np.full(length, default, dtype='<U16')
    for color, regions in color_regions.items():
        for r in regions:
            c_list[r[0]-1: r[1]] = color
    return list(c_list)


class Sample():
    """
    The main RNAvigate object, representing an RNA probing experiment.
    """

    def __init__(self,
                 sample=None,
                 inherit=None,
                 pdb=None,
                 ct=None,
                 compct=None,
                 ss=None,
                 fasta=None,
                 log=None,
                 shapemap=None,
                 dmsmap=None,
                 dancemap=None,
                 rnpmap=None,
                 ringmap=None,
                 shapejump=None,
                 pairmap=None,
                 allcorrs=None,
                 pairprob=None,
                 allpossible=None,
                 dance_prefix=None,
                 annotations=None):
        """Creates a sample object which connects all chemical probing and
        structural data for a single experiment. Contains convenience methods
        to plot, filter, compare and retrieve this data. Every argument is
        optional and defaults to None.

        Args:
            sample (str, optional): Label to be used in plot legends.
            inherit (Sample, optional): inherit all data from this Sample,
                Does NOT make copies: operations on inherit change self and
                vice versa. Saves time on expensive operations and memory on
                large data structures.
            pdb (dict, optional): dictionary containing the following:
                "filepath": "path/to/file", can be .cif or .pdb
                "chain": "A", corresponding to the chain ID within PDB/CIF file
                "fasta": "path/to/file.fa", required if CIF provided, or if
                    PDB file lacks the REFSEQ entry in the header.
            ct (str, optional): path/to/file containing a secondary structure
                .ct, .cte, .nsd, .varna, .xrna, or .json (R2DT)
            compct (str, optional): same as ct above
            ss (str, optional): same as ct above, but should contain drawing
                coordinates.
            fasta (str, optional): path/to/file containing a single sequence
                in fasta format, can be used as a seq_source for various
                arguments below
            log (str, optional): path/to/file containing shapemapper log info
                typical suffix is _shapemapper_log.txt.
            shapemap (str, optional): path/to/file containing shapemapper2
                reactivities, typically the _reactivities.txt file
            dmsmap (str, optional): path/to/file, same as above, but when data
                is parsed, the Norm_profile column is recomputed based on DMS
                conventions
            dancemap (dict, optional): a dictionary containing the following:
                "filepath": "path/to/reactivities.txt" output from dancemapper
                "component": int, the dancemapper component to load
            rnpmap (str, optional): path/to/file, RNPMapper output file
            ringmap (str, optional): path/to/file, RingMapper output file,
                typically _rings.txt
            shapejump (dict, optional): a dictionary containing:
                "filepath": "path/to/file", a ShapeJumper data file
                and one of the following:
                "fasta": "path/to/file", a fasta containing a single sequence
                "sequence": "AUCGCUGA...", a sequence string.
            pairmap (str, optional): "path/to/file"
                PairMapper _pairmap.txt file
            allcorrs (str, optional): "path/to/file"
                PairMapper _allcorrs.txt file
            pairprob (str, optional): "path/to/file"
                ProbabilityPlot text file
            allpossible (str, optional): str can be one of the following:
                data argument used above, retrieves sequence from that data
                a sequence string e.g. "AUCGUCAUGAUGCA"
            dance_prefix (str, optional): "path/to/file_prefix"
                path and prefix for DanceMapper data, automatically locates
                data files and stores DANCE components as a list of Sample
                objects at self.dance
            annotations (dict, optional): A dictionary with the following keys:
                "seq_source": one of the arguments above, in which case the
                    sequence from that data is retrived, or a sequence string
                    e.g. "AUGCGUAC"
                any number of data keyword strings (to access this annotation):
                    a dictionary defining the annotation positions, type, and
                    colors used in plotting functions. Must contain "color"
                    (described below), and an annotations type ("sites",
                    "spans", "primers", "groups", "motif", "orf")
                    (described below). All positions are 1-indexed. Start and
                    end positions are inclusive.
                    "color": a matplotlib color, search "matplotlib specifying
                        color". This color will be used to represent the
                        annotation positions on plots. "groups" includes
                        "colors" in the annotations list and does not
                        require this.
                    "sites": a list of single-nucleotide positions
                    "spans": a list of lists of start and stop positions
                    "primers": similar to spans above, but reverse primer
                        should be in reverse order. e.g. [[1, 20], [300, 280]]
                    "groups": a list of dictionaries containing "sites" and
                        "color" as described above
                    "motif": a string representing a sequence motif of interest
                        using conventional alphabet, e.g. "DRACH"
                        This will create span annotations for every motif match
                    "orfs": True. This will create span annotations for every
                        possible open reading frame
        """
        if dancemap is None:
            dancemap = {"filepath": None}
        if annotations is None:
            annotations = {}
        if shapejump is None:
            shapejump = {"filepath": None}
        if pdb is None:
            pdb = {"filepath": None}

        # for each input
        # [0] filepath
        # [1] instantiation class
        # [2] sequence source
        # [3] kwargs
        self.inputs = {
            "fasta": {
                "filepath": fasta,
                "instantiator": data.Data,
                "seq_source": "self",
                "kwargs": {}},
            "log": {
                "filepath": log,
                "instantiator": data.Log,
                "seq_source": "self",
                "kwargs": {}},
            "shapemap": {
                "filepath": shapemap,
                "instantiator": data.SHAPEMaP,
                "seq_source": "self",
                "kwargs": {}},
            "dmsmap": {
                "filepath": dmsmap,
                "instantiator": data.SHAPEMaP,
                "seq_source": "self",
                "kwargs": {"dms": True}},
            "dancemap": {
                "filepath": dancemap.pop("filepath"),
                "instantiator": data.DanceMaP,
                "seq_source": "self",
                "kwargs": dancemap},
            "rnpmap": {
                "filepath": rnpmap,
                "instantiator": data.RNPMaP,
                "seq_source": "self",
                "kwargs": {}},
            "ringmap": {
                "filepath": ringmap,
                "instantiator": data.RINGMaP,
                "seq_source": "profile",
                "kwargs": {}},
            "pairmap": {
                "filepath": pairmap,
                "instantiator": data.PAIRMaP,
                "seq_source": "profile",
                "kwargs": {}},
            "allcorrs": {
                "filepath": allcorrs,
                "instantiator": data.RINGMaP,
                "seq_source": "profile",
                "kwargs": {}},
            "shapejump": {
                "filepath": shapejump.pop("filepath", None),
                "instantiator": data.SHAPEJuMP,
                "seq_source": "self",
                "kwargs": shapejump},
            "pairprob": {
                "filepath": pairprob,
                "instantiator": data.PairProb,
                "seq_source": "profile",
                "kwargs": {}},
            "ct": {
                "filepath": ct,
                "instantiator": get_ss_class,
                "seq_source": "self",
                "kwargs": {}},
            "compct": {
                "filepath": compct,
                "instantiator": get_ss_class,
                "seq_source": "self",
                "kwargs": {}},
            "ss": {
                "filepath": ss,
                "instantiator": get_ss_class,
                "seq_source": "self",
                "kwargs": {}},
            "pdb": {
                "filepath": pdb.pop("filepath", None),
                "instantiator": data.PDB,
                "seq_source": "self",
                "kwargs": pdb},
            "dance_prefix": {
                "filepath": dance_prefix,
                "instantiator": self.init_dance,
                "seq_source": "self",
                "kwargs": {}},
            "allpossible": {
                "filepath": "",
                "instantiator": data.AllPossible,
                "seq_source": allpossible,
                "kwargs": {}},
        }

        # add annotations to the above inputs
        if isinstance(annotations, str):
            with open(annotations) as f:
                annotations = json.load(f)
        annotation_seq = annotations.pop("seq_source", None)
        self.annotations = list(annotations)
        for name, kwargs in annotations.items():
            if name in self.inputs:
                print(f"Choose a different name: '{name}' is already used.")
                continue
            if "motif" in kwargs:
                instantiator = data.Motif
            elif "orfs" in kwargs:
                instantiator = data.ORFs
            else:
                instantiator = data.Annotation
            self.inputs[name] = {
                "filepath": "",
                "instantiator": instantiator,
                "seq_source": annotation_seq,
                "kwargs": kwargs}

        self.default_profiles = ["shapemap", "dmsmap", "dancemap", "rnpmap"]

        self.sample = sample
        self.parent = None
        self.data = {}  # stores all data objects
        # inherit all data objects from another Sample
        # These data objects are NOT copies: any operation on inherit sample
        # changes this sample and vice versa. Saves time and memory usage.
        if hasattr(inherit, "data") and isinstance(inherit.data, dict):
            self.data |= inherit.data
        # load data
        for name, kwargs in self.inputs.items():
            if None not in [kwargs["filepath"], kwargs["seq_source"]]:
                kwargs.update(kwargs.pop("kwargs", {}))
                self.set_data(name=name, **kwargs)

    def set_data(self, name, filepath=None, instantiator=None, seq_source=None,
                 **kwargs):
        """Convenience function to add data to Sample after initialization.
        Also, sets data to default profile if appropriate and none exists yet.
        Otherwise, equivalent to:
            self.data[name] = instatiator(filepath=filepath,
                sequence=self.get_sequence(seq_source), **kwargs)
        ...or if filename is an existing Data subclass object:
            self.data[name] = filename

        Args:
            name (str): name of data
            filepath (str, optional): "path/to/file"
                Defaults to None.
            instantiator (instantiator of rnavigate.Data, optional):
                Data class instantiator used to create data object.
                Defaults to None.
            seq_source (str, optional): a key of self.data or a sequence string
                Defaults to None.
        """
        # if given previously instantiated data class, attach and end function
        if isinstance(filepath, data.Data):
            self.data[name] = filepath
            if name in self.default_profiles and "profile" not in self.data:
                self.data["profile"] = self.data[name]
            return
        # check for default instantiator
        if instantiator is None:
            instantiator = self.inputs[name]["instantiator"]
        # check for default seq_source
        if seq_source is None:
            seq_source = self.inputs[name]["seq_source"]
        # get sequence in case provided instead of passed to instantiator
        if seq_source != "self":
            sequence = get_sequence(seq_source=seq_source, sample=self)
            sequence = sequence.sequence
            kwargs.update({"sequence": sequence})
        # for lists of files, set_data() for each with a number indicator
        if isinstance(filepath, list):
            for i, path in enumerate(filepath):
                self.set_data(name=f"{name}_{i+1}",
                              filepath=path,
                              instantiator=instantiator,
                              seq_source=seq_source,
                              **kwargs)
        # update inputs
        self.inputs[name] = {
            "filepath": filepath,
            "instantiator": instantiator,
            "seq_source": seq_source,
            "kwargs": kwargs}
        # instantiate the data class and add to self.data dictionary
        self.data[name] = instantiator(filepath=filepath, **kwargs)
        # set as default profile if appropriate
        if name in self.default_profiles and "profile" not in self.data:
            self.data["profile"] = self.data[name]

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
                "ct": [f"{filepath}-{i}.f.ct",  # if using --pk
                       f"{filepath}-{i}.ct"],  # if regular fold used
                "pairprob": f"{filepath}-{i}.dp"
            }
            for key in ["ringmap", "pairmap", "pairprob"]:
                if not os.path.isfile(kwargs[key]):
                    kwargs.pop(key)
            for ct in kwargs["ct"]:
                if os.path.isfile(ct):
                    kwargs["ct"] = ct
            if isinstance(kwargs["ct"], list):
                kwargs.pop("ct")
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

    def get_data(self, key):
        """Internal function that returns a data object associated with key.
        Prevents an error and prints a message if data object not found. It
        will also search the parent of a DANCE model component for structure
        information.

        Args:
            key (str): Descriptor of data type desired. Can be any of the keys
               of Sample.data or "label" or a Data object

        Returns:
            Data object
        """
        if key == "label":
            return self.sample
        elif key is None:
            return None
        if isinstance(key, data.Data):
            return key
        try:
            return self.data[key]
        except (KeyError, TypeError):
            try:
                return self.parent.data[key]
            except (KeyError, AttributeError):
                print(f"{key} data not found in {self.sample}")
                return None

    def get_data_list(self, keys):
        """Calls Sample.get_data() on a list of keys, also accepts "ctcompare",
        which is equivalent to ["ct", "compct"]

        Returns:
            list of data objects
        """
        if isinstance(keys, list):
            return [self.get_data_list(key) for key in keys]
        elif keys == "ctcompare":
            return self.get_data_list(["ct", "compct"])
        else:
            return self.get_data(keys)

    def filter_interactions(self, interactions, fit_to,
                            suppress=False, metric=None, cmap=None,
                            min_max=None, prefiltered=False, **kwargs):
        """Aligns sequence to fit_to, sets properties, filters and aligns data.

        For example, for plotting interaction data containing structure
        cassettes against a CT without, or plotting sequence variants together
        using a best pairwise alignment. If kwargs contains metric, cmap, or
        min_max, these properties of the Interactions object are set. Other
        kwargs are passed to self.data[interactions].filter(). See
        rnavigate.Interactions.filter for more detail.

        Args:
            interactions (Interactions): Interactions object to be filtered
            fit_to (str or Data subclass object): Data object containing a
                sequence for Interactions to be fitted to for plotting purposes
                can either be a Data object or a key of self.data
            metric (str, optional): column of interactions data to be used as
                metric for coloring interactions, "Distance" will compute 3D
                distance in "pdb", defaulting to 2'OH atom. "Distance_DMS" or
                "Distance_[atom id]" will use those atoms to compute distance.
            cmap (str or list, optional): sets the interactions colormap, used
                to color interactions according to metric values. See
                rnavigate.Interactions.cmap for more detailed behavior.
            min_max (list, length 2, optional): cmap will apply to metric
                values between given minimum and maximum, values outside of
                these limits will be given the most extreme color values.
            **kwargs: Other arguments are passed to interactions.filter()
        """
        # check for valid interactions data
        if (interactions is None) or prefiltered:
            return
        elif interactions not in self.data.keys():
            if not suppress:
                print(f"{interactions} not found in sample data")
            return

        interactions = self.data[interactions]
        if not isinstance(interactions, data.Interactions):
            if not suppress:
                print(f"{interactions} is not an Interactions datatype")
            return

        if metric is not None:
            if metric.startswith("Distance"):
                metric = (metric, self.data["pdb"])
            interactions.metric = metric
        else:
            interactions.metric = interactions.default_metric
        if cmap is not None:
            interactions.cmap = cmap
        if min_max is not None:
            interactions.min_max = min_max
        for datatype in ["profile", "ct"]:
            if datatype in kwargs.keys():
                kwargs[datatype] = self.data[kwargs[datatype]]
            elif datatype in self.data.keys():
                kwargs[datatype] = self.data[datatype]
        if not hasattr(fit_to, "sequence"):
            fit_to = self.get_data_list(fit_to)
        interactions.filter(fit_to, **kwargs)

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
                seq_source=fit_to, sample=dance, default='ct')
            dance.data["ringmap"].filter(fit_to=fit_to, ct=dance_ct, **kwargs)
            dance.data["ringmap"].mask_on_ct(ctlist, min_cd=cdfilter)
            dance.data["pairmap"].filter(fit_to=fit_to, ct=dance_ct,
                                         paired_only=True)

###############################################################################
# sample plotting functions
#     plot_shapemapper
###############################################################################

    def plot_shapemapper(self, panels=["profile", "rates", "depth"]):
        """Makes a standard ShapeMapper2 profile plot with 3 panels: Normalized
        Reactivities, modified and untreated mutation rates, and modified and
        untreated read depths.

        Args:
            panels (list, optional): Which of the three panels to include.
                Defaults to ["profile", "rates", "depth"].

        Returns:
            Plot object:
        """
        plot = plots.SM(self.data["profile"].length, panels=panels)
        plot.add_sample(self, profile="profile", label="label")
        plot.set_figure_size()
        return plot

###############################################################################
# accessory functions
#   get_sequence
#   fit_data_list
###############################################################################


def get_sequence(seq_source, sample=None, default=None):
    """Flexible function that returns a Data object containing a sequence.

    Args:
        seq_source (any): Can be a sequence string, a key in sample.data, or a
            Data object with a sequence
        sample (Sample, optional): sample to get retrieve data from.
            Defaults to None.
        default (any, optional): same format as seq_source above. default to
            this if seq_source is None.
            Defaults to None.

    Raises:
        ValueError: If seq_source and default are both None
        ValueError: If Data object could not be retreived based on inputs

    Returns:
        Data object or subclass: contains a sequence and can be used in factory
            methods
    """
    if seq_source is None and default is None:
        raise ValueError("A seq_source must be provided.")
    elif seq_source is None:
        seq_source = default
    if seq_source in sample.data.keys():
        sequence = sample.data[seq_source]
    elif hasattr(seq_source, "sequence"):
        sequence = seq_source
    elif all([nt.upper() in "AUCGT." for nt in seq_source]):
        sequence = data.Data(sequence=seq_source)
    else:
        raise ValueError(f'Cannot find sequence from {seq_source}')
    return sequence


def fit_data_list(sample, data_list, fit_to):
    """Given a sample and list of sample.data keys, Data objects are mapped to
    fit_to

    Args:
        sample (rnavigate.Sample): sample to retrieve data from
        data_list (list): list of sample.data keys or None
        fit_to (rnavigate.Data): Data object with a sequence to fit to
    """
    for data_obj in data_list:
        if data_obj in sample.data.keys():
            sample.data[data_obj].fit_to(fit_to)
        elif isinstance(data_obj, data.Data):
            data_obj.fit_to(fit_to)


###############################################################################
# Plotting functions that accept a list of samples
#   plot_qc_multisample
#   plot_skyline_multisample
#   plot_arcs_multisample
#   plot_ss_multisample
#   plot_mol_multisample
#   plot_heatmap_multisample
#   plot_circle_multisample
#   plot_linreg_multisample
#   plot_roc_multisample
#   plot_disthist_multisample
#   plot_qc
#   plot_skyline
#   plot_arcs
#   plot_ss
#   plot_mol
#   plot_heatmap
#   plot_circle
#   plot_linreg
#   plot_roc
#   plot_disthist
###############################################################################


def plot_qc_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_qc() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_qc starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_qc(*args, **kwargs)


def plot_skyline_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_skyline() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_skyline starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_skyline(*args, **kwargs)


def plot_arcs_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_arcs() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_arcs starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_arcs(*args, **kwargs)


def plot_ss_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_ss() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_ss starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_ss(*args, **kwargs)


def plot_mol_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_mol() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_mol starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_mol(*args, **kwargs)


def plot_heatmap_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_heatmap() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_heatmap starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_heatmap(*args, **kwargs)


def plot_circle_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_circle() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_circle starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_circle(*args, **kwargs)


def plot_linreg_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_linreg() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_linreg starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_linreg(*args, **kwargs)


def plot_roc_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_roc() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_roc starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_roc(*args, **kwargs)


def plot_disthist_multisample(*args, **kwargs):
    """Deprecation Warning: this function name will be removed in favor of
    rnavigate.plot_disthist() starting with version 1.0.0
    """
    print("Deprecation warning, this function will be removed in favor of "
          "rnavigate.plot_disthist starting with version 1.0.0. This function "
          "simply calls the new function, please start using that one.")
    return plot_disthist(*args, **kwargs)


def plot_qc(samples, labels=None, plot_kwargs=None, **kwargs):
    """Makes a multipanel quality control plot displaying mutations per
    molecule, read length distribution, and mutation rate distributions for
    modified and unmodified samples.

    Args:
        samples (list of rnavigate.Sample): samples to plot.
        labels (list of str, optional): same length as samples list. labels to
            to be used on plot legends.
            Defaults to sample attribute of each sample.
        plot_kwargs (dict, optional): kwargs dictionary passed to QC().
            Defaults to {}.
        **kwargs: passed to QC.plot_data

    Returns:
        rnavigate.plots.QC plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    if labels is None:
        labels = ["label"]*len(samples)
    if plot_kwargs is None:
        plot_kwargs = {}
    plot = plots.QC(num_samples=len(samples), **plot_kwargs)
    for sample, label in zip(samples, labels):
        plot.add_sample(sample=sample, log="log", profile="profile",
                        label=label, **kwargs)
    plot.set_figure_size()
    return plot


def plot_skyline(samples, seq_source=None, profile="profile", labels=None,
                 annotations=None, region="all", plot_kwargs=None, **kwargs):
    """Plots multiple per-nucleotide datasets on a single axis.

    Args:
        samples (list of rnavigate.Sample): samples to plot.
        seq_source (str or data object, optional): a key from sample.data, a
            sequence string, or a Data object. All data will be mapped to this
            string using a user-defined or pairwise sequence alignment.
            Defaults to the value of the profile argument below.
        profile (str, optional): per-nucleotide data to retrieve from sample.
            Defaults to "profile".
        labels (list of str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        annotations (list of str, optional): annotations to retrive from
            sample. Defaults to [].
        region (list of int: length 2, optional): start and end positions to
            plot. 1-indexed, inclusive. Defaults to [1, sequence length].
        plot_kwargs (dict, optional): kwargs dictionary passed to Skyline().
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Skyline.plot_data.
            see rnavigate.plots.Skyline.plot_data for more detail.

    Returns:
        rnavigate.plots.Skyline plot: object containing matplotlib figure and
            axes with additional plotting and file saving methods
    """
    if annotations is None:
        annotations = []
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = ["label"]*len(samples)
    sequence = get_sequence(seq_source, samples[0], profile)
    plot = plots.Skyline(num_samples=len(samples),
                         nt_length=sequence.length,
                         region=region,
                         **plot_kwargs)
    for sample, label in zip(samples, labels):
        fit_data_list(sample, annotations + [profile], sequence)
        plot.add_sample(sample, profile=profile, annotations=annotations,
                        label=label, **kwargs)
    plot.set_figure_size()
    return plot


def plot_profile(samples, seq_source=None, profile="profile", labels=None,
                 annotations=None, region="all", plot_kwargs=None, **kwargs):
    """Aligns reactivity profiles by sequence and plots them on seperate axes.

    Args:
        samples (list of rnavigate.Sample): samples to plot.
        seq_source (str or data object, optional): a key from sample.data, a
            sequence string, or a Data object. All data will be mapped to this
            string using a user-defined or pairwise sequence alignment.
            Defaults to the value of the profile argument below.
        profile (str, optional): per-nucleotide data to retrieve from sample.
            Defaults to "profile".
        labels (list of str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        annotations (list of str, optional): annotations to retrive from
            sample. Defaults to [].
        region (list of int: length 2, optional): start and end positions to
            plot. 1-indexed, inclusive. Defaults to [1, sequence length].
        plot_kwargs (dict, optional): kwargs dictionary passed to Profile().
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Profile.plot_data.
            see rnavigate.plots.skyline.Profile.plot_data for more detail.

    Returns:
        rnavigate.plots.Profile plot: object containing matplotlib figure and
            axes with additional plotting and file saving methods
    """
    if annotations is None:
        annotations = []
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = ["label"]*len(samples)
    sequence = get_sequence(seq_source, samples[0], profile)
    plot = plots.skyline.Profile(num_samples=len(samples),
                                 nt_length=sequence.length,
                                 region=region,
                                 **plot_kwargs)
    for sample, label in zip(samples, labels):
        fit_data_list(sample, annotations + [profile], sequence)
        plot.add_sample(sample, seq=sequence, profile=profile,
                        annotations=annotations, label=label, **kwargs)
    plot.set_figure_size()
    return plot


def plot_alignment(data1, data2, labels=None, plot_kwargs=None, **kwargs):
    """Plots the sequence alignment used to compare two data objects.

    Args:
        data1 (tuple (rnavigate.Sample, str)): a sample and data class keyword
        data2 (tuple (rnavigate.Sample, str)): a sample and data class keyword
        labels (list of str, optional): length 2. Labels used in plots
            Defaults to default sample name + data class keyword
        plot_kwargs (dict, optional): kwargs dictionary passed to Alignment().
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to
            Alignment.plot_data. see rnavigate.plots.Alignment.plot_data for
            more detail.

    Returns:
        rnavigate.plots.Alignment plot: object containing matplotlib figure and
            axes with additional plotting and file saving methods
    """
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = [f"{s.sample}: {seq}" for s, seq in [data1, data2]]
    plot = plots.Alignment(num_samples=1, **plot_kwargs)
    plot.add_sample(sample=None,
                    data1=data1[0].data[data1[1]],
                    data2=data2[0].data[data2[1]],
                    label=labels, **kwargs)
    plot.set_figure_size()
    return plot


def plot_arcs(samples, seq_source=None, ct="ct", comp=None, interactions=None,
              interactions_filter=None, interactions2=None,
              interactions2_filter=None, filters=None, profile="profile",
              annotations=[], labels=None, region="all", plot_kwargs=None,
              colorbar=True, **kwargs):
    """Generates a multipanel arc plot displaying combinations of secondary
    structures, per-nucleotide data, inter-nucleotide data, and sequence
    annotations. Each plot may display a unique sample and/or filtering scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        seq_source (str or data object, optional): a key from sample.data, a
            sequence string, or a Data object. All data will be mapped to this
            string using a user-defined or pairwise sequence alignment.
            Defaults to the value of the ct argument below.
        ct (str, optional): a key from sample.data to retreive a secondary
            structure. This will be plotted on the top half of each panel.
            Defaults to "ct".
        comp (str, optional): same as ct. basepairs from ct and comp will be
            plotted on the top half of each panel. Basepairs are colored by
            which structure contains them (shared, ct only, comp only).
            Defaults to None.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to seq_source coordinates,
            filtered using interactions_filter arguments, and displayed on the
            bottom half of each panel.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=seq_source. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        interactions2_filter (dict, optional): same as interactions_filter
            above but applied to interactions2.
            Defaults to {}.
        filters (list of dict, optional): For plotting multiple filtering
            schemes applied to each sample, each in a new panel. Each
            dictionary in this list follows a similar structure as
            interactions_filter, but also requires the key-value pair:
            {"interactions": interactions} where interactions is as described
            above. interactions and interactions_filter arguments above will be
            ignored. interactions2 and interactions2_filter will be displayed
            in every panel.
            Defaults to [].
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are displayed in center of each panel.
            Defaults to "profile".
        annotations (list, optional): a list of keys from sample.data used to
            retrieve sequence annotations. These annotations are displayed in
            the center of each panel.
            Defaults to [].
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        region (list of int: length 2, optional): start and end position of
            seq_source to be plotted. 1-indexed, inclusive.
            Defaults to [0, sequence length].
        plot_kwargs (dict, optional): kwargs passed to AP(). See
            rnavigate.plots.AP for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to AP.plot_data.
            see rnavigate.plots.AP.plot_data for more detail.

    Returns:
        rnavigate.plots.AP plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if interactions2_filter is None:
        interactions2_filter = {}
    if labels is None:
        labels = ["label"]*len(samples)
    if plot_kwargs is None:
        plot_kwargs = {}
    # if filters list given, rows = # samples, columns = # filters
    if ((filters is not None)
            and (len(samples) > 1)
            and ("rows" not in plot_kwargs)
            and ("cols" not in plot_kwargs)):
        plot_kwargs["rows"] = len(samples)
        plot_kwargs["cols"] = len(filters)
    # coerce interactions and interactions_filter into filters format
    elif filters is None:
        filters = [{"interactions": interactions} | interactions_filter]
    # initialize plot using all structure drawings
    num_samples = len(samples) * len(filters)
    seq = get_sequence(seq_source=seq_source, sample=samples[0], default=ct)
    plot = plots.AP(num_samples=num_samples, nt_length=seq.length,
                    region=region, **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    for sample, label in zip(samples, labels):
        fit_data_list(sample, annotations + [ct, comp, profile], seq)
        sample.filter_interactions(interactions=interactions2,
                                   fit_to=seq, **interactions2_filter)
        for filt in filters:
            sample.filter_interactions(fit_to=seq, **filt)
            plot.add_sample(sample=sample, seq=seq, ct=ct, comp=comp,
                            interactions=filt["interactions"],
                            interactions2=interactions2, profile=profile,
                            label=label, annotations=annotations,
                            **kwargs)
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_arcs_compare(samples, seq_source=None, ct="ct", comp=None,
                      interactions=None, interactions_filter=None,
                      interactions2=None, interactions2_filter=None,
                      profile="profile", labels=None, region="all",
                      plot_kwargs=None, colorbar=True, **kwargs):
    """Generates a single arc plot displaying combinations of secondary
    structures, per-nucleotide data, inter-nucleotide data, and sequence
    annotations. The first sample will be on top, the second on the bottom.
    Center shows how these sequences are being aligned.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        seq_source (str or data object, optional): a key from sample.data, a
            sequence string, or a Data object. All data will be mapped to this
            string using a user-defined or pairwise sequence alignment.
            Defaults to the value of the ct argument below.
        ct (str, optional): a key from sample.data to retreive a secondary
            structure. This will be plotted on the top half of each panel.
            Defaults to "ct".
        comp (str, optional): same as ct. basepairs from ct and comp will be
            plotted on the top half of each panel. Basepairs are colored by
            which structure contains them (shared, ct only, comp only).
            Defaults to None.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to seq_source coordinates,
            filtered using interactions_filter arguments, and displayed on the
            bottom half of each panel.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=seq_source. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        interactions2_filter (dict, optional): same as interactions_filter
            above but applied to interactions2.
            Defaults to {}.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are displayed in center of each panel.
            Defaults to "profile".
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        region (list of int: length 2, optional): start and end position of
            seq_source to be plotted. 1-indexed, inclusive.
            Defaults to [0, sequence length].
        plot_kwargs (dict, optional): kwargs passed to AP(). See
            rnavigate.plots.AP for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to AP.plot_data.
            see rnavigate.plots.AP.plot_data for more detail.

    Returns:
        rnavigate.plots.AP plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    seq1 = get_sequence(seq_source=seq_source, sample=samples[0], default=ct)
    seq2 = get_sequence(seq_source=seq_source, sample=samples[1], default=ct)
    _, _, al1, al2 = seq1.get_alignment_map(fit_to=seq2, return_alignment=True)
    seq1_full = al1.replace('-', '.')
    seq2_full = al2.replace('-', '.')
    seq1_full = get_sequence(seq_source=seq1_full,
                             sample=samples[0], default=ct)
    seq2_full = get_sequence(seq_source=seq2_full,
                             sample=samples[0], default=ct)
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if interactions2_filter is None:
        interactions2_filter = {}
    if labels is None:
        labels = ["label"]*len(samples)
    if plot_kwargs is None:
        plot_kwargs = {}
    if region != "all":
        al1 = [i for i, nt in enumerate(seq1_full.sequence) if nt != '.']
        region = [al1[region[0]], al1[region[1]]]
    # coerce interactions and interactions_filter into filters format
    filters = [{"interactions": interactions} | interactions_filter]
    # initialize plot using all structure drawings
    plot = plots.AP(num_samples=1, nt_length=seq1_full.length, region=region,
                    **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    for sample, seq, panel in zip(samples, [seq1_full, seq2_full],
                                  ["top", "bottom"]):
        fit_data_list(sample, [ct, comp, profile], seq)
        sample.filter_interactions(interactions=interactions2,
                                   fit_to=seq, **interactions2_filter)
        for filt in filters:
            sample.filter_interactions(fit_to=seq, **filt)
            plot.add_sample(ax=0, seq=None, sample=sample, annotation_gap=10,
                            ct=ct, comp=comp, ct_panel=panel,
                            interactions=filt["interactions"],
                            interactions_panel=panel,
                            interactions2=interactions2,
                            interactions2_panel=panel,
                            profile=profile, profile_panel=panel, label='',
                            annotations=[], seqbar=False, **kwargs)
    plots.alignment.plot_alignment(plot=plot, ax=plot.axes[0, 0], data1=seq1,
                                   data2=seq2,
                                   label=labels, center=-5, offset=4,
                                   spines_positions={"top": 0, "bottom": -10})
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_ss(samples, ss="ss", profile="profile", annotations=[],
            interactions=None, interactions_filter=None, interactions2=None,
            interactions2_filter=None, filters=None, labels=None,
            plot_kwargs=None, colorbar=True, **kwargs):
    """Generates a multipanel secondary structure drawing with optional
    coloring by per-nucleotide data and display of inter-nucleotide data and/or
    sequence annotations. Each plot may display a unique sample and/or
    inter-nucleotide data filtering scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        ss (str, optional): a key from sample.data used to retrieve a secondary
            structure containing drawing coordinates.
            Defaults to "ss"
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data may be used to color nucleotides in the
            structure drawing.
            Defaults to "profile".
        annotations (list, optional): a list of keys from sample.data used to
            retrieve sequence annotations. These annotations are highlighted on
            the structure drawing.
            Defaults to [].
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to ss sequence coordinates,
            filtered using interactions_filter arguments, and displayed as
            lines connecting nucleotides in the structure drawing.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=ss. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        interactions2_filter (dict, optional): same as interactions_filter
            above but applied to interactions2.
            Defaults to {}.
        filters (list of dict, optional): For plotting multiple filtering
            schemes applied to each sample, each in a new panel. Each
            dictionary in this list follows a similar structure as
            interactions_filter, but also requires the key-value pair:
            {"interactions": interactions} where interactions is as described
            above. interactions and interactions_filter arguments above will be
            ignored. interactions2 and interactions2_filter will be displayed
            in every panel.
            Defaults to [].
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to SS(). See
            rnavigate.plots.SS for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to SS.plot_data.
            see rnavigate.plots.SS.plot_data for more detail.

    Returns:
        rnavigate.plots.SS plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if interactions2_filter is None:
        interactions2_filter = {}
    if labels is None:
        labels = ["label"]*len(samples)
    if plot_kwargs is None:
        plot_kwargs = {}
    # if filters list given, rows = # samples, columns = # filters
    if ((filters is not None)
            and (len(samples) > 1)
            and ("rows" not in plot_kwargs)
            and ("cols" not in plot_kwargs)):
        plot_kwargs["rows"] = len(samples)
        plot_kwargs["cols"] = len(filters)
    # coerce interactions and interactions_filter into filters format
    elif filters is None:
        filters = [{"interactions": interactions} | interactions_filter]
    # initialize plot using all structure drawings
    num_samples = len(samples) * len(filters)
    plot = plots.SS(num_samples=num_samples, **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    for sample, label in zip(samples, labels):
        fit_data_list(sample, annotations + [profile], sample.data[ss])
        sample.filter_interactions(interactions=interactions2, fit_to=ss,
                                   **interactions2_filter)
        for filt in filters:
            sample.filter_interactions(fit_to=ss, **filt)
            plot.add_sample(sample, structure=ss,
                            interactions=filt["interactions"],
                            interactions2=interactions2, profile=profile,
                            annotations=annotations, label=label, **kwargs)
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_mol(samples, structure="pdb", interactions=None,
             interactions_filter=None, filters=None, profile="profile",
             labels=None, show=True, hide_cylinders=False, colorbar=True,
             custom_function=None, plot_kwargs=None, **kwargs):
    """Generates a multipanel interactive 3D molecular rendering of a PDB
    structure. Nucleotides may be colored by per-nucleotide data or custom
    color lists. Inter-nucleotide data may be displayed as cylinders connecting
    atoms or residues. Each plot may display a unique sample and/or filtering
    scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        structure (str, optional): a key from sample.data to retrieve PDB data
            with atomic coordinates.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to structre sequence
            coordinates, filtered using interactions_filter arguments, and
            displayed as cylinders connecting nucleotides in the 3D structure.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=structure. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        filters (list of dict, optional): For plotting multiple filtering
            schemes applied to each sample, each in a new panel. Each
            dictionary in this list follows a similar structure as
            interactions_filter, but also requires the key-value pair:
            {"interactions": interactions} where interactions is as described
            above. interactions and interactions_filter arguments above will be
            ignored.
            Defaults to [].
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data may be used to color nucleotides.
            Defaults to "profile".
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends.
            Defaults to default sample name.
        show (bool, optional): whether to display the interactive rendering.
            Defaults to True
        hide_cylinders (bool, optional): whether to display cylinders
            representing nucleoside orientation. Setting to false will display
            only the backbone as a ribbon.
            Defaults to False.
        plot_kwargs (dict, optional): kwargs passed to Mol(). See
            rnavigate.plots.Mol for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Mol.plot_data.
            see rnavigate.plots.Mol.plot_data for more detail.

    Returns:
        rnavigate.plots.Mol plot: object containing py3dmol viewer with
            additional plotting and file saving methods
    """
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = ["label"]*len(samples)
    # if filters list given, rows = # samples, columns = # filters
    if ((filters is not None)
            and (len(samples) > 1)
            and ("rows" not in plot_kwargs)
            and ("cols" not in plot_kwargs)):
        plot_kwargs["rows"] = len(samples)
        plot_kwargs["cols"] = len(filters)
    # coerce interactions and interactions_filter into filters format
    elif filters is None:
        filters = [{"interactions": interactions} | interactions_filter]
    num_samples = len(samples) * len(filters)
    # initialize plot using 1st 3D structure (applies to all samples)
    plot = plots.Mol(num_samples=num_samples, pdb=samples[0].data[structure],
                     **plot_kwargs)
    # loop through samples and filters, adding each as a new viewer
    for sample, label in zip(samples, labels):
        fit_data_list(sample, [profile], sample.data[structure])
        for filt in filters:
            sample.filter_interactions(fit_to=structure, **filt)
            plot.add_sample(sample=sample,
                            interactions=filt["interactions"],
                            profile=profile, label=label, **kwargs)
    # apply custom function
    if custom_function is not None:
        custom_function(plot)
    # hide nucleotide cylinders in all viewers
    if hide_cylinders:
        plot.hide_cylinders()
    # show viewer grid
    if show:
        plot.view.show()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_heatmap(samples, structure=None, interactions=None,
                 interactions_filter=None, filters=None, labels=None,
                 plot_kwargs=None, **kwargs):
    """Generates a multipanel plot displaying a heatmap of inter-nucleotide
    data (nucleotide resolution of 2D KDE) and/or contour map of structure
    distances. Each plot may display a unique sample and/or filtering scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        structure (str, optional): a key from sample.data to retrieve either
            PDB data with atomic coordinates (contours outline 3D distance) or
            secondary structure data (contours outline contact distance).
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to seq_source coordinates,
            filtered using interactions_filter arguments, and displayed as
            either nucleotide-resolution heatmaps, or 2D kernel density
            estimate.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=seq_source. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        filters (list of dict, optional): For plotting multiple filtering
            schemes applied to each sample, each in a new panel. Each
            dictionary in this list follows a similar structure as
            interactions_filter, but also requires the key-value pair:
            {"interactions": interactions} where interactions is as described
            above. interactions and interactions_filter arguments above will be
            ignored.
            Defaults to [].
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to Heatmap(). See
            rnavigate.plots.Heatmap for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Heatmap.plot_data.
            see rnavigate.plots.Heatmap.plot_data for more detail.

    Returns:
        rnavigate.plots.Heatmap plot: object containing matplotlib figure and
            axes with additional plotting and file saving methods
    """
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = ["label"]*len(samples)
    # if filters list given, rows = # samples, columns = # filters
    if ((filters is not None)
            and (len(samples) > 1)
            and ("rows" not in plot_kwargs)
            and ("cols" not in plot_kwargs)):
        plot_kwargs["rows"] = len(samples)
        plot_kwargs["cols"] = len(filters)
    # coerce interactions and interactions_filter into filters format
    elif filters is None:
        filters = [{"interactions": interactions} | interactions_filter]
    # initialize plot using 1st 3D structure (applies to all samples)
    num_samples = len(samples) * len(filters)
    plot = plots.Heatmap(
        num_samples, samples[0].data[structure], **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    for sample, label in zip(samples, labels):
        for filt in filters:
            sample.filter_interactions(fit_to=structure, **filt)
            plot.add_sample(sample, interactions=filt["interactions"],
                            label=label, **kwargs)
    plot.set_figure_size()
    return plot


def plot_circle(samples, seq_source=None, ct=None, comp=None,
                interactions=None, interactions_filter=None,
                interactions2=None, interactions2_filter=None, filters=None,
                annotations=None, profile="profile", labels=None,
                colorbar=True, plot_kwargs=None, **kwargs):
    """Generates a multipanel circle plot displaying combinations of secondary
    structures, per-nucleotide data, inter-nucleotide data, and sequence
    annotations. Each plot may display a unique sample and/or filtering scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        seq_source (str or data object, optional): a key from sample.data, a
            sequence string, or a Data object. All data will be mapped to this
            sequence using a user-defined or pairwise sequence alignment.
            Defaults to the value of the ct argument below.
        ct (str, optional): a key from sample.data to retreive a secondary
            structure. Basepairs are plotted as grey arcs within the circle.
            Defaults to "ct".
        comp (str, optional): same as ct. basepairs from ct and comp will be
            plotted. Basepairs are colored by which structure contains them
            (shared, ct only, comp only).
            Defaults to None.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to seq_source coordinates,
            filtered using interactions_filter arguments, and displayed as arcs
            within the circle.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=seq_source. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        interactions2_filter (dict, optional): same as interactions_filter
            above but applied to interactions2.
            Defaults to {}.
        filters (list of dict, optional): For plotting multiple filtering
            schemes applied to each sample, each in a new panel. Each
            dictionary in this list follows a similar structure as
            interactions_filter, but also requires the key-value pair:
            {"interactions": interactions} where interactions is as described
            above. interactions and interactions_filter arguments above will be
            ignored. interactions2 and interactions2_filter will be displayed
            in every panel.
            Defaults to [].
        annotations (list, optional): a list of keys from sample.data used to
            retrieve sequence annotations. These annotations are displayed next
            to the sequence, outside of the circle.
            Defaults to [].
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data may be used to color nucleotides.
            Defaults to "profile".
        labels (str, optional): Same length as samples list. Labels to
            be used as titles.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to Circle(). See
            rnavigate.plots.Circle for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Circle.plot_data.
            see rnavigate.plots.Circle.plot_data for more detail.

    Returns:
        rnavigate.plots.Circle: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if interactions2_filter is None:
        interactions2_filter = {}
    if annotations is None:
        annotations = []
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = ["label"]*len(samples)
    # if filters list given, rows = # samples, columns = # filters
    if ((filters is not None)
            and (len(samples) > 1)
            and ("rows" not in plot_kwargs)
            and ("cols" not in plot_kwargs)):
        plot_kwargs["rows"] = len(samples)
        plot_kwargs["cols"] = len(filters)
    # coerce interactions and interactions_filter into filters format
    elif filters is None:
        filters = [{"interactions": interactions} | interactions_filter]
    # initialize plot
    sequence = get_sequence(seq_source, samples[0], profile)
    num_samples = len(samples) * len(filters)
    plot = plots.Circle(num_samples=num_samples,
                        seq_source=sequence, **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    for sample, label in zip(samples, labels):
        fit_data_list(sample=sample, data_list=annotations+[ct, comp, profile],
                      fit_to=sequence)
        sample.filter_interactions(interactions=interactions2, fit_to=sequence,
                                   **interactions2_filter)
        for filt in filters:
            sample.filter_interactions(fit_to=sequence, **filt)
            plot.add_sample(sample, ct=ct, comp=comp,
                            interactions=filt["interactions"],
                            interactions2=interactions2, profile=profile,
                            annotations=annotations, label=label, **kwargs)
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_linreg(samples, seq_source=None, ct="ct", profile="profile",
                annotations=None, colorby=None, labels=None, column=None,
                plot_kwargs=None, **kwargs):
    """Performs linear regression analysis and generates scatter plots of all
    sample-to-sample profile vs. profile comparisons. Colors nucleotides by
    identity or base-pairing status.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list squared.
        ct (str, optional): a key from sample.data to retreive a secondary
            structure. Scatter plot points may be colored by base-pairing
            status in this structure.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are used for the linear regression.
        labels (str, optional): Same length as samples list. Labels to
            be used in titles.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to LinReg(). See
            rnavigate.plots.LinReg for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to LinReg.plot_data.
            see rnavigate.plots.LinReg.plot_data for more detail.

    Returns:
        rnavigate.plots.LinReg: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    if labels is None:
        labels = ["label"] * len(samples)
    if annotations is None:
        annotations = []
    if plot_kwargs is None:
        plot_kwargs = {}
    plot = plots.LinReg(len(samples), **plot_kwargs)
    sequence = get_sequence(seq_source, samples[0], profile)
    for sample, label in zip(samples, labels):
        fit_data_list(sample=sample, data_list=annotations+[ct, profile],
                      fit_to=sequence)
        plot.add_sample(sample, sequence=sequence, ct=ct, profile=profile,
                        label=label, colorby=colorby, column=column,
                        annotations=annotations, **kwargs)
    plot.set_figure_size()
    return plot


def plot_roc(samples, ct="ct", profile="profile", labels=None,
             plot_kwargs=None, **kwargs):
    """Performs receiver operator characteristic analysis (ROC), calculates
    area under ROC curve (AUC), and generates ROC plots to assess how well
    per-nucleotide data predicts base-paired status. Does this for all
    positions as well as positions categorized by nucleotide
    5 plots: All, A, U, C, G

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            All samples are plotted on the same set of axes.
        ct (str, optional): a key from sample.data to retreive a secondary
            structure. Base-pairing status retreived from this data.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are used for the ROC/AUC analysis.
        labels (str, optional): Same length as samples list. Labels to
            be used in legends.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to ROC(). See
            rnavigate.plots.ROC for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to ROC.plot_data.
            see rnavigate.plots.ROC.plot_data for more detail.

    Returns:
        rnavigate.plots.ROC: object containing matplotlib figure and axes with
            additional plotting and file saving methods
    """
    if labels is None:
        labels = ["label"] * len(samples)
    if plot_kwargs is None:
        plot_kwargs = {}
    plot = plots.ROC(len(samples), **plot_kwargs)
    for sample, label in zip(samples, labels):
        plot.add_sample(sample, ct=ct, profile=profile, label=label, **kwargs)
    plot.set_figure_size()
    return plot


def plot_disthist(samples, structure="pdb", interactions=None,
                  interactions_filter=None, bg_interactions=None,
                  bg_interactions_filter=None, filters=None, labels=None,
                  same_axis=False, plot_kwargs=None, **kwargs):
    """Calculates 3D distance of nucleotides in inter-nucleotide data and plots
    the distribution of these distances. Compares this to a 'background'
    distribution consisting of either all pairwise distances in structure, or
    those defined by bg_interactions and bg_interactions_filter

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list unless
            same_axis is set to True.
        structure (str, optional): a key from sample.data to retreive a PDB
            structure with atomic coordinates.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to PDB sequence coordinates,
            filtered using interactions_filter arguments, and used to calculate
            distance distribution histograms.
            Defaults to None.
        interactions_filter (dict, optional): These key-value pairs are passed
            as keyword arguments to sample.filter_interactions along with
            interactions=interactions and fit_to=seq_source. See
            rnavigate.Sample.filter_interactions for more detail.
            Defaults to {}.
        bg_interactions (str, optional): same as interactions above. used to
            calulate 'background' distance distribution histograms.
            Defaults to None.
        bg_interactions_filter (dict, optional): same as interactions_filter
            above but applied to bg_interactions.
            Defaults to {}.
        filters (list of dict, optional): For plotting multiple filtering
            schemes applied to each sample. Each dictionary in this list
            follows a similar structure as interactions_filter, but also
            requires the key-value pair: {"interactions": interactions} where
            interactions is as described above. interactions and
            interactions_filter arguments above will be ignored.
            Defaults to [].
        labels (str, optional): Same length as samples list. Labels to
            be used as titles.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to DistHist(). See
            rnavigate.plots.DistHist for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to DistHist.plot_data
            see rnavigate.plots.DistHist.plot_data for more detail.

    Returns:
        rnavigate.plots.DistHist: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    # use mutable defaults
    if interactions_filter is None:
        interactions_filter = {}
    if bg_interactions_filter is None:
        bg_interactions_filter = {}
    if plot_kwargs is None:
        plot_kwargs = {}
    if labels is None:
        labels = ["label"]*len(samples)
    # if filters list given, rows = # samples, columns = # filters
    if ((filters is not None)
            and (len(samples) > 1)
            and ("rows" not in plot_kwargs)
            and ("cols" not in plot_kwargs)
            and (not same_axis)):
        plot_kwargs["rows"] = len(samples)
        plot_kwargs["cols"] = len(filters)
    # coerce interactions and interactions_filter into filters format
    elif filters is None:
        filters = [{"interactions": interactions} | interactions_filter]
    # initialize plot
    num_samples = len(samples) * len(filters)
    if same_axis:
        plot = plots.DistHist(num_samples=1, **plot_kwargs)
        ax = plot.axes[0, 0]
    else:
        plot = plots.DistHist(num_samples=num_samples, **plot_kwargs)
        ax = None
    # loop through samples and filters, adding each as a new axis
    for sample, label in zip(samples, labels):
        sample.filter_interactions(interactions=bg_interactions,
                                   fit_to=structure, **bg_interactions_filter)
        for filt in filters:
            sample.filter_interactions(fit_to=structure, **filt)
            plot.add_sample(sample, structure=structure,
                            interactions=filt["interactions"],
                            bg_interactions=bg_interactions, label=label,
                            ax=ax, **kwargs)
    plot.set_figure_size()
    return plot
