#!/usr/bin/env python

# general python packages
import os.path
import numpy as np

# modules in RNAvigate
from .data import Data, PDB, Log
from .data import Annotation, Motif, ORFs
from .data import CT, DotBracket, XRNA, VARNA, NSD, CTE, get_ss_class
from .data import Interactions, RINGMaP, PAIRMaP, PairProb, SHAPEJuMP
from .data import Profile, SHAPEMaP, DanceMaP, RNPMaP
from .plots import AP, Circle, DistHist, Heatmap, LinReg, Mol
from .plots import QC, Skyline, SM, SS, ROC
from .analysis import LogCompare, LowSS


def create_code_button():
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
                 dance_prefix=None,
                 sites=None,
                 spans=None,
                 groups=None,
                 primers=None,
                 orfs=None,
                 motif=None):
        """Creates a sample object which connects all chemical probing and
        structural data for a single experiment. Contains convenience methods
        to plot, filter, compare and retrieve this data. Every argument is
        optional and defaults to None.

        Args:
            sample (str, optional): This string will serve as a label for this
                sample within plots.
            pdb (str, optional): A dictionary, which must contain at least
                "filepath": the path to your pdb file
                "chain": the chain ID of your RNA in the pdb file
                Other options include "offset" and "fasta".
            ct (str, optional): Path to a .ct, .dbn, or .bracket structure
                file.
            compct (str, optional): Same as ct above, but serves as a
                comparison structure for plotting.
            ss (str, optional): Path to a .xrna, .varna, .nsd, or .cte
                secondary structure drawing file.
            fasta (str, optional): Path to a fasta file which will provide the
                reference sequence for ShapeJumper Data.
            log (str, optional): Path to a ShapeMapper shapemapper_log.txt.
            shapemap (str, optional): Path to a ShapeMapper profile.txt.
            dmsmap (str, optional): Path to a ShapeMapper profile.txt, which
                will be renormalized using DMS conventions.
            dancemap (str, optional): Use dance_prefix instead.
            rnpmap (str, optional): Path to an RNPMapper .csv file.
            ringmap (str, optional): Path to a RingMapper output file.
            shapejump (str, optional): Path to a ShapeJumper output file.
            pairmap (str, optional): Path to a PairMapper pairmap.txt output
                file.
            allcorrs (str, optional): Path to a PairMapper allcorrs.txt output
                file.
            pairprob (str, optional): Path to an RNAStructure .dp pairing
                probability file.
            dance_prefix (str, optional): Path prefix for DanceMapper output
                files.
        """
        if dancemap is None:
            dancemap = {"filepath": None}
        if sites is None:
            sites = {}
        if spans is None:
            spans = {}
        if groups is None:
            groups = {}
        if primers is None:
            primers = {}
        if orfs is None:
            orfs = {}
        if motif is None:
            motif = {}
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
                "instantiator": Data,
                "seq_source": "self",
                "kwargs": {}},
            "log": {
                "filepath": log,
                "instantiator": Log,
                "seq_source": "self",
                "kwargs": {}},
            "shapemap": {
                "filepath": shapemap,
                "instantiator": SHAPEMaP,
                "seq_source": "self",
                "kwargs": {}},
            "dmsmap": {
                "filepath": dmsmap,
                "instantiator": SHAPEMaP,
                "seq_source": "self",
                "kwargs": {"dms": True}},
            "dancemap": {
                "filepath": dancemap.pop("filepath"),
                "instantiator": DanceMaP,
                "seq_source": "self",
                "kwargs": dancemap},
            "rnpmap": {
                "filepath": rnpmap,
                "instantiator": RNPMaP,
                "seq_source": "self",
                "kwargs": {}},
            "ringmap": {
                "filepath": ringmap,
                "instantiator": RINGMaP,
                "seq_source": "profile",
                "kwargs": {}},
            "pairmap": {
                "filepath": pairmap,
                "instantiator": PAIRMaP,
                "seq_source": "profile",
                "kwargs": {}},
            "allcorrs": {
                "filepath": allcorrs,
                "instantiator": RINGMaP,
                "seq_source": "profile",
                "kwargs": {}},
            "shapejump": {
                "filepath": shapejump.pop("filepath", None),
                "instantiator": SHAPEJuMP,
                "seq_source": "self",
                "kwargs": shapejump},
            "pairprob": {
                "filepath": pairprob,
                "instantiator": PairProb,
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
                "instantiator": PDB,
                "seq_source": "self",
                "kwargs": pdb},
            "dance_prefix": {
                "filepath": dance_prefix,
                "instantiator": self.init_dance,
                "seq_source": "self",
                "kwargs": {}},
            "sites": {
                "filepath": "",
                "instantiator": Annotation,
                "seq_source": sites.pop("seq_source", ""),
                "kwargs": sites},
            "spans": {
                "filepath": "",
                "instantiator": Annotation,
                "seq_source": spans.pop("seq_source", ""),
                "kwargs": spans},
            "groups": {
                "filepath": "",
                "instantiator": Annotation,
                "seq_source": groups.pop("seq_source", ""),
                "kwargs": groups},
            "primers": {
                "filepath": "",
                "instantiator": Annotation,
                "seq_source": primers.pop("seq_source", ""),
                "kwargs": primers},
            "orfs": {
                "filepath": "",
                "instantiator": ORFs,
                "seq_source": orfs.pop("seq_source", ""),
                "kwargs": orfs},
            "motif": {
                "filepath": "",
                "instantiator": Motif,
                "seq_source": motif.pop("seq_source", ""),
                "kwargs": motif},
        }
        self.default_profiles = ["shapemap", "dmsmap", "dancemap", "rnpmap"]

        self.sample = sample
        self.parent = None
        self.data = {}  # stores profile, interactions, and structure objects
        # load data
        for input, kwargs in self.inputs.items():
            if kwargs["filepath"] is not None and kwargs["seq_source"] != "":
                # additional_kwargs = kwargs.pop("kwargs", {})
                kwargs.update(kwargs.pop("kwargs", {}))
                self.set_data(name=input, **kwargs)

    def set_data(self, name, filepath=None, instantiator=None, seq_source=None,
                 **kwargs):
        if instantiator is None:
            instantiator = self.inputs[name]["instantiator"]
        if seq_source is None:
            seq_source = self.inputs[name]["seq_source"]
        if seq_source != "self":
            kwargs.update({"sequence": self.data[seq_source].sequence})
        if isinstance(filepath, list):
            for i, path in enumerate(filepath):
                self.set_data(name=f"{name}_{i+1}",
                              filepath=path,
                              instantiator=instantiator,
                              seq_source=seq_source,
                              **kwargs)
        self.inputs[name] = {
            "filepath": filepath,
            "instantiator": instantiator,
            "seq_source": seq_source,
            "kwargs": kwargs}
        self.data[name] = instantiator(filepath=filepath, **kwargs)
        if name in self.default_profiles and "profile" not in self.data:
            self.data["profile"] = self.data[name]

    def init_dance(self, filepath):
        """Initializes a list of Sample objects which each represent a
        component of the DANCE model, available as Sample.dance.

        Args:
            prefix (str): Path to DanceMapper output file prefixes. Finds
                profiles, rings, pairs, and structure predictions if present.
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
               of Sample.data or "label"

        Returns:
            Data object
        """
        if key == "label":
            return self.sample
        elif key is None:
            return None
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
                            suppress=False, **kwargs):
        """Aligns sequence to fit_to, sets properties, filters and aligns data.

        For example, for plotting interaction data containing structure
        cassettes against a CT without, or plotting sequence variants together
        using a best pairwise alignment. If kwargs contains metric, cmap, or
        min_max, these properties of the Interactions object are set. Other
        kwargs are passed to Interactions.filter(). Refer to that method for
        more detail on filters: help(MaP.Interactions.filter)

        Args:
            interactions (Interactions): Interactions object to be filtered
            fit_to (Data (any)): Data object containing a sequence for
                Interactions to be fitted to for plotting purposes.
            **kwargs: metric, cmap, and min_max. Others passed to
                interactions.filter() For metric="Distance" you can specify an
                atom or reagent: metric="Distance_O2'" or metric="Distance_DMS"
        """
        data = self.data[interactions]
        # check for valid interactions data
        if interactions not in self.data.keys():
            if not suppress:
                print(f"{interactions} not found in sample data")
            return
        if not isinstance(data, Interactions):
            if not suppress:
                print(f"{Interactions} is not an Interactions datatype")
            return

        if "metric" in kwargs.keys():
            metric = kwargs.pop("metric")
            if metric.startswith("Distance"):
                metric = (metric, self.data["pdb"])
            data.metric = metric
        else:
            data.metric = data.default_metric
        if "cmap" in kwargs.keys():
            cmap = kwargs.pop("cmap")
            data.cmap = cmap
        if "min_max" in kwargs.keys():
            min_max = kwargs.pop("min_max")
            data.min_max = min_max
        for datatype in ["profile", "ct"]:
            if datatype in kwargs.keys():
                kwargs[datatype] = self.data[kwargs[datatype]]
            elif datatype in self.data.keys():
                kwargs[datatype] = self.data[datatype]
        data.filter(self.get_data_list(fit_to), **kwargs)

    def dance_filter(self, filterneg=True, cdfilter=15, sigfilter=23,
                     ssfilter=True):
        """Applies a standard filter to plot DANCE rings, pairs, and predicted
        structures together.

        Args:
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
        """
        kwargs = {}
        if filterneg:
            kwargs["Sign"] = 0
        if ssfilter:
            kwargs["ss_only"] = True
        kwargs["Statistic"] = sigfilter
        ctlist = [dance.data["ct"] for dance in self.dance]
        for dance in self.dance:
            dance_ct = dance.data["ct"]
            dance.data["ringmap"].filter(dance_ct, ct=dance_ct, **kwargs)
            dance.data["ringmap"].mask_on_ct(ctlist, min_cd=cdfilter)
            dance.data["pairmap"].filter(dance_ct, ct=dance_ct,
                                         paired_only=True)

###############################################################################
# sample plotting functions
#     plot_qc
#     plot_ss
#     plot_mol
#     plot_heatmap
#     plot_circle
#     plot_disthist
#     plot_skyline
#     plot_arcs
#     plot_shapemapper
#     plot_arcs_multifilter
#     plot_ss_multifilter
#     plot_mol_multifilter
#     plot_circle_multifilter
#     plot_disthist_multifilter
###############################################################################

    def plot_qc(self, **kwargs):
        """Makes a QC plot. See help(MaP.plot_qc_multisample) for more."""
        return plot_qc_multisample([self], **kwargs)

    def plot_ss(self, dance=False, **kwargs):
        if dance:
            self.dance_filter()
            plot = plot_ss_multisample(self.dance, prefiltered=True, **kwargs)
            for i, dance in enumerate(self.dance):
                ax = plot.get_ax(i)
                ax.set_title(
                    f"DANCE component: {i}, Percent: {self.dance_percents[i]}")
            return plot
        return plot_ss_multisample([self], **kwargs)

    def plot_mol(self, dance=False, **kwargs):
        if dance:
            self.dance_filter()
            plot = plot_arcs_multisample(
                self.dance, prefiltered=True, **kwargs)
            for i, dance in enumerate(self.dance):
                ax = plot.get_ax(i)
                ax.set_title(
                    f"DANCE component: {i}, Percent: {self.dance_percents[i]}")
            return plot
        return plot_mol_multisample([self], **kwargs)

    def plot_heatmap(self, **kwargs):
        return plot_heatmap_multisample([self], **kwargs)

    def plot_circle(self, **kwargs):
        return plot_circle_multisample([self], **kwargs)

    def plot_disthist(self, **kwargs):
        return plot_disthist_multisample([self], **kwargs)

    def plot_skyline(self, dance=False, **kwargs):
        if dance:
            plot = plot_skyline_multisample(self.dance, **kwargs)
            plot.axes[0, 0].legend(title="Comp: Percent")
            plot.axes[0, 0].set_title(f"{self.sample}: DANCE Reactivities")
            return plot
        plot = plot_skyline_multisample([self], **kwargs)
        return plot

    def plot_arcs(self, dance=False, **kwargs):
        if dance:
            self.dance_filter()
            plot = plot_arcs_multisample(samples=self.dance,
                                         prefiltered=True, **kwargs)
            return plot
        return plot_arcs_multisample([self], **kwargs)

    def plot_shapemapper(self, plots=["profile", "rates", "depth"]):
        plot = SM(self.data["profile"].length, plots=plots)
        plot.add_sample(self, profile="profile", label="label")
        return plot

    def plot_arcs_multifilter(self, filters, ct="ct", comp=None,
                              interactions2=None, profile="profile",
                              label="label", region="all"):
        """Makes an array of arc plots of different filtered views of data from
        Sample.

        Args:
            filters (dict): Dictionary containing filtering kwargs for each
                plot in this array. Requires "interactions" and "fit_to".
            ct (str, optional): key for data object to use as ct.
                Defaults to "ct".
            comp (str, optional): key for data object to use as comparison ct.
                Defaults to "compct".
            interactions2 (_type_, optional): key for data object to use as
                second interactions.
                Defaults to None.
            profile (str, optional): key for data object to use as profile.
                Defaults to "profile".
            label (str, optional): Defaults to "label", same as Sample.sample.

        Returns:
            rnavigate.plots.AP
        """
        plot = AP(len(filters), self.get_data_list(ct).length, region=region)
        for filter in filters:
            interactions = filter.pop("interactions")
            self.filter_interactions(interactions, ct, **filter)
            plot.add_sample(self, ct=ct, comp=comp, interactions=interactions,
                            interactions2=interactions2,
                            profile=profile, label=label)
        return plot

    def plot_ss_multifilter(self, filters, ss="ss", profile="profile",
                            label="label", interactions2=None, **kwargs):
        plot = SS(len(filters), self.get_data_list(ss))
        pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
        for filter in filters:
            interactions = filter.pop("interactions")
            self.filter_interactions(interactions, ss, **filter)
            plot.add_sample(self, interactions=interactions,
                            interactions2=interactions2, profile=profile,
                            label=label, **pt_kwargs)
        return plot

    def plot_mol_multifilter(self, filters, profile="profile", label="label",
                             show=True, **kwargs):
        plot = Mol(len(filters), self.data["pdb"])
        pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
        for filter in filters:
            interactions = filter.pop("interactions")
            self.filter_interactions(interactions, "pdb", **filter)
            plot.add_sample(self, interactions=interactions, profile=profile,
                            label=label, **pt_kwargs)
        if show:
            plot.view.show()
        return plot

    def plot_circle_multifilter(self, filters, ct=None, comp=None,
                                interactions2=None, profile=None,
                                label="label"):
        plot = Circle(len(filters), self.data["profile"].length)
        for filter in filters:
            interactions = filter.pop("interactions")
            self.filter_interactions(interactions, "profile", **filter)
            plot.add_sample(self, ct=ct, comp=comp, interactions=interactions,
                            interactions2=interactions2,
                            profile=profile, label=label)
        return plot

    def plot_disthist_multifilter(self, filters, structure="pdb",
                                  interactions=None, label="label",
                                  same_axis=True):
        if same_axis:
            plot = DistHist(1)
            ax = plot.axes[0, 0]
        else:
            plot = DistHist(len(filters))
            ax = None
        for filter in filters:
            interactions = filter.pop("interactions")
            self.filter_interactions(interactions, "pdb", **filter)
            plot.add_sample(structure, interactions, label, ax)
        return plot

###############################################################################
# Plotting functions that accept a list of samples
#   extract_passthrough_kwargs
#   plot_qc_multisample
#   plot_skyline_multisample
#   plot_arcs_multisample
#   plot_ss_multisample
#   plot_mol_multisample
#   plot_heatmap_multisample
#   plot_circle_multisample
#   plot_linreg_multisample
#   plot_disthist_multisample
###############################################################################


def extract_passthrough_kwargs(plot, kwargs):
    pt_kwargs = {}
    for kw in plot.pass_through:
        if kw in kwargs.keys():
            pt_kwargs[kw] = kwargs.pop(kw)
    return pt_kwargs


def plot_qc_multisample(samples=[], labels=None, plot_kwargs={}, **kwargs):
    plot = QC(num_samples=len(samples), **plot_kwargs)
    if labels is None:
        labels = ["label"]*len(samples)
    for sample, label in zip(samples, labels):
        plot.add_sample(sample=sample, log="log", profile="profile",
                        label=label, **kwargs)
    return plot


def plot_skyline_multisample(samples, profile="profile", labels=None,
                             annotations=[], region="all", plot_kwargs={},
                             **kwargs):
    plot = Skyline(num_samples=len(samples),
                   nt_length=samples[0].data[profile].length,
                   region=region,
                   **plot_kwargs)
    for sample in samples:
        plot.add_sample(sample, profile=profile, annotations=annotations,
                        label="label", **kwargs)
    return plot


def plot_arcs_multisample(samples, ct="ct", comp=None,
                          interactions=None, interactions_filter={},
                          interactions2=None, interactions2_filter={},
                          profile="profile", annotations=[], labels=None,
                          region="all", plot_kwargs={}, prefiltered=False,
                          **kwargs):
    plot = AP(num_samples=len(samples), nt_length=samples[0].data[ct].length,
              region=region, **plot_kwargs)
    if labels is None:
        labels = ["label"]*len(samples)
    for sample, label in zip(samples, labels):
        if not prefiltered:
            if ct not in ["ss", "ct", "compct"]:
                sample.filter_interactions(interactions=ct, fit_to=ct)
            if interactions is not None:
                sample.filter_interactions(interactions=interactions,
                                           fit_to=ct, **interactions_filter)
            if interactions2 is not None:
                sample.filter_interactions(interactions=interactions2,
                                           fit_to=ct, **interactions2_filter)
        plot.add_sample(sample=sample, ct=ct, comp=comp,
                        interactions=interactions, interactions2=interactions2,
                        profile=profile, label=label, annotations=annotations,
                        **kwargs)
    return plot


def plot_ss_multisample(samples, ss="ss", profile="profile", annotations=[],
                        interactions=None, interactions_filter={},
                        interactions2=None, interactions2_filter={},
                        labels=None, plot_kwargs={}, prefiltered=False,
                        **kwargs):
    plot = SS(len(samples), samples[0].data[ss], **plot_kwargs)
    if labels is None:
        labels = ["label"]*len(samples)
    for sample, label in zip(samples, labels):
        if not prefiltered:
            if interactions is not None:
                sample.filter_interactions(interactions, ss,
                                           **interactions_filter)
            if interactions2 is not None:
                sample.filter_interactions(interactions2, ss,
                                           **interactions2_filter)
        plot.add_sample(sample, interactions=interactions,
                        interactions2=interactions2, profile=profile,
                        annotations=annotations, label=label, **kwargs)
    return plot


def plot_mol_multisample(samples, structure="pdb",
                         interactions=None, interactions_filter={},
                         profile="profile", labels=None, show=True,
                         width=400, height=400, background_alpha=1,
                         hide_cylinders=False,
                         prefiltered=False, **kwargs):
    plot = Mol(len(samples), samples[0].data[structure], width=width,
               height=height, background_alpha=background_alpha)
    if labels is None:
        labels = ["label"]*len(samples)
    for sample, label in zip(samples, labels):
        if not prefiltered:
            if interactions is not None:
                sample.filter_interactions(interactions, structure,
                                           **interactions_filter)
        plot.add_sample(sample, interactions=interactions, profile=profile,
                        label=label, **kwargs)
    if hide_cylinders:
        plot.hide_cylinders()
    if show:
        plot.view.show()
    return plot


def plot_heatmap_multisample(samples, structure=None, interactions=None,
                             interactions_filter={}, label="label", **kwargs):
    plot = Heatmap(len(samples), samples[0].data[structure])
    for sample in samples:
        sample.filter_interactions(interactions, structure,
                                   **interactions_filter)
        plot.add_sample(sample, interactions=interactions,
                        label=label, **kwargs)
    return plot


def plot_circle_multisample(samples, seq_source="profile", ct=None, comp=None,
                            interactions=None, interactions_filter={},
                            interactions2=None, interactions2_filter={},
                            annotations=[], prefiltered=False,
                            profile="profile", label="label", **kwargs):
    plot = Circle(len(samples), samples[0].data[seq_source])
    for sample in samples:
        if not prefiltered:
            if interactions is not None:
                sample.filter_interactions(interactions, seq_source,
                                           **interactions_filter)
            if interactions2 is not None:
                sample.filter_interactions(interactions2, seq_source,
                                           **interactions2_filter)
        plot.add_sample(sample, ct=ct, comp=comp, interactions=interactions,
                        interactions2=interactions2, profile=profile,
                        annotations=annotations, label=label, **kwargs)
    return plot


def plot_linreg_multisample(samples, ct="ct", profile="profile", label="label",
                            **kwargs):
    plot = LinReg(len(samples))
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        plot.add_sample(sample, ct=ct, profile=profile,
                        label=label, **pt_kwargs)
    return plot


def plot_roc_multisample(samples, ct="ct", profile="profile", label="label",
                         **kwargs):
    plot = ROC(len(samples))
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        plot.add_sample(sample, ct=ct, profile=profile,
                        label=label, **pt_kwargs)
    return plot


def plot_disthist_multisample(samples, structure="pdb", interactions=None,
                              interactions_filter={},
                              label="label", same_axis=False, **kwargs):
    if same_axis:
        plot = DistHist(1)
        ax = plot.axes[0, 0]
    else:
        plot = DistHist(len(samples))
        ax = None
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        sample.filter_interactions(interactions=interactions, fit_to="pdb",
                                   **interactions_filter)
        plot.add_sample(sample, structure=structure, interactions=interactions,
                        label=label, ax=ax, **pt_kwargs)
