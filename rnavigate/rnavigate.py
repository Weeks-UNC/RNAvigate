#!/usr/bin/env python

# general python packages
import os.path
import numpy as np

# modules in RNAvigate
from .data import Data, PDB, Log
from .data import Annotation, Motif, ORFs
from .data import CT, DotBracket, XRNA, VARNA, NSD, CTE
from .data import IJ, RINGMaP, PAIRMaP, PairProb, SHAPEJuMP
from .data import Profile, SHAPEMaP, DanceMaP, RNPMaP
from .plots import AP, Circle, DistHist, Heatmap, LinReg, Mol
from .plots import QC, Skyline, SM, SS
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
                 pdb_kwargs={'chain': 'A'},
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
                 orfs=None,
                 motif=None):
        """Creates a sample object which connects all chemical probing and
        structural data for a single experiment. Contains convenience methods
        to plot, filter, compare and retrieve this data. Every argument is
        optional and defaults to None.

        Args:
            sample (str, optional): This string will serve as a label for this
                sample within plots.
            pdb (str, optional): Path to a .pdb atom coordinates file.
            pdb_kwargs (dict, optional): Dictionary containing info needed to
                parse provided .pdb. "chain" is required. "fasta" and "offset"
                may be required if not provided in header.
                Defaults to {'chain':'A'}.
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
            probs (str, optional): Path to an RNAStructure .dp pairing
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
        if orfs is None:
            orfs = {}
        if motif is None:
            motif = {}

        # returns the instantiator using the filepath extension
        def get_ss_class(filepath):
            if filepath is None:
                return None
            extension = filepath.split('.')[-1]
            return {"varna": VARNA,
                    "xrna": XRNA,
                    "nsd": NSD,
                    "cte": CTE,
                    "ct": CT,
                    "dbn": DotBracket,
                    "bracket": DotBracket
                    }[extension]
        # for each input
        # [0] filepath
        # [1] instantiation class
        # [2] sequence source
        # [3] kwargs
        self.inputs = {
            "fasta": [fasta, Data, "self", {}],
            "log": [log, Log, "self", {}],
            "shapemap": [shapemap, SHAPEMaP, "self", {}],
            "dmsmap": [dmsmap, SHAPEMaP, "self", {"dms": True}],
            "dancemap": [dancemap.pop("filepath"), DanceMaP, "self", dancemap],
            "rnpmap": [rnpmap, RNPMaP, "self", {}],
            "ringmap": [ringmap, RINGMaP, "profile", {}],
            "pairmap": [pairmap, PAIRMaP, "profile", {}],
            "allcorrs": [allcorrs, RINGMaP, "profile", {}],
            "shapejump": [shapejump, SHAPEJuMP, "fasta", {}],
            "pairprob": [pairprob, PairProb, "profile", {}],
            "ct": [ct, get_ss_class(ct), "self", {}],
            "compct": [compct, get_ss_class(compct), "self", {}],
            "ss": [ss, get_ss_class(ss), "self", {}],
            "pdb": [pdb, PDB, "self", pdb_kwargs],
            "dance_prefix": [dance_prefix, self.init_dance, "self", {}],
            "sites": ["", Annotation, sites.pop("seq_source", ""), sites],
            "spans": ["", Annotation, spans.pop("seq_source", ""), spans],
            "groups": ["", Annotation, groups.pop("seq_source", ""), groups],
            "orfs": ["", ORFs, orfs.pop("seq_source", ""), orfs],
            "motif": ["", Motif, motif.pop("seq_source", ""), motif],
        }
        self.default_profiles = ["shapemap", "dmsmap", "dancemap", "rnpmap"]

        self.sample = sample
        self.parent = None
        self.data = {}  # stores profile, ij, and structure objects
        # load data
        for input, (path, instantiator, source, kwargs) in self.inputs.items():
            if path is not None and source != "":
                self.set_data(
                    name=input,
                    filepath=path,
                    instantiator=instantiator,
                    seq_source=source,
                    kwargs=kwargs)

    def set_data(self, name, filepath, instantiator, seq_source, kwargs,
                 overwrite=False):
        if not overwrite:
            assert name not in self.data, (f"{name} data already exists. "
                                           "Set overwrite=True to ignore.")
        if seq_source != "self":
            kwargs.update({"sequence": self.data[seq_source].sequence})
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
            for key in ["ringmap", "pairmap"]:
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
#     filter_ij
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
        except (KeyError, TypeError) as _:
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

    def filter_ij(self, ij, fit_to, suppress=False, **kwargs):
        """Aligns sequence to fit_to, sets properties, filters and aligns data.

        For example, for plotting IJ data containing structure cassettes
        against a CT without, or plotting sequence variants together using a
        best pairwise alignment. If kwargs contains metric, cmap, or min_max,
        these properties of the IJ data are set. Other kwargs are passed to
        IJ.filter(). Refer to that method for more detail on filters:
        help(MaP.IJ.filter)

        Args:
            ij (IJ): IJ object to be filtered and fitted
            fit_to (Data (any)): Data object containing a sequence for IJ to be
                fit to for plotting purposes.
            **kwargs: metric, cmap, and min_max. Others passed to ij.filter()
                      For metric="Distance" you can specify an atom or reagent:
                      metric="Distance_O3'" or metric="Distance_DMS"
        """
        # check for valid ij
        if ij not in self.data.keys():
            if not suppress:
                print(f"{ij} not found in sample data")
            return
        if not isinstance(self.data[ij], IJ):
            if not suppress:
                print(f"{ij} is not an IJ datatype")
            return

        if "metric" in kwargs.keys():
            metric = kwargs.pop("metric")
            if metric.startswith("Distance"):
                metric = (metric, self.data["pdb"])
            self.data[ij].metric = metric
        else:
            self.data[ij].metric = self.data[ij].default_metric
        if "cmap" in kwargs.keys():
            cmap = kwargs.pop("cmap")
            self.data[ij].cmap = cmap
        if "min_max" in kwargs.keys():
            min_max = kwargs.pop("min_max")
            self.data[ij].min_max = min_max
        for data in ["profile", "ct"]:
            if data in kwargs.keys():
                continue
            elif data in self.data.keys():
                kwargs[data] = self.data[data]
        self.data[ij].filter(self.get_data_list(fit_to), **kwargs)

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
            dance.data["rings"].filter(dance_ct, ct=dance_ct, **kwargs)
            dance.data["rings"].mask_on_ct(ctlist, cdAbove=cdfilter)
            dance.data["pairs"].filter(dance_ct, ct=dance_ct, paired_only=True)

###############################################################################
# sample plotting functions
#     make_qc
#     make_ss
#     make_mol
#     make_heatmap
#     make_circle
#     make_disthist
#     make_skyline
#     make_ap
#     make_shapemapper
#     make_ap_multifilter
#     make_ss_multifilter
#     make_mol_multifilter
#     make_circle_multifilter
#     make_disthist_multifilter
###############################################################################

    def make_qc(self, **kwargs):
        """Makes a QC plot. See help(MaP.array_qc) for more."""
        return array_qc([self], **kwargs)

    def make_ss(self, dance=False, **kwargs):
        if dance:
            self.dance_filter()
            plot = array_ss(self.dance, prefiltered=True, **kwargs)
            for i, dance in enumerate(self.dance):
                ax = plot.get_ax(i)
                ax.set_title(
                    f"DANCE component: {i}, Percent: {self.dance_percents[i]}")
            return plot
        return array_ss([self], **kwargs)

    def make_mol(self, dance=False, **kwargs):
        if dance:
            self.dance_filter()
            plot = array_ap(self.dance, prefiltered=True, **kwargs)
            for i, dance in enumerate(self.dance):
                ax = plot.get_ax(i)
                ax.set_title(
                    f"DANCE component: {i}, Percent: {self.dance_percents[i]}")
            return plot
        return array_mol([self], **kwargs)

    def make_heatmap(self, **kwargs):
        return array_heatmap([self], **kwargs)

    def make_circle(self, **kwargs):
        return array_circle([self], **kwargs)

    def make_disthist(self, **kwargs):
        return array_disthist([self], **kwargs)

    def make_skyline(self, dance=False, **kwargs):
        if dance:
            plot = array_skyline(self.dance, **kwargs)
            plot.axes[0, 0].legend(title="Comp: Percent")
            plot.axes[0, 0].set_title(f"{self.sample}: DANCE Reactivities")
            return plot
        plot = array_skyline([self], **kwargs)
        return plot

    def make_ap(self, dance=False, **kwargs):
        if dance:
            self.dance_filter()
            plot = array_ap(self.dance, prefiltered=True, **kwargs)
            for i, dance in enumerate(self.dance):
                ax = plot.get_ax(i)
                ax.set_title(
                    f"DANCE component: {i}, Percent: {self.dance_percents[i]}")
            return plot
        return array_ap([self], **kwargs)

    def make_shapemapper(self, plots=["profile", "rates", "depth"]):
        plot = SM(self.data["profile"].length, plots=plots)
        plot.add_sample(self, profile="profile", label="label")
        return plot

    def make_ap_multifilter(self, filters, ct="ct", comp=None, ij2=None,
                            profile="profile", label="label"):
        """Makes an array of arc plots of different filtered views of data from
        Sample.

        Args:
            filters (dict): Dictionary containing filtering kwargs for each
                plot in this array. Requires "ij" and "fit_to".
            ct (str, optional): key for data object to use as ct.
                Defaults to "ct".
            comp (str, optional): key for data object to use as comparison ct.
                Defaults to "compct".
            ij2 (_type_, optional): key for data object to use as second ij.
                Defaults to None.
            profile (str, optional): key for data object to use as profile.
                Defaults to "profile".
            label (str, optional): Defaults to "label", same as Sample.sample.

        Returns:
            rnavigate.plots.AP
        """
        plot = AP(len(filters), self.get_data_list(ct).length)
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, ct, **filter)
            plot.add_sample(self, ct=ct, comp=comp, ij=ij, ij2=ij2,
                            profile=profile, label=label)
        return plot

    def make_ss_multifilter(self, filters, ss="ss", profile="profile",
                            label="label", ij2=None, **kwargs):
        plot = SS(len(filters), self.get_data_list(ss))
        pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, ss, **filter)
            plot.add_sample(self, ij=ij, ij2=ij2, profile=profile,
                            label=label, **pt_kwargs)
        return plot

    def make_mol_multifilter(self, filters, profile="profile", label="label",
                             show=True, **kwargs):
        plot = Mol(len(filters), self.data["pdb"])
        pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "pdb", **filter)
            plot.add_sample(self, ij=ij, profile=profile, label=label,
                            **pt_kwargs)
        if show:
            plot.view.show()
        return plot

    def make_circle_multifilter(self, filters, ct=None, comp=None, ij2=None,
                                profile=None, label="label"):
        plot = Circle(len(filters), self.data["profile"].length)
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "profile", **filter)
            plot.add_sample(self, ct=ct, comp=comp, ij=ij, ij2=ij2,
                            profile=profile, label=label)
        return plot

    def make_disthist_multifilter(self, filters, structure="pdb", ij=None,
                                  label="label", same_axis=True):
        if same_axis:
            plot = DistHist(1)
            ax = plot.axes[0, 0]
        else:
            plot = DistHist(len(filters))
            ax = None
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "pdb", **filter)
            plot.add_sample(structure, ij, label, ax)
        return plot

###############################################################################
# Plotting functions that accept a list of samples
#   extract_passthrough_kwargs
#   array_qc
#   array_skyline
#   array_ap
#   array_ss
#   array_mol
#   array_heatmap
#   array_circle
#   array_linreg
#   array_disthist
###############################################################################


def extract_passthrough_kwargs(plot, kwargs):
    pt_kwargs = {}
    for kw in plot.pass_through:
        if kw in kwargs.keys():
            pt_kwargs[kw] = kwargs.pop(kw)
    return pt_kwargs


def array_qc(samples=[], **kwargs):
    plot = QC(len(samples))
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        plot.add_sample(sample, log="log", profile="profile", label="label",
                        **pt_kwargs)
    return plot


def array_skyline(samples, annotations=[], plot_kwargs={}, **kwargs):
    plot = Skyline(len(samples), samples[0].data["profile"].length,
                   **plot_kwargs)
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        plot.add_sample(sample, profile="profile", annotations=annotations,
                        label="label", **pt_kwargs)
    return plot


def array_ap(samples, ct="ct", comp=None, ij=None, ij2=None, ij2_filter={},
             profile="profile", annotations=[], label="label",
             plot_kwargs={}, prefiltered=False, **kwargs):
    plot = AP(len(samples), samples[0].data[ct].length, **plot_kwargs)
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        if not prefiltered:
            if ct not in ["ss", "ct"]:
                sample.filter_ij(ct, ct)
            if ij is not None:
                sample.filter_ij(ij, ct, **kwargs)
            if ij2 is not None:
                sample.filter_ij(ij2, ct, **ij2_filter)
        plot.add_sample(sample, ct=ct, comp=comp, ij=ij, ij2=ij2,
                        profile=profile, label=label, annotations=annotations,
                        **pt_kwargs)
    return plot


def array_ss(samples, ss="ss", ij=None, ij2=None, ij2_filter={},
             profile="profile", annotations=[], label="label", plot_kwargs={},
             prefiltered=False, **kwargs):
    plot = SS(len(samples), samples[0].data[ss], **plot_kwargs)
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        if not prefiltered:
            if ij is not None:
                sample.filter_ij(ij, ss, **kwargs)
            if ij2 is not None:
                sample.filter_ij(ij2, ss, **ij2_filter)
        plot.add_sample(sample, ij=ij, ij2=ij2, profile=profile,
                        annotations=annotations, label=label, **pt_kwargs)
    return plot


def array_mol(samples, ij=None, profile="profile", label="label", show=True,
              prefiltered=False, **kwargs):
    plot = Mol(len(samples), samples[0].data["pdb"])
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        if not prefiltered:
            if ij is not None:
                sample.filter_ij(ij, "pdb", **kwargs)
        plot.add_sample(sample, ij=ij, profile=profile,
                        label=label, **pt_kwargs)
    if show:
        plot.view.show()
    return plot


def array_heatmap(samples, structure=None, ij=None, label="label", **kwargs):
    plot = Heatmap(len(samples), samples[0].data[structure])
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        sample.filter_ij(ij, structure, **kwargs)
        plot.add_sample(sample, ij=ij, label=label, **pt_kwargs)
    return plot


def array_circle(samples, ct=None, comp=None, ij=None, ij2=None, profile=None,
                 label="label", **kwargs):
    plot = Circle(len(samples), samples[0].data["profile"].length)
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        if ij is not None:
            sample.filter_ij(ij, "profile", **kwargs)
        plot.add_sample(sample, ct=ct, comp=comp, ij=ij, ij2=ij2,
                        profile=profile, label=label, **pt_kwargs)
    return plot


def array_linreg(samples, ct="ct", profile="profile", label="label", **kwargs):
    plot = LinReg(len(samples))
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        plot.add_sample(sample, ct=ct, profile=profile,
                        label=label, **pt_kwargs)
    return plot


def array_disthist(samples, structure="pdb", ij=None, label="label",
                   same_axis=False, **kwargs):
    if same_axis:
        plot = DistHist(1)
        ax = plot.axes[0, 0]
    else:
        plot = DistHist(len(samples))
        ax = None
    pt_kwargs = extract_passthrough_kwargs(plot, kwargs)
    for sample in samples:
        sample.filter_ij(ij, "pdb", **kwargs)
        plot.add_sample(sample, structure=structure, ij=ij,
                        label=label, ax=ax, **pt_kwargs)
