#!/usr/bin/env python

# general python packages
import os.path

# scripts in JNBTools
from .data import *
from .plots import *
from .analysis import *


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


class Sample():

    def __init__(self,
                 sample=None,
                 fasta=None,
                 profile=None,
                 rnp=None,
                 ct=None,
                 compct=None,
                 ss=None,
                 log=None,
                 rings=None,
                 deletions=None,
                 pairs=None,
                 allcorrs=None,
                 pdb=None,
                 pdb_kwargs={'chain': 'A'},
                 probs=None,
                 dance_prefix=None):
        """Creates a sample object which connects all chemical probing and
        structural data for a single experiment. Contains convenience methods
        to plot, filter, compare and retrieve this data.

        Args:
            sample (str, optional): This string will serve as a label for this
                sample within plots.
                Defaults to None.
            fasta (str, optional): Path to a fasta file which will provide the
                reference sequence for ShapeJumper Data.
                Defaults to None.
            profile (str, optional): Path to a ShapeMapper profile.txt output
                file.
                Defaults to None.
            rnp (str, optional): Path to an RNPMapper .csv file.
                Defaults to None.
            ct (str, optional): Path to a .ct, .dbn, or .bracket structure
                file.
                Defaults to None.
            compct (str, optional): Same as ct above, but serves as a
                comparison structure for plotting.
                Defaults to None.
            ss (str, optional): Path to a .xrna, .varna, .nsd, or .cte
                secondary structure drawing file.
                Defaults to None.
            log (str, optional): Path to a ShapeMapper shapemapper_log.txt
                output file.
                Defaults to None.
            rings (str, optional): Path to a RingMapper output file.
                Defaults to None.
            deletions (str, optional): Path to a ShapeJumper output file.
                Defaults to None.
            pairs (str, optional): Path to a PairMapper pairmap.txt output
                file.
                Defaults to None.
            allcorrs (str, optional): Path to a PairMapper allcorrs.txt output
                file.
                Defaults to None.
            pdb (str, optional): Path to a .pdb atom coordinates file.
                Defaults to None.
            pdb_kwargs (dict, optional): Dictionary containing info needed to
                parse provided .pdb. "chain" is required. "fasta" and "offset"
                may be required if not provided in header.
                Defaults to {'chain':'A'}.
            probs (dict, optional): Path to an RNAStructure .dp pairing
                probability file.
                Defaults to None.
            dance_prefix (str, optional): Path prefix for DanceMapper output
                files.
                Defaults to None.
        """
        self.paths = {"fasta": fasta,
                      "profile": profile,
                      "rnp": rnp,
                      "ct": ct,
                      "comptct": compct,
                      "ss": ss,
                      "log": log,
                      "rings": rings,
                      "deletions": deletions,
                      "pairs": pairs,
                      "allcorrs": allcorrs,
                      "pdb": pdb,
                      "probs": probs}
        self.sample = sample
        self.parent = None
        self.data = {}  # stores profile, ij, and structure objects
        if ct is not None:
            filetype = ct.split(".")[-1]
            self.data["ct"] = CT(filetype, ct)
        if compct is not None:
            filetype = ct.split(".")[-1]
            self.data["compct"] = CT(filetype, compct)
        if ss is not None:
            self.data["ss"] = CT("ss", ss)
        if pdb is not None:
            self.data["pdb"] = PDB(pdb, **pdb_kwargs)

        # ShapeMapper downstream analysis requires sequence given from profile
        if profile is not None:
            self.data["profile"] = Profile(profile)
            prof_seq = self.data["profile"].sequence
            self.has_profile = True
        else:
            self.has_profile = False
        no_profile_message = "{} requires a sequence from ShapeMapper profile."
        if rnp is not None:
            self.data["rnp"] = Profile(rnp, 'RNP')
        if log is not None:
            self.data["log"] = Log(log)
        if rings is not None:
            assert self.has_profile, no_profile_message.format("Rings")
            self.data["rings"] = IJ("rings", rings, prof_seq)
        if pairs is not None:
            assert self.has_profile, no_profile_message.format("Pairs")
            self.data["pairs"] = IJ("pairs", pairs, prof_seq)
        if allcorrs is not None:
            assert self.has_profile, no_profile_message.format("allcorrs")
            self.data["allcorrs"] = IJ("rings", allcorrs, prof_seq)
        if probs is not None:
            assert self.has_profile, no_profile_message.format("Probabilities")
            self.data["probs"] = IJ("probs", probs, prof_seq)

        # Deletions requires a reference sequence in fasta format.
        if deletions is not None:
            assert fasta is not None, "Deletions plotting requires fasta"
            self.data["deletions"] = IJ("deletions", deletions, fasta=fasta)
        if dance_prefix is not None:
            self.init_dance(dance_prefix)

    def init_dance(self, prefix):
        """Initializes a list of Sample objects which each represent a
        component of the DANCE model, available as Sample.dance.

        Args:
            prefix (str): Path to DanceMapper output file prefixes. Finds
                profiles, rings, pairs, and structure predictions if present.
        """
        reactivityfile = f"{prefix}-reactivities.txt"
        # read in 2 line header
        with open(reactivityfile) as inf:
            header1 = inf.readline().strip().split()
            header2 = inf.readline().strip().split()
        # number of components
        self.dance_components = int(header1[0])
        # population percentage of each component
        self.dance_percents = header2[1:]
        # dance is a list containing one sample for each component
        self.dance = [Sample() for _ in range(self.dance_components)]
        # build column names for reading in BM file
        for i, sample in enumerate(self.dance):
            sample.sample = f"{self.sample}: {i} - {self.dance_percents[i]}"
            sample.paths = {"profile": reactivityfile,
                            "rings": f"{prefix}-{i}-rings.txt",
                            "pairs": f"{prefix}-{i}-pairmap.txt",
                            "ct": [f"{prefix}-{i}.f.ct",  # if using --pk
                                   f"{prefix}-{i}.ct"]}  # if regular fold used
            # read in "profile" from reactivities
            sample.data["profile"] = Profile(reactivityfile, "dance", i)
            # read in other attributes
            if os.path.isfile(sample.paths["rings"]):
                sample.data["rings"] = IJ("rings", sample.paths["rings"],
                                          sample.data["profile"].sequence)
            if os.path.isfile(sample.paths["pairs"]):
                sample.data["pairs"] = IJ("pairs", sample.paths["pairs"],
                                          sample.data["profile"].sequence)
            # ! possible that these both exist
            for ct_file in sample.paths["ct"]:
                if os.path.isfile(ct_file):
                    sample.data["ct"] = CT("ct", ct_file)
                    sample.paths["ct"] = ct_file
            sample.parent = self

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
        elif key == None:
            return None
        try:
            return self.data[key]
        except (KeyError, TypeError) as e:
            try:
                return self.parent.data[key]
            except (KeyError, AttributeError):
                print(f"{key} data not found in {self.sample}")
                return None

    def get_data_list(self, *keys):
        """Calls Sample.get_data() on a list of keys, also accepts "ctcompare",
        which is equivalent to ["ct", "compct"]

        Returns:
            list of data objects
        """
        data_list = []
        for key in keys:
            if key == "ctcompare":
                data = self.get_data_list("ct", "compct")
                if data[1] is None:
                    data = data[0]
            elif key == "label":
                return self.sample
            elif isinstance(key, list):
                data = self.get_data_list(*key)
            else:
                data = self.get_data(key)
            data_list.append(data)
        return data_list

    def filter_ij(self, ij, fit_to, **kwargs):
        """Fits ij Data object to fit_to Data object's sequence. For example,
        for plotting IJ data containing structure cassettes against a CT
        without, or plotting sequence variants together using a best pairwise
        alignment. Kwargs are passed to IJ.filter(). Refer to that method for
        more detail on filters.

        Args:
            ij (IJ): IJ object to be filtered and fitted
            fit_to (Data (any)): Data object containing a sequence for IJ to be
                fit to for plotting purposes.
        """
        try:
            metric = kwargs.pop("metric")
            if metric == "Distance":
                self.data[ij].set_3d_distances(self.data["pdb"])
            self.data[ij].metric = metric
        except KeyError:
            self.data[ij].metric = self.data[ij].default_metric
        self.data[ij].filter(self.get_data(fit_to),
                             profile=self.get_data("profile"),
                             ct=self.get_data("ct"), **kwargs)

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
#     make_skyline
#     make_shapemapper
#     make_ap
#     make_ap_multifilter
#     make_ss
#     make_ss_multifilter
#     make_mol
#     make_mol_multifilter
#     make_heatmap
#     make_circle
#     make_circle_multifilter
###############################################################################

    def make_qc(self, **kwargs):
        """Makes a plot of ShapeMapper log data and reactivities.

        Returns:
            starmapper.plots.QC
        """
        return array_qc([self], **kwargs)

    def make_skyline(self, dance=False):
        if dance:
            plot = array_skyline(self.dance)
            plot.ax.legend(title="Comp: Percent")
            plot.ax.set_title(f"{self.sample}: DANCE Reactivities")
            return plot
        else:
            plot = array_skyline([self])
            return plot

    def make_shapemapper(self, plots=["profile", "rates", "depth"]):
        plot = SM(self.data["profile"].length, plots=plots)
        plot.add_sample(self, profile="profile", label="label")
        return plot

    def make_ap(self, ct="ct", comp="compct", ij=None, ij2=None,
                profile="profile", label="label", dance=False, **kwargs):
        if dance:
            plot = array_ap(self.dance, ct, comp, "pairs", "rings", profile,
                            label, **kwargs)
            for i, dance in enumerate(self.dance):
                ax = plot.get_ax(i)
                ax.set_title(
                    f"DANCE component: {i}, Percent: {self.dance_percents[i]}")
        return array_ap([self], ct, comp, ij, ij2, profile, label, **kwargs)

    def make_ap_multifilter(self, filters, ct="ct", comp="compct", ij2=None,
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
            starmapper.plots.AP
        """
        plot = AP(len(filters), self.get_data(ct).length)
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, ct, **filter)
            plot.add_sample(self, ct=ct, comp=comp, ij=ij, ij2=ij2,
                            profile=profile, label=label)
        return plot

    def make_ss(self, ss="ss", ij=None, profile="profile", label="label",
                **kwargs):
        return array_ss([self], ss, ij, profile, label, **kwargs)

    def make_ss_multifilter(self, filters, ss="ss", profile="profile",
                            label="label"):
        plot = SS(len(filters), self.get_data(ss))
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, ss, **filter)
            plot.add_sample(self, ij=ij, profile=profile, label=label)
        return plot

    def make_mol(self, ij=None, profile="profile", label="label", show=True,
                 **kwargs):
        return array_mol([self], ij, profile, label, show, **kwargs)

    def make_mol_multifilter(self, filters, profile="profile", label="label",
                             show=True):
        plot = Mol(len(filters), self.data["pdb"])
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "pdb", **filter)
            plot.add_sample(self, ij=ij, profile=profile, label=label)
        if show:
            plot.view.show()
        return plot

    def make_heatmap(self, structure=None, ij=None, levels=None, **kwargs):
        array_heatmap([self], structure, ij, levels, **kwargs)

    def make_circle(self, ct=None, comp=None, ij=None, ij2=None, profile=None,
                    label="label", **kwargs):
        array_circle([self], ct, comp, ij, ij2, profile, label, **kwargs)

    def make_circle_multifilter(self, filters, ct=None, comp=None, ij2=None,
                                profile=None, label="label"):
        plot = Circle(len(filters), self.data["profile"].length)
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "profile", **filter)
            plot.add_sample(self, ct=ct, comp=comp, ij=ij, ij2=ij2,
                            profile=profile, label=label)
        return plot

###############################################################################
# Plotting functions that accept a list of samples
#   array_qc
#   array_skyline
#   array_ap
#   array_ss
#   array_mol
#   array_heatmap
#   array_circle
#   array_linreg
###############################################################################


def array_qc(samples=[], **kwargs):
    plot = QC(len(samples))
    for sample in samples:
        plot.add_sample(sample, log="log", profile="profile", label="label",
                        **kwargs)
    return plot


def array_skyline(samples, **kwargs):
    plot = Skyline(len(samples), samples[0].data["profile"].length)
    for sample in samples:
        plot.add_sample(sample, profile="profile", label="label", **kwargs)
    return plot


def array_ap(samples, ct="ct", comp="compct", ij=None, ij2=None,
             profile="profile", label="label", **kwargs):
    plot = AP(len(samples), samples[0].get_data(ct).length)
    for sample in samples:
        if ij is not None:
            sample.filter_ij(ij, ct, **kwargs)
        plot.add_sample(sample, ct=ct, comp=comp, ij=ij, ij2=ij2,
                        profile=profile, label=label)
    return plot


def array_ss(samples, ss="ss", ij=None, profile="profile", label="label",
             nt_color="profile", markers="o", ** kwargs):
    plot = SS(len(samples), samples[0].data[ss])
    for sample in samples:
        if ij is not None:
            sample.filter_ij(ij, "ss", **kwargs)
        plot.add_sample(sample, ij=ij, profile=profile, label=label,
                        nt_color=nt_color, markers=markers)
    return plot


def array_mol(samples, ij=None, profile="profile", label="label", show=True,
              nt_color="sequence", **kwargs):
    plot = Mol(len(samples), samples[0].data["pdb"])
    for sample in samples:
        if ij is not None:
            sample.filter_ij(ij, "pdb", **kwargs)
        plot.add_sample(sample, ij=ij, profile=profile, label=label,
                        nt_color=nt_color)
    if show:
        plot.view.show()
    else:
        return plot


def array_heatmap(samples, structure=None, ij=None, levels=None, **kwargs):
    plot = Heatmap(len(samples))
    for sample in samples:
        sample.filter_ij(ij, ij, **kwargs)
        plot.add_sample(sample, structure=structure, ij=ij, levels=levels)
    return plot


def array_circle(samples, ct=None, comp=None, ij=None, ij2=None, profile=None,
                 label="label", **kwargs):
    plot = Circle(len(samples), samples[0].data["profile"].length)
    for sample in samples:
        if ij is not None:
            sample.filter_ij(ij, "profile", **kwargs)
        plot.add_sample(sample, ct=ct, comp=comp, ij=ij, ij2=ij2,
                        profile=profile, label=label)
    return plot


def array_linreg(samples, ct="ct", profile="profile", label="label", **kwargs):
    plot = LinReg(len(samples))
    for sample in samples:
        plot.add_sample(sample, ct=ct, profile=profile, label=label, **kwargs)
    return plot
