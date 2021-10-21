#!/usr/bin/env python

# general python packages
import seaborn as sns
import os.path


# scripts in JNBTools
from data import *
from plots2 import *


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


# STYLE SHEET
###############################################################################
sns.set_context("talk")
sns.set_style("ticks")
colors = [
    '#0092edff',  # Blue
    '#ff8300ff',  # Orange
    '#a100ffff',  # Purple
    '#edc600ff',  # Yellow
    '#ff48e9ff',  # Pink
    '#3fd125ff'  # Green
]
sns.set_palette(colors)


class Sample():

    def __init__(self,
                 sample=None,
                 fasta=None,
                 profile=None,
                 ct=None,
                 compct=None,
                 ss=None,
                 log=None,
                 rings=None,
                 deletions=None,
                 pairs=None,
                 pdb=None,
                 pdb_kwargs={},
                 probs=None,
                 dance_prefix=None):
        self.paths = {"fasta": fasta,
                      "profile": profile,
                      "ct": ct,
                      "comptct": compct,
                      "ss": ss,
                      "log": log,
                      "rings": rings,
                      "deletions": deletions,
                      "pairs": pairs,
                      "pdb": pdb,
                      "probs": probs}
        self.sample = sample
        self.parent = None
        self.data = {}  # stores profile, ij, and structure objects
        if ct is not None:
            self.data["ct"] = CT("ct", ct)
        if compct is not None:
            self.data["compct"] = CT("ct", compct)
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
        if log is not None:
            self.data["log"] = Log(log)
        if rings is not None:
            assert self.has_profile, no_profile_message.format("Rings")
            self.data["rings"] = IJ("rings", rings, prof_seq)
        if pairs is not None:
            assert self.has_profile, no_profile_message.format("Pairs")
            self.data["pairs"] = IJ("pairs", pairs, prof_seq)
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
                            "probs": f"{prefix}-{i}.dp",
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
            if os.path.isfile(sample.paths["probs"]):
                sample.data["probs"] = IJ("probs", sample.paths["probs"],
                                          sample.data["profile"].sequence)
            # ! possible that these both exist
            for ct_file in sample.paths["ct"]:
                if os.path.isfile(ct_file):
                    sample.data["ct"] = CT("ct", ct_file)
                    sample.paths["ct"] = ct_file
            sample.parent = self

###############################################################################
# some functions
#     get_data
#     filter_ij
#     filter_dance_rings
###############################################################################

    def get_data(self, key):
        try:
            return self.data[key]
        except KeyError:
            try:
                return self.parent.data[key]
            except (KeyError, AttributeError):
                print(f"{key} data not found in {self.sample}")
                return None

    def get_data_list(self, *keys):
        data_list = []
        for key in keys:
            if key == "ctcompare":
                data = self.get_data_list("ct", "compct")
                if data[1] is None:
                    data = data[0]
            elif isinstance(key, list):
                data = self.get_data_list(*key)
            else:
                data = self.get_data(key)
            data_list.append(data)
        return data_list

    def filter_ij(self, ij, fit_to, **kwargs):
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

    def make_qc(self, **kwargs):
        profiles = [self.data["profile"]]
        logs = [self.data["log"]]
        labels = [self.sample]
        QC(logs, profiles, labels).make_plot(**kwargs)

    def make_skyline(self, dance=False, **kwargs):
        profiles, labels = [], []
        if dance:
            for dance in self.dance:
                profiles.append(dance.data["profile"])
                labels.append(dance.sample)
            kwargs["legend_title"] = "Comp: Percent"
            kwargs["axis_title"] = f"{self.sample}: DANCE Reactivities"
        else:
            profiles.append(self.data["profile"])
            labels.append(self.sample)
        Skyline(profiles, labels).make_plot(**kwargs)

    def make_shapemapper(self, **kwargs):
        SM(self.data["profile"], self.sample).make_plot(**kwargs)

    def make_heatmap(self, heatmap, contour, metric=None):
        heatmap = self.get_data(heatmap)
        contour = self.get_data(contour)
        heatmap.metric = metric
        Heatmap(heatmap, contour).make_plot()

    def make_ap(self, ij, dance=False, **filter_kwargs):
        plot = AP()
        if dance:
            for sample in self.dance:
                sample.filter_ij(ij, "ct", **filter_kwargs)
                plot.add_sample(*sample.get_data_list("ctcompare", ij, "profile"),
                                sample.sample)
        else:
            self.filter_ij(ij, "ct", **filter_kwargs)
            plot.add_sample(*self.get_data_list("ctcompare", ij, "profile"),
                            self.sample)
        return plot.make_plot()

    def make_ap_multifilter(self, filters):
        plot = AP()
        try:
            profile = self.data["profile"]
        except KeyError:
            profile = None
        structure = self.data["ct"]
        if "compct" in self.data.keys():
            structure = [structure, self.data["compct"]]
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "pdb", **filter)
            plot.add_sample(structure, self.data[ij], profile, self.sample)
        return plot.make_plot()

    def make_ss(self, ij, dance=False, label=None, **filter_kwargs):
        plot = SS()
        if label is None:
            label = self.sample
        if dance:
            for sample in self.dance:
                sample.filter_ij(ij, "ss", **filter_kwargs)
                plot.add_sample(sample, "ss", ij, "profile", label)
        else:
            self.filter_ij(ij, "ss", **filter_kwargs)
            plot.add_sample(self, "ss", ij, "profile", label)
        plot.make_plot()

    def make_ss_multifilter(self, filters):
        plot = SS()
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "ss", **filter)
            plot.add_sample(self, "ss", ij, "profile", self.sample)
        return plot.make_plot()

    def make_mol(self, ij, dance=False, **filter_kwargs):
        plot = Mol(self.data["pdb"])
        if dance:
            for sample in self.dance:
                sample.filter_ij(ij, "pdb", **filter_kwargs)
                plot.add_sample(sample.data[ij], sample.data["profile"],
                                sample.sample)
        else:
            try:
                profile = self.data["profile"]
            except KeyError:
                profile = None
            self.filter_ij(ij, "pdb", **filter_kwargs)
            plot.add_sample(self.data[ij], profile, self.sample)
        return plot.make_plot()

    def make_mol_multifilter(self, filters):
        plot = Mol(self.data["pdb"])
        try:
            profile = self.data["profile"]
        except KeyError:
            profile = None
        for filter in filters:
            ij = filter.pop("ij")
            self.filter_ij(ij, "pdb", **filter)
            plot.add_sample(self.data[ij], profile, self.sample)
        return plot.make_plot()

###############################################################################
# Plotting functions that accept a list of samples
#   array_qc
#   array_skyline
#   array_ap
#   array_ss
#   array_3d
###############################################################################


def array_qc(samples=[]):
    logs, profiles, labels = [], [], []
    for sample in samples:
        logs.append(sample.data["log"])
        profiles.append(sample.data["profile"])
        labels.append(sample.sample)
    QC(logs, profiles, labels).make_plot()


def array_skyline(samples):
    profiles, labels = [], []
    for sample in samples:
        profiles.append(sample.data["profile"])
        labels.append(sample.sample)
    Skyline(profiles, labels).make_plot()


def array_ap(samples, ij=None, **kwargs):
    plot = AP()
    for sample in samples:
        structure = sample.get_data("ctcompare")
        sample.filter_ij(ij, "ct", **kwargs)
        plot.add_sample(structure, sample.data[ij],
                        sample.data["profile"], sample.sample)
    plot.make_plot()


def array_ss(samples, ij, **kwargs):
    structures, ijs, profiles, labels = [], [], [], []
    for sample in samples:
        sample.filter_ij(ij, "ss", **kwargs)
        structures.append(sample.data["ss"])
        ijs.append(sample.data[ij])
        profiles.append(sample.data["profile"])
        labels.append(sample.sample)
    SS(structures, labels, ijs, profiles).make_plot()


def array_mol(samples, ij, **kwargs):
    plot = Mol(samples[0].data["pdb"])
    for sample in samples:
        sample.filter_ij(ij, "pdb", **kwargs)
        plot.add_sample(sample.data[ij], sample.data["profile"], sample.sample)
    return plot.make_plot()
