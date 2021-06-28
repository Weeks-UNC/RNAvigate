#!/usr/bin/env python

# general python packages
import matplotlib as mp
import matplotlib.pyplot as plt
from py3Dmol import view
import seaborn as sns
import numpy as np
import os.path


# scripts in JNBTools
from data import *
from plots import *


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

# Hard coded defaults
###############################################################################


def get_nt_color(nt, colors="new"):
    nt_color = {"old": {"A": "#f20000",  # red
                        "U": "#f28f00",  # yellow
                        "G": "#00509d",  # blue
                        "C": "#00c200"},  # green
                "new": {"A": "#366ef0",  # blue
                        "U": "#9bb9ff",  # light blue
                        "G": "#f04c4c",  # red
                        "C": "#ffa77c"}  # light red
                }[colors][nt]
    return nt_color


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

        self.data = {}  # stores profile, ij, and structure objects
        if ct is not None:
            self.data["ct"] = CT("ct", ct)
        if compct is not None:
            self.data["compct"] = CT("ct", compct)
        if ss is not None:
            self.data["ss"] = CT("ss", ss)
        if pdb is not None:
            self.data["pdb"] = PDB(pdb)

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
            sample.sample = f"{i} - {self.dance_percents[i]}"
            sample.paths = {"profile": reactivityfile,
                            "rings": f"{prefix}-{i}-rings.txt",
                            "pairs": f"{prefix}-{i}-pairmap.txt",
                            "ct": [f"{prefix}-{i}.f.ct",  # if using --pk
                                   f"{prefix}-{i}.ct"]}  # if regular fold used
            # read in "profile" from reactivities
            sample.data["profile"] = Profile(reactivityfile, "dance", i)
            # read in other attributes
            if os.path.isfile(sample.paths["rings"]):
                sample.data["rings"] = IJ(sample.paths["rings"], "rings",
                                          sample.data["profile"].sequence)
            if os.path.isfile(sample.paths["pairs"]):
                sample.data["pairs"] = IJ(sample.paths["pairs"], "pairs",
                                          sample.data["profile"].sequence)
            # ! possible that these both exist
            for ct_file in sample.paths["ct"]:
                if os.path.isfile(ct_file):
                    sample.data["ct"] = CT("ct", ct_file)
                    sample.paths["ct"] = ct_file

###############################################################################
# some functions
#     get_data
#     filter_ij
#     filter_dance_rings
###############################################################################

    def get_data(self, key):
        if key == "ctcompare":
            return [self.data["ct"], self.data["compct"]]
        elif key in self.data.keys():
            return self.data[key]
        elif isinstance(key, list):
            return [self.data[k] for k in key]
        else:
            print(f"Key must be one of:\n{self.data.keys()}")

    def filter_ij(self, ij, fit_to, **kwargs):
        try:
            metric = kwargs.pop("metric")
            if metric == "Distance":
                self.data[ij].set_3d_distances(self.data["pdb"])
            self.data[ij].metric = metric
        except KeyError:
            self.data[ij].metric = self.data[ij].default_metric
        self.data[ij].filter(self.data[fit_to], profile=self.data["profile"],
                             ct=self.data["ct"], **kwargs)

    def filter_dance_rings(self, filterneg=True, cdfilter=15, sigfilter=23,
                           ssfilter=True):
        ctlist = [dance.ct for dance in self.dance]
        ringlist = [dance.ij_data["rings"].copy() for dance in self.dance]
        for index, rings in enumerate(ringlist):
            rings[['i_offset', 'j_offset']] = rings[['i', 'j']]
            mask = []
            for _, row in rings.iterrows():
                i, j, G, sign = row[['i', 'j', 'Statistic', '+/-']]
                true_so_far = True
                if ssfilter and true_so_far:
                    true_so_far = (ctlist[index].ct[i-1]
                                   == 0 and ctlist[index].ct[j-1] == 0)
                if filterneg and true_so_far:
                    true_so_far = sign == 1
                if sigfilter is not None and true_so_far:
                    true_so_far = G > sigfilter
                if cdfilter is not None and true_so_far:
                    for ct in ctlist:
                        if not true_so_far:
                            break
                        true_so_far = ct.contactDistance(i, j) >= cdfilter
                mask.append(true_so_far)
            rings['mask'] = mask
            self.dance[index].ij_data["rings"] = rings

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

        def add_sample(sample):
            structure = sample.data["ct"]
            if "compct" in sample.data.keys():
                structure = [structure, sample.data["compct"]]
            sample.filter_ij(ij, "ct", **filter_kwargs)
            plot.add_sample(structure, sample.data[ij],
                            sample.data["profile"], sample.sample)
        if dance:
            for sample in self.dance:
                add_sample(sample)
        else:
            add_sample(self)
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

    def make_ss(self, ij, dance=False, **filter_kwargs):
        def add_sample(sample):
            ss_kwargs["structures"].append(self.data["ss"])
            ss_kwargs["ijs"].append(sample.data[ij])
            ss_kwargs["profiles"].append(sample.data["profile"])
            ss_kwargs["labels"].append(sample.sample)

        ss_kwargs = {"structures": [], "ijs": [], "profiles": [], "labels": []}
        if dance:
            for sample in self.dance:
                sample.filter_ij(ij, "ss", **filter_kwargs)
                add_sample(sample)
        else:
            self.filter_ij(ij, "ss", **filter_kwargs)
            add_sample(self)
        SS(**ss_kwargs).make_plot()

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
    logs, profiles, labels = []
    for sample in samples:
        logs.append(sample.data["logs"])
        profiles.append(sample.data["profile"])
        labels.append(sample.sample)
    QC(logs, profiles, labels).make_plot()


def array_skyline(samples):
    profiles, labels = []
    for sample in samples:
        profiles.append(sample.data["profile"])
        labels.append(sample.sample)
    Skyline(profiles, labels).make_plot()


def array_ap(samples, **kwargs):
    plot = AP()
    for sample in samples:
        structure = sample.data["ct"]
        if "compct" in sample.data.keys():
            structure = [structure, sample.data["compct"]]
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
        profiles.append(sample.data["profiles"])
        labels.append(sample.sample)
    SS(structures, labels, ijs, profiles).make_plot()


def array_mol(samples, ij, **kwargs):
    plot = Mol(samples[0].data["pdb"])
    for sample in samples:
        sample.filter_ij(ij, "pdb", **kwargs)
        plot.add_sample(sample.data[ij], sample.data["profile"], sample.sample)
    return plot.make_plot()
