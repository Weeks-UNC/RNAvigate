"""This module contains example RNAvigate samples for RNaseP data.

All sequences are based on PDB entry 3DHS.

Depending on the data, the sequence may also contain structure cassettes.

There are five samples here:

common_data: data shared between all samples.
    pdb - 3DHS, but with missing loops modelled in.
    shapejump - SHAPE-JuMP data
    ss_ct - The 3DHS base-pairing structure with singlets removed
    ss_pdb - The 3DHS full base-pairing structure, layout reflects 3D structure
    pairprob - pairing probabilities, based on 3DHS + cassettes sequence only

example1, example2, example3, and example4:
    4 SHAPE-MaP experiments with 2 probes under 2 conditions. Each contains:
        all common data keywords
        shapemap - SHAPE-MaP experimental data on 3DHS + cassettes
        ringmap - RING-MaP analysis on SHAPE-MaP data
        pairmap - PAIR-MaP analysis on SHAPE-MaP data

Example usage:

import rnavigate as rnav
from rnavigate.examples.rnasep import samples

plot = rnav.plot_arcs(
    samples=[samples.example1, samples.example2],
    sequence="ct",
    structure="ct",
    interactions="ringmap",
    )

"""

import os
import rnavigate as rnav

rnasep_dir = os.path.dirname(os.path.abspath(__file__))


def abspath(filepath):
    """convert a filepath in this directory to an absolute path"""
    return os.path.join(rnasep_dir, filepath)


common_data = rnav.Sample(
    sample="common data",
    pdb={
        "pdb": abspath("3dhsCrystal_PlusLoops.pdb"),
        "chain": "A"
    },
    shapejump={
        "shapejump": abspath("example-rnasep-deletions.txt"),
        "sequence": abspath("RNaseP-noSC.fasta")
    },
    ss_ct={"ss": abspath("RNaseP.ct")},
    ss_pdb={"ss": abspath("RC_CRYSTAL_STRUCTURE.xrna")},
    pairprob={
        "pairprob": abspath("rnasep.dp"),
        "sequence": abspath("RNaseP-withSC.fasta")}
)

example1 = rnav.Sample(
    sample="Example 1",
    inherit=common_data,
    shapemap={
        "shapemap": abspath("example1_rnasep_profile.txt"),
        "log": abspath("example1_shapemapper_log.txt")},
    pairmap=abspath("example1-rnasep-pairmap.txt"),
    ringmap=abspath("example1-rnasep.corrs"),
)

example2 = rnav.Sample(
    sample="Example 2",
    inherit=common_data,
    shapemap={
        "shapemap": abspath("example2_rnasep_profile.txt"),
        "log": abspath("example2_shapemapper_log.txt")},
    pairmap=abspath("example2-rnasep-pairmap.txt"),
    ringmap=abspath("example2-rnasep.corrs"),
)

example3 = rnav.Sample(
    sample="Example 3",
    inherit=common_data,
    shapemap={
        "shapemap": abspath("example3_rnasep_profile.txt"),
        "log": abspath("example3_shapemapper_log.txt")},
    pairmap=abspath("example3-rnasep-pairmap.txt"),
    ringmap=abspath("example3-rnasep.corrs"),
)

example4 = rnav.Sample(
    sample="Example 4",
    inherit=common_data,
    shapemap={
        "shapemap": abspath("example4_rnasep_profile.txt"),
        "log": abspath("example4_shapemapper_log.txt")},
    pairmap=abspath("example4-rnasep-pairmap.txt"),
    ringmap=abspath("example4-rnasep.corrs"),
)
