"""This module contains example rnav.Samples for different RNAs and data types.

TPP riboswitch in-vitro construct DMS-MaP data:
    All sequences are based on PDB entry 3DHS.
    Except for the PDB file, all other data also contains structure cassettes.
    Data keywords:
        pdb - 2GDI from the PDB
        ss - structure drawing diagram
        dmsmap - DMS-MaP structure probing data
        ringmap - RING-MaP correlations from DMS-MaP experiment
        pairprob - pairing probabilities using DMS-MaP as restraints

RNase P in-vitro construct SHAPE-MaP data
    There are five samples here:
    common reference data shared between all RNase P samples
    variable name (rnasep_common)
    data keywords:
        pdb - 3DHS, but with missing loops modelled in.
        shapejump - SHAPE-JuMP data
        ss_ct - 3DHS + cassettes, base-pairing structure with singlets removed
        ss_pdb - 3DHS, full base-pairing structure diagram
        pairprob - 3DHS + cassettes, pairing probabilities (sequence only)
    4 experimental samples: 2 SHAPE reagents * 2 folding conditions
    variable names (rnasep_1, rnasep_2, rnasep_3, rnasep_4)
    data keywords:
        all rnasep_common data keywords above +
        shapemap -  3DHS + cassettes, SHAPE-MaP experimental data
        ringmap - RING-MaP analysis on SHAPE-MaP data
        pairmap - PAIR-MaP analysis on SHAPE-MaP data

Example usage:

import rnavigate as rnav
from rnavigate.examples import rnasep_1, rnasep_2

plot = rnav.plot_arcs(
    samples=[rnasep_1, rnasep_2],
    sequence="ss_pdb",
    structure="ss_pdb",
    profile="shapemap",
    interactions="ringmap",
    )
"""
from rnavigate.examples.examples import (
    tpp,
    rnasep_1,
    rnasep_2,
    rnasep_3,
    rnasep_4,
    linezolid,
    quinaxoline,
    methyl,
    )

__all__ = [
    "tpp",
    "rnasep_1",
    "rnasep_2",
    "rnasep_3",
    "rnasep_4",
    "linezolid",
    "quinaxoline",
    "methyl",
]
