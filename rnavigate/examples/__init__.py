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

from pathlib import Path
import rnavigate as rnav


example_dir = Path(__file__).parent.absolute()
tpp_dir = example_dir / "tpp"
rnasep_dir = example_dir / "rnasep"
rmrp_dir = example_dir / "rmrp"
rrna_fragmap_dir = example_dir / "rrna_fragmap"

###############################################################################
# TPP Riboswitch
###############################################################################

tpp = rnav.Sample(
    sample="TPP in-vitro DMS-MaP",
    pdb={
        "pdb": str(tpp_dir / "2gdi.pdb"),
        "chain": "X"
    },
    ss=str(tpp_dir / "TPP-2GDI.nsd"),
    dmsmap=str(tpp_dir / "DMS_TPP_profile.txt"),
    ringmap=str(tpp_dir / "DMS_TPP_rings.txt"),
    pairprob=str(tpp_dir / "TPP-dms-bp.dp"),
    )

###############################################################################
# RNase P
###############################################################################

rnasep_common = rnav.Sample(
    sample="common data",
    pdb={
        "pdb": str(rnasep_dir / "3dhsCrystal_PlusLoops.pdb"),
        "chain": "A",
        },
    shapejump={
        "shapejump": str(rnasep_dir / "example-rnasep-deletions.txt"),
        "sequence": str(rnasep_dir / "RNaseP-noSC.fasta"),
        },
    ss_ct={"ss": str(rnasep_dir / "RNaseP.ct")},
    ss_pdb={"ss": str(rnasep_dir / "RC_CRYSTAL_STRUCTURE.xrna")},
    ss_lit={"ss": str(rnasep_dir / "RNaseP-lit-like.nsd")},
    pairprob={
        "pairprob": str(rnasep_dir / "rnasep.dp"),
        "sequence": str(rnasep_dir / "RNaseP-withSC.fasta"),
        },
)

rnasep_1 = rnav.Sample(
    sample="Example 1",
    inherit=rnasep_common,
    shapemap={
        "shapemap": str(rnasep_dir / "example1_rnasep_profile.txt"),
        "log": str(rnasep_dir / "example1_shapemapper_log.txt"),
        },
    pairmap=str(rnasep_dir / "example1-rnasep-pairmap.txt"),
    ringmap=str(rnasep_dir / "example1-rnasep.corrs"),
)

rnasep_2 = rnav.Sample(
    sample="Example 2",
    inherit=rnasep_common,
    shapemap={
        "shapemap": str(rnasep_dir / "example2_rnasep_profile.txt"),
        "log": str(rnasep_dir / "example2_shapemapper_log.txt"),
        },
    pairmap=str(rnasep_dir / "example2-rnasep-pairmap.txt"),
    ringmap=str(rnasep_dir / "example2-rnasep.corrs"),
)

rnasep_3 = rnav.Sample(
    sample="Example 3",
    inherit=rnasep_common,
    shapemap={
        "shapemap": str(rnasep_dir / "example3_rnasep_profile.txt"),
        "log": str(rnasep_dir / "example3_shapemapper_log.txt"),
        },
    pairmap=str(rnasep_dir / "example3-rnasep-pairmap.txt"),
    ringmap=str(rnasep_dir / "example3-rnasep.corrs"),
)

rnasep_4 = rnav.Sample(
    sample="Example 4",
    inherit=rnasep_common,
    shapemap={
        "shapemap": str(rnasep_dir / "example4_rnasep_profile.txt"),
        "log": str(rnasep_dir / "example4_shapemapper_log.txt"),
        },
    pairmap=str(rnasep_dir / "example4-rnasep-pairmap.txt"),
    ringmap=str(rnasep_dir / "example4-rnasep.corrs"),
)

###############################################################################
# RNase MRP
###############################################################################

rmrp = rnav.Sample(
    sample="RMRP in-vitro DMS-MaP",
    rnpmap=str(rmrp_dir / "RMRP-RNPMaP-Example_RESULTS.csv"),
    ss=str(rmrp_dir / "hs-RMRP.nsd"),
)

###############################################################################
# RNase MRP
###############################################################################

common_data = rnav.Sample(
    "6HA1 - LSU",
    sequence=str(rrna_fragmap_dir / "6HA1_LSU.fasta"),
    pdb={"pdb": str(rrna_fragmap_dir / "6HA1_LSU.pdb"),
         "sequence": "sequence",
         "chain": "A"},
    ss=str(rrna_fragmap_dir / "6HA1_LSU.json"),
)
common_data.data["ss"].transform_coordinates(rotate_degrees=180)

quinaxoline = rnav.Sample(
    sample="rRNA quinaxoline",
    inherit=common_data,
    shapemap=str(rrna_fragmap_dir / "2_QN_DMSO_subtracted_LSU_profile.txt"),
)

linezolid = rnav.Sample(
    sample="rRNA linezolid",
    inherit=common_data,
    shapemap=str(rrna_fragmap_dir / "2_ZLD_DMSO_subtracted_LSU_profile.txt"),
)

methyl = rnav.Sample(
    sample="rRNA methyl",
    inherit=common_data,
    shapemap=str(
        rrna_fragmap_dir / "2_Methyl_DMSO_subtracted_LSU_profile.txt"
        ),
)

__all__ = [
    "tpp",
    "rnasep_common",
    "rnasep_1",
    "rnasep_2",
    "rnasep_3",
    "rnasep_4",
    "rmrp",
    "quinaxoline",
    "linezolid",
    "methyl",
]
