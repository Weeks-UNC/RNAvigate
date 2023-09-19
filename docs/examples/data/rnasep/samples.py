import os
import rnavigate as rnav


rnasep_dir = os.path.dirname(os.path.abspath(__file__))


def abspath(filepath):
    """convert a filepath in this directory to an absolute path"""
    return os.path.join(rnasep_dir, filepath)


other_files = {
    "og_pdb": "3dhs_Correct.pdb",
    "ct": "RNaseP_2.ct",
    "dbn": "RNaseP.dbn"}

common_data = rnav.Sample(
    pdb={
        "filepath": abspath("3dhsCrystal_PlusLoops.pdb"),
        "chain": "A"
    },
    shapejump={
        "filepath": abspath("example-rnasep-deletions.txt"),
        "fasta": abspath("RNaseP-noSC.fasta")
    },
    ss=abspath("RC_CRYSTAL_STRUCTURE.xrna"),
    ct=abspath("RNaseP.ct"),
)

example1 = rnav.Sample(
    sample="Example 1",
    inherit=common_data,
    log=abspath("example1_shapemapper_log.txt"),
    shapemap=abspath("example1_rnasep_profile.txt"),
    pairmap=abspath("example1-rnasep-pairmap.txt"),
    ringmap=abspath("example1-rnasep.corrs"),
    dance_prefix=abspath("example1_rnasep"),
    pairprob=abspath("rnasep.dp")
)

example2 = rnav.Sample(
    sample="Example 2",
    inherit=common_data,
    log=abspath("example2_shapemapper_log.txt"),
    shapemap=abspath("example2_rnasep_profile.txt"),
    pairmap=abspath("example2-rnasep-pairmap.txt"),
    ringmap=abspath("example2-rnasep.corrs"),
    dance_prefix=abspath("example2_rnasep"),
    pairprob=abspath("rnasep.dp")
)

example3 = rnav.Sample(
    sample="Example 3",
    inherit=common_data,
    log=abspath("example3_shapemapper_log.txt"),
    shapemap=abspath("example3_rnasep_profile.txt"),
    pairmap=abspath("example3-rnasep-pairmap.txt"),
    ringmap=abspath("example3-rnasep.corrs"),
    dance_prefix=abspath("example3_rnasep"),
    pairprob=abspath("rnasep.dp")
)

example4 = rnav.Sample(
    sample="Example 4",
    inherit=common_data,
    log=abspath("example4_shapemapper_log.txt"),
    shapemap=abspath("example4_rnasep_profile.txt"),
    pairmap=abspath("example4-rnasep-pairmap.txt"),
    ringmap=abspath("example4-rnasep.corrs"),
    dance_prefix=abspath("example4_rnasep"),
    pairprob=abspath("rnasep.dp")
)
