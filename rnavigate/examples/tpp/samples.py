import os
import rnavigate as rnav


rnasep_dir = os.path.dirname(os.path.abspath(__file__))


def abspath(filepath):
    """convert a filepath in this directory to an absolute path"""
    return os.path.join(rnasep_dir, filepath)

example = rnav.Sample(
    sample="Example 1",
    pdb={
        "filepath": abspath("2gdi.pdb"),
        "chain": "X"
    },
    ss=abspath("TPP-2GDI.nsd"),
    shapemap=abspath("DMS_TPP_profile.txt"),
    ringmap=abspath("DMS_TPP_rings.txt"),
    pairprob=abspath("TPP-dms-bp.dp")
    )
