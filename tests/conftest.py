"""Shared pytest fixtures for the RNAvigate test suite.

Matplotlib must switch to the non-interactive Agg backend before any pyplot
import.  Setting it here at module load guarantees that order regardless of
which test file is collected first.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from rnavigate import examples


@pytest.fixture(scope="session")
def tpp():
    """TPP riboswitch DMS-MaP sample (pdb, ss, dmsmap, ringmap, pairprob)."""
    return examples.tpp


@pytest.fixture(scope="session")
def rnasep_1():
    """RNase P SHAPE-MaP sample 1 (shapemap with log, ringmap, pairmap, pdb, ss_*)."""
    return examples.rnasep_1


@pytest.fixture(scope="session")
def rnasep_2():
    """RNase P SHAPE-MaP sample 2 (shapemap with log, ringmap, pairmap, pdb, ss_*)."""
    return examples.rnasep_2


@pytest.fixture(autouse=True)
def close_figures():
    """Close all matplotlib figures after every test to prevent memory leaks."""
    yield
    plt.close("all")
