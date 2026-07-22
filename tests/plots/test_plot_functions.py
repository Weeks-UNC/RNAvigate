"""Smoke tests for all rnavigate plot_*() convenience functions.

Each test calls the function with minimal valid arguments from the bundled
example datasets and asserts:
  1. The call does not raise an exception.
  2. The return value is a rnavigate.plots.Plot instance.

plot_mol is excluded because it requires a Jupyter / py3Dmol rendering
environment that is not available in headless CI.
"""

from __future__ import annotations

import rnavigate as rnav
from rnavigate import plots

# ---------------------------------------------------------------------------
# Per-nucleotide profile plots
# ---------------------------------------------------------------------------


def test_plot_profile(tpp):
    result = rnav.plot_profile([tpp], profile="dmsmap")
    assert isinstance(result, plots.Profile)


def test_plot_skyline(tpp):
    result = rnav.plot_skyline([tpp], profile="dmsmap")
    assert isinstance(result, plots.Skyline)


def test_plot_shapemapper(rnasep_1):
    result = rnav.plot_shapemapper(rnasep_1, profile="shapemap")
    assert isinstance(result, plots.SM)


def test_plot_qc(rnasep_1):
    result = rnav.plot_qc([rnasep_1], profile="shapemap")
    assert isinstance(result, plots.QC)


def test_plot_ntdist(tpp):
    result = rnav.plot_ntdist([tpp], profile="dmsmap")
    assert isinstance(result, plots.NucleotideDistribution)


def test_plot_roc(tpp):
    result = rnav.plot_roc([tpp], structure="ss", profile="dmsmap")
    assert isinstance(result, plots.ROC)


# ---------------------------------------------------------------------------
# Linear regression
# ---------------------------------------------------------------------------


def test_plot_linreg(rnasep_1, rnasep_2):
    result = rnav.plot_linreg([rnasep_1, rnasep_2], profile="shapemap")
    assert isinstance(result, plots.LinReg)


# ---------------------------------------------------------------------------
# Arc plots
# ---------------------------------------------------------------------------


def test_plot_arcs(tpp):
    result = rnav.plot_arcs(
        samples=[tpp],
        sequence="dmsmap",
        structure="ss",
        interactions="ringmap",
        profile="dmsmap",
        profile2="dmsmap",
        profile_scale_factor=[2, 0.5],
    )
    assert isinstance(result, plots.AP)


def test_plot_arcs_compare(rnasep_1, rnasep_2):
    result = rnav.plot_arcs_compare(
        samples=[rnasep_1, rnasep_2], sequence="ss_pdb", structure="ss_pdb"
    )
    assert isinstance(result, plots.AP)


# ---------------------------------------------------------------------------
# Secondary structure
# ---------------------------------------------------------------------------


def test_plot_ss(tpp):
    result = rnav.plot_ss([tpp], structure="ss", profile="dmsmap")
    assert isinstance(result, plots.SS)


# ---------------------------------------------------------------------------
# Circle plot
# ---------------------------------------------------------------------------


def test_plot_circle(tpp):
    result = rnav.plot_circle([tpp], sequence="ss", structure="ss")
    assert isinstance(result, plots.Circle)


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------


def test_plot_heatmap(tpp):
    result = rnav.plot_heatmap(
        [tpp], sequence="ss", structure="ss", interactions="ringmap"
    )
    assert isinstance(result, plots.Heatmap)


# ---------------------------------------------------------------------------
# Distance histogram (contact distances via secondary structure)
# ---------------------------------------------------------------------------


def test_plot_disthist(rnasep_1):
    result = rnav.plot_disthist([rnasep_1], structure="ss_ct", interactions="ringmap")
    assert isinstance(result, plots.DistHist)


# ---------------------------------------------------------------------------
# Alignment visualization
# ---------------------------------------------------------------------------


def test_plot_alignment(rnasep_1, rnasep_2):
    result = rnav.plot_alignment(
        data1=(rnasep_1, "ss_pdb"),
        data2=(rnasep_2, "ss_pdb"),
    )
    assert isinstance(result, plots.Alignment)
