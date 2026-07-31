"""Tests for rnavigate.helper_functions.resolve_data."""

from __future__ import annotations

import rnavigate as rnav
from rnavigate import data


def test_resolve_data_list_of_keywords(tpp):
    """A list of data keywords must resolve element-wise, matching the
    behavior of Sample.get_data, rather than being passed through unresolved.
    """
    sample = rnav.Sample(
        sample="test",
        inherit=tpp,
        motif1={
            "motif": "DRACH",
            "sequence": "dmsmap",
            "name": "motif1",
            "color": "red",
        },
        motif2={
            "motif": "AUCG",
            "sequence": "dmsmap",
            "name": "motif2",
            "color": "blue",
        },
    )

    resolved = rnav.resolve_data(["motif1", "motif2"], data.Annotation, sample=sample)

    assert isinstance(resolved, list)
    assert resolved == [sample.data["motif1"], sample.data["motif2"]]
    assert all(isinstance(each, data.Annotation) for each in resolved)


def test_plot_skyline_with_annotations_list(tpp):
    """annotations=[...] (a list of data keywords) must resolve to Data
    objects end-to-end through a plot_*() call, not remain raw strings.
    """
    sample = rnav.Sample(
        sample="test",
        inherit=tpp,
        motif1={
            "motif": "DRACH",
            "sequence": "dmsmap",
            "name": "motif1",
            "color": "red",
        },
        motif2={
            "motif": "AUCG",
            "sequence": "dmsmap",
            "name": "motif2",
            "color": "blue",
        },
    )

    result = rnav.plot_skyline(
        [sample], profile="dmsmap", annotations=["motif1", "motif2"]
    )

    assert isinstance(result, rnav.plots.Skyline)
