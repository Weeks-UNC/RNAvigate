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


def test_resolve_data_dict_does_not_alias_or_skip_alignment(tpp):
    """Two dict-form interactions arguments pointing at the same data keyword
    (e.g. `interactions={"interactions": "ringmap", "Sign_eq": 1, ...}` and
    `interactions2={"interactions": "ringmap", "Sign_eq": -1, ...}`) must
    resolve to independent objects with independently correct filters, and
    must still be aligned (via `fit_data`) to the plot's target sequence
    rather than silently skipping alignment because the input was a dict.
    """
    sequence = rnav.data.Sequence("A" + tpp.data["ringmap"].sequence)
    alignment = data.SequenceAlignment(tpp.data["ringmap"], sequence)
    original_mask = tpp.data["ringmap"].data["mask"].copy()

    positive = rnav.resolve_data(
        {"interactions": "ringmap", "Sign_eq": 1},
        data.Interactions,
        sample=tpp,
        alignment=alignment,
    )
    negative = rnav.resolve_data(
        {"interactions": "ringmap", "Sign_eq": -1},
        data.Interactions,
        sample=tpp,
        alignment=alignment,
    )

    # independent objects: filtering one must not affect the other
    assert positive is not negative
    assert (positive.data["Sign"] == 1).all()
    assert (negative.data["Sign"] == -1).all()

    # the original sample data must be untouched by either resolution
    assert (tpp.data["ringmap"].data["mask"] == original_mask).all()

    # dict-form input must still be aligned to the target sequence (shifted
    # by 1 nt here), not silently skipped
    assert positive.length == sequence.length
    assert positive.data["i"].min() >= 2


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
