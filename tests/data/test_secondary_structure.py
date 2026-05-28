"""Unit tests for rnavigate.data.SecondaryStructure."""

from __future__ import annotations

import numpy as np

from rnavigate import data


class TestSecondaryStructure:
    def test_has_sequence(self, tpp):
        ss = tpp.data["ss"]
        assert isinstance(ss.sequence, str)
        assert len(ss.sequence) > 0

    def test_is_secondary_structure(self, tpp):
        ss = tpp.data["ss"]
        assert isinstance(ss, data.SecondaryStructure)

    def test_has_pair_nts(self, tpp):
        ss = tpp.data["ss"]
        assert ss.pair_nts is not None
        assert len(ss.pair_nts) == ss.length

    def test_has_drawing_coordinates(self, tpp):
        ss = tpp.data["ss"]
        assert "X_coordinate" in ss.data.columns
        assert "Y_coordinate" in ss.data.columns

    def test_length_matches_sequence(self, tpp):
        ss = tpp.data["ss"]
        assert ss.length == len(ss.sequence)

    def test_pair_nts_dtype_is_numeric(self, tpp):
        ss = tpp.data["ss"]
        assert ss.pair_nts.dtype in (np.int32, np.int64, np.float64, int)

    def test_get_aligned_data_to_self(self, tpp):
        ss = tpp.data["ss"]
        aligned = ss.get_aligned_data(ss.null_alignment)
        assert aligned.length == ss.length

    def test_colors_by_structure(self, tpp):
        ss = tpp.data["ss"]
        colors, colormap = ss.get_colors(source="structure", structure=ss)
        assert isinstance(colors, np.ndarray)
        assert isinstance(colormap, data.ScalarMappable)
        assert len(colors) == ss.length
