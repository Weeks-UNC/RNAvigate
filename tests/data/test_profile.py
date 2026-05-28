"""Unit tests for rnavigate.data.Profile and SHAPEMaP."""

from __future__ import annotations

import numpy as np

from rnavigate import data


class TestProfileFromSample:
    """Tests using the tpp DMS-MaP profile (SHAPEMaP subclass)."""

    def test_has_sequence(self, tpp):
        profile = tpp.data["dmsmap"]
        assert isinstance(profile.sequence, str)
        assert len(profile.sequence) > 0

    def test_length_matches_sequence(self, tpp):
        profile = tpp.data["dmsmap"]
        assert profile.length == len(profile.sequence)

    def test_has_dataframe(self, tpp):
        profile = tpp.data["dmsmap"]
        assert profile.data is not None
        assert len(profile.data) == profile.length

    def test_dataframe_has_nucleotide_column(self, tpp):
        profile = tpp.data["dmsmap"]
        assert "Nucleotide" in profile.data.columns

    def test_is_shapemap_instance(self, tpp):
        profile = tpp.data["dmsmap"]
        assert isinstance(profile, data.SHAPEMaP)

    def test_metric_is_string(self, tpp):
        profile = tpp.data["dmsmap"]
        assert isinstance(profile.metric, str)

    def test_metric_column_exists(self, tpp):
        profile = tpp.data["dmsmap"]
        assert profile.metric in profile.data.columns

    def test_colors_length(self, tpp):
        profile = tpp.data["dmsmap"]
        colors, colormap = profile.get_colors("sequence")
        assert isinstance(colors, np.ndarray)
        assert isinstance(colormap, data.ScalarMappable)
        assert len(colors) == profile.length

    def test_get_aligned_data_to_self(self, tpp):
        profile = tpp.data["dmsmap"]
        aligned = profile.get_aligned_data(profile.null_alignment)
        assert aligned.length == profile.length


class TestProfileNormalization:
    def test_normalize_does_not_crash(self, tpp):
        profile = tpp.data["dmsmap"]
        profile.normalize()

    def test_winsorize_does_not_crash(self, tpp):
        profile = tpp.data["dmsmap"]
        profile.winsorize("Reactivity_profile", 96, 99)
