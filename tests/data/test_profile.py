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


class TestProfileFilter:
    """Exercises Profile.filter() and its constituent mask_on_* methods."""

    def test_filter(self, tpp):
        # start from an independent copy so the mask doesn't leak into the
        # session-scoped `tpp` fixture used by other tests
        profile = tpp.data["dmsmap"].copy()

        # exclude_nts/isolate_nts accept a mix of single positions and
        # (start, end) ranges (tuple or list)
        profile.filter(exclude_nts=[5, (10, 15), [20, 22]])
        excluded = set(profile.data.loc[~profile.data["mask"], "Nucleotide"])
        assert excluded == {5, 10, 11, 12, 13, 14, 15, 20, 21, 22}

        profile.filter(isolate_nts=[1, (3, 4)])
        kept = set(profile.data.loc[profile.data["mask"], "Nucleotide"])
        assert kept == {1, 3, 4}

        # reset_filter=False accumulates filters instead of replacing them
        profile.filter(exclude_nts=[3], reset_filter=False)
        kept = set(profile.data.loc[profile.data["mask"], "Nucleotide"])
        assert kept == {1, 4}

        # sequence= interprets positions in another sequence's coordinate
        # frame and translates them via alignment before masking
        shifted_sequence = data.Sequence(profile.sequence[10:])
        profile.filter(exclude_nts=[1, (2, 5)], sequence=shifted_sequence)
        excluded = set(profile.data.loc[~profile.data["mask"], "Nucleotide"])
        assert excluded == {11, 12, 13, 14, 15}

        # nts filters on nucleotide identity
        profile.filter(nts="AU")
        assert set(profile.data.loc[profile.data["mask"], "Sequence"].str.upper()) <= {
            "A",
            "U",
        }

        # ss_only/ds_only mask by base-pairing status from a structure
        structure = tpp.data["ss"]
        profile.filter(structure=structure, ss_only=True)
        assert profile.data["mask"].sum() > 0
        assert profile.data["mask"].sum() < profile.length

        # column-value keyword filters, e.g. "HQ_profile_ge"
        profile.filter(HQ_profile_ge=0)
        assert (profile.data.loc[profile.data["mask"], "HQ_profile"] >= 0).all()

        # prefiltered=True leaves the existing mask untouched
        profile.data["mask"] = False
        profile.filter(prefiltered=True, exclude_nts=[1])
        assert not profile.data["mask"].any()

        # masked positions are excluded from get_plotting_dataframe()
        profile.reset_mask()
        profile.filter(exclude_nts=[1])
        plotting_dataframe = profile.get_plotting_dataframe()
        assert (
            plotting_dataframe.loc[plotting_dataframe["Nucleotide"] == 1, "Values"]
            .isna()
            .all()
        )
