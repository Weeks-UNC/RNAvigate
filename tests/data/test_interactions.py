"""Unit tests for rnavigate.data.Interactions."""

from __future__ import annotations

from rnavigate import data


class TestInteractions:
    def test_has_sequence(self, tpp):
        ring = tpp.data["ringmap"]
        assert isinstance(ring.sequence, str)
        assert len(ring.sequence) > 0

    def test_is_interactions(self, tpp):
        ring = tpp.data["ringmap"]
        assert isinstance(ring, data.Interactions)

    def test_has_dataframe(self, tpp):
        ring = tpp.data["ringmap"]
        assert ring.data is not None

    def test_dataframe_has_i_j_columns(self, tpp):
        ring = tpp.data["ringmap"]
        assert "i" in ring.data.columns
        assert "j" in ring.data.columns

    def test_length_matches_sequence(self, tpp):
        ring = tpp.data["ringmap"]
        assert ring.length == len(ring.sequence)

    def test_metric_is_string(self, tpp):
        ring = tpp.data["ringmap"]
        assert isinstance(ring.metric, str)

    def test_reset_mask_does_not_crash(self, tpp):
        ring = tpp.data["ringmap"]
        ring.reset_mask()

    def test_get_aligned_data_to_self(self, tpp):
        ring = tpp.data["ringmap"]
        aligned = ring.get_aligned_data(ring.null_alignment)
        assert aligned.length == ring.length

    def test_filter_by_sequence_does_not_crash(self, tpp):
        ring = tpp.data["ringmap"]
        ring.reset_mask()
        ring.mask_on_sequence(compliment_only=True)


class TestPairingProbability:
    def test_has_sequence(self, tpp):
        pp = tpp.data["pairprob"]
        assert isinstance(pp.sequence, str)
        assert len(pp.sequence) > 0

    def test_is_pairing_probability(self, tpp):
        pp = tpp.data["pairprob"]
        assert isinstance(pp, data.PairingProbability)
