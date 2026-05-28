"""Unit tests for rnavigate.data.Sequence."""

from __future__ import annotations

import numpy as np

from rnavigate import data

SIMPLE_SEQ = "AUGCAUGC"


class TestNormalizeSequence:
    def test_t_to_u(self):
        assert data.normalize_sequence("ATGC", t_or_u="U") == "AUGC"

    def test_u_to_t(self):
        assert data.normalize_sequence("AUGC", t_or_u="T") == "ATGC"

    def test_uppercase(self):
        assert data.normalize_sequence("augc") == "AUGC"

    def test_no_conversion(self):
        assert data.normalize_sequence("ATGC", t_or_u=False) == "ATGC"

    def test_accepts_sequence_object(self):
        seq = data.Sequence(SIMPLE_SEQ)
        result = data.normalize_sequence(seq)
        assert result == SIMPLE_SEQ


class TestSequenceFromString:
    def test_sequence_stored(self):
        seq = data.Sequence(SIMPLE_SEQ)
        assert seq.sequence == SIMPLE_SEQ

    def test_length(self):
        seq = data.Sequence(SIMPLE_SEQ)
        assert seq.length == len(SIMPLE_SEQ)

    def test_default_name_is_none(self):
        seq = data.Sequence(SIMPLE_SEQ)
        assert seq.name is None

    def test_name_stored(self):
        seq = data.Sequence(SIMPLE_SEQ, name="test_rna")
        assert seq.name == "test_rna"


class TestSequenceNullAlignment:
    def test_null_alignment_exists(self):
        seq = data.Sequence(SIMPLE_SEQ)
        assert seq.null_alignment is not None

    def test_null_alignment_maps_to_self(self):
        seq = data.Sequence(SIMPLE_SEQ)
        aln = seq.null_alignment
        positions = np.arange(len(SIMPLE_SEQ)) + 1
        assert np.all(aln.map_positions(positions) == positions)


class TestSequenceColors:
    def test_get_colors_by_sequence(self):
        seq = data.Sequence(SIMPLE_SEQ)
        colors, colormap = seq.get_colors("sequence")
        assert isinstance(colors, np.ndarray)
        assert isinstance(colormap, data.ScalarMappable)
        assert len(colors) == len(SIMPLE_SEQ)

    def test_get_colors_by_position(self):
        seq = data.Sequence(SIMPLE_SEQ)
        colors, colormap = seq.get_colors("position")
        assert isinstance(colors, np.ndarray)
        assert isinstance(colormap, data.ScalarMappable)
        assert len(colors) == len(SIMPLE_SEQ)
