"""Smoke tests for rnavigate.analysis classes.

Each test constructs an analysis object from the bundled example datasets and
verifies it does not crash.  Because analysis classes subclass Sample, the
result is also checked to be a Sample instance, confirming it is directly
usable with plot_*() functions.
"""

from __future__ import annotations

import rnavigate as rnav
from rnavigate import analysis


class TestDeltaSHAPE:
    def test_creates_sample(self, rnasep_1, rnasep_2):
        result = analysis.DeltaSHAPE(
            sample1=rnasep_1,
            sample2=rnasep_2,
        )
        assert isinstance(result, rnav.Sample)

    def test_has_deltashape_data(self, rnasep_1, rnasep_2):
        result = analysis.DeltaSHAPE(
            sample1=rnasep_1,
            sample2=rnasep_2,
        )
        assert "deltashape" in result.data


class TestSequenceChecker:
    def test_creates_checker(self, rnasep_1, rnasep_2):
        checker = analysis.SequenceChecker([rnasep_1, rnasep_2])
        assert checker is not None

    def test_has_sequences(self, rnasep_1, rnasep_2):
        checker = analysis.SequenceChecker([rnasep_1, rnasep_2])
        assert hasattr(checker, "sequences")
        assert len(checker.sequences) > 0
