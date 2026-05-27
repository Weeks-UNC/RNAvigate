"""Smoke tests for the top-level rnavigate package."""

from __future__ import annotations


def test_package_imports():
    import rnavigate as rnav

    assert hasattr(rnav, "Sample")
    for fn in [
        "plot_arcs",
        "plot_arcs_compare",
        "plot_circle",
        "plot_disthist",
        "plot_heatmap",
        "plot_linreg",
        "plot_mol",
        "plot_ntdist",
        "plot_profile",
        "plot_qc",
        "plot_roc",
        "plot_shapemapper",
        "plot_skyline",
        "plot_ss",
        "plot_alignment",
        "plot_options",
    ]:
        assert hasattr(rnav, fn), f"missing public function: {fn}"


def test_version():
    import rnavigate as rnav

    assert hasattr(rnav, "__version__")
    assert isinstance(rnav.__version__, str)
    assert len(rnav.__version__) > 0


def test_submodules_importable():
    import rnavigate as rnav

    for mod in ["data", "plots", "analysis", "transcriptomics", "styles"]:
        assert hasattr(rnav, mod), f"missing submodule: {mod}"


def test_examples_module():
    from rnavigate import examples

    assert hasattr(examples, "__getattr__")


def test_data_classes_importable():
    from rnavigate import data

    for cls in [
        "Sequence",
        "Data",
        "Profile",
        "SHAPEMaP",
        "Interactions",
        "SecondaryStructure",
        "Annotation",
        "SequenceAlignment",
    ]:
        assert hasattr(data, cls), f"missing data class: {cls}"
