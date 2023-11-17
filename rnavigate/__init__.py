from rnavigate.rnavigate import Sample
from rnavigate.plotting_functions import (
    plot_alignment,
    plot_arcs,
    plot_arcs_compare,
    plot_circle,
    plot_disthist,
    plot_heatmap,
    plot_linreg,
    plot_mol,
    plot_ntdist,
    plot_profile,
    plot_qc,
    plot_roc,
    plot_shapemapper,
    plot_skyline,
    plot_ss,
    )
from rnavigate import analysis
from rnavigate import data
from rnavigate import plots
from rnavigate import styles
from rnavigate import transcriptomics

__version__ = "1.0.0"
__author__ = "Patrick S. Irving"
__email__ = "psirving@email.unc.edu"

__all__ = [
    "Sample",
    # plotting functions
    "plot_alignment",
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
    # modules
    "analysis",
    "data",
    "plots",
    "styles",
    "transcriptomics",
]
