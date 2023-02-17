Plan for RNAvigate
==================

Coding To-Do List
-----------------

- RNAvigate conda package
- support for undercase nucleotides
- ct.filter()
  - apply various filters
  - fit to other sequences via alignment mapping
- alignment mapping in data objects instead of at point of use
- loading in annotations files
- loading secondary structure files containing more than one structure
- API changes:
  - ct, compct, ss to basepairs or ss or something
  - ct, ij, ij2 to structures, interactions
- better figure scaling using subplotpars -> set_size_inches
- reasonable sizing for all plots
- text scaling factors based on plot size
- find best way to display colormaps for all plots (probably a 2nd fig)
- enable multiple structure inputs for SS plotting
  - scale data points equally for ct-like objects
  - center on 0, 0
  - share x and y axes
  - scale figure size using subplotpars
- implement profiles for interaction data
  - Pairing probability -> per-nucleotide probability or shannon entropy
  - RING-MaP -> RING density
  - use get_data function to return profiles when expected
- get median/average/mode for windows in profile data
- passing override values to init_dance
- better nucleotide colors and display on linear regression plots
- improve labelling for disthist plots
- reading in FORNA JSON: how to distinguish between FORNA and R2DT?
- calling significant sites with log-corrected profile min-diff comparison
- New plots:
  - Stand-alone Profile bar graph
  - Paired/Unpaired KDE

Documentation To-Do List
------------------------

- Reference page for each plot type
  - circleplot, skyline, qc, heatmap, linreg, sm, disthit
- guide for custom use cases
  - loading custom profiles, interactions and annotations
  - plot manipulation with mpl
  - data manipulation with pandas
- good doc strings
