Plan for RNAvigate
==================

Changes checklist
-----------------
- change the code
- change the docstring (as well as higher-order functions)
- add to documentation website

Top priorities (3 maximum)
--------------------------
- completing documentation website and writing doc strings
- find best way to display colormaps for all plots (probably a 2nd fig)
  - store colorbar parameters and only create unique colorbars
- test and check everything before the next major update.

Coding To-Do List
-----------------

- RNAvigate conda package
- clean up namespace: ideally, no redundancy
- support for undercase nucleotides
- alignment mapping in data objects instead of at point of use
- loading in annotations files
- loading secondary structure files containing more than one structure
- API changes:
  - ct, compct, ss to basepairs or ss or something
  - ct, ij, ij2 to structures, interactions
- better figure scaling using subplotpars -> set_size_inches
- reasonable sizing for all plots
- text scaling factors based on plot size
- implement profiles for interaction data
  - Pairing probability -> per-nucleotide probability or shannon entropy
  - RING-MaP -> RING density
  - store as attribute, get attribute at point of use
    - e.g. profile.profile returns profile
- get median/average/mode for windows in profile data
- passing override values to init_dance
- better nucleotide colors and display on linear regression plots
- reading in FORNA JSON: how to distinguish between FORNA and R2DT?
- calling significant sites with log-corrected profile min-diff comparison
- New plots:
  - Stand-alone Profile bar graph
  - Paired/Unpaired KDE
  - disthist as violin plots

Documentation To-Do List
------------------------

- Reference page for each plot type
  - arc plot
  - circle plot
  - skyline
  - qc
  - heatmap
  - linreg
  - sm
  - disthit
- guide for custom use cases
  - loading custom profiles, interactions and annotations
  - plot manipulation with mpl
  - data manipulation with pandas
- good doc strings
  - [X] rnavigate
  - [ ] styles
  - [ ] analysis
    - [x] logcompare
    - [x] auroc
    - [ ] deltashape
    - [x] lowss
  - [ ] data
    - [ ] annotation
    - [ ] ct
    - [ ] data
    - [ ] dp
    - [ ] interactions
    - [ ] log
    - [ ] pdb
    - [ ] profile
  - [ ] plots
    - [ ] arc
    - [ ] circle
    - [ ] disthist
    - [ ] heatmap
    - [ ] linreg
    - [ ] mol
    - [ ] plots
    - [ ] qc
    - [ ] roc
    - [ ] skyline
    - [ ] sm
    - [ ] ss
