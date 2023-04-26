Plan for RNAvigate
==================

Changes workflow
----------------
- make changes to the code
- update the docstring (as well as higher-order functions)
- update the documentation
- add changes to the change log
- push all these changes together

Top priorities (3 maximum)
--------------------------

- completing documentation website and writing doc strings
  - [ ] plots
  - [ ] getting started
  - [ ] analyses
  - [ ] namespace and getting docstrings
- test and check everything before the next major update.

Documentation To-Do List
------------------------

- guide for custom use cases
  - loading custom profiles, interactions and annotations
  - sequence alignments
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
    - [x] annotation
    - [x] ct
    - [x] data
    - [X] interactions
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

Coding To-Do List
-----------------

- check that all functions that compare objects align or check similarity
- unify the scalar -> rgba functionality for all Data objects
- analyses are a sub-class of Sample?
- RNAvigate conda package
- clean up namespace: ideally, no redundancy
- support for undercase nucleotides
- loading secondary structure files containing more than one structure
- combine VARNA, NSD, XRNA, etc classes, rename to SecondaryStructure
- refactor CT to be a subclass of Interactions
- API changes:
  - ct, compct, ss to basepairs or ss or something
  - ct, ij, ij2 to structures, interactions
- reasonable sizing for all plots
- text scaling factors based on plot size (this would be easy with sns.context)
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
  - Paired/Unpaired KDE
  - disthist as violin plots
- refactor code using ProPlot
- refactor circle plot to use polar projection
- sequence + structure alignments
