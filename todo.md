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

### Highest
- text scaling factors based on plot size (this would be easy with context)
- implement the region alignment for all datatypes and plots
- any genomic data
- get Winstons eCLIP -> profile code
- make color scales look good, more consistent.
- all functions that compare data also align sequences (if appropriate)
- normalizing and windowing profiles
  - by nucleotide
  - windowed mean, median, average
  - gaussian smoothing
- loading secondary structure files containing more than one structure
- Class docstrings contain a description, init explains constructor
### Medium
- refactor code using ProPlot
- build heirarchical debugging system?
- RNAvigate conda package
- ss plot scales with figure size and can be scaled down.
- add a printable human-readable identifier for all data:
  - sample + name + datatype + filepath?
### Low
- get median/average/mode for windows in profile data
- better nucleotide colors and display on linear regression plots
- New plots:
  - Paired/Unpaired KDE
  - disthist as violin plots
- refactor circle plot to use polar projection
- secondary structure alignments
### maybe:
- calling significant sites with log-corrected profile min-diff comparison
- passing override values to init_dance
- analyses are a sub-class of Sample?
- refactor CT to be a subclass of Interactions?
- implement profiles for interaction data
  - Pairing probability -> per-nucleotide probability or shannon entropy
  - RING-MaP -> RING density
  - store as attribute, get attribute at point of use
    - e.g. profile.profile returns profile

Finished
--------
- experimental secondary structure based alignments (RNAlign2D)
- get colorbars working again
- all rnavigate imports should be absolute and single level for readability:
  - from rnavigate import namespace
  - from rnavigate.namespace import variable
- annotation groups should act exactly like sites
- simplify arguments:
  - structure(2) <- ss, ct, comp
  - sequence <- fit_to, dataframe, filepath, seq_source
  - interactions(2) <- ij(2)
  - structure2 <- comp
  - input_data <- dataframe, filepath
  - defaults <- column, err_column, color_column, cmap, norm, min_max, etc.
  - read_table_kw <- read_csv_kw, sep
- unify the scalar -> rgba functionality for all Data objects
- changed sequence alignment behavior
  - all data objects have get_aligned_data which returns a matching datatype
    with new positions
- reading in FORNA JSON: how to distinguish between FORNA and R2DT?
- support for undercase nucleotides
- combine VARNA, NSD, XRNA, JSON, CTE, DBN classes, renamed to SecondaryStructure
- ss plots
  - position labels look much nicer
  - axis margins are standardized to 2 data units for all size RNAs
