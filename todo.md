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

Coding To-Do List
-----------------

### Highest
- log file specification is not working for ShapeMapper_v2.2.0
- file name changes:
  - ct.py -> secondary_structure.py
  - pdb.py -> atomic_coordinates.py
  - functions module -> plotting.py
- move plotting functions to plots, unconnected from plot objects
- refactoring analyses
  - [X] fragmapper
  - [x] DeltaSHAPE
  - [ ] AUROC
  - [x] lowSS
  - [ ] log-diff
  - [ ] RNPMapper
- loading secondary structure files containing more than one structure
### Medium
- annotations should be more flexible, be able to take a simple table file
- codon usage bias
- RNAvigate conda package, pip package, or docker image
- ss plot scales with figure size and can be scaled down.
- add a printable human-readable identifier for all data:
  - sample + name + datatype + filepath?
### Low
- make color scales look good, more consistent.
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
- implement profiles for interaction data
  - Pairing probability -> per-nucleotide probability or shannon entropy
  - RING-MaP -> RING density
  - store as attribute, get attribute at point of use
    - e.g. profile.profile returns profile
