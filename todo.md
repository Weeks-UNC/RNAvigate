Plan for RNAvigate
==================
Plot-MaP should be an easy to use, universal tool for plotting Mutational
Profiling (MaP) and Juxtaposed Merged Pairs (JuMP) data and quality control
information. It should be easily extensible to accomodate new plots, data, and
analysis.

Coding To-Do List
-----------------
- calling significant sites with log-corrected profile min-diff comparison.
- better nucleotide colors and display on linear regression plots.
- loading in annotations files
- Reformatting heatmaps to be more familiar.
- deltaSHAPE
- send to file for all plots

Documentation To-Do List
------------------------
- Installation instructions
- One page for each plot type
  - arcplot, circleplot, ss-diagram, 3d mol, skyline, qc, heatmap, linreg, sm
- guide for custom use cases
  = loading profile data, loading ij data, plots with mpl, data with pandas
- good doc strings

Features
--------
A python interface and CLI interface for all of the following:

###Loading MaP/JuMP data and plotting
- Implemented python interface: Y
- Implemented CLI interface: YY
- Not applicable: ---

| plot type | log | profiles | dance | rings | pairs | deletions | array | probs |
|-----------|-----|----------|-------|-------|-------|-----------|-------|-------|
| load data | Y   | Y        | Y     | Y     | Y     | Y         | ----- | Y     |
| QC plots  | Y   | Y        | ---   | ---   | ---   | ---       | Y     | ---   |
| Regression| --- | Y        | Y     | ---   | ---   | ---       | Y     | ---   |
| skyline   | --- | Y        | Y     | ---   | ---   | ---       | Y     | ---   |
| arc plots | --- | Y        | ---   | Y     | Y     | Y         | Y     | Y     |
| secondary | --- | Y        | ---   | Y     | Y     | Y         | Y     | Y     |
| tertiary  | --- | Y        | ---   | Y     | Y     | Y         | Y     | Y     |
| heatmaps  | --- |          | Y     | Y     | Y     | Y         | Y     | Y     |

- Reading in structural information

| ct | xrna | varna | nsd | cte | pdb |
|----|------|-------|-----|-----|-----|
| Y  | Y    | Y     | Y   | Y   | Y   |

- Analyses
  - [ ] RNP-MaP
  - [ ] log(+/-) - k*log(+/-) Normalization
  - [ ] deltaSHAPE
  - [ ] destabilization analysis

CLI examples
------------
Currently, there is no CLI interface, but this is how I imagine it working.
- aliases
  - sl - skyline
  - ap - arcplot
  - ss - secondarystructure
  - 3d - 3d structure
  - qc - muts/mol, read length, and boxplots
  - rnp - RNP-Map
  - delta - deltaSHAPE
  - lr - linear regression
```
rnavigate.py sl --profiles prof1.txt prof2.txt prof3.txt
rnavigate.py sl --dance sample1-reactivities.txt
rnavigate.py sl --jump jump1.file jump2.file jump3.file jump4.file
rnavigate.py ap --ct targetrna.ct
                 --samples sample1 sample2 sample3 sample4
                 --profile-suffix _targetrna_profile.txt
                 --rings-suffix -targetrna.corrs
rnavigate.py ss --ss targetrna.xrna
                 --samples sample1 sample2 sample3 sample4
                 --profile-suffix _targetrna_profile.txt
                 --rings-suffix -targetrna.corrs
```
