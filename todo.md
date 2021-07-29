Master Plan for Plot-MaP
========================
Plot-MaP should be an easy to use, universal tool for plotting Mutational
Profiling (MaP) and Juxtaposed Merged Pairs (JuMP) data and quality control
information. It should be easily extensible to accomodate new plots, data, and
analysis.

To-Do List
----------
### Broken
- windowed ijs that are (-1, -1) plots one line at (-1, 1)
- colorby in ss not being pass correctly
### Objectified
  - color by ij codepth
  - new_plots.qc
### Plots
- 3d/contact distance histograms.
- Minimized log-log comparison of profiles.
- add colorbars and legends
- Plotting annotations
  - vertical bars in ap
  - highlight nts for ss
  - clouds, colors, or labels in 3D figures.
- Reformatting heatmaps to be more familiar.
- deltaSHAPE
- RNP-MaP
- Bin by MutsPerMol
### Documentation
- Keep documentation updated for testers.
- easy and robust installation guidelines
- guide for more custom use cases.
- Class and method descriptions.

Features
--------
A python interface and CLI interface for all of the following:

###Loading MaP/JuMP data and plotting
- Implemented python interface: Y
- Implemented CLI interface: YY
- Beautifully represents data: Done
- Not applicable: ---

| plot type | log | profiles | dance | rings | pairs | deletions | frags | array | probs |
|-----------|-----|----------|-------|-------|-------|-----------|-------|-------|-------|
| load data | Y   | Y        | Y     | Y     | Y     | Y         |       | ----- | Y     |
| QC plots  | Y   | Y        | ---   | ---   | ---   | ---       | ---   | Y     | ---   |
| Regression| --- | Y        | Y     | ---   | ---   | ---       | ---   | Y     | ---   |
| skyline   | --- | Y        | Y     | ---   | ---   | ---       |       | Y     | ---   |
| arc plots | --- | Y        | ---   | Y     | Y     | Y         |       | Y     | Y     |
| secondary | --- | Y        | ---   | Y     | Y     | Y         |       | Y     | Y     |
| tertiary  | --- | Y        | ---   | Y     | Y     | Y         |       | Y     | Y     |
| heatmaps  | --- |          | Y     | Y     | Y     | Y         |       | Y     | Y     |

- Reading in structural information

| ct | compct | xrna | varna | nsd | cte | pdb |
|----|--------|------|-------|-----|-----|-----|
| Y  | Y      | Y    | Y     | Y   | Y   | Y   |

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
plotmapper.py sl --profiles prof1.txt prof2.txt prof3.txt
plotmapper.py sl --dance sample1-reactivities.txt
plotmapper.py sl --jump jump1.file jump2.file jump3.file jump4.file
plotmapper.py ap --ct targetrna.ct
                 --samples sample1 sample2 sample3 sample4
                 --profile-suffix _targetrna_profile.txt
                 --rings-suffix -targetrna.corrs
plotmapper.py ss --ss targetrna.xrna
                 --samples sample1 sample2 sample3 sample4
                 --profile-suffix _targetrna_profile.txt
                 --rings-suffix -targetrna.corrs
```
