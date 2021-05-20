Master Plan for Plot-MaP
========================
Plot-MaP should be an easy to use, universal tool for plotting Mutational
Profiling (MaP) and Juxtaposed Merged Pairs (JuMP) data and quality control
information. It should be easily extensible to accomodate new plots, data, and
analysis.

To-Do List
----------
- [ ] Keep documentation updated for testers.
- [ ] consider moving to kwarg dicts for more customizability
- [ ] implement single sample arrays with filters as kwarg dicts.
- [ ] Implement support for probabilities in ap, ss, maybe 3d.
- [ ] In ss and 3D, add filter for only ij_data that are base-paired.
- [ ] Support for different annotation types
  - [ ] vertical bars in ap
  - [ ] highlight nts for ss
  - [ ] not sure how to show this in 3D figures.
- [ ] 3D-distance histograms of filtered ij data.
- [ ] add colorbar and legends for ap, ss, pdb plots
- [ ] pairwise regression scatterplots with kde along diagonal.
- [ ] support for filtering by multiple cts simultaneously.
- [ ] easy and robust installation guidelines
- [ ] guide for more custom use cases.
- [ ] Class and method descriptions.
- [ ] Reformatting heatmaps to be more familiar.
- [ ] implement 2 sample analyses, rnp, deltaSHAPE, klog normalization

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
| load data | Y   | Y        | Y     | Y     | Y     | Y         |       | ----- |       |
| QC plots  | Y   | Y        | ---   | ---   | ---   | ---       | ---   | Y     |       |
| Regression| --- | Y        |       | ---   | ---   | ---       | ---   |       |       |
| skyline   | --- | Y        | Y     | ---   | ---   | ---       |       | Y     |       |
| arc plots | --- | Y        | ---   | Y     | Y     | Y         |       | Y     |       |
| secondary | --- | Y        | ---   | Y     | Y     | Y         |       | Y     |       |
| tertiary  | --- | Y        | ---   | Y     | Y     | Y         |       |       |       |
| heatmaps  | --- |          |       | Y     | Y     | Y         |       |       |       |

- Reading in structural information

| ct | compct | xrna | varna | nsd | cte | pdb |
|----|--------|------|-------|-----|-----|-----|
| Y  | Y      | Y    | Y     | Y   | Y   | Y   |

- Analyses
  - [ ] RNP-MaP
  - [ ] log(+/-) - k*log(+/-) Normalization
  - [ ] deltaSHAPE

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
