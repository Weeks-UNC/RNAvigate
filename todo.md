Master Plan for Plot-MaP
========================
Plot-MaP should be an easy to use, universal tool for plotting Mutational
Profiling (MaP) and Juxtaposed Merged Pairs (JuMP) data and quality control
information. It should be easily extensible to accomodate new plot types and
data types.

To-Do List
----------
- [ ] add command line interface
- [ ] read deletion files and add to arc plots and secondary structure
- [ ] implement for pdb:
  - [ ] read data
  - [ ] assign clip_pad
  - [ ] compute 3d distances
  - [ ] add contacts to heatmap, arc plot, or secondary structure
  - [ ] make py3dmol view
  - [ ] add data to view using addLine or addCylinder
    - [ ] rings, deletions, maybe pairs
- [ ] Make qc plot for multiple samples less ugly.

Features
--------
A python interface and CLI interface for all of the following:

###Loading MaP/JuMP data and plotting
- Implemented python interface: Y
- Implemented CLI interface: YY
- Beautifully represents data: Done
- Not applicable: ---

| plot type | log | profiles | dance | rings | pairs | deletions | frag-jump | array |
|-----------|-----|----------|-------|-------|-------|-----------|-----------|-------|
| load data | Y   | Y        | Y     | Y     | Y     |           |           | ----- |
| QC plots  | Y   | Y        | ---   | ---   | ---   | ---       | ---       | Y     |
| Regression| --- | Y        |       | ---   | ---   | ---       | ---       |       |
| skyline   | --- | Y        | Y     | ---   | ---   | ---       |           | Y     |
| arc plots | --- | Y        | ---   | Y     | Y     |           |           | Y     |
| secondary | --- | Y        | ---   | Y     | ???   |           |           | Y     |
| tertiary  | --- |          | ---   |       |       |           |           |       |
| heatmaps  | --- |          |       | meh   | meh   |           |           |       |

- Reading in structural information

| ct | compct | xrna | varna | nsd | cte | pdb |
|----|--------|------|-------|-----|-----|-----|
| Y  | Y      | Y    | Y     | Y   | Y   |     |

- Analyses
  - [ ] RNP-MaP
  - [ ] log(+/-) - k*log(+/-)
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
