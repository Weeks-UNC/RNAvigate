Visualizing data
================

This is part 3 in the getting started with RNAvigate guide.

1. [Installing RNAvigate](./installing-rnavigate.md)
2. [Loading Data](./loading-data.md)
3. [Visualizing data](./visualizing-data.md)

With your data loaded using an RNAvigate `Sample` object, you can start
creating visualizations and analyses. This document will show you what types of
plots can be created and where to find more information.

Most visualizations can be created with a simple plotting function.
```python
plot = rnavigate.plot_plottype(samples=[sample, another_sample])
```

`sample` and `another_sample` are sample objects created with
`rnavigate.Sample()`. `plot_plottype` should be replaced with one of the valid
plotting functions, described below.

---

[Alignment plots](plot-options/alignment-plots.md)
--------------------------------------------------

Alignment plots visualize how two sets data will be positionally aligned in
RNAvigate plots. It is a good idea to check the automatic alignment if two
sequences differ significantly. For the most part, this is not necessary for
deletion or point mutants or subsequences.

```python
rnavigate.plot_alignment(
    data1=(sample, "data_keyword"),
    data2=(another_sample, "data_keyword"))
```

[Skyline plots](plot-options/skyline-plots.md)
----------------------------------------------

Skyline plots flexibly display and compare per-nucleotide data sets.

```python
rnavigate.plot_skyline(samples=[sample])
```

[Profile plots](plot-options/profile-plots.md)
----------------------------------------------

Profile plots display per-nucleotide data as colored bar graphs similar to ShapeMapper style bar graphs, but very flexible.

```python
rnav.plot_profile(samples=[sample])
```

[Arc plots](plot-options/arc-plots.md)
--------------------------------------

Arc plots flexibly display inter-nucleotide relationships and secondary
structures as arcs. Per-nucleotide measurements and sequence annotations can
also be displayed without over-crowding.

```python
rnavigate.plot_arcs(samples=[sample])
```

[Arc compare plots](plot-options/arc-compare-plots.md)
------------------------------------------------------

Arc compare plots are the same as arc plots above, but compare two samples
on the same axes.

```python
rnavigate.plot_arcs_compare(samples=[sample, another_sample])
```

[Circle plots](plot-options/circle-plots.md)
--------------------------------------------

Circle plots are similar to arc plots, but display nucleotides in a circle so
that any size RNA fits in a square area.

```python
rnavigate.plot_circle(samples=[sample])
```

[Secondary structure diagrams](plot-options/ss-plots.md)
--------------------------------------------------------
Secondary structure diagrams display per-nucleotide and inter-nucleotide
measurements and sequence annotations on a secondary structure diagram.

```python
rnavigate.plot_ss(samples=[sample])
```

[Interactive 3D molecule renderings](plot-options/mol-plots.md)
---------------------------------------------------------------
These renderings display per-nucleotide and inter-nucleotide data on
3D RNA molecular structures.

```python
rnavigate.plot_mol(samples=[sample])
```

[Distance distribution histograms](plot-options/disthist-plots.md)
------------------------------------------------------------------

Distance distribution histograms display the 3D distance distribution of sets
of inter-nucleotide measurements.

```python
rnavigate.plot_disthist(samples=[sample])
```

[Heatmaps and contour maps](plot-options/heatmap-plots.md)
----------------------------------------------------------

Heatmaps and contour maps are useful for displaying dense inter-nucleotide
measurements or inter-nucleotide data density while highlighting defined
regions such as helices or 3D interactions.

```python
rnavigate.plot_heatmap(samples=[sample])
```

[Linear regression plots](plot-options/linreg-plots.md)
-------------------------------------------------------

Linear regression plots quickly compare multiple sets of per-nucleotide data
and display slope, R^2, and density.

```python
rnavigate.plot_linreg(samples=[sample])
```

[Receiver operator characteristic curves](plot-options/roc-plots.md)
--------------------------------------------------------------------

ROC curves are a way to determine how well per-nucleotide data predict a
classifier, such as base-paired vs. unpaired status.

```python
rnavigate.plot_roc(samples=[sample])
```

[ShapeMapper2 quality control plots](plot-options/qc-plots.md)
--------------------------------------------------------------

`rnavigate.plot_qc(samples=[sample])`

ShapeMapper2 QC plots display useful quality control metrics from ShapeMapper2
analyses.

[ShapeMapper2 profiles](plot-options/sm-plots.md)
-------------------------------------------------

`sample.plot_shapemapper()`, there is no method for multiple samples.

ShapeMapper2 profiles display normalized SHAPE profiles, mutation rates, and/or
read depths in ShapeMapper2's default layout.


Analyses
--------

These analyses perform more in-depth data manipulation in addition to plotting.

* [Windowed AUROC](analysis-options/auroc.md)
* [DeltaSHAPE](analysis-options/deltashape.md)
* [Log profile comparisons](analysis-options/logcompare.md)
* [Low SHAPE, Low Shannon entropy](analysis-options/lowss.md)
