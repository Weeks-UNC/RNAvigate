Visualizing data
================

This is part 3 in the getting started with RNAvigate guide.
1. [Installing RNAvigate](installing-rnavigate.md)
2. [Loading Data](loading-data.md)
3. [Visualizing data](vizualizing-data.md)

Now that you can load your data files into an RNAvigate object, you can start
creating visualizations and analyses.

Visualizations
--------------

Most visualizations can be created with either of two methods, one for a single
sample, and one for any number of samples:

```python
one_sample_plot = sample.plot_plottype()
two_sample_plot = rnavigate.plot_plottype(samples=[sample, another_sample])
```

`sample` and `another_sample` are sample objects created with
`rnavigate.Sample()`. `plottype` should be replaced with one of the valid plot
types provided by RNAvigate. Below are each of the valid visualization
methods, with a brief description and a link to a more detailed guide on how to
use them:

---

`sample.plot_qc()` and `rnavigate.plot_qc(samples=[sample])`

[ShapeMapper2 quality control plots](plot-options/qc-plots.md) display useful
quality control metrics from ShapeMapper2 analyses.

---

`sample.plot_shapemapper()`, there is no method for multiple samples.

[ShapeMapper2 profiles](plot-options/sm-plots.md) display normalized SHAPE
profiles, mutation rates, and/or read depths in ShapeMapper2's default layout.

---

`sample.plot_skyline()` and `rnavigate.plot_skyline(samples=[sample])`

[Skyline plots](plot-options/skyline-plots.md) flexibly display and compare
per-nucleotide data sets.

---

`sample.plot_arcs()` and `rnavigate.plot_arcs(samples=[sample])`

[Arc plots](plot-options/arc-plots.md) flexibly display per-nucleotide and
inter-nucleotide data, secondary structures, and sequence annotations simultaneously.

---

`sample.plot_circle()` and `rnavigate.plot_circle(samples=[sample])`

[Circle plots](plot-options/circle-plots.md) are similar to arc plots, but
display nucleotides in a circle so that any size RNA fits in a square area.

---

`sample.plot_ss()` and `rnavigate.plot_ss(samples=[sample])`

[Secondary structure diagrams](plot-options/ss-plots.md) displays
per-nucleotide and inter-nucleotide data and sequence annotations on a
secondary structure diagram.

---

`sample.plot_mol()` and `rnavigate.plot_mol(samples=[sample])`

[Interactive 3D molecule renderings](plot-options/mol-plots.md) display
per-nucleotide and inter-nucleotide data on 3D renderings of RNA molecules.

---

`sample.plot_disthist()` and `rnavigate.plot_disthist(samples=[sample])`

[Distance distribution histograms](plot-options/disthist-plots.md) display the
3D distance distribution of sets of inter-nucleotide interactions.

---

`sample.plot_heatmap()` and `rnavigate.plot_heatmap(samples=[sample])`

[Heatmap and contour maps](plot-options/heatmap-plots.md) are useful for
displaying dense inter-nucleotide data sets and inter-nucleotide data density
with close-in-space regions highlighted.

---

`sample.plot_linreg()` and `rnavigate.plot_linreg(samples=[sample])`

[Linear regression plots](plot-options/lingreg-plots.md) compare multiple sets
of per-nucleotide data and calculate slope and R^2.

---

`sample.plot_roc()` and `rnavigate.plot_roc(samples=[sample])`

[Receiver operator characteristic curves](plot-options/roc-plots.md) are a way
to determine how well per-nucleotide data predict base-pairing status.

---

Analyses
--------

* [Windowed AUROC](analysis-options/auroc.md)
* [DeltaSHAPE](analysis-options/deltashape.md)
* [Log profile comparisons](analysis-options/logcompare.md)
* [Low SHAPE, Low Shannon entropy](analysis-options/lowss.md)
