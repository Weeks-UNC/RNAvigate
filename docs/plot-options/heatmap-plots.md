Heatmap and contour plots
=========================

Heatmaps and contour plots are a good way to view inter-nucleotide data when
the data are very dense, or when you wish to view inter-nucleotide data
density. The x and y axes are nucleotide positions. XY-positions are colored by
inter-nucleotide data values or density (heatmap) and overlayed with a contour
map. This contour map can show either 3D distance or secondary structure graph
path distance (contact distance). These countours will outline areas where
nucleotide *X* is close-in-space to nucleotide *Y*.

There are two ways to quickly make heatmap and contour plots:

```python
plot1 = sample.plot_heatmap()
plot2 = rnavigate.plot_heatmap(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
single panel heatmap using the data from `sample`. With the second method, you
can also create a multi-panel heatmap by passing multiple samples to the
samples argument. e.g.:

```python
plot3 = rnavigate.plot_heatmap(samples=[sample, another_sample])
```

Using either method, you can explore multiple filtering schemes per sample
using the `filters` argument, explained below. Below are all of the optional
arguments that work with each of the methods above, along with their default
values. `plot4` below would produce exactly the same result as `plot1` and
`plot2`.

```python
plot4 = sample.plot_heatmap(
    structure=None,
    interactions=None,
    interactions_filter={},
    filters=None,
    labels=None,
    levels=None,
    atom="O2'",
    regions=None,
    interpolation=None,
    plot_type="heatmap",
    plot_kwargs={"rows": None, "cols": None, "figsize": None},
)
```

Some of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

---

`structure`

* A sample.data key pointing to either a secondary structure or a 3D structure
  with atomic coordinates.
* Inter-nucleotide data are mapped to this sequence using either a pairwise
  sequence alignment, or an alignment previously defined by the user.

---

`interactions`

* A sample.data key that points to inter-nucleotide data, e.g.:
  `"ringmap"`, `"pairmap"`, `"pairprob"`, `"shapejump"`, etc.
* These data are mapped to `structure`, filtered using the arguments below, then
  plotted as a heatmap.

---

`interactions_filter`

* A dictionary of key-value pairs that specifies how `interactions` are
  filtered and displayed.
* See [interactions guide](../guides/filters.md) for more detail.

---

`filters`

* A list of dictionaries specifying interactions data and filtering schemes.
* Each filtering scheme will be applied to each sample, and plotted on separate
  plots. e.g. 3 samples and 2 filtering schemes produces 6 plots.
* This is an alternative to `interactions` and `interactions_filter`, those
  arguments will be ignored.
* See [interactions guide](../guides/filters.md) for more detail.

---

`labels`

* A list of strings, one for each sample, used as titles of each plot.
* Defaults to using the sample name by retrieving `sample.sample`.

---

`levels`

* A list of distance values around which to draw contours.
* For secondary structures, this defaults to `[5]` for contact distance, and
  for 3D structures, this defaults to `[20]` for angstroms.

---

`atom`

* A string specifying which atom use to calculate atomic distances.
* `"DMS"` specifies `"N1"` for A and C, `"N3"` for U and G.

---

`regions`

* A list of list of tuples. e.g.:
  `regions=[[(1, 10), (45, 55)], [(25, 35), (70, 80)]]`
* This would plot two boxes, instead of a contour plot based on `structure`
  distances.
* In this case, the boxes would be drawn from 1-10 on the x-axis to 45-55 on
  the y-axis and from 25-35 on the x-axis to 70-80 on the y-axis.

---

`plot_type`

* `"heatmap"` or `"kde"`, whether to plot inter-nucleotide data values
  (heatmap) or the kernel density estimate (KDE).

---

`interpolation`

* Image interpolation type passed to `matplotlib.pyplot.imshow()` when
  `plot_type="heatmap"`.
* Usually, `None` is best.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.AP()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"rows"` and `"cols"` specifies the number of axes rows and columns of the
  matplotlib figure.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
