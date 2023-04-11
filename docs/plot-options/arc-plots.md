Arc plots
=========

Arc plots are a flexible plot type for comparing multiple layers of RNA data
simultaneously. Sequence annotations, per-nucleotide data, inter-nucleotide
data, and secondary structures can all be displayed on arc plots.

There are two ways to quickly make arc plots:

```python
plot1 = sample.plot_arcs()
plot2 = rnavigate.plot_arcs(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
single panel arc plot using the data from `sample`. With the second method, you
can also create a multi-panel plot by passing multiple samples to the samples argument. e.g.:

```python
plot3 = rnavigate.plot_arcs(samples=[sample, another_sample])
```

Using either method, you can explore multiple filtering schemes per sample
using the `filters` argument, explained below. Below are all of the optional
arguments that work with each of the methods above, along with their default
values. `plot4` below would produce exactly the same result as `plot1` and
`plot2`.

```python
plot4 = sample.plot_arcs(
    seq_source=None,
    ct="ct",
    comp=None,
    ct_panel="top",
    interactions=None,
    interactions_filter={},
    interactions_panel="bottom",
    filters=None,
    interactions2=None,
    interactions2_filter={},
    interactions2_panel="bottom",
    profile="profile",
    plot_error=True,
    annotations=[],
    annotation_mode="track",
    labels=None,
    title=True,
    region="all",
    colorbar=True,
    seqbar=True,
    plot_kwargs={"rows": None, "cols": None, "figsize": None},
)
```

Many of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

---

`seq_source`

* A sequence string, sample.data key, or Data object.
* Secondary structures, per-nucleotide and inter-nucleotide data, and sequence
  annotations are mapped to this sequence using either a pairwise sequence
  alignment, or an alignment previously defined by the user, prior to plotting.
* If `seq_source` is not provided, it is set to the value of the `ct` argument.

---

`ct` and `comp`

* A sample.data key that points to a secondary structure, e.g.: `"ct"`,
  `"compct"`, `"ss"`, etc.
* If only `ct` is provided, basepairs are drawn as arcs.
* If `comp` is also provided, basepairs from both structures are drawn as arcs,
  which are colored by whether they appear in `ct`, `comp` or both.

---

`ct_panel`

* This can be `"top"` or `"bottom"`. It determines whether to plot `ct` arcs
  above or below the central x-axis.

---

`interactions` and `interactions2`

* These values are sample.data keys that point to inter-nucleotide data, e.g.:
  `"ringmap"`, `"pairmap"`, `"pairprob"`, `"shapejump"`, etc.
* These data are mapped to seq_source, filtered using the arguments below, then
  plotted as arcs.

---

`interactions_filter` and `interactions2_filter`

* A dictionary of key-value pairs that specifies how `interactions` and
  `interactions2` are filtered and displayed
* See [interactions guide](../guides/filters.md) for more detail.

---

`interactions_panel` and `interactions2_panel`

* `"Top"` or `"Bottom"`, where to plot `interactions` and `interactions2` data.

---

`filters`

* A list of dictionaries specifying interactions data and filtering schemes.
* Each filtering scheme will be applied to each sample, and plotted on separate
  plots. e.g. 3 samples and 2 filtering schemes produces 6 plots.
* This is an alternative to `interactions` and `interactions_filter`, those
  arguments will be ignored.
* See [interactions guide](../guides/filters.md) for more detail.

---

`profile`

* A sample.data key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]
* These data are mapped to `seq_source` and plotted as a bar graph along the
  central axis.

---

`plot_error`

* `True` or `False`, whether to plot error bars on `profile` data.

---

`annotations`

* A list of sample.data keys that point to sequence annotations.
* These annotations are mapped to `seq_source`, then plotted along the central
  x-axis.

---

`annotation_mode`

* Either `"track"` or `"vbar"`.
* `"track"` uses markers along the x-axis.
* `"vbar"` uses transparent vertical bars spanning the y-axis.

---

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

---

`title`

* `True` or `False`. Display titles using `labels`.

---

`region`

* A list containing a start and end position (1-indexed, inclusive)
* e.g. `region=[40, 100]` plots nucleotides 40 through 100.
* Defaults to plotting the entire sequence from `seq_source`.

---

`colorbar`

* `True` or `False`. Display the color scale for interactions data.

---

`seqbar`

* `True` or `False`. Display the sequence along the x-axis.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.AP()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"rows"` and `"cols"` specifies the number of axes rows and columns of the
  matplotlib figure.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
