Circle plots
============

Circle plots arrange nucleotides 5' to 3' around a circle and are a useful
visualization for multiple layers of RNA data, especially long RNAs with long
primary distance interactions. Sequence annotations, per-nucleotide data,
inter-nucleotide data, and secondary structures can be displayed.

There are two ways to quickly make circle plots:

```python
plot1 = sample.plot_circle()
plot2 = rnavigate.plot_circle(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
single panel circle plot using the data from `sample`. With the second method,
you can also create a multi-panel plot by passing multiple samples to the
samples argument. e.g.:

```python
plot3 = rnavigate.plot_arcs(samples=[sample, another_sample])
```

Using either method, you can explore multiple filtering schemes per sample
using the `filters` argument, explained below. Below are all of the optional
arguments that work with each of the methods above, along with their default
values. `plot4` below would produce exactly the same result as `plot1` and
`plot2`.

```python
plot4 = sample.plot_circle(
    seq_source=None,
    ct="ct",
    comp=None,
    interactions=None,
    interactions_filter={},
    filters=None,
    interactions2=None,
    interactions2_filter={},
    profile="profile",
    annotations=[],
    labels=None,
    plot_kwargs={"rows": None, "cols": None, "figsize": None},
    colors="sequence",
    apply_color_to="sequence",
    positions=True,
    title=True,
    colorbar=True,
)
```

Many of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

---

`seq_source`

* A sequence string, `sample.data` key, or Data object.
* Secondary structures, per-nucleotide and inter-nucleotide data, and sequence
  annotations are mapped to this sequence using either a pairwise sequence
  alignment or an alignment previously defined by the user, prior to plotting.
* If `seq_source` is not provided, it is set to the value of the `ct` argument.

---

`ct` and `comp`

* A sample.data key that points to a secondary structure, e.g.: `"ct"`,
  `"compct"`, `"ss"`, etc.
* If only `ct` is provided, basepairs are drawn as arcs.
* If `comp` is also provided, basepairs from both structures are drawn as arcs,
  which are colored by whether they appear in `ct`, `comp` or both.

---

`interactions` and `interactions2`

* These values are sample.data keys that point to inter-nucleotide data, e.g.:
  `"ringmap"`, `"pairmap"`, `"pairprob"`, `"shapejump"`, etc.
* These data are mapped to `seq_source`, filtered , then plotted as arcs.

---

`interactions_filter` and `interactions2_filter`

* A dictionary of key-value pairs that specifies how `interactions` and
  `interactions2` are filtered and displayed
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

`profile`

* A `sample.data` key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]
* These data are mapped to `seq_source` and used to color nucleotides.

---

`annotations`

* A list of `sample.data` keys that point to sequence annotations.
* These annotations are mapped to `seq_source`, then used to highlight
  nucleotides.

---

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

---

`title`

* `True` or `False`. Display titles using `labels`.

---

`colorbar`

* `True` or `False`. Display the color scale for interactions data.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.Circle()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"rows"` and `"cols"` specifies the number of axes rows and columns of the
  matplotlib figure.
* `"figsize"` specifies the total size of the matplotlib figure in inches.

---

`colors`

* Can be any valid matplotlib color, a list of colors with length equal
  to the sequence of `ss`, or one of the following:
* `"position"` colors by position using a spectrum.
* `"sequence"` colors by nucleotide.
    * A: red
    * U: light red
    * G: blue
    * C: light blue
* `"profile"` colors using `profile` data.

---

`apply_color_to`

* `"background"` applies `colors` to nucleotide markers.
* `"sequence"` applies `colors` to nucleotide letters.
