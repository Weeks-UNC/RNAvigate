Profile plots
=============

Profile plots are a useful way to compare multiple sets of per-nucleotide data,
similar to [skyline plots](skyline-plots.md). Unlike skyline plots, values are
displayed using a color-coded bar graph, like ShapeMapper plots, with sequences
and sequence annotations along the x-axis. Unlike ShapeMapper plots, this are
more flexible, and color-coding can be defined in the data class.

There is only one method to quickly make profile plots:

```python
plot1 = rnavigate.plot_profile(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, this code will create a single panel profile plot using the data from
`sample`. You can also compare data between multiple samples by passing
multiple samples to the `samples` argument. e.g.:

```python
plot2 = rnavigate.plot_skyline(samples=[sample, another_sample])
```

Below are all of the optional arguments that work with each of the methods
above, along with their default values. `plot3` below produces exactly the
same result as `plot1`.

```python
plot3 = rnav.plot_profile(
    samples=[sample],
    seq_source=None,
    profile="profile",
    column=None,
    error_column=None,
    plot_errors=True,
    annotations=[],
    seqbar=True,
    labels=None,
    region="all",
    plot_kwargs={"figsize": None},
)
```


Some of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

---

`seq_source`

* A sequence string, `sample.data` key, or Data object.
* Per-nucleotide data and sequence annotations are mapped to this sequence
  using either a pairwise sequence alignment or an alignment previously
  defined by the user, prior to plotting.
* If `seq_source` is not provided, it is set to the value of `profile` and
  retrieved from the first sample.

---

`profile`

* A `sample.data` key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  `["shapemap", "dmsmap", "dancemap", "rnpmap"]`
* These data are mapped to `seq_source` and used to plot bars.

---

`column`

* A column name of `sample.data[profile].data`.
* The default is `sample.data[profile].default_column`.

---

`errorbars`

* A column name of `sample.data[profile].data`.
* The default is to use `sample.data[profile].default_err_column`

---

`plot_errors`

* `True` or `False`, whether to plot error bars on bar graphs.

---

`annotations`

* A list of `sample.data` keys that point to sequence annotations.
* These annotations are mapped to `seq_source`, then used to highlight
  nucleotides.

---

`seqbar`

* `True` or `False`, whether to display the sequence along the x-axis.

---

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

---

`region`

* A list containing a start and end positions, 1-indexed, inclusive.
* e.g., `region=[40, 100]` will plot nucleotide positions 40 to 100.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.plots.Profile()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
