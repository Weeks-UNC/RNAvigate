Skyline plots
=============

Skyline plots are useful for comparing multiple sets of per-nucleotide data.
Values are displayed using a stepped line graph (skyline) with
sequence and sequence annotations along the x-axis.

There are two ways to quickly make skyline plots:

```python
plot1 = sample.plot_skyline()
plot2 = rnavigate.plot_skyline(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
single panel skyline plot using the data from `sample`. With the second method,
you can also compare data between multiple samples by passing multiple samples
to the `samples` argument. e.g.:

```python
plot3 = rnavigate.plot_skyline(samples=[sample, another_sample])
```

Below are all of the optional arguments that work with each of the methods
above, along with their default values. `plot4` below produces exactly the
same result as `plot1` and `plot2`.

```python
plot4 = sample.plot_skyline(
    seq_source=None,
    profile="profile",
    columns="Reactivity_profile",
    errorbars=None,
    annotations=[],
    labels=None,
    region="all",
    seqbar=True,
    plot_kwargs={"rows": None, "cols": None, "figsize": None},
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
* These data are mapped to `seq_source` and used to color nucleotides.

---

`columns`

* A list of column names of `sample.data[profile].data`.
* The default is `"Reactivity_profile"` which is valid for `"shapemap"`,
  `"dmsmap"` and `"dancemap"`.
* Run `print(sample.data[profile].data.columns)`, replacing `profile` with the
  intended value of `profile` to see valid column names.

---

`errorbars`

* A column name of `sample.data[profile].data`.
* If given, errors will be plotted based on these values.

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

`region`

* A list containing a start and end positions, 1-indexed, inclusive.
* e.g., `region=[40, 100]` will plot nucleotide positions 40 to 100.

---

`seqbar`

* `True` or `False`, whether to display the sequence along the x-axis.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.Circle()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"rows"` and `"cols"` specifies the number of axes rows and columns of the
  matplotlib figure.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
