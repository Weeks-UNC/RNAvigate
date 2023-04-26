Arc compare plots
=================

Arc compare plots are a flexible plot type for comparing multiple layers of RNA
data between two samples. Sequence annotations, per-nucleotide data,
inter-nucleotide data, and secondary structures can all be displayed.

There is one method to quickly make arc compare plots:

`plot1 = rnavigate.plot_arcs_compare(samples=[sample, another_sample])`

`sample` and `another_sample` here are hypothetical rnavigate.Sample objects
containing data. As written, this code creates a single panel arc compare plot
displaying the `"ct"` secondary structure data from the two samples.

Below are all of the optional arguments that work with this method, along with
their default values. `plot2` below would produce exactly the same result as
`plot1`.

```python
plot2 = rnavigate.plot_arcs_compare(
    samples=[sample, another_sample],
    seq_source=None,
    ct="ct",
    comp=None,
    interactions=None,
    interactions_filter={},
    interactions2=None,
    interactions2_filter={},
    profile="profile",
    plot_error=True,
    annotations=[],
    annotation_mode="track",
    labels=None,
    region="all",
    colorbar=True,
    plot_kwargs={"figsize": None},
)
```

Many of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

---

`samples`

* This argument is required, and is a list of 2 `rnavigate.Sample` objects.
* The first sample is plotted on top, and the second is plotted below.

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

`profile`

* A sample.data key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  `["shapemap", "dmsmap", "dancemap", "rnpmap"]`
* These data are mapped to `seq_source` and plotted as a bar graph along the
  central axis.

---

`plot_error`

* `True` or `False`, whether to plot error bars on `profile` data.

---

`annotations`

* A list of sample.data keys that point to sequence annotations.
* These annotations are mapped to `seq_source`, then plotted as vertical bars.

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

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.AP()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
