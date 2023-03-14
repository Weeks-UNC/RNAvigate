Distance distribution histograms
================================

Distance distribution histograms provide a method of determining how well
inter-nucleotide data and filtering schemes isolate close-in-space contacts.
A histogram of the atomic 3D distances of filtered internucleotide data is
compared to a "background" 3D distance histogram. The background histogram by
default is the distribution of all pairwise 3D distances in the given 3D
structure. It can also be defined by first creating an inter-nucleotide data
set consisting of all nucleotide pairs using the `allpossible` argument with
`rnav.Sample()`, and filtering these.

There are two ways to quickly make distance distribution histograms:

```python
plot1 = sample.plot_disthist()
plot2 = rnavigate.plot_disthist(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
single panel distance distribution histogram using the data from `sample`. With
the second method, you can also create a multi-panel plot by passing multiple samples to the samples argument, or compare multiple samples on a single axis.
e.g.:

```python
plot3 = rnavigate.plot_disthist(samples=[sample, another_sample])
```

Using either method, you can explore multiple filtering schemes per sample
using the `filters` argument, explained below. Below are all of the optional
arguments that work with each of the methods above, along with their default
values. `plot4` below would produce exactly the same result as `plot1` and
`plot2`.

```python
plot4 = sample.plot_disthist(
    structure="pdb",
    interactions=None,
    interactions_filter={},
    filters=None,
    bg_interactions=None,
    bg_interactions_filter={},
    labels=None,
    same_axis=False,
    atom="O2'",
    plot_kwargs={"rows": None, "cols": None, "figsize": None},
)
```

Some of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

`structure`

* A sample.data key that points to a 3D structure with atomic coordinates,
  usually `"pdb"`.
* Inter-nucleotide data are mapped to this sequence using either a pairwise
  sequence alignment, or an alignment previously defined by the user.

`interactions` and `bg_interactions`

* These values are sample.data keys that point to inter-nucleotide data, e.g.:
  `"ringmap"`, `"pairmap"`, `"pairprob"`, `"shapejump"`, `"allpossible"`, etc.
* These data are mapped to seq_source, filtered using the arguments below, used
  to calculate distance distributions.
* `bg_interactions` are used for the background distribution.

`interactions_filter` and `bg_interactions_filter`

* A dictionary of key-value pairs that specifies how `interactions` and
  `bg_interactions` are filtered and displayed.
* See [interactions guide](../filters.md) for more detail.

`filters`

* A list of dictionaries specifying interactions data and filtering schemes.
* Each filtering scheme will be applied to each sample, and plotted on separate
  plots. e.g. 3 samples and 2 filtering schemes produces 6 plots.
* This is an alternative to `interactions` and `interactions_filter`, those
  arguments will be ignored.
* See [interactions guide](../filters.md) for more detail.

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

`same_axis`

* `True` or `False`, whether to plot distributions on the same axis.

`atom`

* A string specifying the atom from which distances are calculated.
* Defaults to `"O2'"`.
* `"DMS"` will use the `"N1"` for A and C, and `"N3"` for
  U and G.

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.DistHist()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"rows"` and `"cols"` specifies the number of axes rows and columns of the
  matplotlib figure.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
