Linear regression plots
=======================

Linear regression plots are useful to determine reproducibility between
replicates, or to quickly quantify the difference between structural states.
RNAvigate creates scatter plots of per-nucleotide values from one sample on the
x-axis for another sample on the y-axis. Slope and R^2 values are displayed.
Nucleotides can be colored by sequence or base-pairing status. A KDE of
paired/unpaired reactivity distributions may also be plotted for each sample.

There are two ways to quickly make linear regression plots:

```python
plot1 = sample.plot_linreg()
plot2 = rnavigate.plot_linreg(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a linear
regression plot and KDE plots using the data from `sample`. With the second
method, you can also cross-compared many samples simultaneously by passing
multiple samples to the `samples` argument. e.g.:

```python
plot3 = rnavigate.plot_linreg(samples=[sample, another_sample])
```

Below are all of the optional arguments that work with each of the methods
above, along with their default values. `plot4` below would produce exactly the
same result as `plot1` and `plot2`.

```python
plot4 = sample.plot_linreg(
    ct="ct",
    profile="profile",
    labels=None,
    colorby="structure",
    plot_kwargs={"figsize": None},
)
```

Some of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

`ct`

* A sample.data key that points to a secondary structure, e.g.: `"ct"`,
  `"compct"`, `"ss"`, etc.
* `ct` is used to classify base-pairing status for KDE plots and optionally to
  color points on the regression scatter plot.

`profile`

* A sample.data key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]
* These data are mapped to each other using a sequence alignment and
  per-nucleotide values are plotted as scatter plot.

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

`colorby`

* `"structure"` or `"sequence"`, how scatter plot points will be colored.

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.AP()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
