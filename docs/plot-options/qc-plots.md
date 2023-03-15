ShapeMapper2 quality control plots
==================================

ShapeMapper2 quality control plots display useful quality control metrics:
mutations per read distributions, read length histograms, and modified vs.
untreated per-nucleotide reactivity distributions as a violin plot. These
metrics are contained in the shapemapper_log.txt files output by ShapeMapper2,
but only if the `--per-read-histogram` flag is used.

First, Mutations per read is an important metric to determine if reads in the
modified sample contain more mutations than reads from the untreated sample.
The shift in distribution will vary depending on many factors, such as read
length, intrinsic reactivity of the RNA and the probe, cell permeability, 
enzyme fidelity, etc. For single-molecule correlated chemical probing
experiments, the shift in distribution should be dramatic, with modified
samples typically containing 4-10 mutations per read. A bimodal distribution
can be indicative of DNA contamination or compartmental difference in reagent
permeability, as sometimes happens with clumped cells.

Read length distributions for amplicon data should be punctate. Large numbers
of shorter reads indicate poor read merging, possibly due to short overlaps or
low NGS base-calling quality scores.

Per-nucleotide reactivity distributions should show an upward shift in modified
mutation rate compared to untreated rate. This will depend on reagent
reactivity, RNA structuredness, enzyme fidelity, etc.

There are two ways to quickly make qc plots:

```python
plot1 = sample.plot_qc()
plot2 = rnavigate.plot_qc(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create the
multi-panel quality control plot using the data from `sample`. With the second
method, you can also compare samples by passing multiple samples to the
`samples` argument. e.g.:

```python
plot3 = rnavigate.plot_qc(samples=[sample, another_sample])
```

Below are all of the optional arguments that work with each of the methods
above, along with their default values. `plot4` below would produce exactly the
same result as `plot1` and `plot2`.

```python
plot4 = sample.plot_qc(
    labels=None,
    plot_kwargs={"figsize": None},
)
```

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.AP()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"figsize"` specifies the total size of the matplotlib figure in inches.
