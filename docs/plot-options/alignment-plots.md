alignment plots
===============

Alignment plots are a quick way to visualize how RNAvigate is positionally
aligning two data sets. The two sequences will be displayed along the x-axis,
with dashes inserted where there are insertions/deletions. A simplified visualization of the alignment appears between the sequences, with grey bars
indicating the sequence agreement and red bars indicating mismatches.

There is one method to quickly display alignment plots:

```python
plot1 = rnavigate.plot_arcs(
    data1=(sample, "data_keyword"),
    data2=(another_sample, "data_keyword"))
```

`sample` and `another_sample` here are hypothetical rnavigate.Sample objects
containing data, and `"data_keyword"` points to data within each sample.

Below are all of the optional arguments that work with this method, along with
their default values. `plot2` below would produce exactly the same result as
`plot1`.

```python
plot2 = rnavigate.plot_alignment(
    data1=(sample, "data_keyword"),
    data2=(another_sample, "data_keyword"),
    labels=None
)
```

---

`data1` and `data2`

* These arguments are required, and each is a list or tuple of an
  `rnavigate.Sample` object and a data keyword.
* The first sample is plotted on top, and the second is plotted below.

---

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.
