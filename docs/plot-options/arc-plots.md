Arc plots
=========

Arc plots are a flexible and useful plot type for displaying and making
comparisons between secondary structures, pairing probabilities, annotations,
and experimentally-measured inter-nucleotide and per-nucleotide data.

Any combination of the sequence, annotations, and per-nucleotide data can be
displayed centrally, and base-pairing, pairing probabilities, and 
inter-nucleotide data can be displayed as semi-circular arcs above or below the
central axis (top or bottom panel).

There are three functions for making arc plots:

```python
plot1 = sample.plot_arcs()
plot2 = sample.plot_arcs_multifilter(filters=filters)
plot3 = rnavigate.plot_arcs_multisample(samples=samples)
```

* `plot1` creates a single plot.
* `plot2` creates multiple plots, displaying interactions data from one
  experimental sample, filtered in different ways.
* `plot3` creates multiple plots, displaying interactions data from multiple
  experimental samples.

Below are all of the optional arguments that work with each of the functions
above, along with their default values.

```python
plot1 = sample.plot_arcs(
    ct="ct",
    comp=None,
    ct_panel="top",
    interactions=None,
    interactions_filter={},
    interactions_panel="bottom"
    interactions2=None,
    interactions2_filter={},
    interactions2_panel="bottom"
    profile="profile",
    annotations=[],
    labels=None,
    region="all",
    prefiltered=False,
    colorbar=True,
    seqbar=True,
    title=True,
    annotation_mode="track")
```

`ct` and `comp`

* These can be any of "ss", "ct", or "compct".
* If `comp` is used, arcs will be colored by whether they appear in
  `ct`, `comp` or both.

`ct_panel`

* This can be "top" or "bottom". It determines where to plot ct data.

`interactions` and `interactions2`

* These can be any interactions datatype: "ringmap", "pairmap", "pairprob",
  "shapejump", etc.
* `_filter` defines the filtering scheme to be applied.
* See [interactions guide](../filters.md) for more info.
* `_panel` can be "top" or "bottom", and determines where to plot
  interactions.

`profile`

* This can be any per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]

`annotations`

* This must be a list of valid annotation datatypes.
* `annotation_mode` can be "track" or "vbar".
* See [annotations guide](../annotations.md) for more info.

`labels`

* A list of strings, one for each sample. Defaults to using the sample name.

`region`

* You can specifiy a start and end position to plot. e.g. [1, 100] plots
  nucleotides 1 through 100, inclusive.
* Defaults to plotting the entire RNA.

`prefiltered`

* True or False. Preserves previous filtering schemes.

`colorbar`

* True or False. Display the color scale for interactions data.

`seqbar`

* True or False. Display the sequence along the x-axis.

`title`

* True or False. Display a title (usually the sample name).
