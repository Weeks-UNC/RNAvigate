Filtering arguments
===================

Any visualization function that can display inter-nucleotide data has arguments
to specify which inter-nucleotide data to use, how to filter them, and how to
display them. These visualizations include arc plots, circle plots, secondary
structure diagrams, 3D molecules, distance distributions, and heatmaps.

To display a single set of filtered inter-nucleotide data, there
are two arguments:

* `interactions` is a key of `sample.data` specifying the inter-nucleotide data
  set. e.g, `interactions="ringmap"` displays RING-MaP data.
* `interactions_filter` is a dictionary that controls how inter-nucleotide data
  are filtered and displayed. e.g,
  `interactions_filter={"cmap": "Greens", "metric": "Zij", "ss_only": True}`
  displays the Zij metric of each RING-MaP correlation on a green color scale,
  only if the correlation connects single-stranded nucleotides.

To display multiple data sets or filtering schemes, there is one argument:

* `filters` is a list of dictionaries specifying the inter-nucleotide data,
  display options, and filtering scheme. e.g.,

```python
filters=[{
  "interactions": "ringmap",
  "cmap": "Greens",
  "metric": "Zij",
  "ss_only": True
  }, {
  "interactions": "shapejump",
  "metric": "Percentile",
  "Percentile_ge": 0.95
  }]
```
* This set of filters would display the same one as the `interactions` and
  `interactions_filter` example above, but would make a second plot with the
  top 5th percentile of a SHAPE-JuMP data set colored by percentile.

Here is every keyword that can be passed to `interactions_filter` or `filters` and it's default value:
```python
interactions_filter = {
  # display arguments
  "cmap": None,
  "metric": None,
  "min_max": None,
  # filtering arguments
  "prefiltered": False,
  # sequence based filters
  "compliments_only": False,
  "nts": None,
  "max_distance": None,
  "min_distance": None,
  "exclude_nts": None,
  "isolate_nts": None,
  # per-nucleotide data filters
  "profile": None,
  "min_profile": None,
  "max_profile": None,
  # secondary structure filters
  "ct": None,
  "min_cd": None,
  "max_cd": None,
  "paired_only": False,
  "ss_only": False,
  "ds_only": False,
  # experimental method to resolve conflicts
  "resolve_conflicts": None,
  # RING-MaP specific filters
  "positive_only": False,
  "negative_only": False,
  # PAIR-MaP specific filters
  "all_pairs": False,
  # additional options to filter by inter-nucleotide data values
  **kwargs
}
```

Below is a breif explanation of all of the available filters. I will use the
short hand *i* and *j* to describe the 5' and 3' ends of an interaction.

Color display arguments
-----------------------

`"metric"`

* The name of a column from the interactions data. Defaults are:
  * `"Statistic"` for RING-MaP data
  * `"Class"` for PAIR-MaP data
  * `"Probability"` for pairing probabilities
  * `"Percentile"` for SHAPE-JuMP data
* The colors of each *i*-to-*j* pair will be taken from this value.
* The column name can also be `"Distance"`, which will be set to the 3D
  distance of the `"O2'"` atoms of *i* and *j* in the given PDB file.
* You can specifiy the atom that distances are calculated from by appending
  `"_atom"` , e.g. `"Distance_O3'"`
* `"Distance_DMS"` will use `"N1"` for A and C, and `"N3"` for U and G.

`"cmap"`

* This can be the name of a colormap, a predefined colormap, a list of colors,
  or a single color.
* Color and colormap names must be valid for matplotlib. e.g., `"red"`,
  `"#f0f0f0"`, `"DodgerBlue"`, `(0, 0.5, 0.5)`, etc. Search the web for
  "matplotlib colormaps" for colormap names.
* Defaults are:
  * `"bwr"` (blue to white to red, continuous) for RING-MaP data
  * `matplotlib.colors.ListedColormap([[0.3, 0.3, 0.3, 0.2], [0.0, 0.0, 0.95, 0.6], [0.12, 0.76, 1.0, 0.6]])` (grey, light blue, dark blue) for PAIR-MaP class.
  * `seaborn.cubehelix_palette(10, 0.7, 0.9, 1.5, 2.5, 1, 0.4, False, True)` for pairing probabilities.
  * `"jet"` (dark blue to light yellow to red full spectrum) for 3D distances.
  * `"YlGnBu"` (yellow to green to blue) for SHAPE-JuMP data.

`"min_max"`

* A list of lower and upper limits to apply the colormap, e.g.,
  `"min_max": [-20, 20]`.
* In this example, values from `"metric"` that are less than or equal to `-20`
  are given the first color value from the colormap. Values greater than or
  equal to `20` are given the last color value. Within these bounds, the colors
  values are evenly disbursed.

Filtering arguments
-------------------

`"prefiltered"`

* `True` or `False`, whether to retain previously applied filters.
* Use this if you are manually filtering data.
* Do not use it to avoid re-typing a long list of filters. Instead, store a
  dictionary in a python variable, and pass the variable name to
  `interactions_filters`

Filtering on sequence
---------------------

`"compliments_only"`

* Set to `True` to only keep an interaction if *i* is the reverse compliment of
  *j*

`"nts"`

* Set to a string of valid nucleotides, e.g., to exclude interactions in which
  *i* or *j* contains a G nucleotide `"nts": "AUC"`

Filtering by position
---------------------

`"max_distance"` and `"min_distance"`

* A minimum and maximum primary sequence distance between *i* and *j*.

`"exclude_nts"` and  `"isolate_nts"`

* A list of nucleotide positions to exclude or to isolate (meaning all other
  positions are excluded).

Filtering on per-nucleotide data
--------------------------------

`"profile"`

* A key of `sample.data`, values to be filtered will be taken from
  `profile.default_column`, usually normalized reactivity.

`"min_profile"` and `"max_profile"`

* The minimum and maximum allowable values of per-nucleotide data at *i* and
  *j*


Filtering on structure
----------------------

`"ct"`

* a key of `sample.data` pointing to a secondary structure. Contact distances
  and base-pairing status for the below filters are taken from this data.

`"min_cd"` and `"max_cd"`

* A minimum or maximum contact distance between *i* and *j*. Contact distance
  is the shortest path distance through the secondary structure graph.

`"paired_only"`

* `True` or `False`, whether to require that *i* and *j* are base-paired to
  each other.

`"ss_only"` and `"ds_only"`

* `True` or `False`, whether to require that *i* and *j* are both
  single-stranded or both double-stranded, respectively.

RING-MaP specific filters
-------------------------

`"positive_only"` and `"negative_only"`

* `True` or `False`, whether to require that RING-MaP correlation direction is
  positive or negative.

PAIR-MaP specific filters
-------------------------

`"all_pairs"`

* `True` or `False`, whether to include complimentary PAIRs that are not
  classified as Primary or Secondary PAIRs.

Filtering on data values
------------------------

`"ColumnName_operator"`

* `"ColumnName"` is a valid column of interactions data
* `"operator"` is one of:
  * `_ge` (greater that or equal to)
  * `_le` (less than or equal to)
  * `_gt` (greater than)
  * `_lt` (less than)
  * `_eq` (equal to)
  * `_ne` (not equal to)
* Interactions data from this column is compared to the given value using the
  given operator. Interactions where that operation is `True` are kept.

Experimental filter
-------------------

`"resolve_conflicts"`

* `True` or `False`, whether to reduce the number of interactions by finding
  the maximum weighted set of non-conflicting interactions (parallel
  correlations do not conflict)
