Filtering arguments
===================

Any plot type that accepts `interactions` also accepts `interactions_filter`,
which controls which interactions get plotted as well as the color gradient used
to plot them. Here is what that looks like:

```python
plot = sample.plot_plottype(
    interactions="datatype",
    interactions_filter={
        "filter_1": some_value,
        "filter_2": some_value})
```

* `plottype`s that accept interactions: `ss`, `arc`, `mol`, `circle`, `disthist`, `heatmap`
* built-in `interactions` datatypes are: `ringmap`, `pairmap`, `allcorrs`, `pairprob`

Here is a reference of all of the available arguments for `interactions_filter`:

Filtering on sequence

* `"compliment_only": True`: only include interaction if *i* is the reverse compliment of *j*
* `"nts": "AUG"`: only include interactions consisting only of *A*s, *U*s, and
  *G*s

Filtering on structure

* `"cdAbove": #`: minumum *i*-to-*j* contact distance in the structure graph must be greater than #
* `"cdBelow": #`: minimum *i*-to-*j* contact distance must be less than #
* `"ss_only": True`: include only if *i* and *j* are both unpaired (only works for 1-nt windows)
* `"ds_only": True`: include only if *i* and *j* are both paired (only works for 1-nt windows)
* `"paired_only": True`: include only if *i* and *j* are base-paired to each other

Filtering by profile

* `"profAbove": #`: include only if Normalized SHAPE reactivity is greater than #
* `"profBelow": #`: include only if Normalized SHAPE reactivity is less than #

Filtering by position

* `"exclude_nts": [#, #, #]`: exclude interaction if *i* or *j* is in this list of nucleotide positions.
* `"isolate_nts": [#, #, #]`: include only if each *i* and *j* is entirely within list of numbers.

Filtering on data values

* `"ColumnName_operator": #`

    * Data from "ColumnName" of data frame must be greater than or equal to #
    * ColumnName is a valid column of the data frame. To see valid column names
      for `"datatype"`, type `sample_name.valid_columns("datatype")`.
    * operator is one of:
        * `_ge` (greater that or equal to)
        * `_le` (less than or equal to)
        * `_gt` (greater than)
        * `_lt` (less than)
        * `_eq` (equal to)
        * `_ne` (not equal to)

Other filters

* `"all_pairs": True`: include all complimentary PAIRs along with primary and secondary PAIRs
* `"positive_only": True`: include only if "+/-" column of RING-MaP file is equal to 1
* `"negative_only": True`: include only if "+/-" column of RING-MaP file is equal to -1
* `"resolve_conflicts": True`: Finds the maximum weighted set of non-conflicting interactions (parallel correlations do not conflict)

Arguments that control the color gradient

* `"metric": "ColumnName"`
    * The colors of each ij pair will be taken from this value.
    * "ColumnName" is a valid column name of the data
    * "ColumnName" can also be "Distance", which will be set to the 3D distance of the "O2`" atoms of *i* and *j* in the given PDB file.
    * You can specifiy the atom that distances are calculated from by appending "_atom_name" , e.g. "Distance_DMS" or "Distance_O3`"
    * For DMS data, use "Distance_DMS" which will use "N1" for A and G, and "N3" for U and C.
* `"cmap": "mpl_cmap_name"`
    * The colormap used to display ij data will be retreived from matplotlib.
    * Can also be a list of matplotlib colors, and a continuous colormap will be created.
    * Can also be a single matplotlib color. All values would be displayed in that color.
    * Special colormaps are used for PAIR-MaP "Class" data and "Distance" data.
* `"min_max": [#, #]`
    * The colormap will start at the first # and end at the second #.
    * Data outside this range will be plotted using the first or last color value.
