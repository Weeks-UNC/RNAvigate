3D plots
========

3D plots are useful to visualize interactions and profile data on 3D models
of RNA. Profile data is used to color nucleotides, and interactions are
plotted as cylinders that connect two nucleotides through space.

Current limitations
-------------------
* Plots don't appear if the first O5' atom is missing.
  * This is an issue in 3Dmol.js and is not likely to be fixed.
* For now, plots can only display one chain at a time.
* Saving as png is easy in Jupyter, but can't currently be done in VS code.

These plots are interactive. Here are the controls:
* scroll up/down to zoom in/out
* click and drag to rotate
* 3rd mouse button and drag to pan

There are three functions for making these plots:

```python
plot1 = sample.plot_mol()
plot2 = sample.plot_mol_multifilter(filters=filters)
plot3 = rnavigate.plot_mol_multisample(samples=samples)
```

Below are all of the optional arguments that work with each of the functions
above, along with their default values.

```python
plot = sample.plot_mol(
    structure="pdb"
    interactions=None,
    interactions_filter={},
    profile="profile",
    labels=None,
    show=True,
    prefiltered=False,
    nt_color="sequence",
    atom="O2'",
    title=True)
```

`structure`

* This must be "pdb" or a user-defined PDB object.

`interactions`

* This can be any interactions datatype: "ringmap", "pairmap", "pairprob",
  "shapejump", etc.
* `_filter` defines the filtering scheme to be applied.
* See [interactions guide](../filters.md) for more info.

`profile`

* This can be any per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]

`labels`

* A list of strings, one for each sample. Defaults to using the sample name.

`title`

* True or False, whether to display a title using the above defined label.

`show`

* True or False, whether to show the plot or not.

`prefiltered`

* True or False. Preserves previous filtering schemes.

`nt_color`

* Can be any valid matplotlib color ("grey" is a good choice) or one of the
  following:
* `"position"` colors by position using a spectrum.
* `"sequence"` colors by nucleotide.
    * A: red
    * U: light red
    * G: blue
    * C: light blue
* `"profile"` assigns colors using profile data.

`atom`

* The atom from which to draw interactions. Defaults to `"O2'"`.
