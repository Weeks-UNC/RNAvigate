Interactive 3D molecular renderings
===================================

Interactive 3D molecular renderings are useful to visualize per-nucleotide and
inter-nucleotide data on 3D models of RNA. Per-nucleotide data may be used to
color nucleotides, and inter-nucleotide data may be drawn as cylinders that
connect nucleotides through space.

These plots are interactive. Here are the controls:
* scroll up/down to zoom in/out
* click and drag to rotate
* 3rd mouse button and drag to pan

There are two ways to quickly render 3D structures:

```python
plot1 = sample.plot_mol()
plot2 = rnavigate.plot_mol(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
single panel rendering using the data from `sample`. With the second method,
you can also create a multi-panel rendering by passing multiple samples to the
samples argument. e.g.:

```python
plot3 = rnavigate.plot_mol(samples=[sample, another_sample])
```

Using either method, you can explore multiple filtering schemes per sample
using the `filters` argument, explained below. Below are all of the optional
arguments that work with each of the methods above, along with their default
values. `plot4` below would produce exactly the same result as `plot1` and
`plot2`.

```python
plot4 = sample.plot_mol(
    structure="pdb",
    interactions=None,
    interactions_filter={},
    filters=None,
    profile="profile",
    nt_color="grey",
    labels=None,
    plot_kwargs={
      "rows": None, "cols": None,
      "width": 400, "height": 400,
      "background_alpha": 1,
      "rotation": {"x":0, "y":0, "z":0}
      "orientation": None},
    hide_cylinders=False,
    title=True,
    get_orientation=False,
    atom="O2'",
    colorbar=True,
    custom_function=None,
    show=True,
)
```

Many of these arguments accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data files are
loaded. To see these keys, run: `print(sample.data.keys())`.

---

`structure`

* a `sample.data` key that points to a 3D structure, usually `"pdb"`.

---

`interactions`

* a `sample.data` key that points to inter-nucleotide data, e.g.:
  `"ringmap"`, `"pairmap"`, `"pairprob"`, `"shapejump"`, etc.
* These data are mapped to `structure`, filtered using the arguments below,
  then plotted as cylinders.

---

`interactions_filter`

* A dictionary of key-value pairs that specifies how `interactions` are
  filtered and displayed.
* See [interactions guide](../guides/filters.md) for more detail.

---

`filters`

* A list of dictionaries specifying interactions data and filtering schemes.
* Each filtering scheme will be applied to each sample, and plotted on separate
  plots. e.g. 3 samples and 2 filtering schemes produces 6 plots.
* This is an alternative to `interactions` and `interactions_filter`, those
  arguments will be ignored.
* See [interactions guide](../guides/filters.md) for more detail.

---

`profile`

* A `sample.data` key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]
* These data are mapped to `structure`, and used to color nucleotides if 
  `nt_color` argument is set to `"profile"`.

---

`nt_color`

* Can be any valid matplotlib color (`"grey"` is a good choice) or one of the
  following:
* `"position"` colors by position using a spectrum.
* `"sequence"` colors by nucleotide.
    * A: red
    * U: light red
    * G: blue
    * C: light blue
* `"profile"` colors using `profile` data.

---

`labels`

* A list of strings, one for each sample, used as titles.
* Defaults to using the sample name by retrieving `sample.sample`.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.Mol()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"rows"` and `"cols"` are integers specifying the number of axes rows and columns of the py3dmol viewer.
* `"width"` and `"height"` specify the size of each panel in pixels.
* `"background_alpha"` specifies the background opacity, `0` is invisible.
* `"rotation"` specifies x, y, and z axis rotation in degrees.
* `"orientation"` specifies an exact orientation of the molecule.

---

`hide_cylinders`

* `True` or `False`, whether to hide cylinders representing nucleotides and
  only show the backbone ribbon.

---

`title`

* `True` or `False`, whether to display titles using `labels`.

---

`get_orientation`

* `True` or `False`, if `True`, when an structure is clicked on, the
  orientation array is displayed. This list of float values can be used with
  the `"orientation"` key of `plot_kwargs` to set the structure to this
  orientation.

---

`atom`

* A string specifying the atom from which to draw interactions.
* `"DMS"` will use the DMS reactive `"N1"` (A, C) and `"N3"` (U, G) atoms.

---

`colorbar`

* `True` or `False`. Display the color scale for interactions data.

---

`custom_function`

* a user defined function to be called to alter molecule appearance.

---

`show`

* `True` or `False`, whether to show the plot or not.
