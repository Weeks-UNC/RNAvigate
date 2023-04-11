Loading custom data
===================

After creating an `rnav.Sample` object and loading data using the standard
method (see [loading data](../loading-data.md)), adding additional data is
straightforward and very flexible using the `rnav.Sample.set_data()` method.

One common example is with annotations. It is common to want multiple
annotations of the same type, but with different colors/names, associated with
a sample object.

Here, we call the `set_data` method with the required keywords and generic
values, these will be discussed further below. `sample` is a hypothetical
`rnav.Sample` object.

```
sample.set_data(
    name="my_custom_data",
    filepath="path/to/data/file.txt",
    instantiator=rnav.data.Data,
    seq_source="self",
    **kwargs)
```

`name`

* This variable is set to an arbitrary string of your choosing.
* This string becomes a key to the `sample.data` and `sample.inputs`
  dictionaries, and can be used with plotting functions that require a key
  of the `sample.data` dictionary.

`filepath`

* This is a string that contains the path to the data file you are loading.
* If a list is passed, then a new data object is created for each file in the
  list. The keys to `sample.data` are then `"name_1"`, `"name_2"`, etc.

`instantiator`

* This is an `rnav.data.Data` class or subclass to be called to instantiate
  (create) the new data object. One of the following (or your own class):
  * `rnav.data.Data`
  * `rnav.data.Profile`
  * `rnav.data.SHAPEMaP`
  * `rnav.data.DanceMaP`
  * `rnav.data.DeltaProfile`
  * `rnav.data.RNPMaP`
  * `rnav.data.CT`
  * `rnav.data.DotBracket`
  * `rnav.data.XRNA`
  * `rnav.data.VARNA`
  * `rnav.data.NSD`
  * `rnav.data.CTE`
  * `rnav.data.Interactions`
  * `rnav.data.RINGMaP`
  * `rnav.data.PAIRMaP`
  * `rnav.data.PairProb`
  * `rnav.data.SHAPEJuMP`
  * `rnav.data.Log`
  * `rnav.data.PDB`
  * `rnav.data.Annotation`
  * `rnav.data.Motif`
  * `rnav.data.ORFs`
* use `help(rnav.data.Data)` to get info for the Data class, or any other class

`seq_source`

* A sequence string, data object, a key of `sample.data`, or `"self"`. This
  sequence becomes the sequence associated with the new data object.
* `"self"` only applies when the file provides the sequence, e.g. ShapeMapper2
  profiles.

`**kwargs`

* Any additional keyword-argument pairs are passed to `instantiator`.
* Use the `help()` function with a data class to see appropriate arguments
