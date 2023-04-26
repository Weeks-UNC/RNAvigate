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
  * `rnav.data.Data`: base classs to handle sequences for all subclasses
    * `rnav.data.Log`: Shapemapper2 log files
    * `rnav.data.Profile`: base class for per-nucleotide data
      * `rnav.data.SHAPEMaP`: Shapemapper2 profile.txt
      * `rnav.data.DanceMaP`: DanceMapper component from reactivities.txt
      * `rnav.data.RNPMaP`: RNPMapper data
      * `rnav.data.DeltaProfile`: differences between profiles
    * `rnav.data.Annotations`: base class for sequence annotations
      * `rnav.data.Motif`: regions matching a sequence motif
      * `rnav.data.ORFs`: regions of open-reading frames
    * `rnav.data.Interactions`: base class for inter-nucleotide data
      * `rnav.data.RINGMaP`: RingMapper rings.txt or PairMapper allcorrs.txt
      * `rnav.data.PAIRMaP`: PairMapper pairmap.txt
      * `rnav.data.SHAPEJuMP`: ShapeJumper output file
      * `rnav.data.PairProb`: pairing probabilities file
      * `rnav.data.AllPossible`: all possible nucleotide pairs in a sequence
    * `rnav.data.CT`: base class for secondary structures and drawings
      * `rnav.data.DotBracket`: secondary structure in .db format
      * `rnav.data.XRNA`: drawing in XRNA format
      * `rnav.data.VARNA`: drawing in VARNA format
      * `rnav.data.R2DT`: drawing in R2DT json format
      * `rnav.data.NSD`: drawing in StructureEditor NSD format
      * `rnav.data.CTE`: drawing in StructureEditor CTE format
    * `rnav.data.PDB`: 3D atomic coordinates in RCSD PDB or CIF format
* use `help(rnav.data.Data)` to get info for the Data class, or any other class

`seq_source`

* A sequence string, data object, a key of `sample.data`, or `"self"`. This
  sequence becomes the sequence associated with the new data object.
* `"self"` only applies when the file provides the sequence, e.g. ShapeMapper2
  profiles.

`**kwargs`

* Any additional keyword-argument pairs are passed to `instantiator`.
* Use the `help()` function with a data class to see appropriate arguments
