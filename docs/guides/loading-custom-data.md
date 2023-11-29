Loading custom data
===================

After creating an `rnav.Sample` object and loading data using the standard
method (see [loading data](../loading-data.md)), adding additional data is
is accomplished using the `rnav.Sample.set_data()` method.

Here, we call the `set_data` method with the required keywords and generic
values, these will be discussed further below. `sample` is a hypothetical
`rnav.Sample` object.

```
my_sample.set_data(
    name="my_custom_data",
    filepath="path/to/data/file.txt",
    data_class="data_keyword",
    **kwargs)
```

`name`

* This variable is set to an arbitrary string of your choosing.
* This string becomes a data keyword and can be provided to plotting functions.
* Data keywords can also be used to access the underlying data:
  - `my_sample.data["my_custom_data"]`

`filepath`

* This is a string that contains the path to the data file you are loading.

``

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
    * `rnav.data.SecondaryStructure`: secondary structures and drawings
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
