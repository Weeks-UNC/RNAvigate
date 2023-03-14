Loading data
============

RNAvigate is built around the `Sample` object, which is used to define the
data files that are associated with a particular experimental sample or a
particular RNA. Creating a `Sample` and providing data files will load each of
your data files and store the data as Python objects. The Data objects allow
access to all of the visualization and analysis tools RNAvigate has to offer.
Here is how to create a sample:

```python
import rnavigate as rnav

my_sample = rnav.Sample(
    sample="My Sample Name",
    datatype="Path/to/my/data.file")
```

`import rnavigate as rnav` imports the RNAvigate module and gives you access to
all of its functionality using the alias `rnav`.

`my_sample` is the python variable that your sample is stored to.

`rnav.Sample()` calls the `Sample` initialization method, creating a new sample.

`sample="My Sample Name"` is always the first argument-value pair, it defines
the label that will be used in legends and figure titles when plotting data
from `my_sample`. It should succinctly describe the experiment and/or the RNA.
e.g., "RNaseP DMS-MaP".

`datatype="Path/to/my/data.file"` loads a data file and adds it to
`my_sample.data`, a Python dictionary. You can access this data using
`my_sample.data["datatype"]`. However, `datatype` is just a placeholder, not a
valid argument. Below are all of the valid arguments and their default values.
Further below are explainations of the values that can be used with each
argument. As support for new data types are added, this list will grow.

```python
my_sample = rnav.Sample(
    sample=None,
    # This is a comment. This section contains general RNA structure data
    inherit=None,
    pdb=None,
    ct=None,
    compct=None,
    ss=None,
    fasta=None,
    # This section contains per-nucleotide experimental measurements.
    log=None,
    shapemap=None,
    dmsmap=None,
    dancemap=None,
    rnpmap=None,
    # This section contains inter-nucleotide experimental measurements or data.
    ringmap=None,
    shapejump=None,
    pairmap=None,
    allcorrs=None,
    pairprob=None,
    allpossible=None,
    # This section contains sequence annotations.
    sites=None,
    spans=None,
    groups=None,
    primers=None,
    orfs=None,
    motif=None,
    # dance_prefix loads all data files from a DANCE model.
    dance_prefix=None,
)
```

RNA structure data
------------------

`inherit`

* Another `rnav.Sample` object.
* The new sample will inherit all of the data associated with the given sample.
* These data objects are shared, so any operation on one sample will affect the
  other.
* This saves on memory and computation time if large data sets are to be shared
  between samples.

`pdb`

* A dictionary containing a PDB or CIF file, a chain ID, and optionally a fasta
  reference sequence matching the PDB model.
* e.g., `pdb={"filepath": "my_rna.cif", "chain": "A", "fasta": "my_rna.fasta"}`

`ct`, `compct`, and `ss`

* A .ct (connection table), .dbn (dot-bracket) file containing secondary
  structure base-pairing information.
* A .cte or .nsd file from StructureEditor
* A .varna file from VARNA
* A .xrna file from XRNA
* A .json file from R2DT
* The default behavior of RNAvigate functions expects `ss` to contain
  structure drawing coordinates (cte, nsd, varna, xrna, or json files).

`fasta`

* A .fasta or .fa file containing a single sequence.
* Not usually necessary, but can be useful if you want to map your data onto a
  different sequence. All data are mapped to a common sequence prior to
  visualization. Some functions allow specifying a custom sequence using an
  optional `seq_source` argument.

Per-nucleotide data
-------------------

`log`

* A shapemapper_log.txt file that is output from ShapeMapper2. This file
  includes important quality control metrics for the experiment if the `--per-read-histograms` flag was used with ShapeMapper2.

`shapemap` or `dmsmap`

* A profile.txt file which is the main data output from ShapeMapper2. It
  includes all of the most important information for a SHAPE or DMS-MaP
  experiment.
* If `dmsmap` is used, the reactivity profile is renormalized according to DMS
  conventions.

`dancemap`

* A dictionary containing a DanceMapper _reactivities.txt file and a model
  component number.
* e.g., `dancemap={"filepath": "my_experiment_reactivities.txt", "component": 0}`

`rnpmap`

* An rnpmap.csv file from RNP-MaP

Inter-nucleotide data
---------------------

`ringmap`

* The only default output file of RingMapper. The extension is specified when
  executing RingMapper, typically rings.txt.

`shapejump`

* A dictionary containing a data file and a fasta file.
* e.g. `shapejump={"filepath":"Path/to/data.file", "fasta":"Path/to/fasta.fa"}`

`pairmap`

* A pairmap.txt file from PairMapper.

`allcorrs`

* An allcorrs.txt file from PairMapper.

`pairprob`

* A .dp plain text file from SuperFold or from RNAStructure's ProbabilityPlot.

`allpossible`

* A key from `sample.data` to take a sequence from.
* This will contain a list of every possible pair of nucleotides.

Sequence annotations
--------------------

`sites`

* A dictionary containing:
  * a sequence source (a key of `my_sample.data`)
  * a list of positions of interest within the sequence, 1-indexed
  * a color value to use with visualizations
```python
    sites={
      "seq_source": "shapemap",
      "annotations": [10, 14, 20],
      "color": "red"}
```

`spans`

* A dictionary containing:
  * a sequence source (a key of `my_sample.data`)
  * a list of regions of interest within the sequence, using start and end
    positions, 1-indexed, inclusive.
  * a color value to use with visualizations
```python
    spans={
        "seq_source": "shapemap",
        "annotations": [[10, 20], [53, 87]],
        "color": "blue"}
```

`groups`

* A dictionary containing:
  * a sequence source (a key of `my_sample.data`)
  * a list of dictionaries containing:
    * nucleotide positions within the sequence to be grouped, 1-indexed
    * a color value to use with visualizations
```python
    groups={
        "seq_source": "shapemap",
        "annotations": [
            {"sites": [10, 14, 20], "color": "red"},
            {"sites": [40, 43, 48, 70], "color": "blue"}
        ]}
```

`primers`

* `primers` acts just like `spans` above, but the start and end positions are
  directional. i.e. for reverse primers start is greater than end.

`orfs`

* `orfs` acts just like `spans`, but the regions are automatically identified.
  Every in-frame start and stop codon is added as a span. Honestly, this is too
  common an occurance for this to be useful.
* A dictionary containing:
  * a sequence source (a key of `my_sample.data`)
  * a color value to use with visualizations
```python
    orfs={
        "seq_source": "shapemap",
        "color": "orange"}
```

`motif`

* `motif` acts just like `spans`, but the regions are automatically identified.
  Every match to the given sequence motif is added as a span. Unlike `orfs`,
  this can be very useful.
* A dictionary containing:
  * a sequence source (a key of `my_sample.data`)
  * a sequence motif to search for using the standard alphabet:
    * A, U, C, G, T matches that nucleotide
    * B = not A
    * D = not C
    * H = not G
    * V = not U or T
    * W = weak (A, U, T)
    * S = strong (C, G)
    * M = amino (A, C)
    * K = ketone (G, U, T)
    * R = purine (A, G)
    * Y = pyrimidine (C, U, T)
    * N = any (A, U, C, G, T)
  * a color value to use with visualizations
```python
    motif={
        "seq_source": "shapemap",
        "motif": "DRACH"
        "color": "orange"}
```

Convenient function for DANCE data
----------------------------------

`dance_prefix`

* The file prefix provided to DanceMapper. RNAvigate will locate each of the
  following for each component *n* of the model if it exists:
    * prefix_reactivities.txt, prefix_*n*\_rings.txt, prefix_*n*_pairmap.txt,
      prefix_allcorrs.txt
* It will also find the outputs of foldClusters, if you gave these outputs the
  same prefix:
    * prefix_*n*.f.ct (if `--pk` was used)
    * prefix_*n*.dp (if `--probs` was used)
    * prefix_*n*.ct (if default folding was performed)
* These data are stored in `my_sample.dance` as a list of `rnav.Sample` objects,
  one for each component of the DANCE model.
