Loading data
============

RNAvigate is built around the `Sample` object, which defines the datasets
that are associated with a particular experimental sample. The first step to
quickly performing vizualization and analysis is to create a `Sample` object.
Here is the basic call signiture:

```python
import rnavigate as rnav

my_sample = rnav.Sample(
    sample="My Sample Name",
    datatype="Path/to/my/data.file")
```

Here `sample="My Sample Name"` defines the label that will be used in legends
and figure titles when plotting data from `my_sample`. It should succinctly
describe the experiment. e.g. "RNaseP DMS-MaP"

The following are all of the datatypes that can be loaded into `my_sample`.

`log`

* A shapemapper_log.txt file that is output from ShapeMapper2. This file
  includes important quality control metrics for the experiment if the `--per-read-histograms` was used.

`shapemap` or `dmsmap`

* A profile.txt file which is the main data output from ShapeMapper2. It
  includes all of the most important information for a SHAPE or DMS-MaP experiment.
* If `dmsmap` is used, the reactivity profile is renormalized according to DMS conventions.

`ringmap`

* The only default output file of RingMapper. The extension is specified when
  executing RingMapper, typically rings.txt.

`pairmap`

* A pairmap.txt file from PairMapper.

`allcorrs`

* An allcorrs.txt file from PairMapper.

`dance_prefix`

* The file prefix provided to DanceMapper. MaP.Sample will locate each of the
  following for each component *n* of the model if it exists:
    * prefix_reactivities.txt, prefix_*n*\_rings.txt, prefix_*n*_pairmap.txt,
      prefix_allcorrs.txt
* It will also find the outputs of foldClusters, if you gave these outputs the
  same prefix:
    * prefix_*n*.f.ct (if `--pk` was used)
    * prefix_*n*.dp (if `--probs` was used)
    * prefix_*n*.ct (if default folding was performed)

`rnpmap`

* An rnpmap.csv file from RNP-MaP

`shapejump`

* Instead of a file path, this argument requires a python dictionary containing
  a data file and a fasta file.
* e.g. `shapejump={"filepath":"Path/to/data.file", "fasta":"Path/to/fasta.fa"}`

`pairprob`

* A .dp file from SuperFold or from ProbabilityPlot.

`ct` and `compct`

  * A .ct (connection table) or .dbn (dot-bracket) containing base-pairing info.

`ss`

* A .cte or .nsd file from StructureEditor
* A .varna file from VARNA
* A .xrna file from XRNA

`pdb`

* Instead of a file path, this argument requires a python dictionary containing
  a pdb file and the chain ID of the molecule of interest.
* e.g. `pdb={"filepath":"Path/to/My_RNA.pdb", "chain": "A"}`
* In addition to `filepath` and `chain`, the dictionary may contain the path to
  a fasta file, and the offset. This is only needed if that information is
  missing from the PDB header.
* Offsets are a weird feature of PDBs. Most of the time there is no offset.
  e.g., residue_id 10 in the atom entries normally corresponds to the 10th
  nucleotide in the sequence. If it corresponds to the 1st nucleotide, the offset
  would be 9.
