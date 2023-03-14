Annotations
===========

Currently, you cannot load annotation data using a file. You have to provide
this type of data as a python object. There are several convenience arguments
to create annotations at Sample creation: `sites`, `spans`, `groups`, `orfs`,
and `motif`.

```python
my_sample = rnavigate.Sample(
    sample="Annotations example",
    fasta="Path/to/fasta.fa",
    spans={"seq_source": "fasta",
           "span_list": [[10, 15], [35, 45]],
           "color": "grey"},
    sites={"seq_source": "fasta",
           "site_list": [12, 27, 48, 92],
           "color": "green"},
    groups={"seq_source": "fasta",
            "groups": [
                {"color": "blue",
                 "sites": [12, 13, 19, 20]},
                {"color": "red",
                 "sites": [45, 47, 63, 64, 65]}
                 ]
            },
    orfs={"seq_source": "fasta",
          "color": "grey"},
    motif={"seq_source": "fasta",
           "motif": "DRACH",
           "color": "yellow"})
```

For all annotation types:

* `"seq_source"` indicates the data from which the sequence will be taken.
* `"color"` indicates the color that will be used on plots to show this data.

`spans`

* These are regions of interest in the RNA sequence.
* `"span_list"` is a list of start and end positions for each span.
* In the example above, spans would be drawn from 10 to 15 and 35 to 45, inclusive.

`sites`

* These are nucleotide positions of interest in the RNA sequence.
* `"site_list"` is a list of nucleotide positions.

`groups`

* These are lists of dictionaries. Each dictionary contains:
    * `"sites"`: a list of related nucleotide positions
    * `"color"`: the color to use for plotting that group

`orfs`

* This will create an annotation that behaves like `spans`, in which each span
  is an open reading frame.

`motif`

* This will also create an annotation that behaves like `spans`, in which each
  span matches a given sequence motif.
* `"motif"` is a sequence motif:
    * A, U, C, G: matches the exact nucleotide.
    * B/V/D/H: not A/U/C/G respectively
    * W: A or U (weak)
    * S: C or G (strong)
    * M: A or C (amino)
    * K: G or U (ketone)
    * R: A or G (purine)
    * Y: C or U (pyrimidine)
    * N: any nucleotide
    * e.g. `"motif": "DRACH"` would match 18 possible 5-mers.
