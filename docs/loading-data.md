Loading data
============

This is part 2 in the getting started with RNAvigate guide.

1. [Installing RNAvigate](./installing-rnavigate.md)
2. [Loading Data](./loading-data.md)
3. [Visualizing data](./visualizing-data.md)

RNAvigate is built around `rnav.Sample`, which defines and organizes the data
that are associated with an experimental sample or RNA. Creating a `Sample` and
providing data files will load and store the data contained within. This step
allows access to all of the visualization and analysis tools RNAvigate offers.

Here is how to create a sample:

```python
import rnavigate as rnav         # Load RNAvigate and give it the alias "rnav"

my_sample = rnav.Sample(         # create a new sample
    sample="My sample name",     # provide a name for plot labels
    data_keyword="my_data.txt")  # load data, data_keywords are explained below
```

Our new variable `my_sample`, now contains the data loaded from `my_data.txt`.
The `data_keyword` used here is not a real data keyword. Data keywords are used
to tell RNAvigate how to parse the data file, and for plotting functions, how
to access the data. `rnav.Sample` can accept any number of data keyword
arguments, thus organizing all of the data related to an experimental sample.

Sections:

- [data keywords](#data-keywords)
- [other optional arguments](#other-optional-arguments)

---

## Data keywords

Data keywords perform two major functions in RNAvigate. When loading data, they
tell RNAvigate how to parse that data. When plotting data or performing
analyses they tell RNAvigate which data to access. I'll use a made-up data
keyword (`keyword`) to illustrate:

`rnav.Sample` initialization:

```python
my_sample = rnav.Sample(
    sample="my sample name",
    keyword="expected_input")
```

`rnav.plot_something` functions:

```python
plot = rnav.plot_skyline(
    samples=[my_sample],
    profile="keyword")
```

Data keywords can either be a standard keyword or an arbitrary one. Arbitrary
data keywords are useful for when multiple data sources of the same kind are
associated with a single sample.

### Standard data keywords

| Sequence     | Annotation  | Profile      | Interactions     | SecondaryStructure | PDB     |
| :----------: | :---------: | :----------: | :--------------: | :----------------: | :-----: |
| [sequence][] | [sites][]   | [profile][]  | [interactions][] | [ss][]             | [pdb][] |
|              | [spans][]   | [shapemap][] | [ringmap][]      |                    |         |
|              | [group][]   | [dancemap][] | [pairmap][]      |                    |         |
|              | [primers][] | [rnpmap][]   | [shapejump][]    |                    |         |
|              | [motif][]   |              | [pairprob][]     |                    |         |
|              | [orfs][]    |              | [allpossible][]  |                    |         |
|              | [domains][] |              |                  |                    |         |

### Arbitrary data keywords

An arbitrary keyword is useful if you are loading 2 or more of the same data
type into a single sample. Arbitrary keywords must follow some simple rules:

1. Cannot already be a data keyword
2. Cannot be `inherit` or `keep_inherited_defaults`
3. Cannot consist only of valid nucleotides: `AUCGTaucgt`
4. Cannot start with a number: `0123456789`
5. Cannot contain any symbols accept underscore: `!@#$%^&`, etc.
6. Cannot contain spaces or tabs

If an arbitrary data keyword is used, a dictionary must be provided, specifying
the standard data keyword to use for parsing inputs, and any other options.

Example:

```python
my_sample = rnav.Sample(
    sample="example",
    standard_keyword="input_file_1.txt",
    arbitrary_keyword={"standard_keyword": "input_file_2.txt"}
)
```

### sequence

an RNA sequence

example uses:

- aligning data between sequences
- all data in RNAvigate is associated with a sequence and can be aligned to
  other data, or vice versa.

input explaination:

- Input should be a fasta file, a sequence string, or another data keyword. If
  another data keyword is provided, the sequence from that data is retrieved.

example inputs:

```python
# fasta file
my_sample = rnav.Sample(
    name="example",
    sequence="path/to/my_sequence.fa",
    )

# sequence string
my_sample = rnav.Sample(
    name="example",
    sequence="AUCAGCGCUAUGACUGCGAUGACUGA",
    )

# data keyword
my_sample = rnav.Sample(
    name="example",
    data_keyword="some_data_with_a_sequence"
    sequence="data_keyword",
    )
```

### profile

per-nucleotide data that does not have a more specific data keyword:

- [shapemap][] for SHAPE-, DMS-, or other MaP method
- [dancemap][] for DanceMapper reactivities
- [rnpmap][] for RNPMapper data
- profile: for everything else

example uses:

- visualizing per-nucleotide data on profile, skyline, or arc plots.
- coloring nucleotides in a secondary structure diagram, circle plot, or 3D
  molecular rendering
- calculating profile-to-profile linear regressions
- Renormalizing per-nucleotide data

input explaination:

- These inputs allow a lot of customization in loading data.
- For a full explaination, see [customizing profiles][]

back to [standard data keywords][]

### shapemap

SHAPE, DMS, or other reagent per-nucleotide reactivities

[ShapeMapper2 software][]

example uses:

- same as [profile][]
- plus: visualizing quality control metrics if a log file is specified

input explaination:

- Input should be a ShapeMapper2 profile.txt file. This file contains the most
  complete per-nucleotide data from a ShapeMapper2 run.

example inputs:

```python
my_sample = rnav.Sample(
    sample="example",
    shapemap="shapemap_profile.txt"
)
```

other optional inputs:

- `"normalize"`:
  - Defaults to not performing any renormalization.
  - `"DMS"` will perform DMS-MaP renormalization
  - `"eDMS"` will perform eDMS-MaP renormalization
  - `"boxplot"` will perform ShapeMapper2 renormalization (with 1 improvement)
  - By default, renormalization is performed on the HQ_profile and HQ_stderr
    columns, and overwrites the Norm_profile and Norm_stderr columns
  - Normalization can also be done after `rnav.Sample` creation.
    - type: `help(rnav.data.Profile.normalize)`
- `"log"` is used to specify a log file. If ShapeMapper2 was run with the
  `--per-read-histograms` flag, this file will contain read length distribution
  and mutations-per-read distribution. This data can then be visualized with
  `rnav.plot_QC`.
- `"metric"`, `"metric_defaults"`, `"sequence"`, and `"read_table_kw"` are
  explained in [customizing profiles][], but are not recommended for standard
  ShapeMapper2 files.

typical optional input example:

```python
my_sample = rnav.Sample(
    sample="example",
    dmsmap={
        "shapemap": "shapemap_profile.txt",
        "normalize": "dms",
    }
)
```

back to [standard data keywords][]

### dancemap

Reactivity profile of a single component of a DanceMapper model

[DanceMapper software][]

example uses:

- same as [profile][]

input explaination:

- Input should be a dictionary containing:
  - `"dancemap"`: the DanceMapper reactivities.txt file
  - `"component"`: which component of the DANCE model to load
- This works best if each component is a seperate `rnav.Sample`

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    dancemap={
        "dancemap": "mydancemap_reactivities.txt",
        "component": 0},
    )
my_sample_2 = rnav.Sample(
    sample="example2",
    dancemap={
        "dancemap": "mydancemap_reactivities.txt",
        "component": 1},
    )
```

other optional inputs:

- `"metric"`, `"metric_defaults"`, `"sequence"`, and `"read_table_kw"` are
  explained in [customizing profiles][], but are not recommended for standard
  DanceMapper files.

back to [standard data keywords][]

### rnpmap

RNP-MaP per-nucleotide reactivities.

[RNPMapper software][]

example uses:

- same as [profile][]

input explaination:

- Input should be the output csv file from RNPMapper

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    rnpmap="myrnpmap_output.csv",
    )
```

other optional inputs:

- `"metric"`, `"metric_defaults"`, `"sequence"`, and `"read_table_kw"` are
  explained in [customizing profiles][], but are not recommended for standard
  RNPMapper files.

back to [standard data keywords][]

### interactions

inter-nucleotide data that does not have a more specific data keyword:

- [ringmap][] for RingMapper correlations
- [pairmap][] for PairMapper correlations
- [shapejump][] for ShapeJump deletion events
- [pairprob][] for pairing probabilities
- [allpossible][] for every possible nucleotide pairing from a sequence
- interactions: for everything else

example uses:

- visualizing interaction networks in arc and circle plots, secondary structure
  diagrams and 3D molecule renderings
- filtering interactions based on many different factors.
  - see [filtering interactions][] guide.
- calculating a distance distribution histogram of a set of interactions

input explaination:

- These inputs allow a lot of customization in loading data.
- For a full explaination, see [customizing interactions][]

back to [standard data keywords][]

### ringmap

single-molecule correlations from a DMS- or eDMS-MaP experiment

[RingMapper software][]

example uses:

- same as [interactions][]

input explaination:

- Input should be the correlations file output from RingMapper.

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    ringmap="myringmap_corrs.txt",
    )
```

other optional inputs:

- `"sequence"` is used to specify the sequence, and accepts the same inputs
  as the [sequence][] keyword
- `"metric"`, `"metric_defaults"`, `"read_table_kw"`, and `"window"` are
  explained in [customizing interactions][], but are generally not recommended
  for RingMapper files.

typical optional input example:

```python
my_sample = rnav.Sample(
    sample="example",
    ringmap={
        "ringmap": "myringmap_corrs.txt"
        "sequence": "my_rna.fa"
        }
    )
```

back to [standard data keywords][]

### pairmap

single-molecule correlations from a DMS- or eDMS-MaP experiment reflective of
base pairing

[PairMapper software][] (part of RingMapper)

example uses:

- same as [interactions][]

input explaination:

- Input should be the pairmap.txt output file from PairMapper.

example inputs:

```python
my_sample = rnav.Sample(
    sample="example",
    pairmap="mydata_pairmap.txt",
    )
```

other optional inputs:

- `"sequence"` is used to specify the sequence, and accepts the same inputs
  as the [sequence][] keyword
- `"metric"`, `"metric_defaults"`, `"read_table_kw"`, and `"window"` are
  explained in [customizing interactions][], but are generally not recommended
  for PairMapper files

typtical optional input example:

```python
my_sample = rnav.Sample(
    sample="example",
    pairmap={
        "pairmap": "mydata_pairmap.txt",
        "sequence": "my_rna.fa",
    }
)

```

back to [standard data keywords][]

### shapejump

ShapeJump inter-nucleotide RT deletion events

[ShapeJumper software][]

example uses:

- same as [interactions][]

input explaination:

- Input should be the deletions.txt output file from ShapeJumper.

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    shapejump="mydata_deletions.txt",
    )
```
other optional inputs:

- `"sequence"` is used to specify the sequence, and accepts the same inputs
  as the [sequence][] keyword
- `"metric"`, `"metric_defaults"`, `"read_table_kw"`, and `"window"` are
  explained in [customizing interactions][], but are generally not recommended
  for ShapeJumper files.

typical optional argument example:

```python
my_sample = rnav.Sample(
    sample="example",
    shapejump={
        "shapejump": "mydata_deletions.txt",
        "sequence": "my_rna.fa",
    }
)
```

back to [standard data keywords][]

### pairprob

inter-nucleotide predicted pairing probabilities

[RNAStructure software][]

example uses:

- same as [interactions][]
- plus: calculating shannon entropy

input explaination:

- Input should be a dotplot plain text file from running RNAstructure
  `partition` followed by `ProbabilityPlot` with `-t` option

```bash
# for example:
partition my_sequence.fa pair_probabilities.dp
ProbabilityPlot -t pair_probabilities.dp pair_probabilities.txt
```

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    pairprob={"pairprob": "pair_probabilities.txt",
              "sequence": "my_sequence.fa"}
    )
```

other optional inputs:

- `"sequence"` is used to specify the sequence, and accepts the same inputs
  as the [sequence][] keyword
- `"metric"`, `"metric_defaults"`, `"read_table_kw"`, and `"window"` are
  explained in [customizing interactions][], but are generally not recommended
  for "" files

typical optional argument example:

```python
my_sample = rnav.Sample(
    sample="example",
    pairprob={
        "pairprob": "pair_probabilities.txt",
        "sequence": "my_rna.fa",
    }
)
```

back to [standard data keywords][]

### allpossible

All possible inter-nucleotide pairings for a given sequence

example uses:

- same as [interactions][]
- plus: calculating the expected distance distribution of a filtering scheme

input explaination:

- This keyword has the same expected inputs as the [sequence][] keyword.
- Note: the size of the data increases with the sequence length squared.

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    allpossible="my_rna.fa"
    )
```

other optional inputs:

- `"window"` is used to specify the window size of the interacting regions.
    - Default: `"window": 1` means nucleotide *i* to nucleotide *j*
    - `"window": 3` means nucleotides *i*:*i*+3 to nucleotides *j*:*j*+3
- `"metric"`, `"metric_defaults"`, and `"read_table_kw"` are explained
  in [customizing interactions][].

typical optional argument example:

```python
my_sample = rnav.Sample(
    sample="example",
    allpossible={
        "allpossible": "my_rna.fa",
        "window": 3,
    }
)
```

back to [standard data keywords][]

### ss

A secondary structure with optional diagram drawing

example uses:

- visualizing base pairs on arc plots, circle plots, and secondary structure
  diagrams
- calculating contact distances
  - the shortest path between nucleotides in a secondary structure graph
- determining how well per-nucleotide data predict base pairing status (AUROC)

input explaination:

- Input should be one of the following formats:
  - secondary structure files (no diagram)
    - connection table (.ct)
    - dotbracket notation (.dot, .dbn, etc.)
  - secondary structure diagram files
    - [StructureEditor][] (.nsd or .cte)
    - [XRNA][] (.xrna)
    - [VARNA][] (.varna)
    - [FORNA][] (.json)
      - click "add molecule" and paste in a dotbracket notation structure.
      - arrange it how you like
      - click the download button in the lower-right, then click "json"
    - [R2DT][] (.json)
      - Type in an RNA sequence, R2DT creates the secondary structure
      - Click on the R2DT paper link to learn more about how it works
      - Once the structure is drawn, click "Edit in XRNA"
      - Arrange it how you like it
      - In the upper-left, type in a file name, choose "json", click "download"
- Note: The file format is determined by the file extension. Since FORNA and R2DT
  both produce json, the extension should be provided.

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    ss="my_structure.ct"
    )
```

other optional inputs:

- `"extension"` is used to specify a file extension.
- `"autoscale"` is used to scale coordinates to look good in RNAvigate plots
- `"structure_number"` is used to specify which structure to load if the file
  contains multiple structures (0-indexed). Default is to load the first
  structure. This currently only works with "ct" files.

typical optional argument examples:

```python
# specify r2dt vs forna json
my_sample = rnav.Sample(
    sample="example",
    ss={
        "ss": "my_rna.json",
        "extension": "r2dt,
    }
)

#specify structure number
my_sample = rnav.Sample(
    sample="example",
    ss={
        "ss": "my_rna.ct",
        "structure_number": 3,
    }
)
```

back to [standard data keywords][]

### pdb

A tertiary structure with atomic coordinates

example uses:

- rendering 3D molecules with data overlayed
- computing 3D distances between nucleotides

input explaination:

- Input should be a dictionary containing these keys:
  - "pdb": a standard [PDB][] file (.pdb or .cif)
  - "chain": A chain ID
  - "sequence": same inputs as [sequence][] keyword
    - This is not needed if a sequence is found in the file header.

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    pdb={
        "pdb": "my_structure.pdb",
        "chain": "X",
        "sequence": "my_rna.fa" # not needed if sequence in pdb header
        }
    )
```

back to [standard data keywords][]

### motif

Annotation of occurances of a sequence motif

example uses:

- highlighting nucleotides in skyline, profile, arc, circle, or secondary
  structure diagram plots
- coloring nucleotides in circle plots, secondary structure diagrams, 3D
  molecule renderings, or linear regression scatter plots

input explaination:

- Input should be a dictionary containing:
  - `"motif"`: a string that uses the nucleotide alphabet
    - e.g.: "DRACH" for potential m6A modification sites

| alphabet   | meaning    | matches    |
| :--------: | :--------: | :--------: |
| A, U, C, G | identity   | A, U, C, G |
| B          | not A      | U/C/G      |
| D          | not C      | A/U/G      |
| H          | not G      | A/U/C      |
| V          | not U      | A/C/G      |
| W          | weak       | A/U        |
| S          | strong     | C/G        |
| M          | amino      | A/C        |
| K          | ketone     | U/G        |
| R          | purine     | A/G        |
| Y          | pyrimidine | U/C        |
| N          | any        | A/U/C/G    |

  - `"sequence"`: same as [sequence][] keyword
  - `"color"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"name"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    m6A={
        "motif": "DRACH",
        "sequence": "my_rna.fa",
        "color": "blue",
        "name": "m6A motif"
        }
    )
```

back to [standard data keywords][]

### orfs

Annotation of open-reading frames

example uses:

- same as [motif][]
- coming soon: displaying amino acid translation and codon usage scores

input explaination:

- Input should be a dictionary containing:
  - `"orfs"`: 
    - `"all"` annotates all open-reading frames
    - `"longest"` annotates only the longest open reading frame
  - `"sequence"`: same as [sequence][] keyword
  - `"color"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"name"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    main_orf={
        "orfs": "longest",
        "sequence": "my_sequence.fa",
        "color": "green",
        "name": "Longest ORF"
        }
    )
```

back to [standard data keywords][]

### spans

Annotation of any regions of interest

example uses:

- same as [motif][]

input explaination:

- input is a dictionary containing
  - `"spans"`: a list of lists of 2 integers. Each inner list specifies a
    start and end position of a span (1-indexed, inclusive)
  - `"sequence"`: same as [sequence][] keyword
  - `"color"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"name"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    regions={
        "spans": [[10, 13], [65, 72]],
        "sequence": "my_sequence.fa",
        "color": "purple",
        "name": "interesting regions"
        }
    )
```

back to [standard data keywords][]

### sites

Annotation of any sites of interest

example uses:

- same as [motif][]

input explaination:

- input is a dictionary containing
  - `"sites"`: a list of nucleotide positions (1-indexed, inclusive)
  - `"sequence"`: same as [sequence][] keyword
  - `"color"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"name"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    m6a_sites={
        "sites": [10, 13, 65, 72],
        "sequence": "my_sequence.fa",
        "color": "purple",
        "name": "m6A sites"
        }
    )
```

back to [standard data keywords][]

### group

Annotation of any group of nucleotides, such as a binding pocket

example uses:

- same as [motif][]

input explaination:

- input is a dictionary containing
  - `"group"`: a list of nucleotide positions (1-indexed, inclusive)
  - `"sequence"`: same as [sequence][] keyword
  - `"color"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"name"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    ligand_pocket={
        "sites": [10, 13, 65, 72],
        "sequence": "my_sequence.fa",
        "color": "purple",
        "name": "ligand-binding pocket"
        }
    )
```

back to [standard data keywords][]

### primers

Annotation of primer binding sites

example uses:

- same as [motif][]

input explaination:

- input is a dictionary containing
  - `"primers"`: a list of lists of 2 integers. Each inner list specifies a
    start and end position of a primer (1-indexed, inclusive). A reverse primer
    is specified by listing the 3' -> 5' start and end, e.g. [300, 278]
  - `"sequence"`: same as [sequence][] keyword
  - `"color"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"name"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    pcr_primers={
        "primers": [[1, 22], [300, 278]],
        "sequence": "my_sequence.fa",
        "color": "purple",
        "name": "primer-binding sites"
        }
    )
```

back to [standard data keywords][]

### domains

Annotation of RNA domains

example uses:

- same as [motif][]
- plus: labelling domains across the x-axis of skyline, profile, and arc plots

input explaination:

- input is a dictionary containing
  - `"domains"`: a list of lists of 2 integers. Each inner list specifies a
    start and end position of a primer (1-indexed, inclusive).
  - `"sequence"`: same as [sequence][] keyword
  - `"colors"`: a valid color or hexcode, e.g. `"blue"`, `"grey"`, or `"#fa4ce2"`
  - `"names"`: an arbitrary name to use on plots

example inputs:

```python
my_sample_1 = rnav.Sample(
    sample="example1",
    mrna_domains={
        "domains": [[1, 62], [63, 205], [206,300]],
        "sequence": "my_rna.fa",
        "colors": ["purple", "green", "orange"],
        "names": ["5'UTR", "CDS", "3'UTR"],
        }
    )
```

back to [standard data keywords][]

---

## Other optional arguments

Besides data keyword arguments, there are two other optional arguments to know.
Here is the call signature for `rnav.Sample`:

```python
rnav.Sample(
    sample,
    inherit=None,
    keep_inherited_defaults=False,
    **data_keywords)
```

### inherit

This value defaults to `None`. If another `rnav.Sample` is provide, the created
sample will inherit all of the data keywords from the provided sample.

These data are shared: any change to one sample applies to the other. This
sharing saves on memory and computation time for large data sets.

Example:

```python
shared_data = rnav.Sample(
    name='shared data',
    keyword1='big_structure.pdb')

sample1 = rnav.Sample(
    name='knockout sample',
    inherit=shared_data,
    keyword2='sample1-data.txt')

sample2 = rnav.Sample(
    name='knockout sample',
    inherit=shared_data,
    keyword2='sample2-data.txt')
```

`sample1` and `sample2` now share the data stored in `keyword1`. Another way to
achieve the same result is to pass the data object directly to a data keyword,
but this can get cumbersome if there is a lot of shared data.

```python
sample1 = rnav.Sample(
    name='knockout sample',
    inherit=shared_data,
    keyword1='big_structure.pdb',
    keyword2='sample1-data.txt')

sample2 = rnav.Sample(
    name='knockout sample',
    inherit=shared_data,
    keyword1=sample1.get_data("keyword1"),
    keyword2='sample2-data.txt')
```

### keep_inherited_defaults

This value defaults to `False`. If set to `True`, default keywords from
`inherit` will be kept. Each sample has default keywords for each of
the The provided value should be a boolean (`True` or `False`)

[customizing profiles][]
[customizing interactions][]
[filtering guide][]

[standard data keywords]: #standard-data-keywords
[sequence]: #sequence
[sites]: #sites
[spans]: #spans
[group]: #group
[primers]: #primers
[motif]: #motif
[orfs]: #orfs
[domains]: #domains
[profile]: #profile
[shapemap]: #shapemap
[dmsmap]: #dmsmap
[dancemap]: #dancemap
[rnpmap]: #rnpmap
[interactions]: #interactions
[ringmap]: #ringmap
[pairmap]: #pairmap
[shapejump]: #shapejump
[pairprob]: #pairprob
[allpossible]: #allpossible
[ss]: #ss
[pdb]: #pdb
[RNPMapper software]: https://github.com/Weeks-UNC/RNP-MaP
[ShapeMapper2 software]: https://github.com/Weeks-UNC/shapemapper2
[DanceMapper software]: https://github.com/MustoeLab/DanceMapper
[RingMapper software]: https://github.com/Weeks-UNC/RingMapper
[PairMapper software]: https://github.com/Weeks-UNC/RingMapper
[RNAStructure software]: https://rna.urmc.rochester.edu/RNAstructure.html
[StructureEditor]: https://rna.urmc.rochester.edu/RNAstructure.html
[XRNA]: http://rna.ucsc.edu/rnacenter/xrna/xrna.html
[VARNA]: https://varna.lisn.upsaclay.fr/
[FORNA]: http://rna.tbi.univie.ac.at/forna/
[R2DT]: https://rnacentral.org/r2dt
