Getting started
===============

This guide is meant to help users get started, whether they are seasoned python
pros or if they have never used python.

If you have not already installed RNAvigate, visit the
[installation guide](installation.md).

Open a Jupyter Notebook
-----------------------

The easiest way is to install Anaconda, and start a notebook using Anaconda
Navigator.

Import RNAvigate
----------------

```python
import rnavigate as MaP
```

Load your data files
--------------------

Data loading is simple. All values are optional. The MaP.Sample object is just
a group of data that provides convenient ways to parse, access, manipulate,
compare, align, filter, and plot different datatypes together.

Anatomy of sample creation:

```python
# Comment: the commas, quotes, and brackets here are strictly enforced.
sample_name = MaP.Sample(
    sample="Short Arbitrary Name",
    datatype_1="path/to/your/data.txt",
    datatype_2="path/to/your/other_data.txt",
    annotation_type={"list": [],
                     "seq_source": "datatype_1"}
)
```

`sample_name` should be a short code that you will use during plotting and
analysis to refer to this sample.

`datatype_1`, `datatype_2`, ..., `datatype_n` should be replaced with any of:

- `ct`, `compct`, or `ss`: a structure file - DB, CT, NSD, CTE, XRNA, or VARNA
- `shapemap`: Shapemapper2 output _profile.txt (recommended) or .map
- `dmsmap`: same file as `shapemap`, but RNAvigate will renormalize the profile
- `log`: Shapemapper2 output _log.txt (only if you used --per-read-histograms)
- `dance_prefix`: prefix to DANCE-MaP files
  - creates a list of new sample objects in self.dance containing
    - `ringmap`, `pairmap`, `profile`, `ct`, `pairprob` if the file is found.
- `rnpmap`: RNP-MaPper.py output .csv file
- `ringmap`: RingMapper output file (file extension provided by user)
- `pairmap`: PairMapper output _pairmap.txt file
- `allcorrs`: PairMapper output _allcorrs.txt file
- `pairprob`: a dotplot file (.bp) containing base-pairing probabilities
  - many RNA structure prediction packages can output this file

`annotation_type` should be replaced with:

- `sites={"site_list": [list], "seq_source": "datatype"}` 
  - `[list]` is a comma-seperated list of sites of interest in square brackets
- `spans={"span_list": [[span_1], ..., [span_n]], "seq_source": "datatype"}`
  - `[span_n]` is a comma-seperated start, stop pair of spans of interest
  - e.g. `[[1, 10], [50, 80]]` would be the regions 1-10 and 50-80.
- `groups={"group_list": [[group_1], ..., [group_n]], "seq_source": "datatype"}`
  - `[group_n]` is a comma-seperated list of related nucleotides
  - e.g. `[[1, 5, 10], [50, 65, 72, 80]]`
- `orfs={"seq_source": "datatype"}`
  - `orfs` behaves like spans, but it generates the span list automatically by
    finding ORFs in the sequence from `datatype`.
- `motif={"motif": "DRACH", "seq_source": "datatype"}`
  - `motif` also behaves like `spans`, but generates the span list by search
    for the motif in the sequence from `datatype`.

`"datatype"` in all of these cases should be one of the datatypes listed above
enclosed in double-quotes. e.g. `orfs={"seq_source": "shapemap"}` will make an
annotation of the ORFs in the sequence contained in the file passed to
`shapemap="path/to/my_sample_profile.txt"`.

Explore your data
-----------------

Anatomy of a simple plot call:

```python
# Single sample
sample_name.plot_plottype(
    data="datatype",
    data_filter={"filter": threshold},
    plot_keyword=value)
```

- `sample_name` is the short code you generated in the sample creation step.
- `plottype` can be one of:
  - `qc`: a shapemapper quality control plot (requires `shapemap` and `log`)
  - `shapemapper`: the classic ShapeMapper profile plot
  - `skyline`: a skyline plot of per-nucleotide data, also known as a step plot
  - `circle`: a circle diagram where nucleotides are arranged in a circle
  - `arcs`: an arc plot diagram where the RNA is the x-axis, with arcs
    representing base-pairs or other interactions
  - `ss`: a secondary structure diagram
  - `mol`: a 3d molecular rendering
  - `heatmap`: both axes are nucleotide positions to plot (x,y) interactions
  - `disthist`: a histogram of the 3D or 2D distances of interactions
- `data`: This differs for each plot type. See [Plots](plots.md) for more info.
- `"datatype"` is one of the datatypes provided to the sample
- `"data_filter"`: This is very flexible. See [Filters](filters.md) for more info.

```python
# Single sample, one plot per filter
sample_name.plot_plottype_multifilter(
    filters=[{"data":dataname,
              "filter_1": threshold},
             {"data": dataname,
              "filter_2": threshold}]
    plot_keyword=value
)
```

```python
# Multiple samples
MaP.plot_plottype_multisample(
    samples=[sample_name_1, sample_name2],
    data=dataname,
    data_filter={"filter_1": threshold,
                 "filter_1": threshold},
    plot_keyword=value)
```
