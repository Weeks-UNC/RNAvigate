Welcome to RNAvigate docs!
==========================

[Github](https://github.com/Weeks-UNC/RNAvigate) | [Publications](publications.md)

- [What is it and who is it for?](#what-is-it-and-who-is-it-for)
- [What problem does it solve?](#what-problem-does-it-solve)
- [Workflow overview](#workflow-overview)
- [Questions, requests, and issues](#questions-requests-and-issues)

If you are already interested in getting started, skip ahead to these guides:

1. [Installing RNAvigate](installing-rnavigate.md)
2. [Loading data files](loading-data.md)
3. [Visualizing data](visualizing_data.ipynb)

## What is it and who is it for?

RNAvigate is a toolset for exploring RNA structure and chemical probing data.
It is built for RNA bench scientists and with the following goals:
1. Define analyses in code so they are easier to iterate, reproduce, and share
2. Be easy to set up and to learn without coding experience
3. Improve workflow and reduce time-consuming tasks
4. Handle a wide variety of input data formats and useful analyses
5. Be extensible to new data formats and analyses

## What problem does it solve?

Many tools and scripts exist to analyze RNA structure data via command line interfaces
(CLIs) or graphical-user interfaces (GUIs). CLIs and GUIs each have trade-offs and are
well suited to particular tasks. However, generally, neither is well suited for
exploratory data analysis.

Interactive programming interfaces, such as Jupyter, provide a better workflow.
Jupyter notebooks are a gold standard for reproducible data exploration, analysis,
figure creation and reporting. However, they require programming fluency.

RNAvigate solves this by providing simple, intuitive commands that handle the minutiae
of parsing and manipulating RNA-centric data and of implementing many of the common
analyses and visualizations utilized by the RNA structure community.

## Workflow overview

RNAvigate separates the performance of three main tasks:

### Curate inputs

Search databases or perform experiments and data pre-processing to
obtain data files. These are the inputs for RNAvigate.

### Load and organize data

The first step in using RNAvigate is to create a **Sample**. Samples contain all of the
data related to a single experimental condition. These data can be:
- sequence annotations
- secondary structures, optionally with drawing coordinates
- tertiary structures with 3D atomic coordinates
- per-nucleotide measurements
- internucleotide measurements

Samples are created by providing `rnav.Sample()` with keyword-argument pairs. Keywords
specify file format, and the arguments are file paths. In this simple example, two
experiments, each with associated SHAPE-MaP reactivities and predicted structures, are
loaded as seperate samples.

```python
import rnavigate as rnav

my_sample = rnav.Sample(
    shapemap="input_file.txt",
    ss="secondary_structure.ct"
    )

my_2nd_sample = rnav.Sample(
    shapemap="input_file_2.txt",
    ss="secondary_structure_2.ct"
    )
```

### Analyze and visualize data

The second step in RNAvigate is to explore the data using the built in analysis and
visualization functions. These functions accept sample names (e.g. `my_sample`) and
data keywords (e.g. `"shapemap"`) to define which data to inspect.

For example, to compare our two data sets using arc plots:

```python
rnav.plot_arcs(
    samples=[my_sample, my_2nd_sample],
    sequence="ss",
    structure="ss",
    profile="shapemap",
    )
```

## Questions, requests, and issues

Use [GitHub issues](https://github.com/Weeks-UNC/RNAvigate/issues), liberally,
to request new features, to report bugs or unexpected behavior, or simply to
ask questions. The future development directions of RNAvigate will be based
largely on this feedback.

In almost all cases, RNAvigate can be easily expanded to accept new data file
formats. Use the link above to request support for a new format or a new visual
language (like the red-orange-black-grey scale for structure probing
reactivity), include an example file and/or example visualization. For more
information on the broad categories of data that are well-suited for inclusion,
see [data types](data-types.md).

## Developers

[Contact me](mailto:psirving@email.unc.edu) if you are interested in
helping to improve RNAvigate or in using RNAvigate in your own projects.

These developer resources are works-in-progress.

- [Change log](dev/changelog.md)
- [Full API](dev/full_api.md)
- [Developer's style guide](dev/style_guide.md)
