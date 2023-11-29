Welcome to RNAvigate docs!
==========================

![demo plots](./images/gallery.png)

[Github](https://github.com/Weeks-UNC/RNAvigate) | [Publications](publications.md)

- [What is it and who is it for?](#what-is-it-and-who-is-it-for)
- [What problem does it solve?](#what-problem-does-it-solve)
- [Workflow overview](#workflow-overview)
- [Questions, requests, and issues](#questions-requests-and-issues)

Getting started

1. [Installing RNAvigate](installing-rnavigate.md)
2. [Loading data files](loading-data.md)
3. [Visualizing data](visualizing_data.ipynb)


## What is it and who is it for?

RNAvigate is an RNA structure and chemical probing data exploration toolset
built for RNA bench scientists. It is designed to be simple to learn, to accept
many input data formats, and to provide useful tools, common analyses, and
visualizations for the RNA structure community.

## What problem does it solve?

Many tools and scripts exist to analyze RNA structure data via the command line
or through GUI interfaces. These types of tools are generally not well suited
for rapid data exploration because of the additional overhead of managing and
organizing intermediate data files, or of clicking through to reproduce an
analysis.

Interactive programming with Jupyter notebooks offers an elegant solution to
this problem and most data scientists utilize these notebooks for
reproducible data exploration, analysis, and figure creation. However, these
notebooks require fluency in data science programming. RNAvigate solves this by
providing a simple to understand interface that implements many of the common
analyses and visualizations utilized by the RNA structure community.

## Workflow overview

1) Search databases or perform experiments and data pre-processing to
   obtain data file inputs for RNAvigate.

2) Create an `rnavigate.Sample` and load data by assigning file names to
   **data keywords**.

```python
import rnavigate as rnav

sample_name = rnav.Sample(
    data_keyword="input_file.txt",
)
```

3) Create `rnavigate.plots` by providing `rnavigate.Sample` names and
   **data keywords** to plotting functions.

```python
rnav.plotting_function(
    samples=[sample_name],
    plotting_keyword="data_keyword",
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
