Welcome to RNAvigate docs!
==========================

**RNA vi**sualization and **g**raphical **a**nalysis **t**ools**e**t

[Github](https://github.com/Weeks-UNC/RNAvigate) | [Publication](publications.md)

RNAvigate is a plotting and analysis framework for diverse RNA structure data.
It provides an easy interface for quickly producing visualizations for
data analysis and exploration within Jupyter Notebooks or python scripts.
RNAvigate is designed to be simple to install and use with minimal knowledge of
python or command line tools.

Getting started
---------------

1. [Installing RNAvigate](installing-rnavigate.md)
2. [Loading data files](loading-data.md)
3. [Visualizing data](visualizing-data.md)

Very brief overview
-------------------

1. open a Jupyter notebook and import RNAvigate

```python
import rnavigate as rnav
```

2. Create samples by providing data file inputs with data keywords

```python
sample_name = rnav.Sample(
    data_keyword="input_file.txt"
)
```

3. Create visualizations by providing plotting functions with sample names and
   data keywords.

```python
rnav.plotting_function(
    samples=[sample_name]
    plotting_keyword="data_keyword"
)
```

Bugs, requests, and questions
-----------------------------

Use [GitHub issues](https://github.com/Weeks-UNC/RNAvigate/issues) to request
new features, to report bugs, and to ask questions. This feedback is extremely
useful for improving RNAvigate, so please don't be shy!

In most cases, RNAvigate can be very efficiently expanded to accept new data
file formats. To request a new format, submit a GitHub issue with an example
file, file format specification, and example visualizations if they exist.
For more information on the broad categories of data that RNAvigate is
well-suited for see [data types](data-types.md).

Developers
----------

Please [contact me](mailto:psirving@email.unc.edu) if you are interested in
helping to improve RNAvigate or in using RNAvigate in your own projects.
