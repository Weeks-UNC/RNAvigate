Welcome to RNAvigate docs!
==========================

**RNA VI**sualization and **G**raphical **A**nalysis **T**ools**E**t

[Github](https://github.com/Weeks-UNC/RNAvigate) | [Publication](publications.md)

RNAvigate is a plotting and analysis framework for diverse RNA structure data.
It provides an easy-to-use interface for exploratory data analysis within
Jupyter Notebooks or python scripts. Knowledge of python or command line tools
is not necessary.

Create a **Sample** by assigning file names to **data keywords**.

```python
import rnavigate as rnav

sample_name = rnav.Sample(
    data_keyword="input_file.txt",
)
```

Create plots by providing **Sample** names and **data keywords** to plotting functions

```python
rnav.plotting_function(
    samples=[sample_name],
    plotting_keyword="data_keyword",
)
```

Getting started guides
----------------------

1. [Installing RNAvigate](installing-rnavigate.md)
2. [Loading data files](loading-data.md)
3. [Visualizing data](visualizing-data.md)


Bugs, requests, and questions
-----------------------------

Use [GitHub issues](https://github.com/Weeks-UNC/RNAvigate/issues) to request
new features, to report bugs, or to ask questions. This feedback is extremely
useful for improving RNAvigate, so please don't be shy!

In most cases, RNAvigate can be very efficiently expanded to accept new data
file formats. To request a new format, submit a GitHub issue with an example
file, additionally a file format specification and/or example visualizations if
they exist. For more information on the broad categories of data that RNAvigate
is well-suited for see [data types](data-types.md).

Developers
----------

[Full API](full_api_index.md)

Please [contact me](mailto:psirving@email.unc.edu) if you are interested in
helping to improve RNAvigate or in using RNAvigate in your own projects.
