Welcome to RNAvigate docs!
==========================

[Github](https://github.com/Weeks-UNC/RNAvigate)

RNA visualization and graphical analysis toolset (RNAvigate) is a plotting
package for diverse RNA structure data in standard formats. It provides
an easy text interface for quickly producing visualizations for data analysis
and exploration within Jupyter Notebooks or python scripts. RNAvigate is
designed to be simple to install and use with minimal knowledge of python or
command line tools.

Check out our [publication](publications.md).

Getting started
---------------

To start exploring your own RNA data sets, follow these three easy steps:
1. [Installing RNAvigate](installing-rnavigate.md)
2. [Loading data files](loading-data.md)
3. [Visualizing data](visualizing-data.md)

Brief overview
--------------

Once you know the data class keywords and plotting functions available,
RNAvigate is very easy to use:
1. open a Jupyter notebook and import RNAvigate
```python
import rnavigate as rnav
```
2. Create samples by providing data file inputs with data class keywords
```python
sample_name = rnav.Sample(
    data_class_keyword="input_file.txt"
)
```
3. Create visualizations by providing plotting functions with sample names and
   data class keywords.
```python
rnav.plotting_function(
    samples=[sample_name]
    plotting_keyword="data_class_keyword"
)
```

Feature requests, bug reporting, and collaborating
--------------------------------------------------

RNAvigate is new and under active development. Please, use
[GitHub issues](https://github.com/Weeks-UNC/RNAvigate/issues) to request new
features or report bugs. There's a chance that your feature already exists, but
is not well documented. In this case, I will reply with an explaination of how
to implement it, and I will prioritize adding this to the documentation for
other users to find.

Reach out to me on [Github](https://github.com/Psirving) if you are interested
in helping improve RNAvigate or implementing new data types, plots, or analyses.
