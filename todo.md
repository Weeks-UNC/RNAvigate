Plan for RNAvigate
==================
Plot-MaP should be an easy to use, universal tool for plotting Mutational
Profiling (MaP) and Juxtaposed Merged Pairs (JuMP) data and quality control
information. It should be easily extensible to accomodate new plots, data, and
analysis.

Coding To-Do List
-----------------

- hide code button not working?
- refactoring:
  - ct, compct, ss to basepairs or ss or something
  - ct, ij, ij2 to structures, interactions
- better figure scaling using subplotpars -> set_size_inches
- reasonable sizing for all plots
- find best way to display colormaps for all plots (probably a 2nd fig)
- enable multiple structure inputs for SS plotting
  - scale data points equally for ct-like objects
  - center on 0, 0
  - share x and y axes
  - scale figure size using subplotpars
- Expose arguments in Mol
  - turn off nt cylinders
  - width, height, background color
- implement profiles for interaction data
  - use get_data function to return profiles when expected
- loading secondary structure files containing more than one structure
- passing override values to init_dance
- calling significant sites with log-corrected profile min-diff comparison
- better nucleotide colors and display on linear regression plots
- loading in annotations files
- send to file for all plots
- improve labelling for disthist plots
- reading in FORNA JSON


Documentation To-Do List
------------------------

- Reference page for each plot type
  - circleplot, skyline, qc, heatmap, linreg, sm, disthit
- guide for custom use cases
  - loading custom profiles, interactions and annotations
  - plot manipulation with mpl
  - data manipulation with pandas
- good doc strings

Analyses
--------

- RNP-MaP
- log(+/-) - k*log(+/-) Normalization
- deltaSHAPE

Command-Line interface
----------------------

Currently, there is no CLI interface, but this is how I imagine it working. I
have to learn how to parse args dynamically, given the plotting function.

```
rnavigate plottype --sample1 --datatype filepath --datatype filepath \
                   --sample2 --datatype filepath --datatype filepath \
                   --filter_arguments \
                   --plot_arguments
```
