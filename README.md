Jupyter Notebook Plotting Tools
==============================================================================
Python plotting functions for making pretty figures for common lab experiments.
If you'd like me to add any functionality to this, let me know by openning up an
issue in the issues tab, or email me: psirving@email.unc.edu

plotmapper.py
-------------
* Click here to see [features, plan, and todo list](todo.md)
* Contains all functionality from the other scripts
* Higher level plotting functions
* Simplier API for loading data and plotting

Examples:
* [single sample plots](JNB-example/plotmapper-example.md)
* [multiple samples plots](JNB-example/plotmapper-multiple-examples.md)

### Old code
The old code that made up this repo is still being stored here, but has been
moved to the old-code directory.

#### plottingTools.py
[examples](old-code/plottingTools-example.md)
* Make Skyline plots for comparing nt-by-nt data.
* Add a sequence bar to your plots.
* Fetch histogram data from shapemapper log files.
* Set figure width based on sequence length.

#### shapemapperPlots.py
[examples](old-code/plottingTools-example.md)
* Make standard shapemapper output plots:
  * Normalized profile
  * Raw mutation rates
  * Read depth and effective read depth

#### JNBarcPlot.py
[examples](old-code/JNBarcPlot-example.md)
* Make plots of correlation data:
  * arcPlot figures
  * dotplot figures
* Get pairmap sensitivity and PPV.

#### secondaryStructure.py
[examples](old-code/secondaryStructure-example.md)
* Make secondary structure graphs from XRNA or Structure Editor
* Add Rings, filtered by statistic or contact distance.
