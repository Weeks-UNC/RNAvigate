Jupyter Notebook Plotting Tools
==============================================================================
Python plotting functions for making pretty figures for common lab experiments.
If you'd like me to add any functions to this, let me know by openning up an
issue in the issues tab, or email me: psirving@email.unc.edu

For more in-depth exploration of the functions and example plots, see the
[wiki](https://github.com/Weeks-UNC/JNBTools/wiki).

TODO list:
* Create function to set ylims for arc plots.
* Create a class structure for arc plots.
* Add handling of deletion files to secondary structure.
* Make default secondary structure plots prettier.
* Add example usage notebook for secondary structures.

plottingTools.py [examples](JNB-example/plottingTools-example.md)
* Make Skyline plots for comparing nt-by-nt data.
* Add a sequence bar to your plots.
* Fetch histogram data from shapemapper log files.
* Set figure width based on sequence length.

shapemapperPlots.py [examples](JNB-example/plottingTools-example.md)
* Make standard shapemapper output plots:
  * Normalized profile
  * Raw mutation rates
  * Read depth and effective read depth

JNBarcPlot.py [examples](JNB-example/JNBarcPlot-example.md)
* Make plots of correlation data:
  * arcPlot figures
  * dotplot figures
* Get pairmap sensitivity and PPV.

secondaryStructure.py [examples](JNB-example/secondaryStructure-example.md)
* Make secondary structure graphs from XRNA or Structure Editor
* Add Rings, filtered by statistic or contact distance.
