Jupyter Notebook Plotting Tools
==============================================================================
Python plotting functions for making pretty figures for common lab experiments.
If you'd like me to add any functionality to this, let me know by openning up an
issue in the issues tab, or email me: psirving@email.unc.edu

plotmapper.py
-------------
* Click here to see [features, plan, and todo list](todo.md)
* High level plotting functions
* Simple API for loading data and plotting

Examples:
* [single sample plots](JNB-example/plotmapper-example.md)
* [multiple samples plots](JNB-example/plotmapper-multiple-examples.md)
* [Arc Plot options](JNB-example/ap_test.md)
* [Secondary Structure options](JNB-example/ss_test.md)
* [3D plot options](JNB-example/3d_test.md)
* [colorbar options](JNB-example/colors_test.md)
* [heatmap options](JNB-example/heatmap_test.md)

Issues:
* PDB format can be picky, and there are mistakes in some files from the PDB.
* When creating many plots, it may boost performance to run `plt.close('all')`.
