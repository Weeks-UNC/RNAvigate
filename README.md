StarMapper.py - Secondary and Tertiary Structural Annotation of RNA
==============================================================================
StarMapper provides a framework for rapidly exploring interrelationships between
chemical probing data, RNA structures, and motif annotations for multiple
experimental samples within a Jupyter Notebook environment.
* Click here to see [features, plan, and todo list](todo.md)

Data types
----------
* ShapeMapper2: profiles.txt, shapemapper_log.txt
* RingMapper: rings.txt
* PairMapper: pairmap.txt, allcorrs.txt
* DanceMapper: profiles, rings, pairs, allcorrs, structure predictions
* ShapeJumper: deletions.txt, .fasta
* Base-pairing: .ct, .db
* Secondary Structure diagrams: .xrna, .varna, .nsd, .cte
* 3-D molecular structures: .pdb
* RNP-MaP
* Forthcoming:
  * annotations
  * .cif and .pdbx formats

Plot types
----------
* skyline profile plots
* profile linear regression plots
* arc plots
* circle plots
* secondary structure diagrams
* interactive 3-D molecules
* heatmaps with structure contour maps
* ShapeMapper profile plots
* ShapeMapper quality control plots
* Forthcoming:
  * 3D distance distributions
  * Annotations on skyline, arc, circle, secondary structure, and 3D structure
  * information funnel diagram

Analyses
--------
* Log Compare à la Steve

Other features
--------------
* Filtering data on many parameters
* Calculating 3D distance and contact distances
* Automatic data alignment for comparing sequence variants

Example Jupyter Notebooks
-------------------------
* [single sample plots](JNB-example/starmapper-example.md)
* [multiple samples plots](JNB-example/starmapper-multiple-examples.md)
* [Arc Plot options](JNB-example/ap_test.md)
* [Secondary Structure options](JNB-example/ss_test.md)
* [3D plot options](JNB-example/3d_test.md)
* [heatmap options](JNB-example/heatmap_test.md)

General notes
-------------
* PDB format can be picky, and there are mistakes in some files from the PDB.
* When creating many plots, it may boost performance to run `plt.close('all')`.
