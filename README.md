RNAvigate (RNA visualization and graphical analysis toolset)
==============================================================================
RNAvigate provides a framework for rapidly exploring interrelationships between
chemical probing data, RNA structures, and motif annotations for multiple
experimental samples within a Jupyter Notebook environment.
* [Full documentation](https://rnavigate.readthedocs.io/en/latest/)
* [RNA 2022 poster](https://rnavigate.readthedocs.io/en/latest/rna2022.html)

Data types
----------
* ShapeMapper2: profiles.txt, shapemapper_log.txt
* RingMapper: rings.txt
* PairMapper: pairmap.txt, allcorrs.txt
* DanceMapper: profiles, rings, pairs, allcorrs, structure predictions
* ShapeJumper: deletions.txt (requires .fasta)
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
* 3D distance distributions
* Annotations on arc and secondary structure
* Forthcoming:
  * Annotations on skyline, circle, and 3D structure
  * information funnel diagram

Analyses
--------
* delta log SHAPE
* low SHAPE, low Shannon entropy

Other features
--------------
* Filtering data on many parameters
* Calculating 3D distance and contact distances
* Automatic data alignment for comparing sequence variants

General notes
-------------
* PDB format can be picky, and there are mistakes in some files from the PDB.
* When creating many plots, it may boost performance to run `plt.close('all')`.
