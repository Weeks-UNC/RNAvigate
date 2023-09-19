RNAvigate (RNA visualization and graphical analysis toolset)
==============================================================================
RNAvigate provides a framework to explore and compare chemical probing data,
RNA structures, and motif annotations between experimental samples.
It is useful for scripting, but is most powerful in a Jupyter Notebook.

RNAvigate uses Python, but it is designed to be very easy to learn for
non-programmers.

For Python developers, the underlying object-oriented interface is powerful and
flexible.

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
* 3-D molecular structures: .pdb or .cif
* RNP-MaP
* annotations

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

Analyses
--------
* DeltaSHAPE
* Frag-MaP
* logSHAPE comparison
* low SHAPE, low Shannon entropy
* Windowed area under ROC curve

Other features
--------------
* Filtering data on many parameters
* Calculating 3D distance and contact distances
* Automatic data alignment for comparing sequence variants and subsequences

General notes
-------------
* PDB format can be very picky.
* When creating many plots, it may boost performance to run `plt.close('all')`.
