Version history
===============

1.0.0-alpha (November 2023)

This version is a large departure from v0.2.0, while trying to maintain a
familiar interface.

If you're coming from version 0.2.0, there are many syntax changes. In order to
preserve the functionality of your old notebooks, do not pull these changes to
your existing RNAvigate folder. Instead, follow the [installation][] guide.

[installation]: ../installing-rnavigate.md

NOTES:

- Argument names are fewer, clearer, and more standardized:

| new names                     | old names                                             |
| :---------------------------: | ----------------------------------------------------- |
| structure or structure2       | ss, ct, comp                                          |
| sequence                      | fit_to, dataframe, filepath, seq_source, fasta        |
| interactions or interactions2 | ij, ij2                                               |
| colors                        | nt_color, color, colormap, etc.                       |
| input_data                    | dataframe, filepath                                   |
| metric_defaults               | column, err_column, color_column, cmap, norm, min_max |

- SequenceChecker analysis adds tools to quickly compare sequences.
- Multiple sequence alignment is partially implemented.
- All arguments are explicitly named for high level functions (no more kwargs)
- Removed these standard data keywords:
  - log: log files can now be passed during ShapeMaP profile creation
  - dmsmap: shapemap keyword accepts normalization="DMS", "eDMS", "boxplot"
  - allcorrs: this was just a duplicate for ringmap
  - groups: now called "group" and annotates 1 group at a time
- data keywords are more flexible:
  - can be arbitrary strings
  - optional arguments can be passed to data class constructor
- Loading a genome fasta and transcriptome annotation (GTF) allows:
  - extracting Transcript sequences given a transcript ID
  - annotating exon junctions and coding sequences
  - extracting transcript profile or annotation from bed and NarrowPeak files
- calculations on profiles:
  - Normalizations: (uses code from RNATools ReactivityProfile)
    - "DMS", "eDMS", "boxplot", "percentile"
    - options for setting bounds and flexible grouping by nucleotide
  - Rolling windows:
    - mean, median, average, custom function
    - setting window size
- SecondaryStructure class (new name) combines all file formats
  - Added parsing for FORNA and R2DT structure diagrams in json format
    - `"extension":"forna"` and `"extension":"r2dt"` differentiate these files
- Profile plots can now have horizontal annotation tracks
- Sequence alignment changes:
  - data objects have get_aligned_data which returns a copy with new positions
  - option to use secondary structure alignments (uses RNAlign2D algorithm)
  - more understandable sequence alignment figures.
  - flipping an alignment (seq1->seq2 becomes seq2->seq1)
- Interactions and Profile have a unified and more flexible coloring interface
  - This also controls colorbar appearance
- all plots changes:
  - Scaling text size on a plot using sns.context (doesn't always look nice)
  - undercase nucleotides in sequences are preserved
  - figure size in inches is more reasonable for exporting images and svgs
  - figure size scaling is standardized and correctly calculated for all plots
  - nt_ticks parameter defines major and minor tick marks
    - always includes the first nucleotide
    - skips over indels in the sequence
- secondary structure diagrams changes:
  - position labels look much nicer
  - axis margins are standardized to 2 data units for all size RNAs
- linear regression plots changes:
  - removed KDE plots
  - simplified and improved readability
  - scale options (linear or log) and regression options (Pearson or Spearman)
- Added plot_ntdist, which plots reactivity distributions by nucleotide.
- The following analyses have been refactored:
  - DeltaSHAPE, fragmapper, lowSS
  - (not stable) AUROC, LogDiff


---

0.2.0 (September 2023)
----------------------

NOTES:

- Got rid of `sample.plot_function()` because it is redundant: use `rnav.plot_function([sample])` instead
- Position labelling on secondary structure is prettier
- New features for linear regression plots:
  - pearson or spear correlations
  - any column of data
  - log or linear regression
- Fragmapper analysis

---

0.1.0 (April 2023)
---------------------
NOTES:

- Beginning of version history
- removed py3dmol version from env.yaml (uses most recent)

New plots:

- ROC plots
- DeltaSHAPE
- Windowed AUROC
- Alignment
- Profiles

New data classes:

- AllPossible: creates Interactions data given a sequence. One interaction for
  every possible nucleotide pairing.

Features:

- Many small bug fixes and aesthetic changes
- PDB data objects:
  - Added support for .cif files
  - improved parsing and cross-indexing, offset argument removed
- Improved data fitting flexibility
  - Added fit_to() method for secondary structures, profiles, and annotations
    objects
  - Added seq_source for arc and circle plots
  - alignment maps for a data object can be predefined
- Added features to retrieve/set orientation of 3D molecule plots
- Secondary structure data and plotting:
  - Added xy-coordinate normalization:
    - median base-pair distance = 1, center of structure = (0, 0)
  - Changed SS interface to allow multiple different structures in 1 figure
- data objects can be passed to Sample() arguments
  - this allows rnav.Sample objects to share a data object
  - cuts down on computation time and memory usage
  - using inherit argument, all data objects from a sample are inherited
- Add set_figure_size to Plot class
  - sets figure size so that axis unit to inches ratio is consistent
