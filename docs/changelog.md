Version history
===============

---

0.1.0 (February 2023)
---------------------
NOTES:

- Beginning of version history

CHANGES:

- Many small bug fixes and aesthetic changes
- removed py3dmol version from env.yaml (uses most recent)
- Added ROC plots, DeltaSHAPE, Windowed AUROC
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
- Added AllPossible interactions object for computing known-truth data
- data objects can be passed to Sample() arguments
  - this allows rnav.Sample objects to share a data object
  - cuts down on computation time and memory usage
  - using inherit argument, all data objects from a sample are inherited
- Add set_figure_size to Plot class
  - sets figure size so that axis unit to inches ratio is consistent
  - not yet implemented in all plots
- Extensive updates to Readthedocs documentation website.

---
