# TODO: RNAvigate

> Development backlog for turning RNAvigate into a well-structured, tested, and
> published conda/pip package. Items are roughly ordered within each section by
> priority. Check off items here AND update `ARCHITECTURE.md` when structural
> changes are made.

---

## In Progress

---

## Bugs

- [ ] log file specification is not working for ShapeMapper v2.2.0
  <!-- The read_log() method in SHAPEMaP (data/profile.py) parses the v2
       histogram format. The ShapeMapper 2.2 log format appears to have changed.
       This needs a real v2.2 log file to inspect and update the parser. -->
- [ ] annotation colorbar not added with `plot_ss`
  <!-- plot_ss calls AP.plot_data() which calls add_colorbar_args(), but the
       Annotation colormap is not being passed through. Likely a missing
       add_colorbar_args() call in plots/ss.py plot_data(). -->

---

## Infrastructure & Packaging

- [ ] Test and improve pytest code coverage.
- [ ] Automate release pipeline
  - version bumping (bump2version or hatch version)
- [ ] Write `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`
  <!-- Blocks accepting external contributions. Should describe: how to set up
       dev environment, code style, branching model, PR process. -->

---

## Refactoring
<!-- improve long-term maintainability. -->

- [ ] Fix circular imports between `plots/__init__.py` and subclasses
- [ ] Fix circular import between `data/data.py` ↔ `data/alignments.py`
- [ ] Split `plotting_functions.py` (1,624 lines) into per-plot modules
  <!-- Each plot_*() function is largely independent. Moving them into
       plots/functions/ or a plots/api/ submodule would make the file
       navigable and allow lazy imports. -->
- [ ] Standardize the `analysis/` module pattern
  <!-- Currently mixes: Sample subclasses (DeltaSHAPE, LogCompare, LowSS,
       Fragmapper, FragmapperReplicates), Profile subclasses (FragMaP,
       DeltaSHAPEProfile, LogProfile), and standalones (WindowedAUROC,
       SequenceChecker). The inconsistency makes the module hard to navigate.
       Recommendation: all analysis entry points should be Sample subclasses;
       data products (FragMaP, etc.) should live in data/. -->
- [ ] Refactor `WindowedAUROC` as a `Sample` subclass
  <!-- There is a TODO comment in analysis/auroc.py noting this. It would make
       AUROC results directly plottable like DeltaSHAPE and LowSS. -->
- [ ] Replace flat `data_keyword_defaults` dict with a registration system
  <!-- Currently adding a new data type requires editing data_loading.py.
       A decorator or register() function would allow data types to self-register
       and would support user-defined or plugin data types. -->
- [ ] Refactor analyses (ongoing)
  - [ ] AUROC
  - [ ] log-diff

---

## New Features

### Data & I/O

- [ ] Tracks that show continuous data (not just per-nucleotide bars/skylines)
- [ ] MSA (Multiple Sequence Alignment) tracks
- [ ] Annotations: support input/output file formats
  - columns: start, end, color, label, track, track height, type
  - make annotations more flexible generally
- [ ] mRNA annotation using width, showing splice sites
- [ ] Region support with sequence alignments
- [ ] Profiles for interaction data
  - pairing probability → per-nucleotide probability or Shannon entropy
  - RING-MaP → RING density per nucleotide
  - store as attribute; retrieve at point of use
    - e.g. `profile.profile` returns profile
- [ ] Codon usage bias for ORF annotations
- [ ] Add printable human-readable identifier for all data objects
  - proposed format: `sample + name + datatype + filepath`
- [ ] Implement SecondaryStructure.get_structure_elements()

### Analysis

- [ ] Calling significant sites with log-corrected profile min-diff comparison
- [ ] Structure compatibility analysis

### Visualization

- [ ] Consistent alpha values for colormaps and colorbars
  <!-- The alpha set on the ScalarMappable in colors.py is not reliably propagated
       to the colorbar artist. Needs investigation in plots/plots.py plot_colorbars(). -->
- [ ] Sequence Analyzer for individual `plot_*` functions
  <!-- Originally listed in triage: a per-plot diagnostic showing which data
       keywords and sequences are present/compatible before plotting. -->
- [ ] `plot_ss` with annotation colorbar (see Bugs above — fix the bug first)

### Interface & Usability

- [ ] Feature to turn off automatic sequence alignments (keep sub-sequences as-is)
- [ ] Validate alignments passed to `get_aligned_data()` with informative errors
- [ ] Simplify `get_aligned_data()` and `AlignmentChain()` interfaces
  <!-- The current API requires understanding the internals of the alignment chain.
       A higher-level convenience method on Sequence or Sample would help. -->

---

## Documentation

- [ ] Numpy-style docstrings across all public classes and functions
  <!-- Coverage is good in data/ but sparse in plots/ and analysis/. -->
- [ ] Guides for custom use cases:
  - [ ] Loading custom annotations
  - [ ] Sequence alignments (manual and automatic)
  - [ ] Plot manipulation with matplotlib
  - [ ] Data manipulation with pandas
- [ ] Keep `ARCHITECTURE.md` current
  <!-- Update this file whenever: a module's responsibilities change, a class is
       added/moved/removed, a new dependency is introduced, or a design pattern
       changes. Treat it as the design contract between contributors. -->
- [ ] Fragmapper and Transcriptome notebooks need fixing.
