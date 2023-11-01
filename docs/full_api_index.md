RNAvigate API
=============

This is not the full documentation, but an index and brief description of all
of the modules, classes, methods and functions available. For every entry here,
the help function is given. Import RNAvigate and run that command to get more
information, including the call signature.

Eventually, I will get the full API documentation on this website in a nicely
formatted way using Sphinx.


Index
-----

[main](#rnavigate-main-directory)

These three submodules, along with the analysis module, are the high-level API.

- [rnavigate](#rnavigate)
- [plotting_functions](#plotting_functions)
- [styles](#styles)

Accessory functions and argument parsing

- [data_loading](#data_loading)
- [helper_functions](#helper_functions)

**[data](#data)**

- [alignments](#alignments)
- [colors](#colors)
- [data](#data)
- [annotation](#annotation)
- [ct](#ct)
- [interactions](#interactions)
- [log](#log)
- [pdb](#pdb)
- [profile](#profile)

**[transcriptomics](#transcriptomics)**

- [transcriptome](#transcriptome)
- [bed](#bed)
- [eclip](#eclip)

**[plots](#plots)**

- [plots](#plots)
- [functions](#functions)
- [alignment](#alignment)
- [arc](#arc)
- [circle](#circle)
- [disthist](#disthist)
- [heatmap](#heatmap)
- [lingreg](#lingreg)
- [mol](#mol)
- [profile](#profile)
- [qc](#qc)
- [roc](#roc)
- [skyline](#skyline)
- [sm](#sm)
- [ss](#ss)

**[analysis](#analysis)**

- [fragmapper](#fragmapper)
- [auroc](#auroc)
- [deltashape](#deltashape)
- [logcompare](#logcompare)
- [lowss](#lowss)
- RNPMapper - coming soon

---

## rnavigate main directory

The main directory of RNAvigate contains the high-level API.

### rnavigate

This submodule contains the Sample class.

**Sample class**: `help(rnav.Sample)`

- The high-level API is centered around the sample class, which parses and
  organizes data related to a single experiment on a single RNA.

*print_data_keywords*: `help(rnav.Sample.print_data_keywords)`

- Print the data keywords associated with this sample, organized by data type.

*set_data*: `help(rnav.Sample.set_data)`

- Create a new data keyword using similar syntax as Sample initialization.

*get_data*: `help(rnav.Sample.get_data)`

- retreive data from data keywords in a flexible way with informative errors.

*filter_interactions*: `help(rnav.Sample.filter_interactions)`

- Applies filters and colorscheme to interactions data.

### plotting_functions

This submodule contains all of the high-level plotting API.

*plot_qc*: `help(rnav.plot_qc)`

- Plots quality control metrics from a ShapeMapper run

*plot_shapemapper*: `help(rnav.plot_shapemapper)`

- Plots the standard ShapeMapper 3-panel figure

*plot_skyline*: `help(rnav.plot_skyline)`

- Plots per-nucleotide data as a stepped line graph for easy comparison

*plot_profile*: `help(rnav.plot_profile)`,

- Plots per-nucleotide data as a colored bar graph

*plot_alignemnt*: `help(rnav.plot_alignemnt)`

- Plots the alignment used to compare datasets with differing sequences

*plot_arcs*: `help(rnav.plot_arcs)`

- Plots inter-nucleotide data or secondary structures as arcs

*plot_arcs_compare*: `help(rnav.plot_arcs_compare)`

- Plots two samples on the same arc plot for easy comparison

*plot_ss*: `help(rnav.plot_ss)`,

- Plots data on a secondary structure diagram

*plot_mol*: `help(rnav.plot_mol)`

- Renders an interactive 3D molecule with overlayed data

*plot_heatmap*: `help(rnav.plot_heatmap)`

- Plots inter-nucleotide data as a heatmap, overlaying structure distances

*plot_circle*: `help(rnav.plot_circle)`

- Plots inter-nucleotide data as arcs along a circle

*plot_linreg*: `help(rnav.plot_linreg)`

- Performs pair-wise linear regressions and plots a scatter plot with results

*plot_roc*: `help(rnav.plot_roc)`,

- Plots Receiver Operator Characteristic and area under the curve

*plot_disthist*: `help(rnav.plot_disthist)`

- Plots 3D or contact distance distributions for inter-nucleotide data

### styles

This submodule contains global style settings and methods for changing them.

*settings* `print(rnav.settings)`

- A dictionary of global style settings.

**Settings class** `help(rnav.Settings)`

- A context class to temporarily change style settings.

*set_defaults* `help(rnav.set_defaults)`

- Sets matplotlib and seaborn defaults: dpi, color palette, context, and style

*get_nt_color* `help(rnav.get_nt_color)`

- Given a single nucleotide, returns the color.

*get_nt_cmap* `help(rnav.get_nt_cmap)`

- Returns the ScalarMappable for nucleotide coloring.

### data_loading

This submodule contains functions for parsing Sample data keyword arguments.

*data_keyword_defaults* `print(rnav.data_loading.data_keyword_defaults)`

- A dictionary of standard data keywords with associated class and arguments

*create_data* `help(rnav.data_loading.create_data)`

- An admittedly convoluted function used to parse Sample data keyword arguments

*get_sequence* `help(rnav.data_loading.get_sequence)`

- Flexibly creates or retrieves an `rnav.data.Sequence`

### helper_functions

This submodule contains the PlottingArgumentParser, which is a standardized
way for plotting functions to interact with Samples.

*fit_data*: `help(rnav.helper_functions.fit_data)`

- Flexibly returns data which has been aligned to a target sequence

*PlottingArgumentParser*: `help(rnav.helper_functions.PlottingArgumentParser)`

- Class for parsing arguments from plotting functions, returns organized,
  aligned data

---

### transcriptomics

#### bed

#### eclip

#### transcriptome

---

## data

This module contains all of the code for aligning, manipulating, and
representing data.

### alignments

This submodule contains all of the code for aligning data between sequences.

*convert_sequence*: `help(rnav.data.convert_sequence)`

- Converts a sequence and dot-bracket notation to pseudo-amino acid sequence
  and vice versa (for secondary structure alignment)

*structure_scoring_dict*: `help(rnav.data.structure_scoring_dict)`

- Dictionary of alignment scores for structural alignments

*set_alignment*: `help(rnav.data.set_alignment)`

- Sets an alignment as the default for comparing the given sequences.

**BaseAlignment class**: `help(rnav.data.BaseAlignment)`

- An abstract base class for other alignment classes

*get_mapping*: `help(rnav.data.BaseAlignment.get_mapping)`

- Returns an indexing map, i.e. `new_index = map[old_index]`

*get_target_sequence*: `help(rnav.data.BaseAlignment.get_target_sequence)`

- Returns sequence 1 string mapped to the indices of the sequence 2

*map_values*: `help(rnav.data.BaseAlignment.map_values)`

- Takes a list (1 per position in sequence 1) and returns the aligned list

*map_indices*: `help(rnav.data.BaseAlignment.map_indices)`

- Takes a list of indices from sequence 1 and returns the aligned indices

*map_positions*: `help(rnav.data.BaseAlignment.map_positions)`

- Takes a list of positions from sequence 1 and returns the aligned positions

*map_dataframe*: `help(rnav.data.BaseAlignment.map_dataframe)`

- Takes a dataframe with indexing columns and returns the aligned dataframe

*map_nucleotide_dataframe*: `help(rnav.data.BaseAlignment.map_nucleotide_dataframe)`

- Takes a per-nucleotide dataframe and returns the aligned dataframe

*get_inverse_alignment*: `help(rnav.data.BaseAlignment.get_inverse_alignment)`

- Reverses the current alignment, i.e. seq1 -> seq2, becomes seq2 -> seq1

**SequenceAlignment class**: `help(rnav.data.SequenceAlignment)`

- A class for performing sequence alignment, uses default alignment if set

*get_inverse_alignment*: `help(rnav.data.SequenceAlignment.get_inverse_alignment)`

- Reverses the current alignment, i.e. seq1 -> seq2, becomes seq2 -> seq1

*get_alignment*: `help(rnav.data.SequenceAlignment.get_alignment)`

- Performs sequence alignment, or retreives a default alignment

*get_mapping*: `help(rnav.data.SequenceAlignment.get_mapping)`

- Returns an indexing map, i.e. `new_index = map[old_index]`

**StructureAlignment class**: `help(rnav.data.StructureAlignment)`

- A class for performing secondary structure alignments

*get_inverse_alignment*: `help(rnav.data.StructureAlignment.get_inverse_alignment)`

- Reverses the current alignment, i.e. seq1 -> seq2, becomes seq2 -> seq1

*get_alignment*: `help(rnav.data.StructureAlignment.get_alignment)`

- Performs a secondary structure alignment using psuedo-amino acid sequences

*get_mapping*: `help(rnav.data.StructureAlignment.get_mapping)`

- Returns an indexing map, i.e. `new_index = map[old_index]`

*set_as_default_alignment*: `help(rnav.data.StructureAlignment.set_as_default_alignment)`

- Sets the structure alignment as the default alignment for these sequences

**AlignmentChain class**: `help(rnav.data.AlignmentChain)`

- Combines a list of alignments serially, i.e. 1->2 and 2->3, becomes 1->3

*get_mapping*: `help(rnav.data.AlignmentChain.get_mapping)`

- Returns an indexing map, i.e. `new_index = map[old_index]`

*get_inverse_alignment*: `help(rnav.data.AlignmentChain.get_inverse_alignment)`

- Reverses the current alignment, i.e. seq1 -> seq2, becomes seq2 -> seq1

### annotation

This submodule contains all of the code for representing and manipulating
sequence annotations.

**Annotation class**: `help(rnav.data.Annotation)`

- A class for storing sites or spans of interest, primer binding sites, and
  nucleotide groups

*from_boolean_array*: `help(rnav.data.Annotation.from_boolean_array)`

- Create an annotation from a list of boolean values

*from_spans*: `help(rnav.data.Annotation.from_spans)`

- Returns a spans dataframe from a list of spans

*from_sites*: `help(rnav.data.Annotation.from_sites)`

- Returns a sites dataframe from a list of sites

*get_aligned_data*: `help(rnav.data.Annotation.get_aligned_data)`

- Returns a new annotation with annotation positions aligned to a new sequence

*get_subsequences*: `help(rnav.data.Annotation.get_subsequences)`

- Returns a list of subsequences for each annotation

**Motif class**: `help(rnav.data.Motif)`

- A class for storing locations of a motif within a sequence

*get_spans_from_motif*: `help(rnav.data.Motif.get_spans_from_motif)`

- Returns a list of spans given a motif

*get_aligned_data*: `help(rnav.data.Motif.get_aligned_data)`

- Returns a new annotation with locations of the motif in a new sequence

**ORFs**: `help(rnav.data.ORFs)`

- A class for storing locations of open reading frames

*get_spans_from_orf*: `help(rnav.data.ORFs.get_spans_from_orf)`

- Returns a list of spans that are potential open reading frames

*get_aligned_data*: `help(rnav.data.ORFs.get_aligned_data)`

- Returns the annotation of open reading frames for a new sequence

*domains*: `help(rnav.data.domains)`

- Returns a list of span annotations used for labelling domains of an RNA

### colors

ScalarMappable
is_equivalent_to
values_to_hexcolors
get_norm
get_cmap

### ct

SecondaryStructure
nts
pair_nts
ycoordinates
xcoordinates
read_ct
read_varna
read_xrna
read_cte
read_nsd
read_dotbracket
read_r2dt
read_forna
write_sto
writeRNAstructureSeq
write_ct
write_cte
get_dotbracket
get_human_dotbracket
break_noncanonical_pairs
break_singleton_pairs
get_pairs
get_paired_nts
get_unpaired_nts
get_junction_nts
add_pairs
from_pairs_list
get_masked_ct
get_nonredundant_ct
get_distance_matrix
contact_distance
get_helices
fill_mismatches
extractPK
compute_ppv_sens
copy
get_aligned_data
get_interactions_df
as_interactions
transform_coordinates

StructureCoordinates
get_center_point
scale
flip
center
rotate

SequenceCircle

### data

Sequence
read_fasta
get_seq_from_dataframe
length
get_aligned_data
get_colors_from_sequence
get_colors_from_positions
get_colors_from_profile
get_colors_from_annotations
get_colors_from_structure
get_colors

Data
add_metric_defaults
metric
error_column
color_column
cmap
colors
read_file

### interactions

Interactions
mask_on_sequence
mask_on_structure
mask_on_profile
mask_on_position
mask_on_distance
mask_on_values
reset_mask
get_aligned_data
update_mask
filter
data_specific_filter
get_ij_colors
get_sorted_data
print_new_file
set_3d_distances
resolve_conflicts

SHAPEJuMP
read_file

RINGMaP
read_file
data_specific_filter
get_sorted_data

PAIRMaP
read_file
data_specific_filter
get_sorted_data

PairingProbability
read_file
get_entropy_profile
data_specific_filter

AllPossible

StructureInteractions

### log

Log
read_log

### pdb

PDB
get_sequence
get_sequence_from_seqres
read_pdb
set_indices
get_pdb_idx
get_seq_idx
is_valid_idx
get_xyz_coord
get_distance
get_distance_matrix

### profile

Profile
recreation_kwargs
get_aligned_data
get_plotting_dataframe
calculate_windows
calculate_gini_index
normalize
normalize_external
norm_boxplot
norm_eDMS
norm_percentiles

SHAPEMaP

DanceMaP
recreation_kwargs
read_file

RNPMaP

DeltaProfile

---

## plots module

### plots

Plot
get_ax
add_colorbar_args
plot_colorbars
get_rows_columns
save
set_figure_size

Colorbar
plot_data
set_figure_size

### functions

plot_interactions_circle
plot_annotation_circle
get_contrasting_colors
adjust_spines
clip_spines
set_x_ticks
box_xtick_labels
plot_sequence_alignment
plot_interactions_arcs
plot_profile_bars
plot_profile_skyline
plot_structure_ss
plot_basepairs_ss
plot_sequence_ss
plot_nucleotides_ss
plot_positions_ss
plot_interactions_ss
plot_annotation_ss
plot_sequence_track
plot_annotation_track
plot_domain_track

### alignment

Alignment
plot_data
set_figure_size

### arc

AP
plot_data
set_axis
set_figure_size

### circle

Circle
set_figure_size
plot_data

### disthist

DistHist
set_figure_size
plot_data
plot_structure_distances
plot_experimental_distances

### heatmap

Heatmap
plot_data
set_figure_size
plot_contour_regions
plot_contour_distances
plot_heatmap_data
plot_kde_data

### lingreg

LinReg
set_figure_size
get_rows_columns
plot_data
finalize
plot_regression
plot_kde

### mol

Mol
get_viewer
plot_data
add_lines
plot_interactions
set_colors
hide_cylinders
save
get_orientation

### profile

Profile
get_rows_columns
set_labels
plot_data
set_figure_size
set_axis

### qc

QC
plot_data
set_figure_size
get_rows_columns
plot_MutsPerMol
plot_ReadLength
make_boxplot

### roc

ROC
plot_data
get_rows_columns
set_figure_size

### skyline

Skyline
plot_data
set_figure_size
get_rows_columns
get_ax
plot_data
set_axis
set_labels

### sm

SM
set_figure_size
plot_data
plot_sm_profile
plot_sm_depth
plot_sm_rates
metric_abbreviate

### ss

SS
set_figure_size
plot_data

---

## analysis module

### fragmapper submodule

Fragmapper
plot_scatter

FragMaP
get_dataframe
calc_zscore
get_annotation

### auroc submodule

WindowedAUROC
plot_auroc

### deltashape submodule

DeltaSHAPE
calculate_deltashape
plot

DeltaSHAPEProfile
calculate_deltashape
get_protections_annotation
get_enhancements_annotation

### logcompare submodule

LogCompare
get_profile_sequence
calc_scale_factor
rescale
load_replicates
make_plots

### lowss submodule

LowSS
reset_lowss
reset_window
plot_lowss
