
RNAvigate full API
==================

This is RNAvigate's (mostly) full API. I apologize that it is a bit ugly at the
moment, and that some docstrings are missing. Good documentation is the current top
priority for RNAvigate's development. I plan to tackle this process in this order:

1. Add tutorials for common tasks
2. Fill in any missing docstrings and update to Numpy style
3. Implement mkdocstrings for automatic documentation
4. Eventually, transition to Sphinx-based documentation

Table of contents:

- rnavigate
    - [Sample](#rnavigatesample)
    - [plot_options](#rnavigateplot_options)
    - [plot_alignment](#rnavigateplot_alignment)
    - [plot_arcs](#rnavigateplot_arcs)
    - [plot_arcs_compare](#rnavigateplot_arcs_compare)
    - [plot_circle](#rnavigateplot_circle)
    - [plot_disthist](#rnavigateplot_disthist)
    - [plot_heatmap](#rnavigateplot_heatmap)
    - [plot_linreg](#rnavigateplot_linreg)
    - [plot_mol](#rnavigateplot_mol)
    - [plot_ntdist](#rnavigateplot_ntdist)
    - [plot_profile](#rnavigateplot_profile)
    - [plot_qc](#rnavigateplot_qc)
    - [plot_roc](#rnavigateplot_roc)
    - [plot_shapemapper](#rnavigateplot_shapemapper)
    - [plot_skyline](#rnavigateplot_skyline)
    - [plot_ss](#rnavigateplot_ss)
    - [analysis](#rnavigateanalysis)
        - [LogCompare](#rnavigateanalysislogcompare)
        - [LowSS](#rnavigateanalysislowss)
        - [DeltaSHAPE](#rnavigateanalysisdeltashape)
        - [DeltaSHAPEProfile](#rnavigateanalysisdeltashapeprofile)
        - [WindowedAUROC](#rnavigateanalysiswindowedauroc)
        - [FragMaP](#rnavigateanalysisfragmap)
        - [Fragmapper](#rnavigateanalysisfragmapper)
        - [SequenceChecker](#rnavigateanalysissequencechecker)
    - [data](#rnavigatedata)
        - [Sequence](#rnavigatedatasequence)
        - [Data](#rnavigatedatadata)
        - [SecondaryStructure](#rnavigatedatasecondarystructure)
        - [StructureCoordinates](#rnavigatedatastructurecoordinates)
        - [SequenceCircle](#rnavigatedatasequencecircle)
        - [set_alignment](#rnavigatedataset_alignment)
        - [set_multiple_sequence_alignment](#rnavigatedataset_multiple_sequence_alignment)
        - [lookup_alignment](#rnavigatedatalookup_alignment)
        - [SequenceAlignment](#rnavigatedatasequencealignment)
        - [AlignmentChain](#rnavigatedataalignmentchain)
        - [StructureAlignment](#rnavigatedatastructurealignment)
        - [ScalarMappable](#rnavigatedatascalarmappable)
        - [Interactions](#rnavigatedatainteractions)
        - [SHAPEJuMP](#rnavigatedatashapejump)
        - [RINGMaP](#rnavigatedataringmap)
        - [PAIRMaP](#rnavigatedatapairmap)
        - [PairingProbability](#rnavigatedatapairingprobability)
        - [AllPossible](#rnavigatedataallpossible)
        - [StructureAsInteractions](#rnavigatedatastructureasinteractions)
        - [StructureCompareMany](#rnavigatedatastructurecomparemany)
        - [StructureCompareTwo](#rnavigatedatastructurecomparetwo)
        - [PDB](#rnavigatedatapdb)
        - [Profile](#rnavigatedataprofile)
        - [SHAPEMaP](#rnavigatedatashapemap)
        - [DanceMaP](#rnavigatedatadancemap)
        - [RNPMaP](#rnavigatedatarnpmap)
        - [DeltaProfile](#rnavigatedatadeltaprofile)
        - [Annotation](#rnavigatedataannotation)
        - [Motif](#rnavigatedatamotif)
        - [ORFs](#rnavigatedataorfs)
        - [domains](#rnavigatedatadomains)
    - [plots](#rnavigateplots)
        - [Alignment](#rnavigateplotsalignment)
        - [AP](#rnavigateplotsap)
        - [Circle](#rnavigateplotscircle)
        - [DistHist](#rnavigateplotsdisthist)
        - [Heatmap](#rnavigateplotsheatmap)
        - [LinReg](#rnavigateplotslinreg)
        - [NucleotideDistribution](#rnavigateplotsnucleotidedistribution)
        - [Mol](#rnavigateplotsmol)
        - [Plot](#rnavigateplotsplot)
        - [ColorBar](#rnavigateplotscolorbar)
        - [Profile](#rnavigateplotsprofile)
        - [QC](#rnavigateplotsqc)
        - [ROC](#rnavigateplotsroc)
        - [Skyline](#rnavigateplotsskyline)
        - [SM](#rnavigateplotssm)
        - [SS](#rnavigateplotsss)
        - [get_contrasting_colors](#rnavigateplotsget_contrasting_colors)
        - [adjust_spines](#rnavigateplotsadjust_spines)
        - [clip_spines](#rnavigateplotsclip_spines)
        - [get_nt_ticks](#rnavigateplotsget_nt_ticks)
        - [set_nt_ticks](#rnavigateplotsset_nt_ticks)
        - [box_xtick_labels](#rnavigateplotsbox_xtick_labels)
        - [plot_interactions_arcs](#rnavigateplotsplot_interactions_arcs)
        - [plot_profile_bars](#rnavigateplotsplot_profile_bars)
        - [plot_profile_skyline](#rnavigateplotsplot_profile_skyline)
        - [plot_sequence_alignment](#rnavigateplotsplot_sequence_alignment)
        - [plot_annotation_track](#rnavigateplotsplot_annotation_track)
        - [plot_domain_track](#rnavigateplotsplot_domain_track)
        - [plot_sequence_track](#rnavigateplotsplot_sequence_track)
        - [plot_annotation_ss](#rnavigateplotsplot_annotation_ss)
        - [plot_basepairs_ss](#rnavigateplotsplot_basepairs_ss)
        - [plot_interactions_ss](#rnavigateplotsplot_interactions_ss)
        - [plot_nucleotides_ss](#rnavigateplotsplot_nucleotides_ss)
        - [plot_positions_ss](#rnavigateplotsplot_positions_ss)
        - [plot_sequence_ss](#rnavigateplotsplot_sequence_ss)
        - [plot_structure_ss](#rnavigateplotsplot_structure_ss)
        - [plot_annotation_circle](#rnavigateplotsplot_annotation_circle)
        - [plot_interactions_circle](#rnavigateplotsplot_interactions_circle)
    - [styles](#rnavigatestyles)
        - [update_copy](#rnavigatestylesupdate_copy)
        - [Settings](#rnavigatestylessettings)
        - [set_defaults](#rnavigatestylesset_defaults)
        - [get_nt_color](#rnavigatestylesget_nt_color)
        - [get_nt_cmap](#rnavigatestylesget_nt_cmap)
        - [apply_style](#rnavigatestylesapply_style)
    - [transcriptomics](#rnavigatetranscriptomics)
        - [BedFile](#rnavigatetranscriptomicsbedfile)
        - [NarrowPeak](#rnavigatetranscriptomicsnarrowpeak)
        - [Transcriptome](#rnavigatetranscriptomicstranscriptome)
        - [Transcript](#rnavigatetranscriptomicstranscript)
        - [eCLIPDatabase](#rnavigatetranscriptomicseclipdatabase)
        - [download_eclip_peaks](#rnavigatetranscriptomicsdownload_eclip_peaks)
## rnavigate.Sample

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Sample in module rnavigate.rnavigate

class Sample(builtins.object)
 |  Sample(sample, inherit=None, keep_inherited_defaults=True, **data_keywords)
 |  
 |  Loads and organizes RNA structural data for use with plotting functions.
 |  
 |  The Sample class stores all of the relevant experimental and computational
 |  structural data for a single RNA experiment. Between samples, common data
 |  types should be given a common data keyword so they can be easily compared.
 |  
 |  Parameters
 |  ----------
 |  sample : str
 |      an arbitrary name. This will be used as a label in plot legends and titles
 |      to differentiate it from other samples
 |  inherit : Sample or list of these, optional
 |      Data keywords and associated data from other samples become data keywords
 |      and associated data of this sample. This does not make additional copies
 |      of the data: i.e. operations that make changes to inherited data change the
 |      original sample, and any other samples that inherited that data. This can
 |      be useful to save time and memory on operations and large data structures
 |      that are shared between samples.
 |  keep_inherited_defaults : bool, default = True
 |      whether to keep inherited default keywords
 |  **data_keywords:
 |      There are many built-in data keywords with different expectations and
 |      behaviors. For a full list with input formats and output behavior, visit:
 |      https://rnavigate.readthedocs.io/en/latest/loading-data/
 |  
 |  Attributes
 |  ----------
 |  sample : str
 |      the name of the sample
 |  inputs : dict
 |      a dictionary of data keywords and their (user-defined) inputs
 |  data : dict
 |      a dictionary of data keywords and their associated data
 |  defaults : dict
 |      a dictionary of data classes and their default data keywords
 |  
 |  Example
 |  -------
 |  >>> sample = rnavigate.Sample(
 |  ...     sample="My sample",
 |  ...     shapemap="path/to/shapmapper_profile.txt",
 |  ...     ss="path/to/structure.ct",
 |  ...     ringmap="path/to/ringmapper_rings.txt",
 |  ...     pdb="path/to/pdb.pdb",
 |  ...     arbitrary_keyword={
 |  ...         "sites": [10, 20, 30],
 |  ...         "name": "sites of interest",
 |  ...         "color": "red",
 |  ...     },
 |  ... )
 |  >>> sample.print_data_keywords()
 |  My sample data keywords:
 |    annotations:
 |      arbitrary_keyword (default)
 |    profiles:
 |      shapemap (default)
 |    structures:
 |      ss (default)
 |    interactions:
 |      ringmap (default)
 |    pdbs:
 |      pdb (default)
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sample, inherit=None, keep_inherited_defaults=True, **data_keywords)
 |      Creates a Sample.
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and applies filters to interactions data.
 |      
 |      Parameters
 |      ----------
 |      interactions : rnavigate.data.Interactions or data keyword string
 |          Interactions object to be filtered. If a string, value is
 |          replaced with self.get_data(interactions)
 |      metric : str, optional
 |          column of interactions data to be used as metric for coloring
 |          interactions.
 |          "Distance" will compute 3D distance in "pdb", defaulting to
 |          2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |          those atoms to compute distance.
 |      cmap (str | list, optional):
 |          sets the interactions colormap, used to color interactions
 |          according to metric values.
 |      normalization (str, optional):
 |          `"norm"`: extreme values in colormap are given to the extreme
 |              values of interactions metric data
 |          `"bins"`: data are colored according to which bin they fall in
 |              `values` defines bins (list, length = 2 less than cmap)
 |          `"min_max"`: extreme values in cmap are given to values beyond
 |              minimum and maximum, defined by `values`
 |      values:
 |          behavior depends on normalization
 |          `"norm"`: values are not needed
 |          `"bins"`: list of floats containing the boundaries between bins
 |              One fewer than the number of categories
 |          `"min_max"`: list of floats containing the minimum and maximum
 |      **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : rnavigate.data.Data or data keyword or list/dict of these
 |          If None, returns None.
 |          If a data keyword, returns associated data from sample
 |          If Data, returns that data.
 |          If a list or dictionary, returns list or dictionary with
 |              data keyword values replaced with associated Data
 |      data_class : rnavigate.data.Data class or subclass, optional
 |          If provided, ensures that returned data is of this type.
 |      
 |      Returns
 |      -------
 |      Same type as data_keyword argument, but data keywords are replaced
 |          with associated data
 |      
 |      Raises
 |      ------
 |      ValueError:
 |          if data is not found in sample
 |      ValueError:
 |          if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Parameters
 |      ----------
 |      inherit : Sample or list of Samples
 |          Other samples from which to inherit data and data keywords
 |      keep_inherited_defaults : bool
 |          Use default values from inherited samples
 |      overwrite : bool
 |          whether to overwrite any existing keywords with inherited keywords
 |  
 |  print_data_keywords(self, return_dict=False)
 |      Print a nicely formatted, organized list of data keywords.
 |      
 |      Returns a dictionary of data keywords, organized by data type, if
 |      return_dict is True.
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          The data keyword to set as the default
 |      overwrite : bool, defaults to ``True``
 |          whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name",
 |              data_keyword=inputs)
 |      
 |      is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name")
 |          my_sample.add_data(
 |              "data_keyword", inputs)
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          a data keyword used to store and/or parse the inputs
 |      inputs : dict or rnavigate.data.Data
 |          a dictionary used to create the data object or a data object itself
 |      overwrite : bool, defaults to False
 |          whether to overwrite a pre-existing data_keyword
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

## rnavigate.plot_options

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_options in module rnavigate.plotting_functions

plot_options(samples)
    Prints a list of plotting functions compatible with a sample or list of samples.
    
    Some plotting functions require specific data classes to be loaded into the
    sample. For plotting multiple samples, data keywords that are not shared, or are
    shared, but are not of the same data class, are considered invalid.
    
    Parameters
    ----------
    samples : rnavigate.Sample or list of rnavigate.Sample
        samples to check for compatible plotting functions
```

## rnavigate.plot_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_alignment in module rnavigate.plotting_functions

plot_alignment(data1, data2, labels=None, plot_kwargs=None)
    Plots the sequence alignment used to compare two sequences
    
    Parameters
    ----------
    data1 : tuple (rnavigate Sample, data keyword string)
        a sample and data keyword to retrieve a sequence
    data2 : tuple (rnavigate Sample, data keyword string)
        another sample and data keyword to retrieve a second sequence
    labels : list of 2 strings, defaults to "sample.sample: data keyword" for each
        Labels used for each sample
    plot_kwargs : dict, defaults to {}
        passed to matplotlib.pyplot.subplots()
    
    Returns
    -------
    rnavigate.plots.Alignment
        the Alignment plot object
```

## rnavigate.plot_arcs

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_arcs in module rnavigate.plotting_functions

plot_arcs(samples, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, annotations=None, domains=None, labels=None, nt_ticks=(20, 5), profile_scale_factor=1, plot_error=False, annotation_mode='track', panels=None, seqbar=True, region='all', colorbars=True, title=True, plot_kwargs=None)
    Plots interactions and/or base-pairs as arcs.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    sequence : data keyword string, data object, or sequence string
        All data are mapped to this sequence before plotting
        If a data keyword string, data from the first sample will be used
    structure : data keyword string or data object, defaults to None
        secondary structure to plot as arcs
    structure2 : data keyword string or data object, defaults to None
        another secondary structure to compare with the first structure
        arcs will be colored depending on which structure they are in
        Defaults to None
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot as arcs, no filtering performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    interactions2 : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot as arcs, no filtering performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
    profile : data or data keyword, defaults to None
        Profile from which values will be plotted
    annotations : list of data keyword strings or data objects, defaults to []
        Annotations used to highlight regions or sites of interest
    domains : data keyword string or data object, defaults to None
        domains to label along x-axis
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    nt_ticks : tuple of two integers, defaults to (20, 5)
        first integer is the gap between major tick marks
        second integer is the gap between minor tick marks
    profile_scale_factor : number, defaults to 1
        small profile values will be hard to see
        large profile values will overwhelm the plot
        e.g. use 1/10 to scale values down 10-fold, use 10 to scale up
    plot_error : bool, defaults to False
        Whether to plot error bars, values are determined by profile.metric
    annotation_mode : "track" or "bars", defaults to "track"
        "track" will highlight annotations along the x-axis
        "bars" will use a vertical transparent bar over the plot
    panels : dict, optional
        a dictionary of whether plot elements are displayed on the "top"
        (above x-axis) or "bottom" (below x-axis)
        Only the values you wish to change from the default are needed
        defaults to {"interactions": "bottom",
                        "interactions2": "bottom",
                        "structure": "top",
                        "profile": "top"}
    seqbar : bool, defaults to ``True``
        whether to display the sequence along the x-axis
    region : list of 2 integers, defaults to [1, length of sequence]
        start and end positions to plot. 1-indexed, inclusive.
    colorbars : bool, defaults to True
        Whether to plot colorbars for all plot elements
    title : bool, defaults to True
        Whether to display titles for each axis
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.AP
        the ArcPlot object
```

## rnavigate.plot_arcs_compare

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_arcs_compare in module rnavigate.plotting_functions

plot_arcs_compare(samples, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, labels=None, profile_scale_factor=1, plot_error=False, region='all', colorbars=True, plot_kwargs=None)
    Generates a single arc plot displaying combinations of secondary
    structures, per-nucleotide data, inter-nucleotide data, and sequence
    annotations. The first sample will be on top, the second on the bottom.
    Center shows how these sequences are being aligned. This view does not
    
    Parameters
    ----------
    samples : list of 2 rnavigate Samples
        samples used to retrieve data
        This plotting function can only compare two samples at a time
    sequence : data keyword string, data object, or sequence string
        All data are mapped to this sequence taken from their respective
        sample before plotting
    structure : data keyword string or data object, defaults to None
        secondary structure to plot as arcs
    structure2 : data keyword string or data object, defaults to None
        another secondary structure to compare with the first structure
        arcs will be colored depending on which structure they are in
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot as arcs, no filtering performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    interactions2 : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot as arcs, no filtering performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
    profile : data keyword string or data object, defaults to None
        Profile from which values will be plotted
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    profile_scale_factor : number, defaults to 1
        small profile values will be hard to see
        large profile values will overwhelm the plot
        e.g. use 1/10 to scale values down 10-fold, use 10 to scale up
    plot_error : bool, defaults to ``False``
        Whether to plot error bars, values are determined by profile.metric
    region : list of 2 integers, defaults to [1, length of sequence]
        start and end positions to plot. 1-indexed, inclusive.
    colorbars : bool, defaults to ``True``
        Whether to plot color scales for all plot elements
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.AP plot
        object containing matplotlib figure and axes with additional plotting
        and file saving methods
```

## rnavigate.plot_circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_circle in module rnavigate.plotting_functions

plot_circle(samples, sequence, structure=None, structure2=None, interactions=None, interactions2=None, annotations=None, profile=None, colors=None, nt_ticks=(20, 5), gap=30, labels=None, colorbars=True, plot_kwargs=None)
    Creates a figure containing a circle plot for each sample given.
    
    Data that can be plotted on circle plots includes annotations (highlights
    regions around the edge.Generates a multipanel secondary structure drawing
    with optional coloring by per-nucleotide data and display of inter-
    nucleotide data and/or sequence annotations. Each plot may display a unique
    sample and/or inter-nucleotide data filtering scheme.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    sequence : data or data keyword
        All data are mapped to this sequence before plotting
    structure : data keyword string or data object, defaults to None
        Structure used to plot base-pairs on circle plot
    structure2 : data keyword str, data obj or list of either, defaults to None
        Structures to compare with Structure. Each base-pair is colored by
        which structure contains it or how many structures contain it.
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot on cirle plot, no filtering
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    interactions2 : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot on circle plot, no filtering
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
    annotations : list of data keyword strings or data objects, defaults to []
        Annotations used to highlight regions or sites of interest
    profile : data keyword string or data object, defaults to None
        Profile used for coloring if "profile" used in colors dictionary
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    colors : dictionary, optional
        a dictionary of element: value pairs that determines how colors
        will be applied to each plot element and if that element is plotted
        only the elements you wish to change need to be included
        Defaults to {"sequence": None,  # sequence not shown
                        "nucleotides": "sequence",
                        "structure": "grey"}
        value options and what the colors represent:
            None: don't plot this elelement
            "sequence": nucleotide identity
            "position": position in sequence
            "annotations": sequence annotations
            "profile": per-nucleotide data from profile
                profile argument must be provided
            "structure": base-pairing status
            matplotlib color: all positions plotted in this color
            array of colors: a different color for each position
                must be the same length as structure
        "sequence" may also use "contrast" which automatically chooses
            white or black, which ever contrasts better with "nucleotide"
            color
    nt_ticks : tuple of two integers, defaults to (20, 5)
        first integer is the gap between major tick marks
        second integer is the gap between minor tick marks
    gap : integer, defaults to 30
        Width of gap between 5' and 3' end in degrees
    colorbars : bool, defaults to True
        Whether to plot color scales for all plot elements
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.Circle
        object containing matplotlib figure and axes with additional plotting and
        file saving methods
```

## rnavigate.plot_disthist

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_disthist in module rnavigate.plotting_functions

plot_disthist(samples, structure, interactions, bg_interactions=None, labels=None, same_axis=False, atom="O2'", rows=None, cols=None, plot_kwargs=None)
    Calculates 3D distance of nucleotides in inter-nucleotide data and plots
    the distribution of these distances. Compares this to a "background"
    distribution consisting of either all pairwise distances in structure, or
    those defined by bg_interactions and bg_interactions_filter
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        Samples from which to retreive data
        There will be one panel for each sample unless same_axis is True
    structure : data keyword string or data object
        secondary structure or 3D structure to calculate inter-nucleotide
        contact distance or 3D distance, respectively
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions used to calculate distance histogram, no filtering
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    bg_interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to calculate background distance histogram, no
            filtering is performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        if not provided, background distance histogram is calculated from
        all pairwise distances in structure
    labels : list of strings, defaults to sample.sample for each sample
        Labels to be used as titles, must be same length as samples list
        Defaults to sample.sample for each sample
    atom : string or dictionary, defaults to "O2'"
        from which atoms to calculate distances
        for DMS reactive atoms (N1 for A and G, N3 for U and C) use "DMS"
        use a dictionary to specify a different atom for each nucleotide
            e.g. "DMS" == {"A": "N1", "G": "N1", "U": "N3", "C": "N3"}
    rows : integer, defaults to None (determined automatically)
        number of rows of plots
    cols : integer, defaults to None (determined automatically)
        number of columns of plots
    plot_kwargs : dictionary, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.DistHist
        object containing matplotlib figure and axes with additional plotting and
        file saving methods
```

## rnavigate.plot_heatmap

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_heatmap in module rnavigate.plotting_functions

plot_heatmap(samples, sequence, structure=None, interactions=None, regions=None, labels=None, levels=None, interpolation='nearest', atom="O2'", plot_type='heatmap', weights=None, rows=None, cols=None, plot_kwargs=None)
    Generates a multipanel plot displaying a heatmap of inter-nucleotide
    data (nucleotide resolution of 2D KDE) and/or contour map of pdb
    distances. Each plot may display a unique sample and/or filtering scheme.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    sequence : data keyword string, data object, or sequence string
        All data are mapped to this sequence before plotting
    structure : data keyword string or data object, defaults to None
        secondary structure or 3D structure used to plot contour lines
        contour lines are drawn according to levels argument
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot as a heatmap, no filtering performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    regions : list of lists of 4 integers, defaults to None (no boxes)
        each inner list defines two regions of the RNA that are interacting
        a box will be drawn around this interaction on the heatmap
        e.g. [[10, 20, 50, 60], [35, 45, 70, 80]] draws 2 boxes
            the first box will connect nucleotides 10-20 and 50-60
            the second box will connect nucleotides 35-45 and 70-80
    labels : list of strings, defaults to sample.sample for each sample
        Labels to be used as titles, must be same length as samples list
    levels : list of floats, defaults to [5] contact distance or [20] 3D distance
        contours are drawn separating nucleotides above and below these
        distances
        if structure argument is a secondary structure
            distance refers to contact distance
        if structure argument is a 3D structure
            distance refers to spatial distance in angstroms
    interpolation : string, defaults to "nearest"
        one of matplotlib's interpolations for heatmap (used with imshow)
        "nearest" works well for shorter RNAs (under 300 nt)
        "none" works well for longer RNAs (over 1200 nt)
    atom : string or dictionary, defaults to "O2'"
        from which atoms to calculate distances
        for DMS reactive atoms (N1 for A and G, N3 for U and C) use "DMS"
        use a dictionary to specify a different atom for each nucleotide
            e.g. "DMS" == {"A": "N1", "G": "N1", "U": "N3", "C": "N3"}
    plot_type : "heatmap" or "kde", defaults to "heatmap"
        how to plot interactions data
        "heatmap" will plot raw data, each interaction is a pixel in a grid
        "kde" will calculate a kernel density estimate and plot 5 levels
    weights : string, defaults to None (no weights)
        weights to be used in kernel density estimation
        must be a column of interactions data
    rows : integer, defaults to None (determined automatically)
        number of rows of plots
    cols : integer, defaults to None (determined automatically)
        number of columns of plots
    plot_kwargs : dictionary, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.Heatmap
        object containing matplotlib figure and axes with additional plotting and
        file saving methods
```

## rnavigate.plot_linreg

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_linreg in module rnavigate.plotting_functions

plot_linreg(samples, profile, sequence=None, structure=None, annotations=None, labels=None, kde=False, scale='linear', regression='pearson', colors='sequence', column=None, region='all', colorbars=True, plot_kwargs=None)
    Performs linear regression analysis and generates scatter plots of all
    sample-to-sample profile vs. profile comparisons. Colors nucleotides by
    identity or base-pairing status.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    profile : data keyword string or data object
        per-nucleotide data to perform linear regression
        all data are mapped to the sequence of the profile data from the
        first sample before plotting, unless sequence is supplied
    sequence : data keyword str, data obj, or sequence str, defaults to None
        a sequence from which to align all profiles
        if a data keyword, uses data from the first sample
    structure : data keyword string or data object, defaults to None
        Structure used for coloring if colors argument is "structure"
    annotations : list of data keyword strings or data objects, defaults to []
        Annotations used for coloring if colors argument is "annotations"
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    kde : bool, defaults to False
        whether to plot kde (density) instead of a scatter plot
    scale : "linear" or "log", defaults to "linear"
        "linear" performs regression on raw values, displays linear units
        "log" performs regression on log10(values), displays log10 units
    regression : "pearson" or "spearman", defaults to "pearson"
        "pearson" calculates Pearson R-squared (standard)
        "spearman" calculates Spearman R-squared (rank-order)
    colors : string or list of colors, defaults to "sequence"
        value options and what the colors represent:
            "sequence": nucleotide identity
            "position": position in sequence
            "annotations": sequence annotations
            "profile": per-nucleotide data from profile
                profile argument must be provided
            "structure": base-pairing status
            matplotlib color: all positions plotted in this color
            array of colors: a different color for each position
                must be the same length as structure
    column : string, defaults to profile.metric
        column name of values from profile to use in regression
    region : list of 2 integers, defaults to [1, length of sequence]
        start and end nucleotide positions to include. 1-indexed, inclusive
    colorbars : bool, defaults to ``True``
        Whether to plot colorbars for scatter plot colors
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.LinReg
        object containing matplotlib figure and axes with additional plotting and
        file saving methods
```

## rnavigate.plot_mol

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_mol in module rnavigate.plotting_functions

plot_mol(samples, structure, profile=None, interactions=None, labels=None, style='cartoon', hide_cylinders=False, colors='grey', atom="O2'", rotation=None, orientation=None, get_orientation=False, title=True, colorbars=True, width=400, height=400, rows=None, cols=None, background_alpha=1, show=True)
    Generates a multipanel interactive 3D molecular rendering of a PDB
    structure. Nucleotides may be colored by per-nucleotide data or custom
    color lists. Inter-nucleotide data may be displayed as cylinders connecting
    atoms or residues. Each plot may display a unique sample and/or filtering
    scheme.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    structure : data keyword string or data object
        3D structure to view as interactive molecule
        All data are mapped to this sequence before plotting
    profile : data keyword string or data object, defaults to None
        Profile used to color nucleotides if colors="profile"
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot on molecule, no filtering performed
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot titles
    style : "cartoon", "cross", "line", "sphere" or "stick", defaults to "cartoon"
        sets the py3Dmol style for drawing the molecule
    hide_cylinders : bool, defaults to False
        whether to hide nucleotide cylinders (only shows backbone ribbon)
    colors : string or list of colors, defaults to "grey"
        "sequence": color by nucleotide identity
        "position": color by position in sequence
        "annotations": color by sequence annotations from `annotations`
        "profile": color by per-nucleotide data from `profile`
        "structure": color by base-pairing status
        matplotlib color: all positions plotted in this color
        array of colors: a different color for each position
            must be the same length as structure
    atom : string or dictionary, defaults to "O2'"
        which atoms to draw interactions between
        for DMS reactive atoms (N1 for A and G, N3 for U and C) use "DMS"
        use a dictionary to specify a different atom for each nucleotide
            e.g. "DMS" == {"A": "N1", "G": "N1", "U": "N3", "C": "N3"}
    rotation : dictionary, defaults to {"x": 0, "y": 0, "z": 0}
        axis-degrees pairs for setting the starting orientation of the
        molecule, only the axes to be rotated are needed
        e.g. {"x": 180} flips the molecule on the x-axis
    orientation : list of 9 floats, defaults to None
        set the precise starting orientation
        see get_orientation for more details
    get_orientation : bool, defaults to False
        allows getting the orientation for use with orientation argument
        all other arguments will be ignored and a larger, single panel view
        window is displayed with no title
            1. adjust the molecule to the desired orientation
            2. click on the molecule to display the orientation vector
            3. copy this orientation vector (manually)
            4. provide this list of values to the orientation argument
    title : bool, defaults to True
        whether to display the title
    colorbars : bool, defaults to True
        Whether to plot color scales for all plot elements
    width : integer, defaults to 400
        width of view window in pixels
    height : integer, defaults to 400
        height of view window in pixels
    rows : integer, defaults to None (set automatically)
        the number of rows in the view window
    cols : integer, defaults to None (set automatically)
        the number of columns in the view window
    background_alpha : float, defaults to 1 (completely opaque)
        the opacity of the view window, must be between 0 and 1
    show : bool, defaults to True
        whether to display the viewer object
    
    Returns
    -------
    rnavigate.plots.Mol:
        object containing py3dmol viewer with additional plotting and file saving
        methods
```

## rnavigate.plot_ntdist

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_ntdist in module rnavigate.plotting_functions

plot_ntdist(samples, profile, labels=None, column=None, plot_kwargs=None)
    Plots the distributions of values at A, U, C, and G.
    
    Calculates the kernel density estimate (KDE) for each nucleobase and plots
    them on one axis per sample.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    profile : data keyword string or data object
        per-nucleotide data to plot per-nt-identity distributions
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    column : string, defaults to profile.metric
        which column of data to use for KDE
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.NucleotideDistribution
        object containing matplotlib figure and axes with additional
        plotting and file saving methods
```

## rnavigate.plot_profile

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_profile in module rnavigate.plotting_functions

plot_profile(samples, profile, sequence=None, annotations=None, domains=None, labels=None, nt_ticks=(20, 5), column=None, plot_error=True, annotations_mode='track', seqbar=True, region='all', colorbars=True, plot_kwargs=None)
    Aligns reactivity profiles by sequence and plots them on seperate axes.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    profile : data keyword string or data object
        Profile from which values will be plotted
    sequence : data keyword str, data obj, or sequence str, defaults to `profile`
        All data are mapped to this sequence before plotting
        If a data keyword, data from the first sample will be used
    annotations : list of data keyword strings or data objects, defaults to []
        Annotations used to highlight regions or sites of interest
    domains : data keyword string or data object, defaults to None
        domains to label along x-axis
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    nt_ticks : tuple of two integers, defaults to (20, 5)
        first integer is the gap between major tick marks
        second integer is the gap between minor tick marks
    column : string, defaults to profile.metric
        column name of values from profile to plot
    plot_error : bool, defaults to True
        Whether to plot error bars, values are determined by profile.metric
    annotations_mode : "track" or "bars", defaults to "track"
        "track" will highlight annotations along the x-axis
        "bars" will use a vertical transparent bar over the plot
    seqbar : bool, defaults to ``True``
        whether to display the sequence along the x-axis
    region : list of 2 integers, defaults to [1, length of sequence]
        start and end positions to plot. 1-indexed, inclusive.
    colorbars : bool, defaults to True
        Whether to plot color scales for per-nucleotide data
    plot_kwargs : dictionary, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.Profile
        the Profile plot object
```

## rnavigate.plot_qc

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_qc in module rnavigate.plotting_functions

plot_qc(samples, profile, labels=None)
    Creates a multipanel quality control plot displaying mutations per
    molecule, read length distribution, and mutation rate distributions for
    modified and unmodified samples.
    
    Parameters
    ----------
    samples : list of rnavigate.Sample
        samples to retrieve data from
    profile : data keyword string or data object
        ShapeMaP or similar data for plotting reactivity distributions
        Must contain data from ShapeMapper log file
    labels : list of str, defaults to sample.sample for each sample in `samples`
        labels to be used on legends, must be same length as samples list
    
    Returns
    -------
    rnavigate.plots.QC
        the quality control plot object
```

## rnavigate.plot_roc

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_roc in module rnavigate.plotting_functions

plot_roc(samples, structure, profile, labels=None, nts='AUCG', plot_kwargs=None)
    Performs receiver operator characteristic analysis (ROC), calculates
    area under ROC curve (AUC), and generates ROC plots to assess how well
    per-nucleotide data predicts base-paired status. Does this for all
    positions as well as positions categorized by nucleotide
    5 plots: All, A, U, C, G
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    structure : data keyword string or data object
        secondary structure to use as classifier (paired or unpaired)
        profile data for each sample is first aligned to this structure
    profile : data keyword string or data object
        per-nucleotide data to perform ROC analysis
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    nts : string, defaults to "AUCG"
        which nucleotides to plot nucleotide-type ROC plots
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.ROC
        object containing matplotlib figure and axes with additional plotting and
        file saving methods
```

## rnavigate.plot_shapemapper

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_shapemapper in module rnavigate.plotting_functions

plot_shapemapper(sample, profile, label=None, panels=None)
    Makes a standard ShapeMapper2 profile plot with 3 panels: Normalized
    Reactivities, mutation rates, and read depths.
    
    Parameters
    ----------
    sample : rnavigate Sample
        The sample from which data profile and label will be retreived
    profile : data keyword string or data object
        ShapeMaP or similar data for plotting profiles
    label : str, defaults to sample.sample
        A label to use as the title of the figure
    panels : list of str, defaults to ["profile", "rates", "depth"]
        Which panels to include: options are "profile", "rates", and "depth"
    
    Returns
    -------
    rnavigate.plots.SM
        the ShapeMapper2 plot object
```

## rnavigate.plot_skyline

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_skyline in module rnavigate.plotting_functions

plot_skyline(samples, profile, sequence=None, annotations=None, domains=None, labels=None, nt_ticks=(20, 5), columns=None, errors=None, annotations_mode='track', seqbar=True, region='all', plot_kwargs=None)
    Plots multiple per-nucleotide datasets on a single axis.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    profile : data keyword string or data object
        Profile from which values will be plotted
    sequence : data keyword str, data obj, or sequence str, defaults to `profile`
        All data are mapped to this sequence before plotting
        If a data keyword, data from the first sample will be used
    annotations : list of data keyword strings or data objects, defaults to []
        Annotations used to highlight regions or sites of interest
    domains : data keyword string or data object, defaults to None
        domains to label along x-axis
    labels : list of str, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
    nt_ticks : tuple of two integers, defaults to (20, 5)
        first integer is the gap between major tick marks
        second integer is the gap between minor tick marks
    columns : string or list of strings, defaults to profile.metric
        columns names of values from profile to plot
    errors : string or list of strings, defaults to None (no error bars)
        column names of error values for plotting error bars
    annotations_mode : "track" or "bars", defaults to "track"
        "track" will highlight annotations along the x-axis
        "bars" will use a vertical transparent bar over the plot
    seqbar : bool, defaults to ``True``
        whether to display the sequence along the x-axis
    region : list of 2 integers, defaults to [1, length of sequence]
        start and end positions to plot. 1-indexed, inclusive.
    plot_kwargs : dictionary, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.Skyline
        the skyline plot object
```

## rnavigate.plot_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_ss in module rnavigate.plotting_functions

plot_ss(samples, structure, profile=None, annotations=None, interactions=None, interactions2=None, labels=None, colors=None, nt_ticks=None, bp_style='dotted', colorbars=True, plot_kwargs=None)
    Generates a multipanel secondary structure drawing with optional
    coloring by per-nucleotide data and display of inter-nucleotide data and/or
    sequence annotations. Each plot may display a unique sample and/or
    inter-nucleotide data filtering scheme.
    
    Parameters
    ----------
    samples : list of rnavigate Samples
        samples used to retrieve data
    structure : data keyword string or data object
        secondary structure to plot as arcs
        All data are mapped to this sequence before plotting
    profile : data keyword string or data object, defaults to None
        Profile used for coloring if "profile" used in colors dictionary
    annotations : list of data keyword strings or data objects, defaults to []
        Annotations used to highlight regions or sites of interest
    interactions : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot on secondary structure, no filtering
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
        format 3 (list of format 2 dictionaries)
            This format allows multiple filtering schemes to be applied,
            each will be plotted on a seperate axis
    interactions2 : one of the formats below, defaults to None
        format 1 (data or data keyword)
            Interactions to plot on secondary structure, no filtering
        format 2 (dictionary)
            e.g. {"interactions": format 1}
            additional filtering options can be added to the dictionary
    labels : list of strings, defaults to sample.sample for each sample
        list containing Labels to be used in plot legends
        Defaults to sample.sample for each sample
    colors : dictionary, optional
        a dictionary of element: value pairs that determines how colors
        will be applied to each plot element and if that element is plotted
        only the elements you wish to change need to be included
        value options and what the colors represent:
            None: don"t plot this elelement
            "sequence": nucleotide identity
            "position": position in sequence
            "annotations": sequence annotations
            "profile": per-nucleotide data from profile
                profile argument must be provided
            "structure": base-pairing status
            matplotlib color: all positions plotted in this color
            array of colors: a different color for each position
                must be the same length as structure
        "sequence" may also use "contrast" which automatically chooses
            white or black, which ever contrasts better with "nucleotide"
            color
        Defaults to {"sequence": None,
                        "nucleotides": "sequence",
                        "structure": "grey",
                        "basepairs": "grey"}
    nt_ticks : integer, defaults to None (no labels)
        gap between major tick marks
    bp_style : "dotted", "line", or "conventional", defaults to "dotted"
        "dotted" plots basepairs as a dotted line
        "line" plots basepairs as a solid line
        "conventional" plots basepairs using Leontis-Westhof conventions
            for canonical and wobble pairs ("G-A" plotted as solid dot)
    colorbars : bool, defaults to True
        Whether to plot color scales for all plot elements
    plot_kwargs : dict, defaults to {}
        Keyword-arguments passed to matplotlib.pyplot.subplots
    
    Returns
    -------
    rnavigate.plots.SS plot
        object containing matplotlib figure and axes with additional plotting and
        file saving methods
```

## rnavigate.analysis

### rnavigate.analysis.LogCompare

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class LogCompare in module rnavigate.analysis.logcompare

class LogCompare(builtins.object)
 |  LogCompare(samples1, samples2, name1, name2, data='profile', region='all')
 |  
 |  Compares 2 experimental samples, given replicates of each sample.
 |  
 |  Algorithm
 |  ---------
 |  1. Calculate the log10(modified/untreated) rate for each replicate.
 |  2. Scale these values to minimize the median of the absolute difference
 |  between samples.
 |  3. Calculate the standard error in these values for each replicate.
 |  4. Calculate the difference between samples.
 |  5. Calculate z-scores between samples.
 |  6. Plot the results in two panels:
 |      (1) the scaled log10(modified/untreated) rate for each sample with
 |      error bars, and
 |      (2) the difference between samples, colored by z-score.
 |  
 |  Methods
 |  -------
 |  
 |      __init__: computes log10(modified/untreated) rates, rescales the data,
 |          then calls make_plot()
 |      get_profile_sequence: gets log10(m/u) rate and sequence from sample
 |      calc_scale_factor: gets scale factor given two profiles
 |      rescale: rescales a profile to minimize difference to another profile
 |      load_replicates: calculates average and standard error of replicates
 |      make_plots: displays the two panels described above.
 |  
 |  Attributes:
 |      data (str): a key of sample.data to retrieve per-nucleotide data
 |      groups (dict): a dictionary containing the following key-value pairs:
 |          1: a dictionary containing these key-value pairs:
 |              self.data: averaged scaled log10(m/u) across replicates
 |              "stderr": the standard errors across replicates
 |              "stacked": 2d array containing each scaled log10(m/u) array
 |              "seq": the sequence string
 |          2: same as 1 above, for the second sample
 |  
 |  Methods defined here:
 |  
 |  __init__(self, samples1, samples2, name1, name2, data='profile', region='all')
 |      Takes replicates of two samples for comparison. Replicates are
 |      required. Calculates the log division profile
 |      (log10(modified/untreated)) and minimizes the median of the absolute
 |      difference between these. Finally, creates a plot.
 |      
 |      Args:
 |          samples1 (list of Sample objects): Replicates of the first sample
 |          samples2 (list of Sample objects): Replicates of the second sample
 |          name1 (string): name of first sample
 |          name2 (string): name of second sample
 |          data (str, optional): Datatype to compare. Defaults to "profile".
 |  
 |  calc_scale_factor(self, profile, target_profile)
 |      Finds the scale factor (x) that minimizes the median value of
 |      |profile + x - target_profile|.
 |      
 |      Args:
 |          profile (np.array): log10 profile
 |          target_profile (np.array): 2nd log10 profile
 |  
 |  get_profile_sequence(self, sample)
 |      retrieves log10(Modified_rate/Untreated_rate) and sequence from
 |      sample.data[self.data].
 |      
 |      Args:
 |          sample (rnavigate.Sample): an rnavigate sample
 |      
 |      Returns:
 |          np.array, np.array: log10 profile and sequence
 |  
 |  load_replicates(self, *samples, group)
 |      calculates average and standard error for a group of replicates,
 |      stores these values, along with sequence and a np.array containing
 |      all profiles in self.groups[group].
 |      
 |      Args:
 |          *samples (list of rnavigate.Sample): replicates to load
 |          group (int): self.groups key to access replicate data
 |  
 |  make_plots(self)
 |      Visualize this analysis.
 |  
 |  rescale(self, profile, target_profile)
 |      scales profile to minimize difference to target_profile using
 |      calc_scale_factor
 |      
 |      Args:
 |          profile (np.array): log10 profile to scale
 |          target_profile (np.array): 2nd log10 profile
 |      
 |      Returns:
 |          np.array: scaled profile
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.LowSS

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class LowSS in module rnavigate.analysis.lowss

class LowSS(rnavigate.rnavigate.Sample)
 |  LowSS(sample, window=55, shapemap='shapemap', pairprob='pairprob', structure='ss')
 |  
 |  Creates a new RNAvigate Sample which computes and displays Low SHAPE,
 |  low Shannon entropy regions (LowSS) given a sample containing SHAPE
 |  reactivities, pairing probabilities, and MFE structure.
 |  
 |  Methods
 |  -------
 |  __init__: performs the analysis
 |  plot_lowss: displays the result and returns plot object
 |  
 |  Attributes
 |  ----------
 |  sample : str
 |      the new label for this Sample's data on plots
 |  parent : rnavigate.Sample
 |      the sample from which data is retrieved
 |  window : int
 |      size of the windows, must be odd
 |  median_shape : float
 |      global median SHAPE reactivity
 |  median_entropy : float
 |      global median Shannon entropy
 |  data : dictionary
 |      dictionary of data keyword: Data objects, keys are:
 |          "structure" (rnav.data.SecondaryStructure)
 |              copy of provided MFE structure
 |          "shapemap" (rnav.data.SHAPEMaP)
 |              copy of provided SHAPE-MaP data aligned to "structure"
 |          "pairprob" (rnav.data.PairingProbability)
 |              copy of pairing probabilities aligned to "structure"
 |          "entropies" (rnav.data.Profile)
 |              Profile of Shannon entropies calculated from "pairprob"
 |          "lowSS" (rnav.data.Annotations)
 |              annotations defining low SHAPE, low Shannon entropy regions
 |  
 |  Method resolution order:
 |      LowSS
 |      rnavigate.rnavigate.Sample
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sample, window=55, shapemap='shapemap', pairprob='pairprob', structure='ss')
 |      Perform Low SHAPE, low Shannon entropy analysis on the given sample.
 |      
 |      Sample must contain:
 |          1) reactivities
 |          2) MFE structure
 |          3) pairing probabilities.
 |      
 |      Parameters
 |      ----------
 |      sample : rnavigate.Sample
 |          sample with shapemap, pairing probabilities and MFE structure
 |      window : int, default=55
 |          window size for calculating median SHAPE and Shannon entropy
 |      shapemap : str, default="shapemap"
 |          data keyword that points to SHAPE-MaP data
 |      pairprob : str, default="pairprob"
 |          data keyword that points to pairing probabilities data
 |      structure : str, default="ss"
 |          data keyword that points to MFE structure data
 |  
 |  plot_lowss(self, region=None, colorbars=True)
 |      Visualize LowSS analysis over the given region.
 |      
 |      Parameters
 |      ----------
 |      region : integer or list of 2 integers, default=None (entire sequence)
 |          If list: lowSS start and end positions to plot.
 |          If integer: region number, +/- 150 nts are shown.
 |      colorbars : bool, default=True
 |          whether to plot colorbars for pairing probability
 |      
 |      Returns
 |      -------
 |      rnavigate.plots.AP
 |          LowSS visualization
 |  
 |  reset_lowss(self, maximum_shape=None, maximum_entropy=0.08)
 |      Generates an annotation of lowSS regions. Stored as self.lowSS
 |      
 |      Parameters
 |      ----------
 |      maximum_shape : float, default=None (median SHAPE reactivity)
 |          maximum normalized SHAPE reactivity to be called lowSS.
 |      maximum_entropy : float, default=0.08
 |          maximum shannon entropy to be called lowSS.
 |  
 |  reset_window(self, window=None)
 |      Resets the window size and recalculates windowed SHAPE reactivities
 |      and shannon entropies and lowSS region annotations.
 |      
 |      Parameters
 |      ----------
 |      window : int, default=None (self.window)
 |          window size for calculating median SHAPE and Shannon entropy, must be odd
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.rnavigate.Sample:
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and applies filters to interactions data.
 |      
 |      Parameters
 |      ----------
 |      interactions : rnavigate.data.Interactions or data keyword string
 |          Interactions object to be filtered. If a string, value is
 |          replaced with self.get_data(interactions)
 |      metric : str, optional
 |          column of interactions data to be used as metric for coloring
 |          interactions.
 |          "Distance" will compute 3D distance in "pdb", defaulting to
 |          2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |          those atoms to compute distance.
 |      cmap (str | list, optional):
 |          sets the interactions colormap, used to color interactions
 |          according to metric values.
 |      normalization (str, optional):
 |          `"norm"`: extreme values in colormap are given to the extreme
 |              values of interactions metric data
 |          `"bins"`: data are colored according to which bin they fall in
 |              `values` defines bins (list, length = 2 less than cmap)
 |          `"min_max"`: extreme values in cmap are given to values beyond
 |              minimum and maximum, defined by `values`
 |      values:
 |          behavior depends on normalization
 |          `"norm"`: values are not needed
 |          `"bins"`: list of floats containing the boundaries between bins
 |              One fewer than the number of categories
 |          `"min_max"`: list of floats containing the minimum and maximum
 |      **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : rnavigate.data.Data or data keyword or list/dict of these
 |          If None, returns None.
 |          If a data keyword, returns associated data from sample
 |          If Data, returns that data.
 |          If a list or dictionary, returns list or dictionary with
 |              data keyword values replaced with associated Data
 |      data_class : rnavigate.data.Data class or subclass, optional
 |          If provided, ensures that returned data is of this type.
 |      
 |      Returns
 |      -------
 |      Same type as data_keyword argument, but data keywords are replaced
 |          with associated data
 |      
 |      Raises
 |      ------
 |      ValueError:
 |          if data is not found in sample
 |      ValueError:
 |          if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Parameters
 |      ----------
 |      inherit : Sample or list of Samples
 |          Other samples from which to inherit data and data keywords
 |      keep_inherited_defaults : bool
 |          Use default values from inherited samples
 |      overwrite : bool
 |          whether to overwrite any existing keywords with inherited keywords
 |  
 |  print_data_keywords(self, return_dict=False)
 |      Print a nicely formatted, organized list of data keywords.
 |      
 |      Returns a dictionary of data keywords, organized by data type, if
 |      return_dict is True.
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          The data keyword to set as the default
 |      overwrite : bool, defaults to ``True``
 |          whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name",
 |              data_keyword=inputs)
 |      
 |      is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name")
 |          my_sample.add_data(
 |              "data_keyword", inputs)
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          a data keyword used to store and/or parse the inputs
 |      inputs : dict or rnavigate.data.Data
 |          a dictionary used to create the data object or a data object itself
 |      overwrite : bool, defaults to False
 |          whether to overwrite a pre-existing data_keyword
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.rnavigate.Sample:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.DeltaSHAPE

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class DeltaSHAPE in module rnavigate.analysis.deltashape

class DeltaSHAPE(rnavigate.rnavigate.Sample)
 |  DeltaSHAPE(sample1, sample2, profile='shapemap', smoothing_window=3, zf_coeff=1.96, ss_thresh=1, site_window=3, site_nts=2)
 |  
 |  Detects meaningful differences in chemical probing reactivity
 |  
 |  References
 |  ----------
 |  doi:10.1021/acs.biochem.5b00977
 |  
 |  Algorithm
 |  ---------
 |  1. Extract SHAPE-MaP sequence, normalized profile, and normalized
 |      standard error from given samples
 |  2. Calculated smoothed profiles (mean) and propagate standard errors
 |      over rolling windows
 |  3. Subtract raw and smoothed normalized profiles and propogate errors
 |  4. Calculate Z-factors for smoothed data. This is the magnitude of the
 |      difference relative to the standard error
 |  5. Calculate Z-scores for smoothed data. This is the magnitude of the
 |      difference in standard deviations from the mean difference
 |  6. Call sites. Called sites must have # nucleotides that pass Z-factor
 |      and Z-score thresholds per window.
 |  
 |  Smoothing window size, Z factor threshold, Z score threshold, site-calling
 |  window size and minimum nucleotides per site can be specified.
 |  
 |  Method resolution order:
 |      DeltaSHAPE
 |      rnavigate.rnavigate.Sample
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sample1, sample2, profile='shapemap', smoothing_window=3, zf_coeff=1.96, ss_thresh=1, site_window=3, site_nts=2)
 |      Performs DeltaSHAPE analysis between samples 1 and 2
 |      
 |      Parameters
 |      ----------
 |      sample1 : rnavigate.Sample
 |          First sample to compare
 |      sample2 : rnavigate.Sample
 |          Second sample to compare
 |      profile : str, default="shapemap"
 |          Data keyword pointing to SHAPE-MaP data in samples 1 and 2
 |      smoothing_window : int, default=3
 |          Size of windows for data smoothing
 |      zf_coeff : float, default=1.96
 |          Sites must have a difference more than zf_coeff standard errors
 |      ss_thresh : int, default=1
 |          Sites must have a difference that is ss_thresh standard
 |          deviations from the mean difference
 |      site_window : int, default=3
 |          Number of nucleotides to include when calling sites
 |      site_nts : int, default=2
 |          Number of nts within site_window that must pass thresholds
 |  
 |  calculate_deltashape(self, smoothing_window=3, zf_coeff=1.96, ss_thresh=1, site_window=2, site_nts=3)
 |      Calculate or recalculate deltaSHAPE profile and called sites
 |      
 |      Parameters
 |      ----------
 |      smoothing_window : int, default=3
 |          Size of windows for data smoothing
 |      zf_coeff : float, default=1.96
 |          Sites must have a difference more than zf_coeff standard errors
 |      ss_thresh : int, default=1
 |          Sites must have a difference that is ss_thresh standard
 |          deviations from the mean difference
 |      site_window : int, default=3
 |          Number of nucleotides to include when calling sites
 |      site_nts : int, default=2
 |          Number of nts within site_window that must pass thresholds
 |  
 |  plot(self, region='all')
 |      Plot the deltaSHAPE result
 |      
 |      Parameters
 |      ----------
 |      region : list of 2 integers, default="all"
 |          start and end positions to plot
 |      
 |      Returns
 |      -------
 |      rnav.plots.Profile
 |          The plot object
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.rnavigate.Sample:
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and applies filters to interactions data.
 |      
 |      Parameters
 |      ----------
 |      interactions : rnavigate.data.Interactions or data keyword string
 |          Interactions object to be filtered. If a string, value is
 |          replaced with self.get_data(interactions)
 |      metric : str, optional
 |          column of interactions data to be used as metric for coloring
 |          interactions.
 |          "Distance" will compute 3D distance in "pdb", defaulting to
 |          2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |          those atoms to compute distance.
 |      cmap (str | list, optional):
 |          sets the interactions colormap, used to color interactions
 |          according to metric values.
 |      normalization (str, optional):
 |          `"norm"`: extreme values in colormap are given to the extreme
 |              values of interactions metric data
 |          `"bins"`: data are colored according to which bin they fall in
 |              `values` defines bins (list, length = 2 less than cmap)
 |          `"min_max"`: extreme values in cmap are given to values beyond
 |              minimum and maximum, defined by `values`
 |      values:
 |          behavior depends on normalization
 |          `"norm"`: values are not needed
 |          `"bins"`: list of floats containing the boundaries between bins
 |              One fewer than the number of categories
 |          `"min_max"`: list of floats containing the minimum and maximum
 |      **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : rnavigate.data.Data or data keyword or list/dict of these
 |          If None, returns None.
 |          If a data keyword, returns associated data from sample
 |          If Data, returns that data.
 |          If a list or dictionary, returns list or dictionary with
 |              data keyword values replaced with associated Data
 |      data_class : rnavigate.data.Data class or subclass, optional
 |          If provided, ensures that returned data is of this type.
 |      
 |      Returns
 |      -------
 |      Same type as data_keyword argument, but data keywords are replaced
 |          with associated data
 |      
 |      Raises
 |      ------
 |      ValueError:
 |          if data is not found in sample
 |      ValueError:
 |          if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Parameters
 |      ----------
 |      inherit : Sample or list of Samples
 |          Other samples from which to inherit data and data keywords
 |      keep_inherited_defaults : bool
 |          Use default values from inherited samples
 |      overwrite : bool
 |          whether to overwrite any existing keywords with inherited keywords
 |  
 |  print_data_keywords(self, return_dict=False)
 |      Print a nicely formatted, organized list of data keywords.
 |      
 |      Returns a dictionary of data keywords, organized by data type, if
 |      return_dict is True.
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          The data keyword to set as the default
 |      overwrite : bool, defaults to ``True``
 |          whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name",
 |              data_keyword=inputs)
 |      
 |      is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name")
 |          my_sample.add_data(
 |              "data_keyword", inputs)
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          a data keyword used to store and/or parse the inputs
 |      inputs : dict or rnavigate.data.Data
 |          a dictionary used to create the data object or a data object itself
 |      overwrite : bool, defaults to False
 |          whether to overwrite a pre-existing data_keyword
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.rnavigate.Sample:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.DeltaSHAPEProfile

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class DeltaSHAPEProfile in module rnavigate.analysis.deltashape

class DeltaSHAPEProfile(rnavigate.data.profile.Profile)
 |  DeltaSHAPEProfile(input_data, metric='Smooth_diff', metric_defaults=None, sequence=None, name=None, **kwargs)
 |  
 |  Profile data class for performing deltaSHAPE analysis
 |  
 |  Method resolution order:
 |      DeltaSHAPEProfile
 |      rnavigate.data.profile.Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, metric='Smooth_diff', metric_defaults=None, sequence=None, name=None, **kwargs)
 |      Create the deltaSHAPE Profile
 |      
 |      Parameters
 |      ----------
 |      input_data : tuple of rnavigate.Profile or pd.DataFrame
 |          If tuple of Profiles, the unified Dataframe will be created
 |      metric : str, default="Smooth_diff"
 |          The metric to use for the profile
 |      metric_defaults : dict, default=None
 |          Default settings for the metric
 |      sequence : str, default=None
 |          The sequence of the profile
 |      name : str, default=None
 |          The name of the profile
 |      **kwargs
 |          Additional keyword arguments to pass to data.Profile
 |  
 |  calculate_deltashape(self, smoothing_window=3, zf_coeff=1.96, ss_thresh=1, site_window=3, site_nts=2)
 |      Calculate the deltaSHAPE profile metrics
 |      
 |      Parameters
 |      ----------
 |      smoothing_window : int, default=3
 |          Size of windows for data smoothing
 |      zf_coeff : float, default=1.96
 |          Sites must have a difference more than zf_coeff standard errors
 |      ss_thresh : int, default=1
 |          Sites must have a difference that is ss_thresh standard
 |          deviations from the mean difference
 |      site_window : int, default=3
 |          Number of nucleotides to include when calling sites
 |      site_nts : int, default=2
 |          Number of nts within site_window that must pass thresholds
 |  
 |  get_enhancements_annotation(self)
 |      Get an annotations object for the significant enhancements
 |  
 |  get_protections_annotation(self)
 |      Get an annotations object for the significant protections
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.profile.Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from rnavigate.data.profile.Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.profile.Profile:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.WindowedAUROC

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class WindowedAUROC in module rnavigate.analysis.auroc

class WindowedAUROC(builtins.object)
 |  WindowedAUROC(sample, window=81, profile='default_profile', structure='default_structure')
 |  
 |  Compute and display windowed AUROC analysis.
 |  
 |  This analysis computes the ROC curve over a sliding window for the
 |  performance of per-nucleotide data (usually SHAPE-MaP or DMS-MaP Normalized
 |  reactivity) in predicting the base-pairing status of each nucleotide. The
 |  area under this curve (AUROC) is displayed compared to the median across
 |  the RNA. Below, an arc plot displays the secondary structure and
 |  per-nucleotide profile.
 |  
 |   AUROC values (should) range from 0.5 (no predictive power) to 1.0
 |  (perfect predictive power). A value of 0.5 indicates that the reactivity
 |  profile does not fit the structure prediction well. These regions are good
 |  candidates for further investigation with ensemble deconvolution.
 |  
 |  References
 |  ----------
 |  Lan, T.C.T., Allan, M.F., Malsick, L.E. et al. Secondary structural
 |      ensembles of the SARS-CoV-2 RNA genome in infected cells. Nat Commun
 |      13, 1128 (2022). https://doi.org/10.1038/s41467-022-28603-2
 |  
 |  Methods
 |  -------
 |  __init__: Computes the AUROC array and AUROC median.
 |  plot_auroc: Displays the AUROC analysis over the given region.
 |      Returns Plot object
 |  
 |  Attributes
 |  ----------
 |  sample : rnavigate.Sample
 |      sample to retrieve profile and secondary structure
 |  structure : str
 |      Data keyword of sample pointing to secondary structure
 |      e.g. sample.data[structure]
 |  profile : str
 |      Data keyword of sample pointing to profile
 |      e.g. sample.data[profile]
 |  sequence : the sequence string of sample.data[structure]
 |  window: the size of the windows
 |  nt_length: the length of sequence string
 |  auroc: the auroc numpy array, length = nt_length, padded with np.nan
 |  median_auroc: the median of the auroc array
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sample, window=81, profile='default_profile', structure='default_structure')
 |      Compute the AUROC for all windows. AUROC is a measure of how well a
 |      reactivity profile predicts paired vs. unpaired nucleotide status.
 |      
 |      Parameters
 |      ----------
 |      sample : rnav.Sample
 |          Your rnavigate sample
 |      window : int, optional
 |          number of nucleotides to include in window
 |          Defaults to 81.
 |      profile (str, optional): data keyword of provided sample pointing
 |          to a profile.
 |          Defaults to "default_profile"
 |      structure (str, optional): data keyword of provided sample pointing
 |          to a secondary structure.
 |          Defaults to "default_structure"
 |  
 |  plot_auroc(self, region=None)
 |      Plot the result of the windowed AUROC analysis, with arc plot of
 |      structure and reactivity profile.
 |      
 |      Args:
 |          region (list of int: length 2, optional): Start and end nucleotide
 |              positions to plot. Defaults to [1, RNA length].
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.FragMaP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class FragMaP in module rnavigate.analysis.fragmapper

class FragMaP(rnavigate.data.profile.Profile)
 |  FragMaP(input_data, parameters, metric='Fragmap_profile', metric_defaults=None, read_table_kw=None, sequence=None, name=None)
 |  
 |  Method resolution order:
 |      FragMaP
 |      rnavigate.data.profile.Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, parameters, metric='Fragmap_profile', metric_defaults=None, read_table_kw=None, sequence=None, name=None)
 |      Initialize the Profile object.
 |  
 |  calc_zscore(self, valid, dataframe, incolumn: str, outcolumn: str, base: list) -> None
 |  
 |  get_annotation(self)
 |  
 |  get_dataframe(self, profile1, profile2, mutation_rate_threshold, depth_threshold, p_significant, ss_threshold, correction_method)
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.profile.Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from rnavigate.data.profile.Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.Fragmapper

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Fragmapper in module rnavigate.analysis.fragmapper

class Fragmapper(rnavigate.rnavigate.Sample)
 |  Fragmapper(sample1, sample2, parameters=None, profile='shapemap')
 |  
 |  Method resolution order:
 |      Fragmapper
 |      rnavigate.rnavigate.Sample
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sample1, sample2, parameters=None, profile='shapemap')
 |      Creates a Sample.
 |  
 |  plot_scatter(self, column='Modified_rate')
 |      Generates scatter plots useful for fragmapper quality control.
 |      
 |      Args:
 |          column (str, optional):
 |              Dataframe column containing data to plot (must be avalible for
 |              the sample and control).
 |              Defaults to "Modified_rate".
 |      
 |      Returns:
 |          (matplotlib figure, matplotlib axis)
 |              Scatter plot with control values on the x-axis, sample values
 |              on the y-axis, and each point representing a nucleotide not
 |              filtered out in the fragmapper pipeline.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.rnavigate.Sample:
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and applies filters to interactions data.
 |      
 |      Parameters
 |      ----------
 |      interactions : rnavigate.data.Interactions or data keyword string
 |          Interactions object to be filtered. If a string, value is
 |          replaced with self.get_data(interactions)
 |      metric : str, optional
 |          column of interactions data to be used as metric for coloring
 |          interactions.
 |          "Distance" will compute 3D distance in "pdb", defaulting to
 |          2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |          those atoms to compute distance.
 |      cmap (str | list, optional):
 |          sets the interactions colormap, used to color interactions
 |          according to metric values.
 |      normalization (str, optional):
 |          `"norm"`: extreme values in colormap are given to the extreme
 |              values of interactions metric data
 |          `"bins"`: data are colored according to which bin they fall in
 |              `values` defines bins (list, length = 2 less than cmap)
 |          `"min_max"`: extreme values in cmap are given to values beyond
 |              minimum and maximum, defined by `values`
 |      values:
 |          behavior depends on normalization
 |          `"norm"`: values are not needed
 |          `"bins"`: list of floats containing the boundaries between bins
 |              One fewer than the number of categories
 |          `"min_max"`: list of floats containing the minimum and maximum
 |      **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : rnavigate.data.Data or data keyword or list/dict of these
 |          If None, returns None.
 |          If a data keyword, returns associated data from sample
 |          If Data, returns that data.
 |          If a list or dictionary, returns list or dictionary with
 |              data keyword values replaced with associated Data
 |      data_class : rnavigate.data.Data class or subclass, optional
 |          If provided, ensures that returned data is of this type.
 |      
 |      Returns
 |      -------
 |      Same type as data_keyword argument, but data keywords are replaced
 |          with associated data
 |      
 |      Raises
 |      ------
 |      ValueError:
 |          if data is not found in sample
 |      ValueError:
 |          if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Parameters
 |      ----------
 |      inherit : Sample or list of Samples
 |          Other samples from which to inherit data and data keywords
 |      keep_inherited_defaults : bool
 |          Use default values from inherited samples
 |      overwrite : bool
 |          whether to overwrite any existing keywords with inherited keywords
 |  
 |  print_data_keywords(self, return_dict=False)
 |      Print a nicely formatted, organized list of data keywords.
 |      
 |      Returns a dictionary of data keywords, organized by data type, if
 |      return_dict is True.
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          The data keyword to set as the default
 |      overwrite : bool, defaults to ``True``
 |          whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name",
 |              data_keyword=inputs)
 |      
 |      is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample="name")
 |          my_sample.add_data(
 |              "data_keyword", inputs)
 |      
 |      Parameters
 |      ----------
 |      data_keyword : str
 |          a data keyword used to store and/or parse the inputs
 |      inputs : dict or rnavigate.data.Data
 |          a dictionary used to create the data object or a data object itself
 |      overwrite : bool, defaults to False
 |          whether to overwrite a pre-existing data_keyword
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.rnavigate.Sample:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.analysis.SequenceChecker

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SequenceChecker in module rnavigate.analysis.check_sequence

class SequenceChecker(builtins.object)
 |  SequenceChecker(samples)
 |  
 |  Check the sequences stored in a list of samples.
 |  
 |  Attributes
 |  ----------
 |  samples : list
 |      samples in which to check sequences
 |  sequences : list
 |      all unique sequence strings stored in the list of samples. These are converted
 |      to an all uppercase RNA alphabet.
 |  keywords : list
 |      all unique data keywords stored in the list of samples.
 |  which_sequences : Pandas.DataFrame
 |      each row is a sample, keyword, and index of self.sequences
 |  
 |  Methods defined here:
 |  
 |  __init__(self, samples)
 |      Creates an instance of SequenceChecker given a list of samples
 |      
 |      Parameters
 |      ----------
 |      samples : list of rnav.Sample
 |          samples for which to compare data keywords and sequences.
 |  
 |  get_keywords(self)
 |      A list of all unique data keywords across samples.
 |  
 |  get_sequences(self)
 |      A list of all unique sequences (uppercase RNA) across samples.
 |  
 |  get_which_sequences(self)
 |      A DataFrame of sequence IDs (integers) for each data keyword.
 |  
 |  print_alignments(self, print_format='long', which='all')
 |      Print alignments in the given format for sequence IDs provided.
 |      
 |      Parameters
 |      ----------
 |      print_format : string, defaults to "long"
 |          What format to print the alignments in:
 |          "cigar" prints the cigar string
 |          "short" prints the numbers of mismatches and indels
 |          "long" prints the location and nucleotide identity of all
 |              mismatches, insertions and deletions.
 |      which : tuple of two of integers, defaults to "all" (every pairwise comparison)
 |          two sequence IDs to compare.
 |  
 |  print_mulitple_sequence_alignment(self, base_sequence)
 |      Print the multiple sequence alignment with nice formatting.
 |      
 |      Parameters
 |      ----------
 |      base_sequence : string
 |          a sequence string that represents the longest common sequence.
 |          Usually, this is the return value from:
 |              rnav.data.set_multiple_sequence_alignment()
 |  
 |  print_which_sequences(self)
 |      Print sequence ID (integer) for each data keyword and sample.
 |  
 |  reset(self)
 |      Reset keywords and sequences from sample list in case of changes.
 |  
 |  write_fasta(self, filename, which='all')
 |      Write all unique sequences to a fasta file.
 |      
 |      This is very useful for using external multiple sequence aligners such
 |      as ClustalOmega.
 |          1) go to https://www.ebi.ac.uk/Tools/msa/clustalo/
 |          2) upload new fasta file
 |          3) under STEP 2 output format, select Pearson/FASTA
 |          4) click 'Submit'
 |          5) wait for your alignment to finish
 |          6) download the alignment fasta file
 |          7) use rnav.data.set_multiple_sequence_alignment()
 |      
 |      Parameters
 |      ----------
 |      filename : string
 |          path to a new file to which fasta entries are written
 |      which : list of integers, defaults to "all" (every sequence)
 |          Sequence IDs to write to file.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

## rnavigate.data

### rnavigate.data.Sequence

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Sequence in module rnavigate.data.data

class Sequence(builtins.object)
 |  Sequence(input_data, name=None)
 |  
 |  A class for storing and manipulating RNA sequences.
 |  
 |  Parameters
 |  ----------
 |  sequence : string or pandas.DataFrame
 |      sequence string, fasta file, or a Pandas dataframe containing a
 |      "Sequence" column
 |  name : string, optional
 |      The name of the sequence, defaults to None
 |  
 |  Attributes
 |  ----------
 |  sequence : string
 |      The sequence string
 |  name : string
 |      The name of the sequence
 |  other_info : dict
 |      A dictionary of additional information about the sequence
 |  null_alignment : SequenceAlignment
 |      An alignment of the sequence to itself
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, name=None)
 |      Initialize the Sequence object.
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_aligned_data(self, alignment)
 |      Get a copy of the sequence positionally aligned to another sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          the alignment to use
 |      
 |      Returns
 |      -------
 |      aligned_sequence : rnavigate.data.Sequence
 |          the aligned sequence
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.Data

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Data in module rnavigate.data.data

class Data(Sequence)
 |  Data(input_data, sequence, metric, metric_defaults, read_table_kw=None, name=None)
 |  
 |  The base class for RNAvigate Profile and Interactions classes.
 |  
 |  Parameters
 |  ----------
 |  input_data : pandas.DataFrame or str
 |      a pandas dataframe or path to a data file
 |  sequence : string or rnavigate.data.Sequence
 |      the sequence to use for the data
 |  metric : string or dict
 |      the column of the dataframe to use as the default metric to visualize
 |  metric_defaults : dict
 |      a dictionary of metric defaults
 |  read_table_kw : dict, optional
 |      kwargs dictionary passed to pd.read_table
 |  name : string, optional
 |      the name of the data, defaults to None
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      the data table
 |  filepath : string
 |      the path to the data file
 |  sequence : string or rnavigate.data.Sequence
 |      the sequence to use for the data
 |  metric : string or dict
 |      the column of the dataframe to use as the metric to visualize
 |  metric_defaults : dict
 |      A dictionary of metric values and default settings for visualization
 |  default_metric : string
 |      the default metric to use for visualization
 |  
 |  Method resolution order:
 |      Data
 |      Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric, metric_defaults, read_table_kw=None, name=None)
 |      Initialize the Data object.
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_aligned_data(self, alignment)
 |      Get a copy of the sequence positionally aligned to another sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          the alignment to use
 |      
 |      Returns
 |      -------
 |      aligned_sequence : rnavigate.data.Sequence
 |          the aligned sequence
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.SecondaryStructure

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SecondaryStructure in module rnavigate.data.secondary_structure

class SecondaryStructure(rnavigate.data.data.Sequence)
 |  SecondaryStructure(input_data, extension=None, autoscale=True, name=None, **kwargs)
 |  
 |  Base class for secondary structures.
 |  
 |  Parameters
 |  ----------
 |  input_data : str or pandas.DataFrame
 |      A dataframe or filepath containing a secondary structure
 |      DataFrame should contain these columns:
 |          ["Nucleotide", "Sequence", "Pair"]
 |      "Pair" column must be redundant.
 |      Filepath parsing is determined by file extension:
 |          varna, xrna, nsd, cte, ct, dbn, bracket, json (R2DT), forna
 |  extension : str, optional
 |      The file extension of the input_data file. If not provided, the
 |      extension will be inferred from the input_data filepath.
 |  autoscale : bool, optional
 |      Whether to automatically scale the x and y coordinates. Defaults to True.
 |  name : str, optional
 |      The name of the RNA sequence. Defaults to None.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      DataFrame storing base-pairs
 |  filepath : str
 |      The path to the input file, if provided, otherwise "dataframe"
 |  sequence : str
 |      The RNA sequence
 |  nts : numpy.array
 |      The "Nucleotide" column of data
 |  pair_nts : numpy.array
 |      The "Pair" column of data
 |  header : str
 |      Header information from CT file
 |  xcoordinates : numpy.array
 |      The "X_coordinate" column of data
 |  ycoordinates : numpy.array
 |      The "X_coordinate" column of data
 |  distance_matrix : numpy.array
 |      The contact distance matrix of the RNA structure
 |  
 |  Method resolution order:
 |      SecondaryStructure
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, extension=None, autoscale=True, name=None, **kwargs)
 |      Creates a SecondaryStructure object from a given file or dataframe.
 |  
 |  __str__(self)
 |      print the filepath and length of the RNA
 |  
 |  add_pairs(self, pairs, break_conflicting_pairs=False)
 |      Add base pairs to current secondary structure.
 |      
 |      Parameters
 |      ----------
 |      pairs : list
 |          1-indexed list of paired residues. e.g. [(1, 20), (2, 19)]
 |      break_conflicting_pairs : bool, defaults to False
 |          Whether to break existing pairs if there is a conflict
 |  
 |  as_interactions(self, structure2=None)
 |      Returns rnavigate.Interactions representation of this, or more, structures.
 |      
 |      Parameters
 |      ----------
 |      structure2 : SecondaryStructure or list of these, defaults to None
 |          If provided, basepairs from all structures are included and labeled by
 |          which structures contain them and how many structures contain them.
 |  
 |  break_noncanonical_pairs(self)
 |      Removes non-canonical basepairs from the secondary structure.
 |      
 |      WARNING: this deletes information.
 |  
 |  break_pairs_nts(self, nt_positions)
 |      break base pairs at the given list of positions.
 |      
 |      WARNING: this deletes information.
 |      
 |      Parameters
 |      ----------
 |      nt_positions : list of int
 |          1-indexed positions to break pairs
 |  
 |  break_pairs_region(self, start, end, break_crossing=True, inverse=False)
 |      Removes pairs from the specified region (1-indexed, inclusive).
 |      
 |      WARNING: this deletes information
 |      
 |      Parameters
 |      ----------
 |      start : int
 |          start position (1-indexed, inclusive)
 |      end : int
 |          end position (1-indexed, inclusive)
 |      break_crossing : bool, defaults to True
 |          Whether to keep pairs that cross over the specified region
 |      inverse : bool, defaults to False
 |          Invert the behavior, i.e. remove pairs that are not in this region
 |  
 |  break_singleton_pairs(self)
 |      Removes singleton basepairs from the secondary structure.
 |      
 |      WARNING: This deletes information.
 |  
 |  compute_ppv_sens(self, structure2, exact=True)
 |      Compute the PPV and sensitivity between this and another structure.
 |      
 |      True and False are determined from this structure.
 |      Positive and Negative are determined from structure2.
 |      
 |      PPV = TP / (TP + FP)
 |      Sensitivity = TP / (TP + FN)
 |      
 |      Parameters
 |      ----------
 |      structure2 : SecondaryStructure
 |          The SecondaryStructure to compare to.
 |      exact : bool, defaults to True
 |          True requires BPs to be exactly correct.
 |          False allows +/-1 bp slippage.
 |      
 |      Returns
 |      -------
 |      float
 |          sensitivity
 |      float
 |          PPV
 |      2-tuple of floats
 |          (TP, TP+FP, TP+FN)
 |  
 |  contact_distance(self, i, j)
 |      Returns the contact distance between positions i and j
 |  
 |  copy(self)
 |  
 |  fill_mismatches(self, mismatch=1)
 |      Adds base pairs to fill 1,1 and optionally 2,2 mismatches.
 |      
 |      Parameters
 |      ----------
 |      mismatch : int, defaults to 1
 |          1 will fill only 1,1 mismatches
 |          2 will fill 1,1 and 2,2 mismatches
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new SecondaryStructure object matching the alignment target.
 |      
 |      Parameters
 |      ----------
 |      alignment : data.Alignment
 |          An alignment object used to map values
 |  
 |  get_distance_matrix(self, recalculate=False)
 |      Get a matrix of pair-wise shortest path distances through the structure.
 |      
 |      This function uses a BFS algorithm. The structure is represented as a complete
 |      graph with nucleotides as vertices and base-pairs and backbone as edges. All
 |      edges are length 1. Matrix is stored as an attribute for future use.
 |      
 |      If the attribute is set (not None) and recalculate is False, the attribute
 |      will be returned.
 |      
 |      Based on Tom's contact_distance, but expanded to return the pairwise matrix.
 |      New contact_distance method added to return the distance between two positions.
 |      
 |      Parameters
 |      ----------
 |      recalculate : bool, defaults to False
 |          Set to True to recalculate the matrix even if the attribute is set.
 |  
 |  get_dotbracket(self)
 |      Get a dotbracket notation string representing the secondary structure.
 |      
 |      Pseudoknot levels:
 |          1: ()
 |          2: []
 |          3: {}
 |          4: <>
 |          5: Aa
 |          6: Bb
 |          7: Cc
 |          etc...
 |      
 |      Returns
 |      -------
 |      str
 |          A dot-bracket representation of the secondary structure
 |  
 |  get_helices(self, fill_mismatches=True, split_bulge=True, keep_singles=False)
 |      Get a dictionary of helices from the secondary structure.
 |      
 |      Keys are equivalent to list indices. Values are lists of paired
 |      nucleotides (1-indexed) in that helix. e.g. {0:[(1,50),(2,49),(3,48)}
 |      
 |      Parameters
 |      ----------
 |      fill_mismatches : bool, defaults to True
 |          Whether 1-1 and 2-2 bulges are replaced with base pairs
 |      split_bulge : bool, defaults to True
 |          Whether to split helices on bulges
 |      keep_singles : bool, defaults to False
 |          Whether to return helices that contain only 1 base-pair
 |      
 |      Returns
 |      -------
 |      dict
 |          A dictionary of helices
 |  
 |  get_human_dotbracket(self)
 |      Get a human-readable dotbracket string representing the secondary structure.
 |      
 |      This is an experimental format designed to be more human readable, i.e. no
 |      counting of brackets required.
 |      
 |      1)  Letters, instead of brackets, are used to denote nested base pairs.
 |      2)  Each helix is assigned a letter, which is incremented one letter
 |          alphabetically from the nearest enclosing stem.
 |      3)  Non-nested helices (pseudoknots) are assigned canonical brackets.
 |      
 |      From this canonical dbn string:
 |          how many bases are in the base stem?
 |          how many nested helices are there?
 |          ((((....(((.[[..)))))(((...(((..]].))))))))
 |      Same question, new format:
 |          AABB....CCC.[[..cccbbBBB...CCC..]].cccbbbaa
 |      Read this as:
 |          ((_______________________________________)) (level 1 = A)
 |            ((_______________))(((______________)))   (level 2 = B)
 |                  (((_____)))        (((_____)))      (level 3 = C)
 |                      [[__________________]]          (pseudoknot = [])
 |      
 |      Pseudoknot levels:
 |          1: Aa, Bb, Cc, etc.
 |          2: [], 3: {}, 4: <>
 |  
 |  get_interactions_df(self)
 |      Returns a DataFrame of i, j basepairs.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A DataFrame with columns:
 |              i: the 5' (1-indexed) position of the base pair
 |              j: the 3' (1-indexed) position of the base pair
 |              Structure: always 1
 |  
 |  get_junction_nts(self)
 |      Get a list of junction nucleotides (paired, but at the end of a chain).
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions of junction nucleotides
 |  
 |  get_nonredundant_ct(self)
 |      Returns the ct attribute in a non-redundant form.
 |      
 |      Only returns pairs in which i < j
 |      For example:
 |          self.ct[i-1] == j
 |          self.ct[j-1] == i
 |          BUT
 |          self.get_nonredundant_ct()[j-1] == 0
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          A non-redundant array of base pairs
 |  
 |  get_paired_nts(self)
 |      Get a list of residues that are paired.
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions of paired nucleotides
 |  
 |  get_pairs(self)
 |      Get a non-redundant list of base pairs i < j as a array of tuples.
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions. e.g., [(1, 50), (2, 49), ...]
 |  
 |  get_pseudoknots(self, fill_mismatches=True)
 |      Get the pk1 and pk2 pairs from the secondary structure.
 |      
 |      Ignores single base pairs. PK1 is defined as the helix crossing the
 |      most other bps. If there is a tie, the most 5' helix is called pk1
 |      returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
 |      
 |      Parameters
 |      ----------
 |      fill_mismatches : bool, defaults to True
 |          Whether 1-1 and 2-2 bulges are replaced with base pairs
 |      
 |      Returns
 |      -------
 |      list of 2 lists of 2-tuples
 |          A list of base pairs for pk1 and pk2
 |  
 |  get_structure_elements(self)
 |      This code is not yet implemented.
 |      
 |      Returns a string with a character for each nucleotide, indicating
 |      what kind of structure element it is a part of.
 |      
 |      Characters:
 |          Dangling Ends (E)
 |          Stems (S)
 |          Hairpin Loops (H)
 |          Bulges (B)
 |          Internal Loops (I)
 |          MultiLoops (M)
 |          External Loops (X)
 |          Pseudoknot (P)
 |  
 |  get_unpaired_nts(self)
 |      Get a list of residues that are unpaired.
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions of unpaired nucleotides
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Normalize the sequence attribute (fix case and/or U <-> T).
 |  
 |  read_ct(self, structure_number=0)
 |      Loads secondary structure information from a given ct file.
 |      
 |      Requires a properly formatted header.
 |      
 |      Parameters
 |      ----------
 |      structure_number : int, defaults to 0
 |          0-indexed structure number to load from the ct file.
 |  
 |  read_cte(self)
 |      Generates SecondaryStructure object data from a CTE file
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_dotbracket(self)
 |      Generates SecondaryStructure object data from a dot-bracket file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_forna(self)
 |      Generates SecondaryStructure object data from a FORNA JSON file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_nsd(self)
 |      Generates SecondaryStructure object data from an NSD file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_r2dt(self)
 |      Generates SecondaryStructure object data from an R2DT JSON file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_varna(self)
 |      Generates SecondaryStructure object data from a VARNA file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_xrna(self)
 |      Generates SecondaryStructure object data from an XRNA file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  transform_coordinates(self, flip=None, scale=None, center=None, rotate_degrees=None)
 |      Perform transformations on X and Y structure coordinates.
 |      
 |      To acheive vertical and horizontal flip together, rotate 180 degrees.
 |      
 |      Parameters
 |      ----------
 |      flip : str, optional
 |          "horizontal" or "vertical"
 |      scale : float, optional
 |          new median distance of basepairs
 |      center : tuple of floats, optional
 |          new center x and y coordinate
 |      rotate_degrees : float, optional
 |          number of degrees to rotate structure
 |  
 |  write_ct(self, out_file)
 |      Write structure to a ct file.
 |  
 |  write_cte(self, out_file)
 |      Write structure to CTE format for Structure Editor.
 |  
 |  write_dbn(self, out_file, rna_name)
 |      Write structure to dot-bracket file.
 |  
 |  write_sto(self, out_file, name='seq')
 |      Write structure to Stockholm (STO) file to use in infernal searches.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_pairs_list(input_data, sequence) from builtins.type
 |      Creates a SecondaryStructure from a list of pairs and a sequence.
 |      
 |      Parameters
 |      ----------
 |      input_data : list
 |          1-indexed list of base pairs. e.g. [(1, 20), (2, 19)]
 |      sequence : str
 |          The RNA sequence. e.g., "AUCGUGUCAUGCUA"
 |  
 |  from_sequence(input_data) from builtins.type
 |      Creates a SecondaryStructure from a sequence string.
 |      
 |      This structure is initialized with no base pairs. If base pairs are
 |      needed, use SecondaryStructure.from_pairs_list().
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  nts
 |  
 |  pair_nts
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  xcoordinates
 |  
 |  ycoordinates
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.StructureCoordinates

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class StructureCoordinates in module rnavigate.data.secondary_structure

class StructureCoordinates(builtins.object)
 |  StructureCoordinates(x, y, pairs=None)
 |  
 |  Helper class to perform structure coordinate transformations
 |  
 |  Parameters
 |  ----------
 |  x : numpy.array
 |      x coordinates
 |  y : numpy.array
 |      y coordinates
 |  pairs : list of pairs, optional
 |      list of base-paired positions
 |      required if scaling coordinates
 |  
 |  Methods defined here:
 |  
 |  __init__(self, x, y, pairs=None)
 |      initialize structure coordinates.
 |  
 |  center(self, x=0, y=0)
 |      Center structure on the given x, y coordinate
 |      
 |      Parameters
 |      ----------
 |      x : int, defaults to 0
 |          x coordinate of structure center
 |      y : int, defaults to 0
 |          y coordinate of structure center
 |  
 |  flip(self, horizontal=True)
 |      Flip structure vertically or horizontally.
 |      
 |      Parameters
 |      ----------
 |      horizontal : bool, defaults to True
 |          whether to flip structure horizontally, otherwise vertically
 |  
 |  get_center_point(self)
 |      Get the x, y coordinates for the center of structure.
 |      
 |      Returns
 |      -------
 |      float
 |          x coordinate of structure center
 |      float
 |          y coordinate of structure center
 |  
 |  rotate(self, degrees)
 |      Rotate structure on current center point.
 |      
 |      Parameters
 |      ----------
 |      degrees : float
 |          number of degrees to rotate structure
 |  
 |  scale(self, median_bp_distance=1.0)
 |      Scale structure such that median base-pair distance is constant.
 |      
 |      Parameters
 |      ----------
 |      median_bp_distance : float, defaults to 1.0
 |          New median distance between all base-paired nucleotides.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.SequenceCircle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SequenceCircle in module rnavigate.data.secondary_structure

class SequenceCircle(SecondaryStructure)
 |  SequenceCircle(input_data, gap=30, name=None, **kwargs)
 |  
 |  A circular SecondaryStructure-like representation of RNA sequence.
 |  
 |  Method resolution order:
 |      SequenceCircle
 |      SecondaryStructure
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, gap=30, name=None, **kwargs)
 |      Creates a SecondaryStructure object from a given file or dataframe.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from SecondaryStructure:
 |  
 |  __str__(self)
 |      print the filepath and length of the RNA
 |  
 |  add_pairs(self, pairs, break_conflicting_pairs=False)
 |      Add base pairs to current secondary structure.
 |      
 |      Parameters
 |      ----------
 |      pairs : list
 |          1-indexed list of paired residues. e.g. [(1, 20), (2, 19)]
 |      break_conflicting_pairs : bool, defaults to False
 |          Whether to break existing pairs if there is a conflict
 |  
 |  as_interactions(self, structure2=None)
 |      Returns rnavigate.Interactions representation of this, or more, structures.
 |      
 |      Parameters
 |      ----------
 |      structure2 : SecondaryStructure or list of these, defaults to None
 |          If provided, basepairs from all structures are included and labeled by
 |          which structures contain them and how many structures contain them.
 |  
 |  break_noncanonical_pairs(self)
 |      Removes non-canonical basepairs from the secondary structure.
 |      
 |      WARNING: this deletes information.
 |  
 |  break_pairs_nts(self, nt_positions)
 |      break base pairs at the given list of positions.
 |      
 |      WARNING: this deletes information.
 |      
 |      Parameters
 |      ----------
 |      nt_positions : list of int
 |          1-indexed positions to break pairs
 |  
 |  break_pairs_region(self, start, end, break_crossing=True, inverse=False)
 |      Removes pairs from the specified region (1-indexed, inclusive).
 |      
 |      WARNING: this deletes information
 |      
 |      Parameters
 |      ----------
 |      start : int
 |          start position (1-indexed, inclusive)
 |      end : int
 |          end position (1-indexed, inclusive)
 |      break_crossing : bool, defaults to True
 |          Whether to keep pairs that cross over the specified region
 |      inverse : bool, defaults to False
 |          Invert the behavior, i.e. remove pairs that are not in this region
 |  
 |  break_singleton_pairs(self)
 |      Removes singleton basepairs from the secondary structure.
 |      
 |      WARNING: This deletes information.
 |  
 |  compute_ppv_sens(self, structure2, exact=True)
 |      Compute the PPV and sensitivity between this and another structure.
 |      
 |      True and False are determined from this structure.
 |      Positive and Negative are determined from structure2.
 |      
 |      PPV = TP / (TP + FP)
 |      Sensitivity = TP / (TP + FN)
 |      
 |      Parameters
 |      ----------
 |      structure2 : SecondaryStructure
 |          The SecondaryStructure to compare to.
 |      exact : bool, defaults to True
 |          True requires BPs to be exactly correct.
 |          False allows +/-1 bp slippage.
 |      
 |      Returns
 |      -------
 |      float
 |          sensitivity
 |      float
 |          PPV
 |      2-tuple of floats
 |          (TP, TP+FP, TP+FN)
 |  
 |  contact_distance(self, i, j)
 |      Returns the contact distance between positions i and j
 |  
 |  copy(self)
 |  
 |  fill_mismatches(self, mismatch=1)
 |      Adds base pairs to fill 1,1 and optionally 2,2 mismatches.
 |      
 |      Parameters
 |      ----------
 |      mismatch : int, defaults to 1
 |          1 will fill only 1,1 mismatches
 |          2 will fill 1,1 and 2,2 mismatches
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new SecondaryStructure object matching the alignment target.
 |      
 |      Parameters
 |      ----------
 |      alignment : data.Alignment
 |          An alignment object used to map values
 |  
 |  get_distance_matrix(self, recalculate=False)
 |      Get a matrix of pair-wise shortest path distances through the structure.
 |      
 |      This function uses a BFS algorithm. The structure is represented as a complete
 |      graph with nucleotides as vertices and base-pairs and backbone as edges. All
 |      edges are length 1. Matrix is stored as an attribute for future use.
 |      
 |      If the attribute is set (not None) and recalculate is False, the attribute
 |      will be returned.
 |      
 |      Based on Tom's contact_distance, but expanded to return the pairwise matrix.
 |      New contact_distance method added to return the distance between two positions.
 |      
 |      Parameters
 |      ----------
 |      recalculate : bool, defaults to False
 |          Set to True to recalculate the matrix even if the attribute is set.
 |  
 |  get_dotbracket(self)
 |      Get a dotbracket notation string representing the secondary structure.
 |      
 |      Pseudoknot levels:
 |          1: ()
 |          2: []
 |          3: {}
 |          4: <>
 |          5: Aa
 |          6: Bb
 |          7: Cc
 |          etc...
 |      
 |      Returns
 |      -------
 |      str
 |          A dot-bracket representation of the secondary structure
 |  
 |  get_helices(self, fill_mismatches=True, split_bulge=True, keep_singles=False)
 |      Get a dictionary of helices from the secondary structure.
 |      
 |      Keys are equivalent to list indices. Values are lists of paired
 |      nucleotides (1-indexed) in that helix. e.g. {0:[(1,50),(2,49),(3,48)}
 |      
 |      Parameters
 |      ----------
 |      fill_mismatches : bool, defaults to True
 |          Whether 1-1 and 2-2 bulges are replaced with base pairs
 |      split_bulge : bool, defaults to True
 |          Whether to split helices on bulges
 |      keep_singles : bool, defaults to False
 |          Whether to return helices that contain only 1 base-pair
 |      
 |      Returns
 |      -------
 |      dict
 |          A dictionary of helices
 |  
 |  get_human_dotbracket(self)
 |      Get a human-readable dotbracket string representing the secondary structure.
 |      
 |      This is an experimental format designed to be more human readable, i.e. no
 |      counting of brackets required.
 |      
 |      1)  Letters, instead of brackets, are used to denote nested base pairs.
 |      2)  Each helix is assigned a letter, which is incremented one letter
 |          alphabetically from the nearest enclosing stem.
 |      3)  Non-nested helices (pseudoknots) are assigned canonical brackets.
 |      
 |      From this canonical dbn string:
 |          how many bases are in the base stem?
 |          how many nested helices are there?
 |          ((((....(((.[[..)))))(((...(((..]].))))))))
 |      Same question, new format:
 |          AABB....CCC.[[..cccbbBBB...CCC..]].cccbbbaa
 |      Read this as:
 |          ((_______________________________________)) (level 1 = A)
 |            ((_______________))(((______________)))   (level 2 = B)
 |                  (((_____)))        (((_____)))      (level 3 = C)
 |                      [[__________________]]          (pseudoknot = [])
 |      
 |      Pseudoknot levels:
 |          1: Aa, Bb, Cc, etc.
 |          2: [], 3: {}, 4: <>
 |  
 |  get_interactions_df(self)
 |      Returns a DataFrame of i, j basepairs.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A DataFrame with columns:
 |              i: the 5' (1-indexed) position of the base pair
 |              j: the 3' (1-indexed) position of the base pair
 |              Structure: always 1
 |  
 |  get_junction_nts(self)
 |      Get a list of junction nucleotides (paired, but at the end of a chain).
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions of junction nucleotides
 |  
 |  get_nonredundant_ct(self)
 |      Returns the ct attribute in a non-redundant form.
 |      
 |      Only returns pairs in which i < j
 |      For example:
 |          self.ct[i-1] == j
 |          self.ct[j-1] == i
 |          BUT
 |          self.get_nonredundant_ct()[j-1] == 0
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          A non-redundant array of base pairs
 |  
 |  get_paired_nts(self)
 |      Get a list of residues that are paired.
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions of paired nucleotides
 |  
 |  get_pairs(self)
 |      Get a non-redundant list of base pairs i < j as a array of tuples.
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions. e.g., [(1, 50), (2, 49), ...]
 |  
 |  get_pseudoknots(self, fill_mismatches=True)
 |      Get the pk1 and pk2 pairs from the secondary structure.
 |      
 |      Ignores single base pairs. PK1 is defined as the helix crossing the
 |      most other bps. If there is a tie, the most 5' helix is called pk1
 |      returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
 |      
 |      Parameters
 |      ----------
 |      fill_mismatches : bool, defaults to True
 |          Whether 1-1 and 2-2 bulges are replaced with base pairs
 |      
 |      Returns
 |      -------
 |      list of 2 lists of 2-tuples
 |          A list of base pairs for pk1 and pk2
 |  
 |  get_structure_elements(self)
 |      This code is not yet implemented.
 |      
 |      Returns a string with a character for each nucleotide, indicating
 |      what kind of structure element it is a part of.
 |      
 |      Characters:
 |          Dangling Ends (E)
 |          Stems (S)
 |          Hairpin Loops (H)
 |          Bulges (B)
 |          Internal Loops (I)
 |          MultiLoops (M)
 |          External Loops (X)
 |          Pseudoknot (P)
 |  
 |  get_unpaired_nts(self)
 |      Get a list of residues that are unpaired.
 |      
 |      Returns
 |      -------
 |      list
 |          A list of 1-indexed positions of unpaired nucleotides
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Normalize the sequence attribute (fix case and/or U <-> T).
 |  
 |  read_ct(self, structure_number=0)
 |      Loads secondary structure information from a given ct file.
 |      
 |      Requires a properly formatted header.
 |      
 |      Parameters
 |      ----------
 |      structure_number : int, defaults to 0
 |          0-indexed structure number to load from the ct file.
 |  
 |  read_cte(self)
 |      Generates SecondaryStructure object data from a CTE file
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_dotbracket(self)
 |      Generates SecondaryStructure object data from a dot-bracket file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_forna(self)
 |      Generates SecondaryStructure object data from a FORNA JSON file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_nsd(self)
 |      Generates SecondaryStructure object data from an NSD file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_r2dt(self)
 |      Generates SecondaryStructure object data from an R2DT JSON file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_varna(self)
 |      Generates SecondaryStructure object data from a VARNA file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  read_xrna(self)
 |      Generates SecondaryStructure object data from an XRNA file.
 |      
 |      Resulting SecondaryStructure object will include nucleotide x and y
 |      coordinates and is compatible with plot_ss.
 |  
 |  transform_coordinates(self, flip=None, scale=None, center=None, rotate_degrees=None)
 |      Perform transformations on X and Y structure coordinates.
 |      
 |      To acheive vertical and horizontal flip together, rotate 180 degrees.
 |      
 |      Parameters
 |      ----------
 |      flip : str, optional
 |          "horizontal" or "vertical"
 |      scale : float, optional
 |          new median distance of basepairs
 |      center : tuple of floats, optional
 |          new center x and y coordinate
 |      rotate_degrees : float, optional
 |          number of degrees to rotate structure
 |  
 |  write_ct(self, out_file)
 |      Write structure to a ct file.
 |  
 |  write_cte(self, out_file)
 |      Write structure to CTE format for Structure Editor.
 |  
 |  write_dbn(self, out_file, rna_name)
 |      Write structure to dot-bracket file.
 |  
 |  write_sto(self, out_file, name='seq')
 |      Write structure to Stockholm (STO) file to use in infernal searches.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from SecondaryStructure:
 |  
 |  from_pairs_list(input_data, sequence) from builtins.type
 |      Creates a SecondaryStructure from a list of pairs and a sequence.
 |      
 |      Parameters
 |      ----------
 |      input_data : list
 |          1-indexed list of base pairs. e.g. [(1, 20), (2, 19)]
 |      sequence : str
 |          The RNA sequence. e.g., "AUCGUGUCAUGCUA"
 |  
 |  from_sequence(input_data) from builtins.type
 |      Creates a SecondaryStructure from a sequence string.
 |      
 |      This structure is initialized with no base pairs. If base pairs are
 |      needed, use SecondaryStructure.from_pairs_list().
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from SecondaryStructure:
 |  
 |  nts
 |  
 |  pair_nts
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from SecondaryStructure:
 |  
 |  xcoordinates
 |  
 |  ycoordinates
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.set_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function set_alignment in module rnavigate.data.alignments

set_alignment(sequence1, sequence2, alignment1, alignment2, t_or_u='U')
    Add an alignment to be used as the default between two sequences.
    
    When objects with these sequences are aligned for visualization, RNAvigate
    uses this alignment instead of an automated pairwise sequence alignment.
    Alignment 1 and 2 must have matching lengths.
    alignment(1,2) and sequence(1,2) must differ only by dashes "-".
    
    e.g.:
        sequence1 ="AAGCUUCGGUACAUGCAAGAUGUAC"
        sequence2 ="AUCGAUCGAGCUGCUGUGUACGUAC"
        alignment1="AAGCUUCG---------GUACAUGCAAGAUGUAC"
        alignment2="AUCGAUCGAGCUGCUGUGUAC---------GUAC"
                     |mm|   | indel |    | indel |
    
    Parameters
    ----------
    sequence1 : string
        the first sequence
    sequence2 : string
        the second sequence
    alignment1 : string
        first sequence, plus dashes "-" indicating indels
    alignment2 : string
        second sequence, plus dashes "-" indicating indels
    t_or_u : "T", "U", or False
        "T" converts "U"s to "T"s
```

### rnavigate.data.set_multiple_sequence_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function set_multiple_sequence_alignment in module rnavigate.data.alignments

set_multiple_sequence_alignment(fasta, set_pairwise=False)
    Set alignments from a multiple sequence alignment Pearson fasta file.
    
    Sets alignments to a base sequence, then returns the base sequence to be
    when a multiple sequence alignment plot is desired. Also sets all pairwise
    alignments, if desired. When setting pairwise alignments, dashes that are
    shared between pairwise sequences are removed first.
    
    Parameters
    ----------
    fasta : string
        location of Pearson fasta file
    set_pairwise : bool, defaults to False
        whether to set every pairwise alignment as well as the multiple
        sequence alignment.
```

### rnavigate.data.lookup_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function lookup_alignment in module rnavigate.data.alignments

lookup_alignment(sequence1, sequence2, t_or_u='U')
    look up a previously set alignment in the _alignments_cache
    
    Parameters
    ----------
    sequence1 : string
        The first sequence to align
    sequence2 : string
        The second sequence to be aligned to
    t_or_u : "T", "U", or False, defaults to "U"
        "T" converts "U"s to "T"s
        "U" converts "U"s to "T"s
        False does nothing
    
    Returns
    -------
    dictionary, if an alignment is found, otherwise None
        {"seqA": sequence1 with gap characters representing alignment,
         "seqB": sequence2 with gap characters representing alignment}
```

### rnavigate.data.SequenceAlignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SequenceAlignment in module rnavigate.data.alignments

class SequenceAlignment(BaseAlignment)
 |  SequenceAlignment(sequence1, sequence2, align_kwargs=None, full=False, use_previous=True)
 |  
 |  The most useful feature of RNAvigate. Maps positions from one sequence
 |  to a totally different sequence using user-defined pairwise alignment or
 |  automatic pairwise alignment.
 |  
 |  Parameters
 |  ----------
 |  sequence1 : string
 |      the sequence to be aligned
 |  sequence2 : string
 |      the sequence to align to
 |  align_kwargs : dict, defaults to None
 |      a dictionary of arguments to pass to pairwise2.align.globalms
 |  full : bool, defaults to False
 |      whether to keep unmapped starting sequence positions.
 |  use_previous : bool, defaults to True
 |      whether to use previously set alignments
 |  
 |  Attributes
 |  ----------
 |  sequence1 : str
 |      the sequence to be aligned
 |  sequence2 : str
 |      the sequence to align to
 |  alignment1 : str
 |      the alignment string matching sequence1 to sequence2
 |  alignment2 : str
 |      the alignment string matching sequence2 to sequence1
 |  starting_sequence : str
 |      sequence1
 |  target_sequence : str
 |      sequence2 if full is False, else alignment2
 |  mapping : numpy.array
 |      the alignment map array.
 |      index of starting_sequence is mapping[index] of target_sequence
 |  
 |  Method resolution order:
 |      SequenceAlignment
 |      BaseAlignment
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sequence1, sequence2, align_kwargs=None, full=False, use_previous=True)
 |      Creates an alignment from sequence1 to sequence2.
 |  
 |  __repr__(self)
 |      a nice text only representation of an alignment
 |  
 |  get_alignment(self)
 |      Gets an alignment that has either been user-defined or previously
 |      calculated or produces a new pairwise alignment between two sequences.
 |      
 |      Returns
 |      -------
 |      alignment1, alignment2 : tuple of 2 str
 |          the alignment strings matching sequence1 and sequence2, respectively.
 |  
 |  get_inverse_alignment(self)
 |      Gets an alignment that maps from sequence2 to sequence1.
 |  
 |  get_mapping(self)
 |      Calculates a mapping from starting sequence to target sequence.
 |      
 |      Returns
 |      -------
 |      mapping : numpy.array
 |          an array that maps to an index of target sequence.
 |          index of starting_sequence is mapping[index] of target_sequence
 |  
 |  print(self, print_format='full')
 |      Print the alignment in a human-readable format.
 |      
 |      Parameters
 |      ----------
 |      print_format : "full", "cigar", "long" or "short", defaults to "full"
 |          how to format the alignment.
 |          "full": the full length alignment with changes labeled "X"
 |          "cigar": the CIGAR string
 |          "long": locations and sequences of each change
 |          "short": total number of matches, mismatches, and indels
 |  
 |  print_all_changes(self)
 |      Print location and sequence of all changes.
 |  
 |  print_cigar(self)
 |      Print the CIGAR string
 |  
 |  print_number_of_changes(self)
 |      Print the total numbers of matches, mismatches, and indels.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from BaseAlignment:
 |  
 |  get_target_sequence(self)
 |      Gets the portion of starting sequence that fits the alignment
 |  
 |  map_dataframe(self, dataframe, position_columns)
 |      Takes a dataframe and maps position columns to target sequence.
 |      
 |      Rows with unmapped positions are dropped.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          a dataframe with position columns
 |      position_columns : list of str
 |          a list of columns containing positions to map
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a new dataframe (copy) with position columns mapped or dropped
 |  
 |  map_indices(self, indices, keep_minus_one=True)
 |      Takes a list of indices (0-index) and maps them to target sequence
 |      
 |      Parameters
 |      ----------
 |      indices : int or list of int
 |          a single or list of integer indices
 |      keep_minus_one : bool, defaults to True
 |          whether to keep unmapped starting sequence indices (-1) in the
 |          returned array.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          the equivalent indices in target sequence
 |  
 |  map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide', sequence_column='Sequence')
 |      Takes a per-nt dataframe and map it to the target sequence.
 |      
 |      Dataframe must have 1 row per nucleotide in starting sequence,
 |      with a position column and a sequence column. Dataframe is
 |      mapped to have the same format, but for target sequence nucleotides and
 |      positions.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          a per-nucleotide dataframe
 |      position_column : string, defaults to "Nucleotide"
 |          name of the position column.
 |      sequence_column : string, defaults to "Sequence"
 |          name of the sequence column.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a new dataframe (copy) mapped to target sequence.
 |          Unmapped starting sequence positions are dropped and unmapped
 |          target sequence positions are filled.
 |  
 |  map_positions(self, positions, keep_zero=True)
 |      Takes a list of positions (1-index) and maps them to target sequence
 |      
 |      Parameters
 |      ----------
 |      positions : int or list of int
 |          a single or list of integer positions
 |      keep_zero : bool, defaults to True
 |          whether to keep unmapped starting sequence positions (0) in the
 |          returned array.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          the equivalent positions in target sequence
 |  
 |  map_values(self, values, fill=nan)
 |      Takes an array of length equal to starting sequence and maps them to
 |      target sequence, unmapped positions in starting sequence are dropped
 |      and unmapped positions in target sequence are filled with fill value.
 |      
 |      Parameters
 |      ----------
 |      values : iterable
 |          values to map to target sequence.
 |      fill : any, defaults to np.nan
 |          a value for unmapped positions in target sequence.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          an array of values equal in length to target sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from BaseAlignment:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.AlignmentChain

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class AlignmentChain in module rnavigate.data.alignments

class AlignmentChain(BaseAlignment)
 |  AlignmentChain(*alignments)
 |  
 |  Combines a list of alignments into one.
 |  
 |  Parameters
 |  ----------
 |  alignments : list of Alignment objects
 |      the alignments to chain together
 |  
 |  Attributes
 |  ----------
 |  alignments : list
 |      the constituent alignments
 |  starting_sequence : str
 |      starting sequence of alignments[0]
 |  target_sequence : str
 |      target sequence of alignments[-1]
 |  mapping : numpy.array
 |      an array which maps from `starting_sequence` to `target_sequence`.
 |      index of starting_sequence is mapping[index] of target sequence
 |  
 |  Method resolution order:
 |      AlignmentChain
 |      BaseAlignment
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, *alignments)
 |      Creates a single alignment from multiple alignments.
 |  
 |  get_inverse_alignment(self)
 |      Alignments require a method to get the inverted alignment
 |  
 |  get_mapping(self)
 |      combines mappings from each alignment.
 |      
 |      Returns
 |      -------
 |      mapping : numpy.array
 |          mapping from initial starting sequence to final target sequence
 |          index of starting_sequence is mapping[index] of target sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from BaseAlignment:
 |  
 |  get_target_sequence(self)
 |      Gets the portion of starting sequence that fits the alignment
 |  
 |  map_dataframe(self, dataframe, position_columns)
 |      Takes a dataframe and maps position columns to target sequence.
 |      
 |      Rows with unmapped positions are dropped.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          a dataframe with position columns
 |      position_columns : list of str
 |          a list of columns containing positions to map
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a new dataframe (copy) with position columns mapped or dropped
 |  
 |  map_indices(self, indices, keep_minus_one=True)
 |      Takes a list of indices (0-index) and maps them to target sequence
 |      
 |      Parameters
 |      ----------
 |      indices : int or list of int
 |          a single or list of integer indices
 |      keep_minus_one : bool, defaults to True
 |          whether to keep unmapped starting sequence indices (-1) in the
 |          returned array.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          the equivalent indices in target sequence
 |  
 |  map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide', sequence_column='Sequence')
 |      Takes a per-nt dataframe and map it to the target sequence.
 |      
 |      Dataframe must have 1 row per nucleotide in starting sequence,
 |      with a position column and a sequence column. Dataframe is
 |      mapped to have the same format, but for target sequence nucleotides and
 |      positions.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          a per-nucleotide dataframe
 |      position_column : string, defaults to "Nucleotide"
 |          name of the position column.
 |      sequence_column : string, defaults to "Sequence"
 |          name of the sequence column.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a new dataframe (copy) mapped to target sequence.
 |          Unmapped starting sequence positions are dropped and unmapped
 |          target sequence positions are filled.
 |  
 |  map_positions(self, positions, keep_zero=True)
 |      Takes a list of positions (1-index) and maps them to target sequence
 |      
 |      Parameters
 |      ----------
 |      positions : int or list of int
 |          a single or list of integer positions
 |      keep_zero : bool, defaults to True
 |          whether to keep unmapped starting sequence positions (0) in the
 |          returned array.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          the equivalent positions in target sequence
 |  
 |  map_values(self, values, fill=nan)
 |      Takes an array of length equal to starting sequence and maps them to
 |      target sequence, unmapped positions in starting sequence are dropped
 |      and unmapped positions in target sequence are filled with fill value.
 |      
 |      Parameters
 |      ----------
 |      values : iterable
 |          values to map to target sequence.
 |      fill : any, defaults to np.nan
 |          a value for unmapped positions in target sequence.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          an array of values equal in length to target sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from BaseAlignment:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.StructureAlignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class StructureAlignment in module rnavigate.data.alignments

class StructureAlignment(BaseAlignment)
 |  StructureAlignment(sequence1, sequence2, structure1=None, structure2=None, full=False)
 |  
 |  Experimental secondary structure alignment based on RNAlign2D algorithm
 |  (https://doi.org/10.1186/s12859-021-04426-8)
 |  
 |  Parameters
 |  ----------
 |  sequence1 : string
 |      the sequence to be aligned
 |  sequence2 : string
 |      the sequence to align to
 |  structure1 : string, defaults to None
 |      the secondary structure of sequence1
 |  structure2 : string, defaults to None
 |      the secondary structure of sequence2
 |  full : bool, defaults to False
 |      whether to align to full length of sequence2 or just mapped length
 |  
 |  Attributes
 |  ----------
 |  sequence1 : str
 |      the sequence to be aligned
 |  sequence2 : str
 |      the sequence to align to
 |  structure1 : str
 |      the secondary structure of sequence1
 |  structure2 : str
 |      the secondary structure of sequence2
 |  alignment1 : str
 |      the alignment string matching sequence1 to sequence2
 |  alignment2 : str
 |      the alignment string matching sequence2 to sequence1
 |  starting_sequence : str
 |      sequence1
 |  target_sequence : str
 |      sequence2 if full is False, else alignment2
 |  mapping : numpy.array
 |      the alignment map array.
 |      index of starting_sequence is mapping[index] of target_sequence
 |  
 |  Method resolution order:
 |      StructureAlignment
 |      BaseAlignment
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sequence1, sequence2, structure1=None, structure2=None, full=False)
 |      Creates an alignment from structure1 to structure2.
 |  
 |  get_alignment(self)
 |      Aligns pseudo-amino-acid sequences according to RNAlign2D rules.
 |      
 |      Returns
 |      -------
 |      alignment1, alignment2 : tuple of 2 str
 |          the alignment strings matching sequence1 and sequence2, respectively.
 |  
 |  get_inverse_alignment(self)
 |      Gets an alignment that maps from sequence2 to sequence1.
 |  
 |  get_mapping(self)
 |      Calculates a mapping from starting sequence to target sequence.
 |      
 |      Returns
 |      -------
 |      mapping : numpy.array
 |          an array which maps an indices to the target sequence.
 |          starting_sequence[idx] == target_sequence[self.mapping[idx]]
 |  
 |  set_as_default_alignment(self)
 |      Set this as the default alignment between sequence1 and sequence2.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from BaseAlignment:
 |  
 |  get_target_sequence(self)
 |      Gets the portion of starting sequence that fits the alignment
 |  
 |  map_dataframe(self, dataframe, position_columns)
 |      Takes a dataframe and maps position columns to target sequence.
 |      
 |      Rows with unmapped positions are dropped.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          a dataframe with position columns
 |      position_columns : list of str
 |          a list of columns containing positions to map
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a new dataframe (copy) with position columns mapped or dropped
 |  
 |  map_indices(self, indices, keep_minus_one=True)
 |      Takes a list of indices (0-index) and maps them to target sequence
 |      
 |      Parameters
 |      ----------
 |      indices : int or list of int
 |          a single or list of integer indices
 |      keep_minus_one : bool, defaults to True
 |          whether to keep unmapped starting sequence indices (-1) in the
 |          returned array.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          the equivalent indices in target sequence
 |  
 |  map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide', sequence_column='Sequence')
 |      Takes a per-nt dataframe and map it to the target sequence.
 |      
 |      Dataframe must have 1 row per nucleotide in starting sequence,
 |      with a position column and a sequence column. Dataframe is
 |      mapped to have the same format, but for target sequence nucleotides and
 |      positions.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          a per-nucleotide dataframe
 |      position_column : string, defaults to "Nucleotide"
 |          name of the position column.
 |      sequence_column : string, defaults to "Sequence"
 |          name of the sequence column.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a new dataframe (copy) mapped to target sequence.
 |          Unmapped starting sequence positions are dropped and unmapped
 |          target sequence positions are filled.
 |  
 |  map_positions(self, positions, keep_zero=True)
 |      Takes a list of positions (1-index) and maps them to target sequence
 |      
 |      Parameters
 |      ----------
 |      positions : int or list of int
 |          a single or list of integer positions
 |      keep_zero : bool, defaults to True
 |          whether to keep unmapped starting sequence positions (0) in the
 |          returned array.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          the equivalent positions in target sequence
 |  
 |  map_values(self, values, fill=nan)
 |      Takes an array of length equal to starting sequence and maps them to
 |      target sequence, unmapped positions in starting sequence are dropped
 |      and unmapped positions in target sequence are filled with fill value.
 |      
 |      Parameters
 |      ----------
 |      values : iterable
 |          values to map to target sequence.
 |      fill : any, defaults to np.nan
 |          a value for unmapped positions in target sequence.
 |      
 |      Returns
 |      -------
 |      numpy.array
 |          an array of values equal in length to target sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from BaseAlignment:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.ScalarMappable

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class ScalarMappable in module rnavigate.data.colors

class ScalarMappable(matplotlib.cm.ScalarMappable)
 |  ScalarMappable(cmap, normalization, values, title='', tick_labels=None, **cbar_args)
 |  
 |  Used to map scalar values to a color and to create a colorbar plot.
 |  
 |  Parameters
 |  ----------
 |  cmap : str, tuple, float, or list
 |      A valid mpl color, list of valid colors or a valid colormap name
 |  normalization : "min_max", "0_1", "none", or "bins"
 |      The type of normalization to use when mapping values to colors
 |  values : list
 |      The values to use when normalizing the data
 |  title : str, defaults to ""
 |      The title of the colorbar.
 |  tick_labels : list, defaults to None
 |      The labels to use for the colorbar ticks. If None, values are
 |      determined automatically.
 |  **cbar_args : dict
 |      Additional arguments to pass to the colorbar function
 |  
 |  Attributes
 |  ----------
 |  rnav_norm : str
 |      The type of normalization to use when mapping values to colors
 |  rnav_vals : list
 |      The values to use when normalizing the data
 |  rnav_cmap : list
 |      The colors to use when mapping values to colors
 |  cbar_args : dict
 |      Additional arguments to pass to the colorbar function
 |  tick_labels : list
 |      The labels to use for the colorbar ticks. If None, values are
 |      determined automatically.
 |  title : str
 |      The title of the colorbar.
 |  
 |  Method resolution order:
 |      ScalarMappable
 |      matplotlib.cm.ScalarMappable
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, cmap, normalization, values, title='', tick_labels=None, **cbar_args)
 |      Initialize the ScalarMappable object
 |  
 |  get_cmap(self, cmap)
 |      Converts a cmap specification to a matplotlib colormap object.
 |      
 |      Parameters
 |      ----------
 |      cmap : string, tuple, float, or list
 |          A valid mpl color, list of valid colors or a valid colormap name
 |      
 |      Returns
 |      -------
 |      matplotlib colormap
 |          a colormap matching the input
 |  
 |  get_norm(self, normalization, values, cmap)
 |      Given a normalization type and values, return a normalization object.
 |      
 |      Parameters
 |      ----------
 |      normalization : "min_max", "0_1", "none", or "bins"
 |          The type of normalization to use when mapping values to colors
 |      values : list
 |          The values to use when normalizing the data
 |      cmap : matplotlib colormap
 |          The colormap to use when normalizing the data
 |      
 |      Returns
 |      -------
 |      matplotlib.colors normalization object
 |          Used to normalize data before mapping to colors
 |  
 |  is_equivalent_to(self, cmap2)
 |      Check if two ScalarMappable objects are equivalent.
 |      
 |      Parameters
 |      ----------
 |      cmap2 : ScalarMappable
 |          The ScalarMappable object to compare to
 |      
 |      Returns
 |      -------
 |      bool
 |          True if the two ScalarMappable objects are equivalent, False
 |          otherwise
 |  
 |  values_to_hexcolors(self, values, alpha=1.0)
 |      Map values to colors and return a list of hex colors.
 |      
 |      Parameters
 |      ----------
 |      values : list
 |          The values to map to colors
 |      alpha : float, defaults to 1.0
 |          The alpha value to use for the colors
 |      
 |      Returns
 |      -------
 |      list of strings
 |          A list of hex colors
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from matplotlib.cm.ScalarMappable:
 |  
 |  add_checker(self, checker)
 |      [*Deprecated*] 
 |      
 |      Notes
 |      -----
 |      .. deprecated:: 3.3
 |         \
 |  
 |  autoscale(self)
 |      Autoscale the scalar limits on the norm instance using the
 |      current array
 |  
 |  autoscale_None(self)
 |      Autoscale the scalar limits on the norm instance using the
 |      current array, changing only limits that are None
 |  
 |  changed(self)
 |      Call this whenever the mappable is changed to notify all the
 |      callbackSM listeners to the 'changed' signal.
 |  
 |  check_update(self, checker)
 |      [*Deprecated*] 
 |      
 |      Notes
 |      -----
 |      .. deprecated:: 3.3
 |         \
 |  
 |  get_alpha(self)
 |      Returns
 |      -------
 |      float
 |          Always returns 1.
 |  
 |  get_array(self)
 |      Return the data array.
 |  
 |  get_clim(self)
 |      Return the values (min, max) that are mapped to the colormap limits.
 |  
 |  set_array(self, A)
 |      Set the image array from numpy array *A*.
 |      
 |      Parameters
 |      ----------
 |      A : ndarray
 |  
 |  set_clim(self, vmin=None, vmax=None)
 |      Set the norm limits for image scaling.
 |      
 |      Parameters
 |      ----------
 |      vmin, vmax : float
 |           The limits.
 |      
 |           The limits may also be passed as a tuple (*vmin*, *vmax*) as a
 |           single positional argument.
 |      
 |           .. ACCEPTS: (vmin: float, vmax: float)
 |  
 |  set_cmap(self, cmap)
 |      Set the colormap for luminance data.
 |      
 |      Parameters
 |      ----------
 |      cmap : `.Colormap` or str or None
 |  
 |  set_norm(self, norm)
 |      Set the normalization instance.
 |      
 |      Parameters
 |      ----------
 |      norm : `.Normalize` or None
 |      
 |      Notes
 |      -----
 |      If there are any colorbars using the mappable for this norm, setting
 |      the norm of the mappable will reset the norm, locator, and formatters
 |      on the colorbar to default.
 |  
 |  to_rgba(self, x, alpha=None, bytes=False, norm=True)
 |      Return a normalized rgba array corresponding to *x*.
 |      
 |      In the normal case, *x* is a 1-D or 2-D sequence of scalars, and
 |      the corresponding ndarray of rgba values will be returned,
 |      based on the norm and colormap set for this ScalarMappable.
 |      
 |      There is one special case, for handling images that are already
 |      rgb or rgba, such as might have been read from an image file.
 |      If *x* is an ndarray with 3 dimensions,
 |      and the last dimension is either 3 or 4, then it will be
 |      treated as an rgb or rgba array, and no mapping will be done.
 |      The array can be uint8, or it can be floating point with
 |      values in the 0-1 range; otherwise a ValueError will be raised.
 |      If it is a masked array, the mask will be ignored.
 |      If the last dimension is 3, the *alpha* kwarg (defaulting to 1)
 |      will be used to fill in the transparency.  If the last dimension
 |      is 4, the *alpha* kwarg is ignored; it does not
 |      replace the pre-existing alpha.  A ValueError will be raised
 |      if the third dimension is other than 3 or 4.
 |      
 |      In either case, if *bytes* is *False* (default), the rgba
 |      array will be floats in the 0-1 range; if it is *True*,
 |      the returned rgba array will be uint8 in the 0 to 255 range.
 |      
 |      If norm is False, no normalization of the input data is
 |      performed, and it is assumed to be in the range (0-1).
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from matplotlib.cm.ScalarMappable:
 |  
 |  update_dict
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from matplotlib.cm.ScalarMappable:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.Interactions

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Interactions in module rnavigate.data.interactions

class Interactions(rnavigate.data.data.Data)
 |  Interactions(input_data, sequence, metric, metric_defaults, read_table_kw=None, window=1, name=None)
 |  
 |  A class for storing and manipulating interactions data.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing interactions data.
 |      If dataframe, the dataframe containing interactions data. The dataframe
 |      must contain columns "i", "j", and self.metric. Dataframe may also
 |      include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the interactions data.
 |  metric : string
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap)
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int
 |      The window size used to generate the interactions data.
 |  name : str
 |      The name of the data object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The interactions data.
 |  window : int
 |      The window size that is being represented by i-j pairs.
 |  
 |  Method resolution order:
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric, metric_defaults, read_table_kw=None, window=1, name=None)
 |      Initializes the Interactions object.
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Does nothing for the base Interactions class, can be overwritten in
 |      subclasses.
 |      
 |      Returns:
 |          dict: dictionary of keyword argument pairs
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.SHAPEJuMP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SHAPEJuMP in module rnavigate.data.interactions

class SHAPEJuMP(Interactions)
 |  SHAPEJuMP(input_data, sequence=None, metric='Percentile', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |  
 |  A class for storing and manipulating SHAPEJuMP data.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing SHAPEJuMP data.
 |      If dataframe, the dataframe containing SHAPEJuMP data. The dataframe
 |      must contain columns "i", "j", "Metric" (JuMP rate) and "Percentile"
 |      (percentile ranking). Dataframe may also include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the SHAPEJuMP data.
 |  metric : string, defaults to "Percentile"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap)
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int
 |      The window size used to generate the SHAPEJuMP data.
 |  name : str
 |      A name for the interactions object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The SHAPEJuMP data.
 |  
 |  Method resolution order:
 |      SHAPEJuMP
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence=None, metric='Percentile', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |      Constructs an Interactions object from SHAPEJuMP data
 |  
 |  read_file(self, input_data, read_table_kw=None)
 |      Parses a deletions.txt file and stores it as a dataframe.
 |      
 |      Also calculates a "Percentile" column.
 |      
 |      Parameters
 |      ----------
 |      input_data : str
 |          path to deletions.txt file
 |      read_table_kw : dict, defaults to {}
 |          kwargs passed to pandas.read_table().
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          the SHAPEJuMP data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Does nothing for the base Interactions class, can be overwritten in
 |      subclasses.
 |      
 |      Returns:
 |          dict: dictionary of keyword argument pairs
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.RINGMaP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class RINGMaP in module rnavigate.data.interactions

class RINGMaP(Interactions)
 |  RINGMaP(input_data, sequence=None, metric='Statistic', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |  
 |  A class for storing and manipulating RINGMaP data.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing RINGMaP data.
 |      If dataframe, the dataframe containing RINGMaP data. The dataframe
 |      must contain columns "i", "j", "Statistic", and "Zij". Dataframe may
 |      also include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the RINGMaP data.
 |  metric : string, defaults to "Statistic"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap)
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the RINGMaP data. If an input file is
 |      provided, this value is overwritten by the value in the header.
 |  name : str, optional
 |      A name for the interactions object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The RINGMaP data.
 |  
 |  Method resolution order:
 |      RINGMaP
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence=None, metric='Statistic', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |      Constructs an Interactions object from RINGMaP data
 |  
 |  data_specific_filter(self, positive_only=False, negative_only=False, **kwargs)
 |      Adds filters for "Sign" column to parent filter() function
 |      
 |      Parameters
 |      ----------
 |      positive_only : bool, defaults to False
 |          If True, only keep positive correlations.
 |      negative_only : bool, defaults to False
 |          If True, only keep negative correlations.
 |      
 |      Returns
 |      -------
 |      kwargs : dict
 |          any additional keyword-argument pairs are returned
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_sorted_data(self)
 |      Sorts on the product of self.metric and "Sign" columns.
 |      
 |      Except when self.metric is "Distance".
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by (self.metric * "Sign") columns
 |  
 |  read_file(self, filepath, read_table_kw=None)
 |      Parses a RINGMaP correlations file and stores data as a dataframe.
 |      
 |      Also sets self.window (usually 1, from header).
 |      
 |      Parameters
 |      ----------
 |      filepath : str
 |          path to correlations file.
 |      read_table_kw : dict, defaults to {}
 |          kwargs passed to pandas.read_table().
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          the RINGMaP data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.PAIRMaP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class PAIRMaP in module rnavigate.data.interactions

class PAIRMaP(RINGMaP)
 |  PAIRMaP(input_data, sequence=None, metric='Class', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |  
 |  A class for storing and manipulating PAIRMaP data.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing PAIRMaP data.
 |      If dataframe, the dataframe containing PAIRMaP data. The dataframe
 |      must contain columns "i", "j", "Statistic", and "Class". Dataframe may
 |      also include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the PAIRMaP data.
 |  metric : string, defaults to "Class"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap)
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the PAIRMaP data. If an input file is
 |      provided, this value is overwritten by the value in the header.
 |  name : str, optional
 |      A name for the interactions object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The PAIRMaP data.
 |  
 |  Method resolution order:
 |      PAIRMaP
 |      RINGMaP
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence=None, metric='Class', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |      Constructs an Interactions object from PAIRMaP data
 |  
 |  data_specific_filter(self, all_pairs=False, **kwargs)
 |      Used by Interactions.filter(). By default, non-primary and
 |      -secondary pairs are removed. all_pairs=True changes this behavior.
 |      
 |      Parameters
 |      ----------
 |      all_pairs : bool, defaults to False
 |          whether to include all PAIRs.
 |      
 |      Returns
 |      -------
 |      kwargs : dict
 |          any additional keyword-argument pairs are returned
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_sorted_data(self)
 |      Same as parent function, unless metric is set to "Class", in which
 |      case ij pairs are returned in a different order.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  read_file(self, filepath, read_table_kw=None)
 |      Parses a pairmap.txt file and stores data as a dataframe
 |      
 |      Sets self.window (usually 3), from header.
 |      
 |      Parameters
 |      ----------
 |      filepath : str
 |          path to pairmap.txt file
 |      read_table_kw : dict, defaults to None
 |          This argument is ignored.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.PairingProbability

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class PairingProbability in module rnavigate.data.interactions

class PairingProbability(Interactions)
 |  PairingProbability(input_data, sequence=None, metric='Probability', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |  
 |  A class for storing and manipulating pairing probability data.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing pairing probability data.
 |      If dataframe, the dataframe containing pairing probability data. The
 |      dataframe must contain columns "i", "j", "Probability", and "log10p".
 |      Dataframe may also include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the pairing probability data.
 |  metric : string, defaults to "Probability"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the pairing probability data.
 |  name : str, optional
 |      A name for the PairingProbability object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The pairing probability data.
 |  
 |  Method resolution order:
 |      PairingProbability
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence=None, metric='Probability', metric_defaults=None, read_table_kw=None, window=1, name=None)
 |      Constructs an Interactions object from pairing probability data
 |  
 |  data_specific_filter(self, **kwargs)
 |      By default, interactions with probabilities less than 0.03 are removed.
 |      
 |      Returns
 |      -------
 |      kwargs : dict
 |          any additional keyword-argument pairs are returned
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_entropy_profile(self, print_out=False, save_file=None)
 |      Calculates per-nucleotide Shannon entropy from pairing probabilities.
 |      
 |      Parameters
 |      ----------
 |      print_out : bool, defaults to False
 |          If True, entropy values are printed to console.
 |      save_file : str, defaults to None
 |          If not None, entropy values are saved to this file.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Profile
 |          a Profile object containing the entropy data
 |  
 |  read_file(self, filepath, read_table_kw=None)
 |      Parses a pairing probability file and stores data as a dataframe.
 |      
 |      Parameters
 |      ----------
 |      filepath : str
 |          path to pairing probability file
 |      read_table_kw : dict, defaults to None
 |          This argument is ignored.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          the pairing probability data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.AllPossible

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class AllPossible in module rnavigate.data.interactions

class AllPossible(Interactions)
 |  AllPossible(sequence, metric='data', input_data=None, metric_defaults=None, read_table_kw=None, window=1, name=None)
 |  
 |  A class for storing and manipulating all possible interactions.
 |  
 |  Parameters
 |  ----------
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the pairing probability data.
 |  metric : string, defaults to "Probability"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the pairing probability data.
 |  name : str, optional
 |      A name for the PairingProbability object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The pairing probability data.
 |  
 |  Method resolution order:
 |      AllPossible
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sequence, metric='data', input_data=None, metric_defaults=None, read_table_kw=None, window=1, name=None)
 |      Constructs an Interactions object from pairing probability data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Does nothing for the base Interactions class, can be overwritten in
 |      subclasses.
 |      
 |      Returns:
 |          dict: dictionary of keyword argument pairs
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.StructureAsInteractions

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class StructureAsInteractions in module rnavigate.data.interactions

class StructureAsInteractions(Interactions)
 |  StructureAsInteractions(input_data, sequence, metric=None, metric_defaults=None, window=1, name=None)
 |  
 |  A class for storing and manipulating structure data.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing structure data.
 |      If dataframe, the dataframe containing structure data. The dataframe
 |      must contain columns "i", "j", and "Structure". Dataframe may also
 |      include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the structure data.
 |  metric : string, defaults to "Structure"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the structure data.
 |  name : str, optional
 |      A name for the StructureAsInteractions object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The structure data.
 |  
 |  Method resolution order:
 |      StructureAsInteractions
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric=None, metric_defaults=None, window=1, name=None)
 |      Constructs an Interactions object from structure data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Does nothing for the base Interactions class, can be overwritten in
 |      subclasses.
 |      
 |      Returns:
 |          dict: dictionary of keyword argument pairs
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.StructureCompareMany

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class StructureCompareMany in module rnavigate.data.interactions

class StructureCompareMany(Interactions)
 |  StructureCompareMany(input_data, sequence, metric=None, metric_defaults=None, window=1, name=None)
 |  
 |  A class for storing and manipulating a comparison of many structures.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing structure data.
 |      If dataframe, the dataframe containing structure data. The dataframe
 |      must contain columns "i", "j", and "Structure". Dataframe may also
 |      include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the structure data.
 |  metric : string, defaults to "Structure"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the structure data.
 |  name : str, optional
 |      A name for the StructureAsInteractions object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The structure data.
 |  
 |  Method resolution order:
 |      StructureCompareMany
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric=None, metric_defaults=None, window=1, name=None)
 |      Constructs an Interactions object from structure data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Does nothing for the base Interactions class, can be overwritten in
 |      subclasses.
 |      
 |      Returns:
 |          dict: dictionary of keyword argument pairs
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.StructureCompareTwo

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class StructureCompareTwo in module rnavigate.data.interactions

class StructureCompareTwo(Interactions)
 |  StructureCompareTwo(input_data, sequence, metric=None, metric_defaults=None, window=1, name=None)
 |  
 |  A class for storing and manipulating a comparison of two structures.
 |  
 |  Parameters
 |  ----------
 |  input_data : string or pandas.DataFrame
 |      If string, a path to a file containing structure data.
 |      If dataframe, the dataframe containing structure data. The dataframe
 |      must contain columns "i", "j", and "Structure". Dataframe may also
 |      include other columns.
 |  sequence : string or rnavigate.data.Sequence
 |      The sequence string corresponding to the structure data.
 |  metric : string, defaults to "Structure"
 |      The column name to use for visualization.
 |  metric_defaults : dict
 |      Keys are metric names and values are dictionaries of metric-specific defaults.
 |      These defaults include:
 |          "metric_column" : string
 |              the column name to use for visualization
 |          "cmap" : string or matplotlib.colors.Colormap
 |              the colormap to use for visualization
 |          "normalization" : "min_max", "0_1", "none", or "bins"
 |              The type of normalization to use when mapping values to colors
 |          "values" : list of float
 |              The values to used with normalization of the data
 |          "title" : string
 |              the title to use for colorbars
 |          "extend" : "min", "max", "both", or "neither"
 |              Which ends to extend when drawing the colorbar.
 |          "tick_labels" : list of string
 |  read_table_kw : dict, optional
 |      kwargs passed to pandas.read_table() when reading input_data.
 |  window : int, defaults to 1
 |      The window size used to generate the structure data.
 |  name : str, optional
 |      A name for the StructureAsInteractions object.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The structure data.
 |  
 |  Method resolution order:
 |      StructureCompareTwo
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric=None, metric_defaults=None, window=1, name=None)
 |      Constructs an Interactions object from structure data
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a copy of the interactions, optionally with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      apply_filter : bool, defaults to False
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          A copy of the interactions.
 |  
 |  count_filter(self, **kwargs)
 |      Counts the number of interactions that pass the given filters.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Does nothing for the base Interactions class, can be overwritten in
 |      subclasses.
 |      
 |      Returns:
 |          dict: dictionary of keyword argument pairs
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Parameters
 |      ----------
 |      prefiltered : bool, defaults to False
 |          If True, the mask is not updated.
 |      reset_filter : bool, defaults to True
 |          If True, the mask is reset before applying filters.
 |      structure : rnavigate.data.SecondaryStructure, defaults to None
 |          The structure to use for filtering.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      profile : rnavigate.data.Profile, defaults to None
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      compliments_only : bool, defaults to False
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      max_distance : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_distance : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      exclude_nts : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate_nts : list of int, defaults to None
 |          A list of positions to isolate.
 |      resolve_conflicts : str, defaults to None
 |          If not None, conflicting windows are resolved using the Maximal
 |          Weighted Independent Set. The weights are taken from the metric
 |          value. The graph is first broken into components to speed up the
 |          identification of the MWIS. Then the mask is updated to only
 |          include the MWIS.
 |      **kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Returns a copy mapped to a new sequence with masked rows removed.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use for mapping the interactions.
 |      apply_filter : bool, defaults to True
 |          If True, masked rows ("mask" == False) are dropped.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Interactions
 |          Interactions mapped to a new sequence.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions.
 |      
 |      i and j are the 5' and 3' ends of each interaction, and colors is the color
 |      to use for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then mapped to
 |      a color using self.cmap.
 |      
 |      Returns
 |      -------
 |      i : list
 |          5' ends of each interaction
 |      j : list
 |          3' ends of each interaction
 |      colors : list
 |          colors to use for each interaction
 |  
 |  get_sorted_data(self)
 |      Returns a copy of the data sorted by self.metric.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          a copy of the data sorted by self.metric
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their distance in sequence space.
 |      
 |      Parameters
 |      ----------
 |      max_dist : int, defaults to None
 |          The maximum distance to allow. If None, no maximum distance is set.
 |      min_dist : int, defaults to None
 |          The minimum distance to allow. If None, no minimum distance is set.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Mask interactions based on their i and j positions.
 |      
 |      Parameters
 |      ----------
 |      exclude : list of int, defaults to None
 |          A list of positions to exclude.
 |      isolate : list of int, defaults to None
 |          A list of positions to isolate.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide measurements.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          The profile to use for masking.
 |      min_profile : float, defaults to None
 |          The minimum profile value to allow.
 |      max_profile : float, defaults to None
 |          The maximum profile value to allow.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence.
 |      
 |      Parameters
 |      ----------
 |      compliment_only : bool, defaults to None
 |          If True, only keep interactions where i and j are complimentary
 |          nucleotides.
 |      nts : str, defaults to None
 |          If compliment_only is False, only keep interactions where i and j
 |          are in nts.
 |      
 |      Returns
 |      -------
 |      numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Masks interactions based on a secondary structure.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          The secondary structure to use for masking.
 |      min_cd : int, defaults to None
 |          The minimum contact distance to allow.
 |      max_cd : int, defaults to None
 |          The maximum contact distance to allow.
 |      ss_only : bool, defaults to False
 |          If True, only keep interactions between single-stranded nucleotides.
 |      ds_only : bool, defaults to False
 |          If True, only keep interactions between double-stranded nucleotides.
 |      paired_only : bool, defaults to False
 |          If True, only keep interactions that are paired in the structure.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions based on values in self.data.
 |      
 |      Parameters
 |      ----------
 |      kwargs : dict
 |          Each keyword should have the format "column_operator" where column
 |          is a valid column name of the dataframe and operator is one of:
 |              "ge": greater than or equal to
 |              "le": less than or equal to
 |              "gt": greater than
 |              "lt": less than
 |              "eq": equal to
 |              "ne": not equal to
 |          The values given to these keywords are then used in the comparison
 |          and False comparisons are filtered out. e.g.:
 |              self.mask_on_values(Statistic_ge=23) evaluates to:
 |              self.update_mask(self.data["Statistic"] >= 23)
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  print_new_file(self, outfile=None)
 |      Create a new file with mapped and filtered interactions.
 |      
 |      Parameters
 |      ----------
 |      outfile : str, defaults to None
 |          path to an output file. If None, file string is printed to console.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Uses an experimental method to resolve conflicts.
 |      
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS. This method is computationally
 |      expensive for large or dense datasets.
 |      
 |      Parameters
 |      ----------
 |      metric : str, defaults to None
 |          The metric to use for weighting the graph. If None, self.metric is used.
 |      
 |      Returns
 |      -------
 |      mask : numpy array
 |          a boolean array of the same length as self.data
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Calculates the distance between atoms in i and j in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for calculating distances
 |      atom : str
 |          atom id to use for calculating distances
 |  
 |  update_mask(self, mask)
 |      Updates the mask by ANDing the current mask with the given mask.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.PDB

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class PDB in module rnavigate.data.pdb

class PDB(rnavigate.data.data.Sequence)
 |  PDB(input_data, chain, sequence=None, name=None)
 |  
 |  A class to represent RNA tertiary structures with atomic coordinates.
 |  
 |  This data can be used to filter interactions by 3D distance, and to visualize
 |  profile and interactions data on interactive 3D structures.
 |  
 |  Parameters
 |  ----------
 |  input_data : str
 |      path to a PDB or CIF file
 |  chain : str
 |      chain identifier of RNA of interest
 |  sequence : rnavigate.Sequence or str, optional
 |      A sequence to use as the reference sequence.
 |      This is required if the sequence cannot be found in the header
 |      Defaults to None.
 |  name : str, optional
 |      A name for the data set. Defaults to None.
 |  
 |  Attributes
 |  ----------
 |  sequence : str
 |      The RNA sequence
 |  length : int
 |      The length of the RNA sequence
 |  name : str
 |      A name for the data set
 |  path : str
 |      The path to the PDB or CIF file
 |  chain : str
 |      The chain identifier of the RNA of interest
 |  offset : int
 |      The offset between the sequence positions and the PDB residue indices
 |  pdb : Bio.PDB.Structure.Structure
 |      The PDB structure
 |  pdb_idx : np.array
 |      The PDB indices of the RNA
 |  pdb_seq : np.array
 |      The PDB sequence of the RNA
 |  distance_matrix : dict
 |      A dictionary of distance matrices for each atom type
 |  
 |  Method resolution order:
 |      PDB
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, chain, sequence=None, name=None)
 |      Construct PDB object based on an input PDB or CIF file.
 |  
 |  get_distance(self, i, j, atom="O2'")
 |      Get the distance between given atom in nucleotides i and j (1-indexed).
 |      
 |      Parameters
 |      ----------
 |      i : int
 |          The first nucleotide
 |      j : int
 |          The second nucleotide
 |      atom : string or dict, defaults to "O2'"
 |          The atom to use for distance calculations. If a string, the same atom
 |          will be used for all residues. If a dict, the atom will be chosen based
 |          on the nucleotide type. If "DMS", the N1 atom will be used for A and G,
 |          and the N3 atom will be used for U and C.
 |      
 |      Returns
 |      -------
 |      distance : float
 |          The distance between the atoms
 |  
 |  get_distance_matrix(self, atom="O2'")
 |      Get the pairwise atomic distance matrix for all residues.
 |      
 |      Parameters
 |      ----------
 |      atom : string or dict, defaults to "O2'"
 |          The atom to use for distance calculations. If a string, the same atom
 |          will be used for all residues. If a dict, the atom will be chosen based
 |          on the nucleotide type. If "DMS", the N1 atom will be used for A and G,
 |          and the N3 atom will be used for U and C.
 |      
 |      Returns
 |      -------
 |      matrix : NxN numpy.ndarray
 |          A 2D array of pairwise distances. N is the length of the RNA.
 |  
 |  get_pdb_idx(self, seq_idx)
 |      Return the PDB index given the sequence index (0-indexed).
 |  
 |  get_seq_idx(self, pdb_idx)
 |      Return the sequence index given the PDB index.
 |  
 |  get_sequence(self, pdb)
 |      Find the sequence in the provided CIF or PDB file.
 |      
 |      Parameters
 |      ----------
 |      pdb : str
 |          path to a PDB or CIF file
 |      
 |      Returns
 |      -------
 |      sequence : string
 |          The RNA sequence
 |  
 |  get_sequence_from_seqres(self, seqres)
 |      Used by get_sequence to parse the SEQRES entries.
 |      
 |      Parameters
 |      ----------
 |      seqres : list
 |          A list of SEQRES entries for the RNA chain of interest
 |      
 |      Returns
 |      -------
 |      sequence : string
 |          The RNA sequence
 |  
 |  get_xyz_coord(self, nt, atom)
 |      Return the x, y, and z coordinates for a given residue and atom.
 |      
 |      Parameters
 |      ----------
 |      nt : int
 |          The nucleotide of interest (1-indexed)
 |      atom : string or dict, defaults to "O2'"
 |          The atom to use for distance calculations. If a string, the same atom
 |          will be used for all residues. If a dict, the atom will be chosen based
 |          on the nucleotide type. If "DMS", the N1 atom will be used for A and G,
 |          and the N3 atom will be used for U and C.
 |      
 |      Returns
 |      -------
 |      xyz : list
 |          A list of x, y, and z coordinates
 |  
 |  is_valid_idx(self, pdb_idx=None, seq_idx=None)
 |      Determines if a PDB or sequence index is in the PDB structure.
 |      
 |      Parameters
 |      ----------
 |      pdb_idx : int, optional
 |          A PDB index (1-indexed). Defaults to None.
 |      seq_idx : int, optional
 |          A sequence index (1-indexed). Defaults to None.
 |      
 |      Returns
 |      -------
 |      bool
 |          True if the index is in the PDB structure, False otherwise.
 |  
 |  read_pdb(self, pdb)
 |      Read a PDB or CIF file into the data structure.
 |      
 |      Parameters
 |      ----------
 |      pdb : str
 |          path to a PDB or CIF file
 |  
 |  set_indices(self)
 |      Uses self.data and self.sequence to set self.offset
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_aligned_data(self, alignment)
 |      Get a copy of the sequence positionally aligned to another sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          the alignment to use
 |      
 |      Returns
 |      -------
 |      aligned_sequence : rnavigate.data.Sequence
 |          the aligned sequence
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.Profile

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Profile in module rnavigate.data.profile

class Profile(rnavigate.data.data.Data)
 |  Profile(input_data, metric='default', metric_defaults=None, read_table_kw=None, sequence=None, name=None)
 |  
 |  A class to represent per-nucleotide data.
 |  
 |  Parameters
 |  ----------
 |  input_data : str or pandas.DataFrame
 |      path to a csv or tab file or a pandas DataFrame
 |      Table must be 1 row for each nucleotide in the sequence.
 |      table columns must contain these columns:
 |          A nucleotide position column labelled "Nucleotide"
 |          A sequence column labelled "Sequence" with 1 of (A, C, G, U, T) per row
 |              These will be added to the table if `sequence` is provided.
 |          A data measurement column labelled "Profile" with a float or integer
 |              Label may be another name if specified in `metric_defaults`
 |          Optionally: A measurement error column.
 |              Label must be specified in `metric_defaults`
 |          Other columns may be present, and set up using `metric_defaults`.
 |              See `metric_defaults` for more information.
 |  read_table_kw : dict, optional
 |      Keyword arguments to pass to pandas.read_table.
 |      Defaults to None.
 |  sequence : rnavigate.Sequence or str, optional
 |      A sequence to use as the reference sequence.
 |      This is required if `input_data` does not contain a "Sequence" column.
 |      Defaults to None.
 |  metric : str, defaults to "default"
 |      The name of the set of value-to-color options to use.
 |      "default" specifies:
 |          "Profile" column is used
 |          No error rates are present
 |          Values are normalized to the range [0, 1]
 |          Values are mapped to colors using the "viridis" colormap
 |      "Distance" specifies:
 |          (3-D) "Distance" column is used
 |          No error rates are present
 |          Values in the range [5, 50] are normalized to the range [0, 1]
 |          Values are mapped to colors using the "cool" colormap
 |      Other options may be defined in `metric_defaults`.
 |  metric_defaults : dict, optional
 |      Keys are metric names, to be used with `metric`.
 |      Values are dictionaries of plotting parameters:
 |          "metric_column" : str
 |              The name of the column to use as the metric.
 |              Plots and analyses that use per-nucleotide data will use this column.
 |              If "color_column" is not provided, this column also defines colors.
 |          "error_column" : str or None
 |              The name of the column to use as the error.
 |              If None, no error bars are plotted.
 |          "color_column" : str or None
 |              The name of the column to use for coloring.
 |              If None, colors are defined by "metric_column".
 |          "cmap" : str or list
 |              The name of the colormap to use.
 |              If a list, the list of colors to use.
 |          "normalization" : str
 |              The type of normalization to use.
 |              In order to be used with colormaps, values are normalized to either
 |              be integers for categorical colormaps, or floats in the range [0, 1]
 |              for continuous colormaps.
 |              "none" : no normalization is performed
 |              "min_max" : values are scaled to floats in the range [0, 1] based on
 |                  the upper and lower bounds defined in "values"
 |              "0_1" : values are scaled to floats in the range [0, 1] based on
 |                  the minimum and maximum values in the data
 |              "bins" : values are scaled an integer based on bins defined by the
 |                  list of bounds defined in "values"
 |              "percentiles" : values are scaled to floats in the range [0, 1]
 |                  based on upper and lower percentile bounds defined by "values"
 |          "values" : list or None
 |              The values to use when normalizing the data.
 |              if "normalization" is "min_max", this should be a list of two values
 |                  defining the upper and lower bounds.
 |              if "normalization" is "bins", this should be a list of values
 |                  of length 1 less than the length of cmap.
 |                  example: [5, 10, 20] defines 4 bins:
 |                      (-infinity, 5), [5, 10), [10, 20), [20, infinity)
 |              if "normalization" is "percentiles", this should be a list of two
 |                  values defining the upper and lower percentile bounds.
 |              if "normalization" is "0_1" or "none", this should be None.
 |          "title" : str, defaults to ""
 |              The title of the colorbar.
 |          "ticks" : list, defaults to None
 |              The tick locations to use for the colorbar. If None, values are
 |              determined automatically.
 |          "tick_labels" : list, defaults to None
 |              The labels to use for the colorbar ticks. If None, values are
 |              determined automatically from "ticks".
 |          "extend" : "neither", "both", "min", or "max", defaults to "neither"
 |              Which ends of the colorbar to extend (places an arrow head).
 |      Defaults to None.
 |  name : str, optional
 |      A name for the data set. Defaults to None.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The data table
 |  
 |  Method resolution order:
 |      Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, metric='default', metric_defaults=None, read_table_kw=None, sequence=None, name=None)
 |      Initialize the Profile object.
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.SHAPEMaP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SHAPEMaP in module rnavigate.data.profile

class SHAPEMaP(Profile)
 |  SHAPEMaP(input_data, normalize=None, read_table_kw=None, sequence=None, metric='Norm_profile', metric_defaults=None, log=None, name=None)
 |  
 |  A class to represent per-nucleotide SHAPE-MaP data.
 |  
 |  Parameters
 |  ----------
 |  input_data : str or pandas.DataFrame
 |      path to a ShapeMapper2 profile.txt or .map file or a pandas DataFrame
 |  normalize : "DMS", "eDMS", "boxplot", "percentiles", or None, defaults to None
 |      The normalization method to use.
 |      "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |          scales the median of 90th to 95th percentiles to 1
 |          As and Cs are normalized seperately from Us and Gs
 |      "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |          Applies the new eDMS-MaP normalization.
 |          Each nucleotide is normalized seperately.
 |      "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |          removes outliers (> 1.5 iqr) and scales median to 1
 |          scales nucleotides together unless specified with nt_groups
 |      "percentiles" uses self.norm_percentile and nt_groups=['AUCG']
 |          scales the median of 90th to 95th percentiles to 1
 |          scales nucleotides together unless specified with nt_groups
 |      Defaults to None: no normalization is performed
 |  read_table_kw : dict, optional
 |      Keyword arguments to pass to pandas.read_table. These are not necessary for
 |      profile.txt and .map files.
 |      Defaults to None.
 |  sequence : rnavigate.Sequence or str, optional
 |      A sequence to use as the reference sequence. This is not necessary for
 |      profile.txt and .map files.
 |      Defaults to None.
 |  metric : str, defaults to "Norm_profile"
 |      The name of the set of value-to-color options to use.
 |      "Norm_profile" specifies:
 |          "Norm_profile" column is used
 |          "Norm_stderr" column is used for error bars
 |          Values are normalized to bins:
 |              (-inf, -0.4), [-0.4, 0.4), [0.4, 0.85), [0.85, 2), [2, inf)
 |          Bins are mapped to "grey", "black", "orange", "red", "red"
 |      Other options may be defined in `metric_defaults`.
 |  metric_defaults : dict, optional
 |      Keys are metric names, to be used with `metric`.
 |      Values are dictionaries of plotting parameters:
 |          "metric_column" : str
 |              The name of the column to use as the metric.
 |              Plots and analyses that use per-nucleotide data will use this column.
 |              If "color_column" is not provided, this column also defines colors.
 |          "error_column" : str or None
 |              The name of the column to use as the error.
 |              If None, no error bars are plotted.
 |          "color_column" : str or None
 |              The name of the column to use for coloring.
 |              If None, colors are defined by "metric_column".
 |          "cmap" : str or list
 |              The name of the colormap to use.
 |              If a list, the list of colors to use.
 |          "normalization" : str
 |              The type of normalization to use.
 |              In order to be used with colormaps, values are normalized to either
 |              be integers for categorical colormaps, or floats in the range [0, 1]
 |              for continuous colormaps.
 |              "none" : no normalization is performed
 |              "min_max" : values are scaled to floats in the range [0, 1] based on
 |                  the upper and lower bounds defined in "values"
 |              "0_1" : values are scaled to floats in the range [0, 1] based on
 |                  the minimum and maximum values in the data
 |              "bins" : values are scaled an integer based on bins defined by the
 |                  list of bounds defined in "values"
 |              "percentiles" : values are scaled to floats in the range [0, 1]
 |                  based on upper and lower percentile bounds defined by "values"
 |          "values" : list or None
 |              The values to use when normalizing the data.
 |              if "normalization" is "min_max", this should be a list of two values
 |                  defining the upper and lower bounds.
 |              if "normalization" is "bins", this should be a list of values
 |                  of length 1 less than the length of cmap.
 |                  example: [5, 10, 20] defines 4 bins:
 |                      (-infinity, 5), [5, 10), [10, 20), [20, infinity)
 |              if "normalization" is "percentiles", this should be a list of two
 |                  values defining the upper and lower percentile bounds.
 |              if "normalization" is "0_1" or "none", this should be None.
 |          "title" : str, defaults to ""
 |              The title of the colorbar.
 |          "ticks" : list, defaults to None
 |              The tick locations to use for the colorbar. If None, values are
 |              determined automatically.
 |          "tick_labels" : list, defaults to None
 |              The labels to use for the colorbar ticks. If None, values are
 |              determined automatically from "ticks".
 |          "extend" : "neither", "both", "min", or "max", defaults to "neither"
 |              Which ends of the colorbar to extend (places an arrow head).
 |      Defaults to None.
 |  log : str, optional
 |      Path to a ShapeMapper v2 shapemap_log.txt file with mutations-per-molecule
 |      and read-length histograms. These will be present if the --per-read-histogram
 |      flag was used when running ShapeMapper v2.
 |      Currently, this is not working with ShapeMapper v2.2 files.
 |      Defaults to None.
 |  name : str, optional
 |      A name for the data set. Defaults to None.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      The data table
 |  
 |  Method resolution order:
 |      SHAPEMaP
 |      Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, normalize=None, read_table_kw=None, sequence=None, metric='Norm_profile', metric_defaults=None, log=None, name=None)
 |      Initialize the SHAPEMaP object.
 |  
 |  read_log(self, log)
 |      Read the ShapeMapper log file.
 |      
 |      Parameters
 |      ----------
 |      log : str
 |          Path to a ShapeMapper v2 shapemap_log.txt file with
 |          mutations-per-molecule and read-length histograms.
 |      
 |      Returns
 |      -------
 |      read_lengths : pandas.DataFrame
 |          A dataframe with the columns "Read_length", "Modified_read_length",
 |          and "Untreated_read_length".
 |      mutations_per_molecule : pandas.DataFrame
 |          A dataframe with the columns "Mutation_count",
 |          "Modified_mutations_per_molecule", and
 |          "Untreated_mutations_per_molecule".
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_rnaframework(input_data, normalize=None) from builtins.type
 |      Construct a SHAPEMaP object from an RNAFramework output file.
 |      
 |      Parameters
 |      ----------
 |      input_data : str
 |          path to an RNAFramework .xml reactivities file
 |      normalize : "DMS", "eDMS", "boxplot", "percentiles", or None, defaults to None
 |          The normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentiles" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to None: no normalization is performed
 |      
 |      Returns
 |      -------
 |      SHAPEMaP
 |          A SHAPEMaP object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Profile:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.DanceMaP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class DanceMaP in module rnavigate.data.profile

class DanceMaP(SHAPEMaP)
 |  DanceMaP(input_data, component, read_table_kw=None, sequence=None, metric='Norm_profile', metric_defaults=None, name=None)
 |  
 |  A class to represent per-nucleotide DanceMaP data.
 |  
 |  Parameters
 |  ----------
 |  input_data : str or pandas.DataFrame
 |      path to a DanceMapper reactivities.txt file or a pandas DataFrame
 |  component : int
 |      Which component of the DanceMapper ensemble to read in (0-indexed).
 |  read_table_kw : dict, optional
 |      Keyword arguments to pass to pandas.read_table. These are not necessary for
 |      reactivities.txt files.
 |      Defaults to None.
 |  sequence : rnavigate.Sequence or str, optional
 |      A sequence to use as the reference sequence. This is not necessary for
 |      reactivities.txt files.
 |      Defaults to None.
 |  metric : str, defaults to "Norm_profile"
 |      The name of the set of value-to-color options to use.
 |  
 |  Method resolution order:
 |      DanceMaP
 |      SHAPEMaP
 |      Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, component, read_table_kw=None, sequence=None, metric='Norm_profile', metric_defaults=None, name=None)
 |      Initialize the SHAPEMaP object.
 |  
 |  read_file(self, input_data, read_table_kw={})
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from SHAPEMaP:
 |  
 |  read_log(self, log)
 |      Read the ShapeMapper log file.
 |      
 |      Parameters
 |      ----------
 |      log : str
 |          Path to a ShapeMapper v2 shapemap_log.txt file with
 |          mutations-per-molecule and read-length histograms.
 |      
 |      Returns
 |      -------
 |      read_lengths : pandas.DataFrame
 |          A dataframe with the columns "Read_length", "Modified_read_length",
 |          and "Untreated_read_length".
 |      mutations_per_molecule : pandas.DataFrame
 |          A dataframe with the columns "Mutation_count",
 |          "Modified_mutations_per_molecule", and
 |          "Untreated_mutations_per_molecule".
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from SHAPEMaP:
 |  
 |  from_rnaframework(input_data, normalize=None) from builtins.type
 |      Construct a SHAPEMaP object from an RNAFramework output file.
 |      
 |      Parameters
 |      ----------
 |      input_data : str
 |          path to an RNAFramework .xml reactivities file
 |      normalize : "DMS", "eDMS", "boxplot", "percentiles", or None, defaults to None
 |          The normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentiles" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to None: no normalization is performed
 |      
 |      Returns
 |      -------
 |      SHAPEMaP
 |          A SHAPEMaP object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.RNPMaP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class RNPMaP in module rnavigate.data.profile

class RNPMaP(Profile)
 |  RNPMaP(input_data, read_table_kw=None, sequence=None, metric='NormedP', metric_defaults=None, name=None)
 |  
 |  Represents per-nucleotide RNPMaP data.
 |  
 |  Parameters
 |  ----------
 |  input_data : str or pandas.DataFrame
 |      path to an RNAModMapper reactivities.txt file or a pandas DataFrame
 |  read_table_kw : dict, optional
 |      Keyword arguments to pass to pandas.read_table. These are not necessary for
 |      reactivities.txt files.
 |      Defaults to None.
 |  sequence : rnavigate.Sequence or str, optional
 |      A sequence to use as the reference sequence. This is not necessary for
 |      reactivities.txt files.
 |      Defaults to None.
 |  metric : str, defaults to "NormedP"
 |      The name of the set of value-to-color options to use.
 |  metric_defaults : dict, optional
 |      Keys are metric names, to be used with `metric`.
 |      Values are dictionaries of plotting parameters.
 |      Defaults to None.
 |  name : str, optional
 |      A name for the data set. Defaults to None.
 |  
 |  Method resolution order:
 |      RNPMaP
 |      Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, read_table_kw=None, sequence=None, metric='NormedP', metric_defaults=None, name=None)
 |      Initialize the Profile object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Profile:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.DeltaProfile

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class DeltaProfile in module rnavigate.data.profile

class DeltaProfile(Profile)
 |  DeltaProfile(profile1, profile2, metric=None, metric_defaults=None, name=None)
 |  
 |  A class to represent the difference between two profiles.
 |  
 |  Parameters
 |  ----------
 |  profile1 : Profile
 |      The first profile to compare.
 |  profile2 : Profile
 |      The second profile to compare.
 |  metric : str, optional
 |      The name of the metric to use.
 |      Defaults to the metric of profile1.
 |  metric_defaults : dict, optional
 |      Keys are metric names, to be used with `metric`.
 |      Values are dictionaries of plotting parameters.
 |      Defaults to None.
 |  name : str, optional
 |      A name for the data set. Defaults to None.
 |  
 |  Method resolution order:
 |      DeltaProfile
 |      Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, profile1, profile2, metric=None, metric_defaults=None, name=None)
 |      Initialize the Profile object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of data.
 |      
 |      Result is stored in a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Parameters
 |      ----------
 |      column : str
 |          name of column to perform operation on
 |      window : int
 |          window size, must be an odd number
 |      method : string or function, defaults to "median"
 |          operation to perform over windows.
 |          if string, must be "median", "mean", "minimum", or "maximum"
 |          if function, must take a 1D numpy array as input and return a scalar
 |      new_name : str, defaults to f"{method}_{window}_nt"
 |          name of new column for stored result.
 |      minimum_points : int, defaults to value of `window`
 |          minimum number of points within each window.
 |      mask_na : bool, defaults to True
 |          whether to mask the result of the operation where the original
 |          column has a nan value.
 |  
 |  copy(self)
 |      Returns a copy of the Profile.
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new Profile object with the data aligned to a sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.SequenceAlignment
 |          The alignment to use to map rows of self.data to a new sequence.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A new Profile object with the data aligned to the sequence in the
 |          alignment.
 |  
 |  get_plotting_dataframe(self)
 |      Returns a dataframe with the data to be plotted.
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |          A dataframe with the columns "Nucleotide", "Values", "Errors", and
 |          "Colors".
 |  
 |  norm_boxplot(self, values)
 |      removes outliers (> 1.5 * IQR) and scales the mean to 1.
 |      
 |      NOTE: This method varies slightly from normalization method used in the
 |      SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and then
 |      uses these values when computing iqr and 90th percentile. Including
 |      these values can skew these result. This method excludes such nan
 |      values. Other elements are the same.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Calculates norm factors following eDMS pernt scheme in ShapeMapper 2.2
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates factors to scale the median between percentile bounds to 1.
 |      
 |      Parameters
 |      ----------
 |      values : 1D numpy array
 |          values to normalize
 |      lower_bound : int or float, optional
 |          percentile of lower bound, Defaults to 90
 |      upper_bound : int or float, optional
 |          percentile of upper bound, Defaults to 99
 |      median_or_mean : string, optional
 |          whether to use the median or mean of the values between the bounds.
 |      
 |      Returns
 |      -------
 |      (float, float)
 |          scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Parameters
 |      ----------
 |      profile_column : string, defaults to self.metric
 |          column name of values to normalize
 |      new_profile : string, defaults to "Norm_profile"
 |          column name of new normalized values
 |      error_column : string, defaults to self.error_column
 |          column name of error values to propagate
 |      new_error : string, defaults to "Norm_error"
 |          column name of new propagated error values
 |      norm_method : string, defaults to "boxplot"
 |          normalization method to use.
 |          "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |              scales the median of 90th to 95th percentiles to 1
 |              As and Cs are normalized seperately from Us and Gs
 |          "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |              Applies the new eDMS-MaP normalization.
 |              Each nucleotide is normalized seperately.
 |          "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |              removes outliers (> 1.5 iqr) and scales median to 1
 |              scales nucleotides together unless specified with nt_groups
 |          "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |              scales the median of 90th to 95th percentiles to 1
 |              scales nucleotides together unless specified with nt_groups
 |          Defaults to "boxplot": the default normalization of ShapeMapper
 |      nt_groups : list of strings, defaults to None
 |          A list of nucleotides to group
 |          e.g. ['AUCG'] groups all nts together
 |                  ['AC', 'UG'] groups As with Cs and Us with Gs
 |                  ['A', 'C', 'U', 'G'] scales each nt seperately
 |          Default depends on norm_method
 |      profile_factors : dictionary, defaults to None
 |          a scaling factor (float) for each nucleotide. keys must be:
 |              'A', 'C', 'U', 'G'
 |          Note: using this argument overrides any calculation of scaling
 |          Defaults to None
 |      **norm_kwargs
 |          these are passed to the norm_method function
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Parameters
 |      ----------
 |      profiles : list of rnavigate.data.Profile
 |          a list of other profiles used to compute scaling factors
 |      
 |      Returns
 |      -------
 |      profile_factors : dict
 |          the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Changes the values in self.data["Sequence"] to the normalized sequence.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T" or "U", Defaults to "U".
 |          Whether to replace T with U or U with T.
 |      uppercase : bool, Defaults to True.
 |          Whether to convert the sequence to uppercase.
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Parameters
 |      ----------
 |      column : string
 |          the column of data to be winsorized
 |      lower_bound : Number or None, defaults to None
 |          Data below this value is set to this value.
 |          If None, no lower bound is applied.
 |      upper_bound : Number or None, defaults to None
 |          Data above this value is set to this value.
 |          If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |      Construct a Profile object from an array of values.
 |      
 |      Parameters
 |      ----------
 |      input_data : list or np.array
 |          A list or array of values to use as the metric.
 |      sequence : str
 |          The RNA sequence.
 |      **kwargs
 |          Additional keyword arguments to pass to the Profile constructor.
 |      
 |      Returns
 |      -------
 |      Profile
 |          A Profile object with the provided values.
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Profile:
 |  
 |  recreation_kwargs
 |      A dictionary of keyword arguments to pass when recreating the object.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |      Add metric defaults to self.metric_defaults
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Parameters
 |      ----------
 |      filepath : string
 |          path to data file containing interactions
 |      read_table_kw : dict
 |          kwargs dictionary passed to pd.read_table
 |      
 |      Returns
 |      -------
 |      dataframe : pandas.DataFrame
 |          the data table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |      Get the colormap to use for colorbars and to retrieve colors.
 |  
 |  color_column
 |      Get the column of the dataframe to use as the color for visualization.
 |  
 |  colors
 |      Get one matplotlib color-like value for each nucleotide in self.sequence.
 |  
 |  error_column
 |      Get the column of the dataframe to use as the error for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |      Get the column of the dataframe to use as the metric for visualization.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.Annotation

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Annotation in module rnavigate.data.annotation

class Annotation(rnavigate.data.data.Sequence)
 |  Annotation(input_data, annotation_type, sequence, name=None, color='blue')
 |  
 |  Basic annotation class to store 1D features of an RNA sequence
 |  
 |  Each feature type must be a seperate instance. Feature types include:
 |      a group of separted nucleotides (e.g. binding pocket)
 |      regions of interest (e.g. coding sequence, Alu elements)
 |      sites of interest (e.g. m6A locations)
 |      primer binding sites.
 |  
 |  Parameters
 |  ----------
 |  input_data : list
 |      List will be treated according to `annotation_type` argument.
 |      Expected behaviors for each value of `annotation_type`:
 |      "sites" or "group": 1-indexed location of sites of interest
 |          example: [1, 10, 20, 30] is four sites, 1, 10, 20, and 30
 |      "spans": 1-indexed, inclusive locations of spans of interest
 |          example: [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30
 |      "primers": Similar to spans, but 5'/3' direction is preserved.
 |          example: [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20
 |  annotation_type : "group", "sites", "spans", or "primers"
 |      The type of annotation.
 |  sequence : str or pandas.DataFrame
 |      Nucleotide sequence, path to fasta file, or dataframe containing a
 |      "Sequence" column.
 |  name : str, defaults to None
 |      Name of annotation.
 |  color : matplotlib color-like, defaults to "blue"
 |      Color to be used for displaying this annotation on plots.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      Stores the list of sites or regions
 |  name : str
 |      The label for this annotation for use on plots
 |  color : valid matplotlib color
 |      Color to represent annotation on plots
 |  sequence : str
 |      The reference sequence string
 |  
 |  Method resolution order:
 |      Annotation
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __getitem__(self, idx)
 |  
 |  __init__(self, input_data, annotation_type, sequence, name=None, color='blue')
 |      Create an Annotation.
 |  
 |  __iter__(self)
 |  
 |  __len__(self)
 |  
 |  from_sites(self, sites)
 |      Create the self.data dataframe from a list of sites.
 |  
 |  from_spans(self, spans)
 |      Create the self.data dataframe from a list of spans.
 |  
 |  get_aligned_data(self, alignment)
 |      Aligns this Annotation to a new sequence and returns a copy.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          Alignment object used to align to a new sequence.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Annotation
 |          A new Annotation with the same name, color, and annotation
 |          type, but with the input data aligned to the target sequence.
 |  
 |  get_sites(self)
 |      Returns a list of nucleotide positions included in this annotation.
 |      
 |      Returns
 |      -------
 |      sites : tuple
 |          a list of nucleotide positions
 |  
 |  get_subsequences(self, buffer=0)
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_boolean_array(values, sequence, annotation_type, name, color='blue', window=1) from builtins.type
 |      Create an Annotation from an array of boolean values.
 |      
 |      True values are used to create the Annotation.
 |      
 |      Parameters
 |      ----------
 |      values : list of True or False
 |          the boolean array
 |      sequence : string or rnav.data.Sequence
 |          the sequence of the Annotation
 |      annotation_type : "spans", "sites", "primers", or "group"
 |          the type of the new annotation
 |          If "spans" or "primers", adjacent True values, or values within
 |          window are collapse to a region.
 |      name : string
 |          a name for labelling the annotation.
 |      color : string, defaults to "blue"
 |          a color for plotting the annotation
 |      window : integer, defaults to 1
 |          a window around True values to include in the annotation.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Annotation
 |          the new Annotation
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.Motif

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Motif in module rnavigate.data.annotation

class Motif(Annotation)
 |  Motif(input_data, sequence, name=None, color='blue')
 |  
 |  Automatically annotates the occurances of a sequence motif as spans.
 |  
 |  Parameters
 |  ----------
 |  input_data : str
 |      sequence motif to search for.
 |      Uses conventional nucleotide codes.
 |      e.g. "DRACH" = [AGTU] [AG] A C [ATUC]
 |  sequence : str or pandas.DataFrame
 |      Nucleotide sequence, path to fasta file, or dataframe containing a
 |      "Sequence" column.
 |  name : str, defaults to None
 |      Name of annotation.
 |  color : matplotlib color-like, defaults to "blue"
 |      Color to be used for displaying this annotation on plots.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      Stores the list of regions that match the motif
 |  name : str
 |      The label for this annotation for use on plots
 |  color : valid matplotlib color
 |      Color to represent annotation on plots
 |  sequence : str
 |      The reference sequence string
 |  
 |  Method resolution order:
 |      Motif
 |      Annotation
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, name=None, color='blue')
 |      Creates a Motif annotation
 |  
 |  get_aligned_data(self, alignment)
 |      Searches the new sequence for the motif and returns a new Motif annotation.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          Alignment object used to align to a new sequence.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Motif
 |          A new Motif with the same name, color, and motif
 |          but with the input data aligned to the target sequence.
 |  
 |  get_spans_from_motif(self, sequence, motif)
 |      Returns a list of spans for each location of motif found within sequence.
 |      
 |      Parameters
 |      ----------
 |      sequence : string
 |          sequence to be searched
 |      motif : string
 |          sequence motif to searched for.
 |      
 |      Returns
 |      -------
 |      spans : list of lists
 |          list of [start, end] positions of each motif occurance
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Annotation:
 |  
 |  __getitem__(self, idx)
 |  
 |  __iter__(self)
 |  
 |  __len__(self)
 |  
 |  from_sites(self, sites)
 |      Create the self.data dataframe from a list of sites.
 |  
 |  from_spans(self, spans)
 |      Create the self.data dataframe from a list of spans.
 |  
 |  get_sites(self)
 |      Returns a list of nucleotide positions included in this annotation.
 |      
 |      Returns
 |      -------
 |      sites : tuple
 |          a list of nucleotide positions
 |  
 |  get_subsequences(self, buffer=0)
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Annotation:
 |  
 |  from_boolean_array(values, sequence, annotation_type, name, color='blue', window=1) from builtins.type
 |      Create an Annotation from an array of boolean values.
 |      
 |      True values are used to create the Annotation.
 |      
 |      Parameters
 |      ----------
 |      values : list of True or False
 |          the boolean array
 |      sequence : string or rnav.data.Sequence
 |          the sequence of the Annotation
 |      annotation_type : "spans", "sites", "primers", or "group"
 |          the type of the new annotation
 |          If "spans" or "primers", adjacent True values, or values within
 |          window are collapse to a region.
 |      name : string
 |          a name for labelling the annotation.
 |      color : string, defaults to "blue"
 |          a color for plotting the annotation
 |      window : integer, defaults to 1
 |          a window around True values to include in the annotation.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Annotation
 |          the new Annotation
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.ORFs

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class ORFs in module rnavigate.data.annotation

class ORFs(Annotation)
 |  ORFs(input_data, name=None, sequence=None, color='blue')
 |  
 |  Automatically annotations occurances of open-reading frames as spans.
 |  
 |  Parameters
 |  ----------
 |  input_data : "longest" or "all"
 |      which ORFs to annotate. "longest" annotates the longest ORF. "all"
 |      annotates all potential ORFs.
 |  sequence : str or pandas.DataFrame
 |      Nucleotide sequence, path to fasta file, or dataframe containing a
 |      "Sequence" column.
 |  name : str, defaults to None
 |      Name of annotation.
 |  color : matplotlib color-like, defaults to "blue"
 |      Color to be used for displaying this annotation on plots.
 |  
 |  Attributes
 |  ----------
 |  data : pandas.DataFrame
 |      Stores the list of regions that match the motif
 |  name : str
 |      The label for this annotation for use on plots
 |  color : valid matplotlib color
 |      Color to represent annotation on plots
 |  sequence : str
 |      The reference sequence string
 |  
 |  Method resolution order:
 |      ORFs
 |      Annotation
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, name=None, sequence=None, color='blue')
 |      Creates an ORF annotation
 |  
 |  get_aligned_data(self, alignment)
 |      Searches the new sequence for ORFs and returns a new ORF annotation.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          Alignment object used to align to a new sequence.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.ORFs
 |          A new ORFs annotation with the same name, color, and input_data
 |          but with the input data aligned to the target sequence.
 |  
 |  get_spans_from_orf(self, sequence, which='all')
 |      Given a sequence string, returns spans for specified ORFs
 |      
 |      Parameters
 |      ----------
 |      sequence : string
 |          RNA nucleotide sequence
 |      which : "longest" or "all", defaults to "all"
 |          "all" returns all spans, "longest" returns the longest span
 |      
 |      Returns
 |      -------
 |      list of tuples
 |          (start, end) position of each ORF 1-indexed, inclusive
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Annotation:
 |  
 |  __getitem__(self, idx)
 |  
 |  __iter__(self)
 |  
 |  __len__(self)
 |  
 |  from_sites(self, sites)
 |      Create the self.data dataframe from a list of sites.
 |  
 |  from_spans(self, spans)
 |      Create the self.data dataframe from a list of spans.
 |  
 |  get_sites(self)
 |      Returns a list of nucleotide positions included in this annotation.
 |      
 |      Returns
 |      -------
 |      sites : tuple
 |          a list of nucleotide positions
 |  
 |  get_subsequences(self, buffer=0)
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Annotation:
 |  
 |  from_boolean_array(values, sequence, annotation_type, name, color='blue', window=1) from builtins.type
 |      Create an Annotation from an array of boolean values.
 |      
 |      True values are used to create the Annotation.
 |      
 |      Parameters
 |      ----------
 |      values : list of True or False
 |          the boolean array
 |      sequence : string or rnav.data.Sequence
 |          the sequence of the Annotation
 |      annotation_type : "spans", "sites", "primers", or "group"
 |          the type of the new annotation
 |          If "spans" or "primers", adjacent True values, or values within
 |          window are collapse to a region.
 |      name : string
 |          a name for labelling the annotation.
 |      color : string, defaults to "blue"
 |          a color for plotting the annotation
 |      window : integer, defaults to 1
 |          a window around True values to include in the annotation.
 |      
 |      Returns
 |      -------
 |      rnavigate.data.Annotation
 |          the new Annotation
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.data.domains

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function domains in module rnavigate.data.annotation

domains(input_data, names, colors, sequence)
    Create a list of Annotations from a list of spans.
    
    Currently, domains functionality in RNAvigate just uses a list of spans. In the
    future, this should be a dedicated class. Generally, domains should cover an entire
    sequence without overlap, but this is not enforced.
    e.g. [[1, 100], [101, 200]] for a 200 nt sequence.
    
    Parameters
    ----------
    input_data : list of lists
        list of spans for each domain
    names : list of strings
        list of names for each domain
    colors : list of valid matplotlib colors
        list of colors for each domain
    sequence : string
        sequence to be annotated
    
    Returns
    -------
    list of rnavigate.data.Annotation
        list of Annotations
```

## rnavigate.plots

### rnavigate.plots.Alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Alignment in module rnavigate.plots.alignment

class Alignment(rnavigate.plots.plots.Plot)
 |  Alignment(num_samples, rows=None, cols=1, **kwargs)
 |  
 |  Class for plotting sequence alignments.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      The number of samples to plot. Always 2.
 |  rows : int, optional
 |      The number of rows of plots. Always 1.
 |  cols : int, optional
 |      The number of columns of plots. Always 1.
 |  **kwargs : dict
 |      Additional keyword arguments to pass to plots.Plot.
 |  
 |  Attributes
 |  ----------
 |  region : tuple
 |      The region of the alignment to plot. Defaults to (1, len(alignment1))
 |  fig : matplotlib.figure.Figure
 |      The figure containing the plot.
 |  axes : numpy.array
 |      A 1x1 array of the axes containing the plot.
 |  
 |  Method resolution order:
 |      Alignment
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, rows=None, cols=1, **kwargs)
 |      Create a new Alignment object.
 |  
 |  plot_data(self, alignment, label, ax=None)
 |      Add the alignment to the next (or given) axes.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.Alignment
 |          The alignment to plot.
 |      label : str
 |          The label for the alignment.
 |      ax : matplotlib.axes.Axes, optional
 |          The axes containing the plot.
 |          Defaults to None.
 |  
 |  set_figure_size(self, height_ax_rel=0.03, width_ax_rel=0.03, width_ax_in=None, height_ax_in=None, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set the figure size for the plot.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          The relative height of each axes.
 |          Defaults to 0.03.
 |      width_ax_rel : float, optional
 |          The relative width of each axes.
 |          Defaults to 0.03.
 |      width_ax_in : float, optional
 |          The width of each axes in inches.
 |          Defaults to None.
 |      height_ax_in : float, optional
 |          The height of each axes in inches.
 |          Defaults to None.
 |      height_gap_in : float, optional
 |          The height of the gap between axes in inches.
 |          Defaults to 1.
 |      width_gap_in : float, optional
 |          The width of the gap between axes in inches.
 |          Defaults to 0.5.
 |      top_in : float, optional
 |          The top margin in inches.
 |          Defaults to 1.
 |      bottom_in : float, optional
 |          The bottom margin in inches.
 |          Defaults to 0.5.
 |      left_in : float, optional
 |          The left margin in inches.
 |          Defaults to 0.5.
 |      right_in : float, optional
 |          The right margin in inches.
 |          Defaults to 0.5.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.AP

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class AP in module rnavigate.plots.arc

class AP(rnavigate.plots.plots.Plot)
 |  AP(num_samples, nt_length, region='all', track_labels=True, **kwargs)
 |  
 |  Class for plotting arc plots
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot
 |  nt_length : int
 |      Length of the sequence
 |  region : tuple of 2 integers, optional
 |      starting and ending positions of the region to plot.
 |      Default is "all", which plots the entire sequence.
 |  track_labels : bool, optional
 |      Whether to plot track labels. Default is True.
 |  **kwargs
 |      Additional keyword arguments are passed to plots.Plot
 |  
 |  Attributes
 |  ----------
 |  nt_length : int
 |      Length of the sequence
 |  region : tuple of 2 integers
 |      starting and ending positions of the region to plot.
 |  track_labels : bool
 |      Whether to plot track labels.
 |  fig : matplotlib.figure.Figure
 |      Figure object containing the plot
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects containing the plots
 |  i : int
 |      Index of the current plot
 |  
 |  Method resolution order:
 |      AP
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, nt_length, region='all', track_labels=True, **kwargs)
 |      Initialize AP object.
 |  
 |  plot_data(self, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, annotations=None, domains=None, label='', ax=None, seqbar=True, title=True, panels=None, annotation_mode='track', track_height=None, profile_scale_factor=1, plot_error=False, nt_ticks=(20, 5))
 |      Add data to the next (or specified) plot axes.
 |      
 |      This function assumes data has already been aligned to a common sequence.
 |      rnavigate.plot_ functions can be used to automatically align data.
 |      
 |      Parameters
 |      ----------
 |      sequence : rnavigate.data.Sequence
 |          Sequence object containing the sequence to plot
 |      structure : rnavigate.data.Structure, optional
 |          Structure object containing a structure to plot
 |      structure2 : rnavigate.data.Structure, optional
 |          Structure object containing a structure to compare to the first
 |      interactions : rnavigate.data.Interactions, optional
 |          Interactions object containing inter-nucleotide data to plot
 |      interactions2 : rnavigate.data.Interactions, optional
 |          Interactions object containing other inter-nucleotide data to plot
 |      profile : rnavigate.data.Profile, optional
 |          Profile object containing per-nucleotide data to plot
 |      annotations : list of rnavigate.data.Annotation, optional
 |          List of Annotation objects containing annotations to plot
 |      domains : list of rnavigate.data.Spans, optional
 |          List of Spans objects containing domains to plot
 |      label : str, defaults to ""
 |          Label for the title of the plot.
 |      ax : matplotlib.axes.Axes, optional
 |          Axes object to plot on. If None, the next axes in the figure will be used.
 |      seqbar : bool, Defaults to True
 |          Whether to plot the sequence track.
 |      title : bool, defaults to True
 |          Whether to show the title.
 |      panels : dict, optional
 |          Dictionary of panels to plot, with keys being the panel name and values
 |          being the panel location. Default is {"interactions": "bottom",
 |          "interactions2": "bottom", "structure": "top", "profile": "top"}
 |      annotation_mode : "track" or "bar", defaults to "track"
 |          Mode for plotting annotations.
 |      track_height : int, optional
 |          Height of the track. If None, the height is automatically determined.
 |      profile_scale_factor : float, defaults to 1
 |          Scale factor for the profile track.
 |      plot_error : bool, defaults to False
 |          Whether to plot the error bars for the profile track.
 |      nt_ticks : tuple of 2 ints, optional
 |          Major and minor tick spacing for the nucleotide axis. Default is (20, 5).
 |  
 |  set_axis(self, ax, sequence, track_height=0, nt_ticks=(20, 5), max_height=300, yticks=None, ylabels=None)
 |      Set up the plotting axis settings for an aesthetic arc plot.
 |      
 |      Sets the following properties of the given axis:
 |      1. spine positions
 |      2. x-axis and y-axis limits
 |      3. x-axis tick labels and positions according to `sequence`
 |      4. background boxes for x-axis tick labels
 |      5. y-axis tick labels and positions according to `track_height` and `ylabels`
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to set up.
 |      sequence : str
 |          Sequence to plot.
 |      track_height : int, optional
 |          Height of the track. If None, the height is automatically determined.
 |      nt_ticks : tuple of 2 ints, optional
 |          Major and minor tick spacing for the nucleotide axis. Default is (20, 5).
 |      max_height : int, optional
 |          Maximum height of the plot. Default is 300.
 |      yticks : list of ints, optional
 |          List of ytick positions. If None, the yticks are automatically determined.
 |      ylabels : list of str, optional
 |          List of ytick labels. If None, the ylabels are automatically determined.
 |  
 |  set_figure_size(self, height_ax_rel=0.03, width_ax_rel=0.03, width_ax_in=None, height_ax_in=None, height_gap_in=0.5, width_gap_in=0.5, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Set the figure size for an arc plot.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, Default is 0.03.
 |          Relative height of each axes.
 |      width_ax_rel : float, Default is 0.03.
 |          Relative width of each axes.
 |      width_ax_in : float, optional
 |          Absolute width of each axes in inches. If None, the width is automatically
 |          determined.
 |      height_ax_in : float, optional
 |          Absolute height of each axes in inches. If None, the height is
 |          automatically determined.
 |      height_gap_in : float, Default is 0.5.
 |          Absolute height of the gap between axes in inches.
 |      width_gap_in : float, Default is 0.5.
 |          Absolute width of the gap between axes in inches.
 |      top_in : float, Default is 1.
 |          Absolute top margin in inches.
 |      bottom_in : float, Default is 1.
 |          Absolute bottom margin in inches.
 |      left_in : float, Default is 1.
 |          Absolute left margin in inches.
 |      right_in : float, Default is 1.
 |          Absolute right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.Circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Circle in module rnavigate.plots.circle

class Circle(rnavigate.plots.plots.Plot)
 |  Circle(num_samples, **kwargs)
 |  
 |  Create a circle plot.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  **kwargs
 |      Keyword arguments passed to :class:`rnavigate.plots.Plot`.
 |  
 |  Attributes
 |  ----------
 |  zorder : dict
 |      Dictionary of zorder values for each plot element.
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  ax : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of current axes object. Increments with each call to
 |      :meth:`rnavigate.plots.Circle.plot_data`.
 |  
 |  Method resolution order:
 |      Circle
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **kwargs)
 |      Initialize Circle object.
 |  
 |  plot_data(self, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, annotations=None, label=None, colors=None, gap=30, nt_ticks=(20, 5))
 |      Plot data on the current (or specified) axes.
 |      
 |      Parameters
 |      ----------
 |      sequence : rnavigate.data.Sequence
 |          Sequence object to plot around the circle.
 |      structure : rnavigate.data.Structure, defaults to None
 |          Structure object to plot as lines within the circle.
 |      structure2 : rnavigate.data.Structure, defaults to None
 |          Structure object to compare to `structure`.
 |      interactions : rnavigate.data.Interactions, defaults to None
 |          Interactions object to plot as lines within the circle.
 |      interactions2 : rnavigate.data.Interactions, defaults to None
 |          Interactions object to plot as lines within the circle.
 |      profile : rnavigate.data.Profile, defaults to None
 |          Profile object used to color the sequence.
 |      annotations : list of rnavigate.data.Annotation, defaults to None
 |          List of Annotation objects to highlight regions around the circle.
 |      label : str, defaults to None
 |          Label for the plot title.
 |      colors : dict, defaults to None
 |          Dictionary of colors for each plot element.
 |      gap : float, defaults to 30
 |          Gap between the start and end of the sequence in degrees.
 |      nt_ticks : tuple of int, defaults to (20, 5)
 |          Gap between major and minor nucleotide ticks in degrees.
 |  
 |  set_axis(self, ax, label, seq_circle, gap, nt_ticks)
 |      Set axis limits and ticks.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to modify.
 |      label : str
 |          Label for the plot title.
 |      seq_circle : rnavigate.data.SequenceCircle
 |          SequenceCircle object.
 |      gap : float
 |          Gap between the start and end of the sequence in degrees.
 |      nt_ticks : tuple of 2 integers
 |          Gap between major and minor nucleotide ticks in degrees.
 |  
 |  set_figure_size(self, height_ax_rel=0.035, width_ax_rel=0.035, width_ax_in=None, height_ax_in=None, height_gap_in=1, width_gap_in=1, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, defaults to 0.035
 |          Height of axes relative to the y-axis limits.
 |      width_ax_rel : float, defaults to 0.035
 |          Width of axes relative to the x-axis limits.
 |      width_ax_in : float, defaults to None (overridden by width_ax_rel)
 |          Width of axes in inches.
 |      height_ax_in : float, defaults to None (overridden by height_ax_rel)
 |          Height of axes in inches.
 |      height_gap_in : float, defaults to 1
 |          Height of gap between axes in inches.
 |      width_gap_in : float, defaults to 1
 |          Width of gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Height of top margin in inches.
 |      bottom_in : float, defaults to 0.5
 |          Height of bottom margin in inches.
 |      left_in : float, defaults to 0.5
 |          Width of left margin in inches.
 |      right_in : float, defaults to 0.5
 |          Width of right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.DistHist

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class DistHist in module rnavigate.plots.disthist

class DistHist(rnavigate.plots.plots.Plot)
 |  DistHist(num_samples, **plot_kwargs)
 |  
 |  Create a distance histogram plot.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  **kwargs
 |      Keyword arguments passed to :class:`rnavigate.plots.Plot`.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  axes2 : dict
 |      Dictionary of twin axes objects. Keys are axes objects from
 |      :attr:`rnavigate.plots.DistHist.axes`.
 |  i : int
 |      Index of current axes object. Increments with each call to
 |      :meth:`rnavigate.plots.DistHist.plot_data`.
 |  
 |  Method resolution order:
 |      DistHist
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **plot_kwargs)
 |      Initialize DistHist object.
 |  
 |  plot_data(self, structure, interactions, bg_interactions, label, atom="O2'", ax=None)
 |      Plot data on the current (or specified) axes object.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
 |          Structure object to compute contact distances or 3D distances.
 |      interactions : rnavigate.data.Interactions
 |          Filtered Interactions object to to visualize pairwise distances.
 |      bg_interactions : rnavigate.data.Interactions
 |          Filtered Interactions object to to visualize pairwise background distances.
 |      label : str
 |          Label for the current axes object.
 |      atom : str, defaults to "O2'"
 |          Atom to compute distances from.
 |      ax : matplotlib.axes.Axes, defaults to None
 |          Axes object to plot on. If None, use the current axes object.
 |  
 |  plot_experimental_distances(self, ax, structure, interactions, atom, histtype='bar')
 |      Plot pairwise distances from the interactions object.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to plot on.
 |      structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
 |          Structure object to compute contact distances or 3D distances.
 |      interactions : rnavigate.data.Interactions
 |          Filtered Interactions object to to visualize pairwise distances.
 |      atom : str
 |          Atom to compute distances from.
 |      histtype : "bar" or "step", defaults to "bar"
 |          Type of histogram to plot.
 |  
 |  plot_structure_distances(self, ax, structure, atom)
 |      Plot all distances in the structure.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to plot on.
 |      structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
 |          Structure object to compute contact distances or 3D distances.
 |      atom : str
 |          Atom to compute distances from.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=1, width_gap_in=0.4, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Set figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, defaults to None
 |          Height of axes objects relative to y-axis limits.
 |      width_ax_rel : float, defaults to None
 |          Width of axes objects relative to x-axis limits.
 |      width_ax_in : float, defaults to 2
 |          Width of axes objects in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of axes objects in inches.
 |      height_gap_in : float, defaults to 1
 |          Height of gap between axes objects in inches.
 |      width_gap_in : float, defaults to 0.4
 |          Width of gap between axes objects in inches.
 |      top_in : float, defaults to 1
 |          Height of top margin in inches.
 |      bottom_in : float, defaults to 1
 |          Height of bottom margin in inches.
 |      left_in : float, defaults to 1
 |          Width of left margin in inches.
 |      right_in : float, defaults to 1
 |          Width of right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.Heatmap

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Heatmap in module rnavigate.plots.heatmap

class Heatmap(rnavigate.plots.plots.Plot)
 |  Heatmap(num_samples, structure, **plot_kwargs)
 |  
 |  Create a heatmap plot.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
 |      Structure object.
 |  **kwargs
 |      Keyword arguments passed to :class:`rnavigate.plots.Plot`.
 |  
 |  Attributes
 |  ----------
 |  structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
 |      Structure object.
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of current axes object. Increments with each call to
 |      :meth:`rnavigate.plots.Heatmap.plot_data`.
 |  
 |  Method resolution order:
 |      Heatmap
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, structure, **plot_kwargs)
 |      Initialize Heatmap object.
 |  
 |  plot_contour_distances(self, ax, levels, atom)
 |      Plot contour distances on the current (or specified) axes.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to plot on.
 |      levels : list of float, defaults to None
 |          Levels for contour lines. If None, levels are set automatically.
 |          [5] for SecondaryStructure, [20] for PDB.
 |      atom : str, defaults to "O2'"
 |          Atom to use for calculating distances.
 |  
 |  plot_contour_regions(self, ax, interactions, regions)
 |      Plot contour regions on the current (or specified) axes.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to plot on.
 |      interactions : rnavigate.data.Interactions
 |          Interactions object containing inter-nucleotide data to plot as heatmap.
 |      regions : list of tuple of int
 |          Regions to plot as contour lines.
 |  
 |  plot_data(self, interactions, label, levels=None, regions=None, interpolation=None, atom="O2'", plot_type='heatmap', weights=None)
 |      Plot heatmap and contour data on the current (or specified) axes.
 |      
 |      Parameters
 |      ----------
 |      interactions : rnavigate.data.Interactions
 |          Interactions object containing inter-nucleotide data to plot as heatmap.
 |      label : str
 |          Label for plot title.
 |      levels : list of float, defaults to None
 |          Levels for contour lines. If None, levels are set automatically.
 |          [5] for SecondaryStructure, [20] for PDB.
 |      regions : list of tuple of int, defaults to None
 |          Regions to plot as contour lines. If None, distances are plotted.
 |      interpolation : str, defaults to None
 |          Interpolation method for heatmap. If None, no interpolation is used.
 |      atom : str, defaults to "O2'"
 |          Atom to use for calculating distances.
 |      plot_type : "heatmap" or "kde", defaults to "heatmap"
 |          Type of plot to create.
 |      weights : array-like, defaults to None
 |          A column name of `interactions` to use as weights for KDE plot.
 |  
 |  plot_heatmap_data(self, ax, interactions, interpolation)
 |      Plot heatmap data on the specified axes.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to plot on.
 |      interactions : rnavigate.data.Interactions
 |          Interactions object containing inter-nucleotide data to plot as heatmap.
 |      interpolation : str, defaults to None
 |          Interpolation method for heatmap. If None, no interpolation is used.
 |  
 |  plot_kde_data(self, ax, interactions, weights=None, **kwargs)
 |      Plot KDE data on the specified axes.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object to plot on.
 |      interactions : rnavigate.data.Interactions
 |          Interactions object containing inter-nucleotide data to plot as KDE.
 |      weights : array-like, defaults to None
 |          A column name of `interactions` to use as weights for KDE plot.
 |      **kwargs
 |          Keyword arguments passed to :func:`seaborn.kdeplot`.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, defaults to None
 |          Height of each axes object relative to y-axis limits.
 |      width_ax_rel : float, defaults to None
 |          Width of each axes object relative to x-axis limits.
 |      width_ax_in : float, defaults to 2
 |          Width of each axes object in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of each axes object in inches.
 |      height_gap_in : float, defaults to 1
 |          Height of gap between axes objects in inches.
 |      width_gap_in : float, defaults to 0.5
 |          Width of gap between axes objects in inches.
 |      top_in : float, defaults to 1
 |          Top margin in inches.
 |      bottom_in : float, defaults to 0.5
 |          Bottom margin in inches.
 |      left_in : float, defaults to 0.5
 |          Left margin in inches.
 |      right_in : float, defaults to 0.5
 |          Right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.LinReg

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class LinReg in module rnavigate.plots.linreg

class LinReg(rnavigate.plots.plots.Plot)
 |  LinReg(num_samples, scale='linear', regression='pearson', kde=False, region='all')
 |  
 |  Plot a linear regression scatter plot, pairwise, between profiles.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  scale : str, optional
 |      Scale of the plot, either 'linear' or 'log'. The default is 'linear'.
 |  regression : 'pearson' or 'spearman', Defaults to 'pearson'
 |      Type of regression to perform.
 |  kde : bool, optional
 |      Whether to plot a kernel density estimate instead of a scatter plot.
 |      The default is False.
 |  region : 'all' or tuple of 2 integers, defaults to 'all'
 |      Start and end positions of the region to plot. If 'all', plot the
 |      entire profile.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure of the plot.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Axes of the plot.
 |  length : int
 |      Number of samples to plot.
 |  lims : list of 2 floats
 |      Limits of the plot applied to all x and y axes.
 |  profiles : list of numpy.ndarray
 |      Each sample's per-nucleotide values over the region and column of interest.
 |  colors : list of numpy.ndarray
 |      A color for each nucleotide in the region applied to the scatter plot.
 |  labels : list of str
 |      A label for each sample.
 |  scale : 'linear' or 'log'
 |      Scale of the plot axes.
 |  kde : bool
 |      Whether to plot a kernel density estimate instead of a scatter plot.
 |  regression : function
 |      Regression function to use.
 |  region : tuple of 2 integers
 |      Start and end positions of the region to plot.
 |  
 |  Method resolution order:
 |      LinReg
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, scale='linear', regression='pearson', kde=False, region='all')
 |      Initialize the linear regression plot.
 |  
 |  finalize(self)
 |      Finalize the plot by formatting the axes and adding labels.
 |  
 |  plot_data(self, structure, profile, annotations, label, column=None, colors='sequence')
 |      Add profile data. If all samples have been added, plot the regression.
 |      
 |      Parameters
 |      ----------
 |      structure : Structure
 |          Structure object.
 |      profile : Profile
 |          Profile object.
 |      annotations : Annotations
 |          Annotations object.
 |      label : str
 |          Label for the sample.
 |      column : str, optional
 |          Column of the profile to plot. The default is None.
 |      colors : str, optional
 |          Color scheme to use. The default is 'sequence'.
 |  
 |  plot_regression(self, i, j)
 |      Plot a linear regression between two profiles.
 |      
 |      Profiles must already have been added using `plot_data`.
 |      
 |      Parameters
 |      ----------
 |      i : int
 |          Index of the first profile.
 |      j : int
 |          Index of the second profile.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=0.3, width_gap_in=0.3, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, defaults to None
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, defaults to None
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, defaults to 2
 |          Width of the axes in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 0.3
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 0.3
 |          Width of the gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Top margin of the figure in inches.
 |      bottom_in : float, defaults to 0.5
 |          Bottom margin of the figure in inches.
 |      left_in : float, defaults to 0.5
 |          Left margin of the figure in inches.
 |      right_in : float, defaults to 0.5
 |          Right margin of the figure in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.NucleotideDistribution

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class NucleotideDistribution in module rnavigate.plots.ntdist

class NucleotideDistribution(rnavigate.plots.plots.Plot)
 |  NucleotideDistribution(num_samples, **plot_kwargs)
 |  
 |  Plot nucleotide distribution of a profile.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  sharex : bool, optional
 |      Whether to share the x-axis between plots.
 |  cols : int, optional
 |      Number of columns in the plot.
 |  **plot_kwargs
 |      Keyword arguments passed to the plot function.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      NucleotideDistribution
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **plot_kwargs)
 |      Initialize the plot.
 |  
 |  plot_data(self, profile, label, column=None, normalize=None, ax=None)
 |      Plot data to the current (or specified) axes.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      label : str
 |          label for the y-axis.
 |      column : str, optional
 |          Column of `profile.data` to plot.
 |      normalize : dict, optional
 |          Keyword arguments passed to `profile.normalize`.
 |      ax : matplotlib.axes.Axes, optional
 |          Axes object to plot to. If not specified, the current axes is used.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=0.2, width_gap_in=0.4, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, optional
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, defaults to 2
 |          Width of the axes in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 0.2
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 0.4
 |          Width of the gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Top margin in inches.
 |      bottom_in : float, defaults to 1
 |          Bottom margin in inches.
 |      left_in : float, defaults to 1
 |          Left margin in inches.
 |      right_in : float, defaults to 1
 |          Right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.Mol

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Mol in module rnavigate.plots.mol

class Mol(rnavigate.plots.plots.Plot)
 |  Mol(num_samples, pdb, width=400, height=400, background_alpha=1, rotation=None, orientation=None, style='cartoon', rows=None, cols=None)
 |  
 |  Method resolution order:
 |      Mol
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, pdb, width=400, height=400, background_alpha=1, rotation=None, orientation=None, style='cartoon', rows=None, cols=None)
 |      Create a 3Dmol.js viewer with a grid of subviewers
 |      
 |      Parameters
 |      ----------
 |      num_samples : int
 |          Number of subviewers to create
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for plotting
 |      width : int, defaults to 400
 |          Width of each subviewer in pixels
 |      height : int, defaults to 400
 |          Height of each subviewer in pixels
 |      background_alpha : float, defaults to 1
 |          Alpha value for the background color
 |      rotation : dict, defaults to None
 |          Dictionary of rotation angles for the viewer
 |      orientation : list of floats
 |          List of floats defining the orientation of the viewer
 |      style : str, defaults to "cartoon"
 |          Style of the viewer
 |      rows : int, defaults to None
 |          Number of rows to use in the grid of subviewers. If None, the
 |          number of rows is determined automatically.
 |      cols : int, defaults to None
 |          Number of columns to use in the grid of subviewers. If None, the
 |          number of columns is determined automatically.
 |      
 |      Attributes
 |      ----------
 |      view : py3Dmol.view
 |          3Dmol.js viewer object
 |      i : int
 |          Index of the current subviewer
 |      rows : int
 |          Number of rows in the grid of subviewers
 |      columns : int
 |          Number of columns in the grid of subviewers
 |      colorbars : list of matplotlib.colorbar.ColorbarBase
 |          List of colorbars to be added to the plot
 |      style : str
 |          Style of the viewer
 |      pdb : rnavigate.pdb.PDB
 |          PDB object to use for plotting
 |      length : int
 |          Number of samples
 |  
 |  add_lines(self, i, j, color, viewer, atom)
 |      Add lines between nucleotides i and j
 |      
 |      Parameters
 |      ----------
 |      i : int
 |          Index of the first nucleotide
 |      j : int
 |          Index of the second nucleotide
 |      color : str
 |          Color to use for the line
 |      viewer : tuple of ints
 |          Tuple of ints defining the subviewer to plot on
 |      atom : str
 |          Atom to use for plotting interactions
 |  
 |  get_orientation(self)
 |      Adds a clickable event to the viewer to display the orientation vector.
 |  
 |  get_viewer(self, i=None)
 |      Get the subviewer at index i
 |  
 |  hide_cylinders(self)
 |      Hide the cylinders that represent nucleotides.
 |  
 |  plot_data(self, interactions=None, profile=None, label=None, colors='grey', atom="O2'", title=True, get_orientation=False, viewer=None)
 |      Plot data on the current subviewer.
 |      
 |      Parameters
 |      ----------
 |      interactions : rnavigate.interactions.Interactions
 |          Interactions object to plot as lines on the 3d structure
 |      profile : rnavigate.profile.Profile
 |          Profile object to use as nucleotide colors
 |      label : str, defaults to None
 |          Label to use as a title on the subviewer
 |      colors : str, defaults to "grey"
 |          Color scheme to use for nucleotides
 |      atom : str, defaults to "O2'"
 |          Atom to use for plotting interactions
 |      title : bool, defaults to True
 |          Whether or not to add a title to the subviewer
 |      get_orientation : bool, defaults to False
 |          Whether or not to get the orientation of the subviewer
 |          This will display the orientation as a label on the subviewer when the
 |          structure is clicked.
 |      viewer : tuple of ints, defaults to None
 |          Tuple of ints defining the subviewer to plot on
 |  
 |  plot_interactions(self, viewer, interactions, atom)
 |      Plot interactions on the current subviewer
 |      
 |      Parameters
 |      ----------
 |      viewer : tuple of ints
 |          Tuple of ints defining the subviewer to plot on
 |      interactions : rnavigate.interactions.Interactions
 |          Interactions object to plot as lines on the 3d structure
 |      atom : str
 |          Atom to use for plotting interactions
 |  
 |  save(self)
 |      Display the current orientation of the viewer as a png image.
 |      
 |      Notes
 |      -----
 |      This method must be run in a new cell after the viewer has been
 |      instantiated. The resulting png image will be saveable as a png file by
 |      clicking and dragging the image to your desktop.
 |  
 |  set_colors(self, viewer, profile, colors)
 |      Set the colors of the nucleotides on the current subviewer
 |      
 |      Parameters
 |      ----------
 |      viewer : tuple of ints
 |          Tuple of ints defining the subviewer to plot on
 |      profile : rnavigate.profile.Profile
 |          Profile object to use as nucleotide colors
 |      colors : str
 |          Color scheme to use for nucleotides
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=None, height_ax_in=None, height_gap_in=None, width_gap_in=None, top_in=None, bottom_in=None, left_in=None, right_in=None)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, optional
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, optional
 |          Width of the axes in inches.
 |      height_ax_in : float, optional
 |          Height of the axes in inches.
 |      height_gap_in : float, optional
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, optional
 |          Width of the gap between axes in inches.
 |      top_in : float, optional
 |          Top margin in inches.
 |      bottom_in : float, optional
 |          Bottom margin in inches.
 |      left_in : float, optional
 |          Left margin in inches.
 |      right_in : float, optional
 |          Right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.Plot

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Plot in module rnavigate.plots.plots

class Plot(abc.ABC)
 |  Plot(num_samples, rows=None, cols=None, **kwargs)
 |  
 |  Abstract base class for plots.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  rows : int, optional
 |      Number of rows in the plot.
 |  cols : int, optional
 |      Number of columns in the plot.
 |  **kwargs
 |      Keyword arguments passed to matplotlib.pyplot.subplots.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, rows=None, cols=None, **kwargs)
 |      Initialize the plot.
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  plot_data(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=None, height_ax_in=None, height_gap_in=None, width_gap_in=None, top_in=None, bottom_in=None, left_in=None, right_in=None)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, optional
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, optional
 |          Width of the axes in inches.
 |      height_ax_in : float, optional
 |          Height of the axes in inches.
 |      height_gap_in : float, optional
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, optional
 |          Width of the gap between axes in inches.
 |      top_in : float, optional
 |          Top margin in inches.
 |      bottom_in : float, optional
 |          Bottom margin in inches.
 |      left_in : float, optional
 |          Left margin in inches.
 |      right_in : float, optional
 |          Right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset({'plot_data'})
```

### rnavigate.plots.ColorBar

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class ColorBar in module rnavigate.plots.plots

class ColorBar(Plot)
 |  ColorBar(num_samples, rows=None, cols=None, **kwargs)
 |  
 |  Plot a colorbar.
 |  
 |  Parameters
 |  ----------
 |  rows : int, optional
 |      Number of rows in the plot.
 |  **kwargs
 |      Keyword arguments passed to matplotlib.pyplot.subplots.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      ColorBar
 |      Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  plot_data(self, colorbar)
 |      Add a colorbar to the current axes.
 |      
 |      Parameters
 |      ----------
 |      colorbar : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=3, height_ax_in=0.1, height_gap_in=0.75, width_gap_in=0.5, top_in=None, bottom_in=None, left_in=None, right_in=None)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, optional
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, defaults to 3
 |          Width of the axes in inches.
 |      height_ax_in : float, defaults to 0.1
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 0.75
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 0.5
 |          Width of the gap between axes in inches.
 |      top_in : float, optional
 |          Top margin in inches.
 |      bottom_in : float, optional
 |          Bottom margin in inches.
 |      left_in : float, optional
 |          Left margin in inches.
 |      right_in : float, optional
 |          Right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Plot:
 |  
 |  __init__(self, num_samples, rows=None, cols=None, **kwargs)
 |      Initialize the plot.
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.Profile

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Profile in module rnavigate.plots.profile

class Profile(rnavigate.plots.plots.Plot)
 |  Profile(num_samples, nt_length, region='all', **kwargs)
 |  
 |  Plot per-nucleotide measurements as colored bars.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  nt_length : int
 |      Length of the nucleotide sequence.
 |  region : tuple of int, optional
 |      Region of the nucleotide sequence to plot.
 |  **kwargs
 |      Keyword arguments passed to the plot function.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      Profile
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, nt_length, region='all', **kwargs)
 |      Initialize the plot.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot.
 |      cols : int, optional
 |          Number of columns in the plot. This is ignored and set to 1.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_data(self, profile, annotations, domains, label, plot_error=True, column=None, seqbar=True, annotations_mode='track', nt_ticks=(20, 5))
 |      Plot data to the current (or specified) axes.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      annotations : list of rnavigate.data.Annotation
 |          List of annotation objects to plot along the sequence.
 |      domains : list of rnavigate.data.Annotation, optional
 |          List of domains to plot along the sequence.
 |      label : str
 |          label for the y-axis.
 |      plot_error : bool, optional
 |          Whether to plot the error bars.
 |      column : str, optional
 |          Column of `profile` to plot.
 |      seqbar : bool, optional
 |          Whether to plot the sequence track.
 |      annotations_mode : "track" or "bar", defaults to "track"
 |          Mode of the annotations track.
 |      nt_ticks : tuple of int, optional
 |          Major and minor tick interval for the nucleotide axis.
 |  
 |  set_axis(self, ax, sequence, nt_ticks=(20, 5))
 |      Set up axis properties for aesthetics.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      sequence : str
 |          Nucleotide sequence.
 |      nt_ticks : tuple of int, optional
 |          Major and minor tick interval for the nucleotide axis.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=0.03, width_ax_in=None, height_ax_in=2, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Set the size of the figure.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, defaults to 0.03
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, optional
 |          Width of the axes in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 1
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 0.5
 |          Width of the gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Height of the top margin in inches.
 |      bottom_in : float, defaults to 1
 |          Height of the bottom margin in inches.
 |      left_in : float, defaults to 1
 |          Width of the left margin in inches.
 |      right_in : float, defaults to 1
 |          Width of the right margin in inches.
 |  
 |  set_labels(self, ax, axis_title='Reactivity Profile', xlabel='Nucleotide Position', ylabel='Reactivity')
 |      Set the labels of the plot.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      axis_title : str, optional
 |          Title of the axis.
 |      xlabel : str, optional
 |          Label of the x-axis.
 |      ylabel : str, optional
 |          Label of the y-axis.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.QC

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class QC in module rnavigate.plots.qc

class QC(rnavigate.plots.plots.Plot)
 |  QC(num_samples)
 |  
 |  Plot QC data from log files.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  ax_muts_unt : matplotlib.axes.Axes
 |      Axes object for the mutations per molecule plot of untreated samples.
 |  ax_muts_mod : matplotlib.axes.Axes
 |      Axes object for the mutations per molecule plot of modified samples.
 |  ax_read_unt : matplotlib.axes.Axes
 |      Axes object for the read length plot of untreated samples.
 |  ax_read_mod : matplotlib.axes.Axes
 |      Axes object for the read length plot of modified samples.
 |  ax_boxplot : matplotlib.axes.Axes
 |      Axes object for the boxplot of mutation rates.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      QC
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples)
 |      Initialize the plot.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows. This is ignored.
 |      cols : int, optional
 |          Number of columns. This is ignored.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows. 1 if there is only one sample, 2 otherwise.
 |      cols : int
 |          Number of columns. 3 if there is only one sample, 4 otherwise.
 |  
 |  make_boxplot(self, labels)
 |      Make the boxplot of mutation rates.
 |      
 |      Parameters
 |      ----------
 |      labels : list of str
 |          Labels for the samples.
 |  
 |  plot_data(self, profile, label)
 |      Plot the data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      label : str
 |          Label for the sample.
 |  
 |  plot_mutations_per_molecule(self, profile, label, upper_limit=12)
 |      Plot the mutations per molecule.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      label : str
 |          Label for the sample.
 |      upper_limit : int, optional
 |          Upper limit of the x-axis.
 |  
 |  plot_read_lengths(self, profile, label, upper_limit=12)
 |      Plot the read lengths.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      label : str
 |          Label for the sample.
 |      upper_limit : int, optional
 |          Upper limit of the x-axis.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=1, width_gap_in=1, top_in=1, bottom_in=0.5, left_in=1, right_in=0.5)
 |      Set the size of the figure.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, optional
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, optional
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, defaults to 2
 |          Width of the axes in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 1
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 1
 |          Width of the gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Height of the top margin in inches.
 |      bottom_in : float, defaults to 0.5
 |          Height of the bottom margin in inches.
 |      left_in : float, defaults to 0.5
 |          Width of the left margin in inches.
 |      right_in : float, defaults to 0.5
 |          Width of the right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.ROC

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class ROC in module rnavigate.plots.roc

class ROC(rnavigate.plots.plots.Plot)
 |  ROC(num_samples, **kwargs)
 |  
 |  Plot ROC curves.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  **kwargs
 |      Keyword arguments passed to `rnavigate.plots.Plot`.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  a_ax : matplotlib.axes.Axes
 |      Axes object for the AUC plot of A nucleotides.
 |  u_ax : matplotlib.axes.Axes
 |      Axes object for the AUC plot of U nucleotides.
 |  g_ax : matplotlib.axes.Axes
 |      Axes object for the AUC plot of G nucleotides.
 |  c_ax : matplotlib.axes.Axes
 |      Axes object for the AUC plot of C nucleotides.
 |  main_ax : matplotlib.axes.Axes
 |      Axes object for the main ROC plot.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      ROC
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **kwargs)
 |      Initialize the plot.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns.
 |      
 |      Parameters
 |      ----------
 |      rows : int
 |          Number of rows. This is ignored.
 |      cols : int
 |          Number of columns. This is ignored.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows. This is always 2.
 |      cols : int
 |          Number of columns. This is always 4.
 |  
 |  plot_data(self, structure, profile, label, nts='AUCG')
 |      Plot the data.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.structure.Structure
 |          Structure object.
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      label : str
 |          Sample name.
 |      nts : str, defaults to "AUCG"
 |          Which nucleotides to plot.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=None, width_ax_in=1.5, height_ax_in=1.5, height_gap_in=0.3, width_gap_in=0.2, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float
 |          Width of the axes in inches.
 |      height_ax_in : float
 |          Height of the axes in inches.
 |      height_gap_in : float
 |          Height of the gap between axes in inches.
 |      width_gap_in : float
 |          Width of the gap between axes in inches.
 |      top_in : float
 |          Height of the top margin in inches.
 |      bottom_in : float
 |          Height of the bottom margin in inches.
 |      left_in : float
 |          Width of the left margin in inches.
 |      right_in : float
 |          Width of the right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.Skyline

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Skyline in module rnavigate.plots.skyline

class Skyline(rnavigate.plots.plots.Plot)
 |  Skyline(num_samples, nt_length, region='all', **kwargs)
 |  
 |  Plot per-nucleotide measurements as stepped line graphs (skyline plots).
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  nt_length : int
 |      Length of the nucleotide sequence.
 |  region : tuple of int, defaults to "all" (entire sequence)
 |      start and end position of the region to plot. If "all", plot the entire
 |      sequence.
 |  **kwargs
 |      Keyword arguments passed to `rnavigate.plots.Plot`.
 |  
 |  Attributes
 |  ----------
 |  nt_length : int
 |      Length of the nucleotide sequence.
 |  region : tuple of int
 |      start and end position of the region to plot.
 |  track_height : float
 |      Height of the tracks in the plot.
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  ax : matplotlib.axes.Axes
 |      Axes object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      Skyline
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, nt_length, region='all', **kwargs)
 |      Initialize the plot.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object. This is ignored.
 |      
 |      Returns
 |      -------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns.
 |      
 |      Parameters
 |      ----------
 |      rows : int
 |          Number of rows. This is ignored.
 |      cols : int
 |          Number of columns. This is ignored.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows. This is always 1.
 |      cols : int
 |          Number of columns. This is always 1.
 |  
 |  plot_data(self, profile, annotations=None, domains=None, label=None, columns=None, seqbar=True, errors=None, annotations_mode='track', nt_ticks=(20, 5))
 |      Add data to the axes.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      annotations : list of rnavigate.annotation.Annotation, optional
 |          List of annotation objects.
 |      domains : list of rnavigate.domain.Domain, optional
 |          List of domain objects.
 |      label : str, optional
 |          Sample name.
 |      columns : str or list of str, optional
 |          Which columns to plot. If None, plot the metric column.
 |      seqbar : bool, defaults to True
 |          Whether to plot a sequence bar.
 |      errors : str, optional
 |          Which error columns to plot. If None, do not plot errors.
 |      annotations_mode : "track" or "bar", defaults to "track"
 |          Whether to plot annotations as a track or as vertical bars.
 |      nt_ticks : tuple of int, defaults to (20, 5)
 |          Major and minor tick frequency for nucleotide positions.
 |  
 |  set_axis(self, ax, sequence, nt_ticks)
 |      Set the axis limits and ticks.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      sequence : rnavigate.data.Sequence
 |          The sequence on which position labels are based. Dashes are ignored.
 |      nt_ticks : tuple of int
 |          Major and minor tick frequency for nucleotide positions.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=0.03, width_ax_in=None, height_ax_in=2, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float
 |          Width of the axes in inches.
 |      height_ax_in : float
 |          Height of the axes in inches.
 |      height_gap_in : float
 |          Height of the gap between axes in inches.
 |      width_gap_in : float
 |          Width of the gap between axes in inches.
 |      top_in : float
 |          Height of the top margin in inches.
 |      bottom_in : float
 |          Height of the bottom margin in inches.
 |      left_in : float
 |          Width of the left margin in inches.
 |      right_in : float
 |          Width of the right margin in inches.
 |  
 |  set_labels(self, ax, axis_title='Raw Reactivity Profile', legend_title='Samples', xlabel='Nucleotide Position', ylabel='Profile')
 |      Set the axis labels and legend.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      axis_title : str, optional
 |          Title of the axes.
 |      legend_title : str, optional
 |          Title of the legend.
 |      xlabel : str, optional
 |          Label of the x-axis.
 |      ylabel : str, optional
 |          Label of the y-axis.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.SM

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SM in module rnavigate.plots.sm

class SM(rnavigate.plots.plots.Plot)
 |  SM(nt_length, region=None, panels=['profile', 'rates', 'depth'])
 |  
 |  Plot classic ShapeMapper plots.
 |  
 |  Parameters
 |  ----------
 |  nt_length : int
 |      Length of the nucleotide sequence.
 |  region : tuple of int, defaults to "all" (entire sequence)
 |      start and end position of the region to plot. If "all", plot the entire
 |      sequence.
 |  panels : list of str, defaults to ["profile", "rates", "depth"]
 |      Which panels to plot. Options are "profile", "rates", and "depth".
 |  
 |  Attributes
 |  ----------
 |  nt_length : int
 |      Length of the nucleotide sequence.
 |  region : tuple of int
 |      start and end position of the region to plot.
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      SM
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, nt_length, region=None, panels=['profile', 'rates', 'depth'])
 |      Initialize the plot.
 |  
 |  metric_abbreviate(self, num)
 |      takes a large number and applies an appropriate abbreviation
 |      
 |      Parameters
 |      ----------
 |          num : int
 |              number to be abbreviated
 |      
 |      Returns
 |      -------
 |          str
 |              abbreviated number
 |  
 |  plot_data(self, profile, label)
 |      Plot the data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |      label : str
 |          Label for the plot.
 |  
 |  plot_sm_depth(self, ax, profile)
 |      Plots classic ShapeMapper read depth on the given axis.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |  
 |  plot_sm_profile(self, ax, profile)
 |      Plots classic ShapeMapper profile on the given axis.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |  
 |  plot_sm_rates(self, ax, profile)
 |      Plots classic ShapeMapper mutation rates on the given axis.
 |      
 |      Parameters
 |      ----------
 |      ax : matplotlib.axes.Axes
 |          Axes object.
 |      profile : rnavigate.profile.Profile
 |          Profile object.
 |  
 |  set_figure_size(self, height_ax_rel=None, width_ax_rel=0.03, width_ax_in=None, height_ax_in=2, height_gap_in=0.5, width_gap_in=0.5, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, defaults to None
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, defaults to 0.03
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float, defaults to None
 |          Width of the axes in inches.
 |      height_ax_in : float, defaults to 2
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 0.5
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 0.5
 |          Width of the gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Height of the top margin in inches.
 |      bottom_in : float, defaults to 0.5
 |          Height of the bottom margin in inches.
 |      left_in : float, defaults to 0.5
 |          Width of the left margin in inches.
 |      right_in : float, defaults to 0.5
 |          Width of the right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.SS

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class SS in module rnavigate.plots.ss

class SS(rnavigate.plots.plots.Plot)
 |  SS(num_samples, **kwargs)
 |  
 |  Plot secondary structure diagrams.
 |  
 |  Parameters
 |  ----------
 |  num_samples : int
 |      Number of samples to plot.
 |  **kwargs
 |      Keyword arguments passed to `rnavigate.plots.Plot`.
 |  
 |  Attributes
 |  ----------
 |  fig : matplotlib.figure.Figure
 |      Figure object.
 |  axes : numpy.ndarray of matplotlib.axes.Axes
 |      Array of axes objects.
 |  xlims : list of float
 |      x limits of the plot.
 |  ylims : list of float
 |      y limits of the plot.
 |  i : int
 |      Index of the current plot.
 |  
 |  Method resolution order:
 |      SS
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **kwargs)
 |      Initialize the plot.
 |  
 |  plot_data(self, structure, interactions=None, interactions2=None, profile=None, annotations=None, label='', colors=None, nt_ticks=None, bp_style='dotted')
 |      Plot the data on the current axes.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          Structure object with diagram drawing coordinates.
 |      interactions : rnavigate.data.Interactions, optional
 |          Interactions object to plot as lines between nucleotides.
 |      interactions2 : rnavigate.data.Interactions, optional
 |          Interactions object to plot as lines between nucleotides.
 |      profile : rnavigate.data.Profile, optional
 |          Profile object used to color nucleotides.
 |      annotations : list of rnavigate.data.Annotation, optional
 |          Annotation objects to highlight regions or nucleotides of interest.
 |      label : str, optional
 |          Label for the plot title.
 |      colors : dict, optional
 |          Dictionary of colors for each plot element. Keys are "sequence",
 |          "nucleotides", "structure", and "basepairs". Values are either
 |          matplotlib colors or strings specifying the color scheme.
 |      nt_ticks : int, optional
 |          Number of nucleotides between tick marks.
 |      bp_style : "dotted", "solid", or "conventional", defaults to "dotted"
 |          Style of base pair lines.
 |  
 |  set_figure_size(self, height_ax_rel=0.2, width_ax_rel=0.2, width_ax_in=None, height_ax_in=None, height_gap_in=0.5, width_gap_in=0.2, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Set the figure size.
 |      
 |      Parameters
 |      ----------
 |      height_ax_rel : float, defaults to 0.2
 |          Height of the axes relative to the y-axis limits.
 |      width_ax_rel : float, defaults to 0.2
 |          Width of the axes relative to the x-axis limits.
 |      width_ax_in : float
 |          Width of the axes in inches.
 |      height_ax_in : float
 |          Height of the axes in inches.
 |      height_gap_in : float, defaults to 0.5
 |          Height of the gap between axes in inches.
 |      width_gap_in : float, defaults to 0.2
 |          Width of the gap between axes in inches.
 |      top_in : float, defaults to 1
 |          Height of the top margin in inches.
 |      bottom_in : float, defaults to 0.5
 |          Height of the bottom margin in inches.
 |      left_in : float, defaults to 0.5
 |          Width of the left margin in inches.
 |      right_in : float, defaults to 0.5
 |          Width of the right margin in inches.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __abstractmethods__ = frozenset()
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.plots.plots.Plot:
 |  
 |  add_colorbar_args(self, cmap)
 |      Add colorbar arguments to the plot.
 |      
 |      Parameters
 |      ----------
 |      cmap : rnavigate.data.ScalarMappable
 |          Colormap object.
 |  
 |  get_ax(self, i=None)
 |      Get the current axes object.
 |      
 |      Parameters
 |      ----------
 |      i : int, optional
 |          Index of the axes object to return. If None, return the current
 |          axes object.
 |      
 |      Returns
 |      -------
 |      matplotlib.axes.Axes
 |          Axes object.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |      Get the number of rows and columns in the plot.
 |      
 |      Parameters
 |      ----------
 |      rows : int, optional
 |          Number of rows in the plot. If None, the number of rows is
 |          determined automatically.
 |      cols : int, optional
 |          Number of columns in the plot. If None, the number of columns is
 |          determined automatically.
 |      
 |      Returns
 |      -------
 |      rows : int
 |          Number of rows in the plot.
 |      cols : int
 |          Number of columns in the plot.
 |  
 |  plot_colorbars(self)
 |      Plot colorbars.
 |      
 |      Returns
 |      -------
 |      ColorBar
 |          ColorBar plot object.
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Parameters
 |      ----------
 |          filename : string
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.plots.plots.Plot:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.plots.get_contrasting_colors

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_contrasting_colors in module rnavigate.plots.functions.functions

get_contrasting_colors(colors)
    Get contrasting colors for a list of colors.
    
    Returns a list of "k" (black) or "w" (white) for each color in the input list,
    which ever contrasts better with the input color.
    
    Parameters
    ----------
    colors : list of str
        List of colors to get contrasting colors for.
    
    Returns
    -------
    list of str
        List of "k" or "w" for each color in the input list.
```

### rnavigate.plots.adjust_spines

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function adjust_spines in module rnavigate.plots.functions.functions

adjust_spines(ax, spines)
    Places the given spines on the given axis, removes others.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to adjust the spines of.
    spines : list of str
        The spines to adjust. Valid options are "left", "right", "top", and "bottom".
```

### rnavigate.plots.clip_spines

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function clip_spines in module rnavigate.plots.functions.functions

clip_spines(ax, spines)
    Clips the given spines to the range of the ticks.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to clip the spines of.
    spines : list of str
        The spines to clip. Valid options are "left", "right", "top", and "bottom".
```

### rnavigate.plots.get_nt_ticks

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_nt_ticks in module rnavigate.plots.functions.functions

get_nt_ticks(sequence, region, gap)
    Get the ticks and labels for a region of a sequence.
    
    Dashes are ignored when counting nucleotide positions.
    
    Parameters
    ----------
    sequence : str or rnavigate.data.Sequence
        The sequence to get ticks and labels for.
    region : tuple of int
        The region of the sequence to get ticks and labels for.
    gap : int
        The gap between major ticks.
    
    Returns
    -------
    ticks, labels : tuple of list of int
        The ticks and labels for the given region of the sequence.
```

### rnavigate.plots.set_nt_ticks

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function set_nt_ticks in module rnavigate.plots.functions.functions

set_nt_ticks(ax, sequence, region, major, minor)
    Set the ticks and labels for a region of a sequence.
    
    Dashes are ignored when counting nucleotide positions.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to set the ticks and labels of.
    sequence : str or rnavigate.data.Sequence
        The sequence to set ticks and labels for.
    region : tuple of int
        The region of the sequence to set ticks and labels for.
    major : int
        The gap between major ticks.
    minor : int
        The gap between minor ticks.
```

### rnavigate.plots.box_xtick_labels

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function box_xtick_labels in module rnavigate.plots.functions.functions

box_xtick_labels(ax)
    Place a transparent box behind the xtick labels of the provided axis.
```

### rnavigate.plots.plot_interactions_arcs

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_interactions_arcs in module rnavigate.plots.functions.functions

plot_interactions_arcs(ax, interactions, panel, yvalue=0, region='all')
    Plot interactions as arcs.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the interactions on.
    interactions : rnavigate.data.Interactions
        The interactions to plot.
    panel : "top" or "bottom"
        The panel to plot the interactions on.
    yvalue : float, optional
        The y-value at which to plot the interactions.
    region : tuple of int, optional
        The region of the sequence to plot interactions for.
```

### rnavigate.plots.plot_profile_bars

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_profile_bars in module rnavigate.plots.functions.functions

plot_profile_bars(ax, profile, scale_factor=1, plot_error=True, bottom=0, region='all')
    Plot per-nucleotide data as colored bars.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the profile on.
    profile : rnavigate.data.Profile
        The profile to plot.
    scale_factor : float, optional
        The factor by which to multiply per-nucleotide values.
    plot_error : bool, optional
        Whether to plot error bars.
    bottom : float, optional
        The y-value at which to plot the bars.
    region : tuple of int, optional
        The region of the sequence to plot bars for.
```

### rnavigate.plots.plot_profile_skyline

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_profile_skyline in module rnavigate.plots.functions.functions

plot_profile_skyline(ax, profile, label, columns, errors)
    Plot per-nucleotide data as a skyline plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the profile on.
    profile : rnavigate.data.Profile
        The profile to plot.
    label : str
        The label to use for the plot legend.
    columns : list of str
        The columns of the profile to plot.
    errors : list of str
        The columns of the profile to use for error bars.
```

### rnavigate.plots.plot_sequence_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_sequence_alignment in module rnavigate.plots.functions.functions

plot_sequence_alignment(ax, alignment, labels, top=5, bottom=-5, ytrans='data')
    Plot a sequence alignment.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the alignment on.
    alignment : rnavigate.data.Alignment
        The alignment to plot.
    labels : tuple of str
        The labels for the two sequences in the alignment.
    top : int, optional
        The y-value to plot the top sequence at.
    bottom : int, optional
        The y-value to plot the bottom sequence at.
    ytrans : str, optional
        The transform to use for the y-values. Valid options are "data" and "axes".
```

### rnavigate.plots.plot_annotation_track

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_annotation_track in module rnavigate.plots.functions.tracks

plot_annotation_track(ax, annotation, yvalue, height, mode, region='all', ytrans='data')
    Plot an annotation track along the x-axis of a plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the annotation track on.
    annotation : rnavigate.data.Annotation
        The annotation to plot.
    yvalue : float
        The y-value of the annotation track.
    height : float
        The height of the annotation track.
    mode : "track" or "bar"
        The annotation mode.
    region : list of 2 int, defaults to "all"
        Start and end positions of the region to plot. If "all", plot the entire
        sequence.
    ytrans : "data" or "axes", defaults to "data"
        The y-axis coordinate system.
```

### rnavigate.plots.plot_domain_track

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_domain_track in module rnavigate.plots.functions.tracks

plot_domain_track(ax, spans, yvalue, height, region='all', ytrans='data')
    Plot a domain track along the x-axis of a plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the domain track on.
    spans : rnavigate.data.Spans
        The spans to plot.
    yvalue : float
        The y-value of the domain track.
    height : float
        The height of the domain track.
    region : list of 2 int, defaults to "all"
        Start and end positions of the region to plot. If "all", plot the entire
        sequence.
    ytrans : "data" or "axes", defaults to "data"
        The y-axis coordinate system.
```

### rnavigate.plots.plot_sequence_track

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_sequence_track in module rnavigate.plots.functions.tracks

plot_sequence_track(ax, sequence, yvalue=-0.05, height=0.05, ytrans='data', verticalalignment='bottom', region='all')
    Plot a sequence track along the x-axis of a plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the sequence track on.
    sequence : str
        The sequence to plot.
    yvalue : float, defaults to -0.05
        The y-value of the sequence track.
    height : float, defaults to 0.05
        The height of the sequence track.
    ytrans : "data" or "axes", defaults to "data"
        The y-axis coordinate system.
    verticalalignment : "top" or "bottom", defaults to "bottom"
        The vertical alignment of the sequence track.
    region : list of 2 int, defaults to "all"
        Start and end positions of the region to plot. If "all", plot the entire
        sequence.
```

### rnavigate.plots.plot_annotation_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_annotation_ss in module rnavigate.plots.functions.ss

plot_annotation_ss(ax, structure, annotation)
    Highlight regions or nucleotides of interest on a secondary structure diagram.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates.
    annotation : rnavigate.data.Annotation
        Annotation to plot.
```

### rnavigate.plots.plot_basepairs_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_basepairs_ss in module rnavigate.plots.functions.ss

plot_basepairs_ss(ax, structure, bp_style)
    Plot the basepairs of a secondary structure diagram.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    bp_style : "conventional", "dotted" or "line"
        Style of basepairs to plot.
```

### rnavigate.plots.plot_interactions_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_interactions_ss in module rnavigate.plots.functions.ss

plot_interactions_ss(ax, structure, interactions)
    Plot the interactions as lines over a secondary structure diagram.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates.
    interactions : rnavigate.data.Interactions
        Interactions to plot.
```

### rnavigate.plots.plot_nucleotides_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_nucleotides_ss in module rnavigate.plots.functions.ss

plot_nucleotides_ss(ax, structure, colors)
    Plot the nucleotides of a secondary structure diagram.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    colors : list
        List of colors to use for each nucleotide in the structure.
```

### rnavigate.plots.plot_positions_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_positions_ss in module rnavigate.plots.functions.ss

plot_positions_ss(ax, structure, xticks=20)
    Plot the positions of a secondary structure diagram.
    
    Label locations are chosen from a point on a circle around each position that is
    the furthest from any other nucleotides. This sometimes causes tick marks and
    labels to overlap with other plot elements.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    xticks : int
        Spacing between position labels.
```

### rnavigate.plots.plot_sequence_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_sequence_ss in module rnavigate.plots.functions.ss

plot_sequence_ss(ax, structure, colors)
    Plot the sequence of a secondary structure diagram.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    colors : list
        List of colors to use for each nucleotide in the structure.
```

### rnavigate.plots.plot_structure_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_structure_ss in module rnavigate.plots.functions.ss

plot_structure_ss(ax, structure, colors)
    Plot the structure of a secondary structure diagram.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    colors : list
        List of colors to use for each nucleotide in the structure.
```

### rnavigate.plots.plot_annotation_circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_annotation_circle in module rnavigate.plots.functions.circle

plot_annotation_circle(ax, seq_circle, annotation, offset=1)
    Plot annotations on a circle plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    seq_circle : rnavigate.data.SequenceCircle
        The sequence circle object containing nucleotide positions.
    annotation : rnavigate.data.Annotation
        The annotation to be plotted.
    offset : float, optional
        The offset from the circle circumference to plot the annotation.
```

### rnavigate.plots.plot_interactions_circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_interactions_circle in module rnavigate.plots.functions.circle

plot_interactions_circle(ax, seq_circle, interactions)
    Plot interactions on a circle plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    seq_circle : rnavigate.data.SequenceCircle
        The sequence circle object containing nucleotide positions.
    interactions : rnavigate.data.Interactions
        The interactions to be plotted as arcs between nucleotides.
```

## rnavigate.styles

### rnavigate.styles.update_copy

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function update_copy in module rnavigate.styles

update_copy(original_settings, user_settings)
    Recursively updates and returns a copy of og settings with new settings applied.
    
    Parameters
    ----------
    original_settings (dict)
        a default settings dictionary, usually rnav.settings
    user_settings (dict)
        a dictionary with only the fields from the original_settings that
        are to be changed
    
    Returns
    -------
    settings dict
        the original_settings dictionary with the new_settings dictionary
        values recursively applied
```

### rnavigate.styles.Settings

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Settings in module rnavigate.styles

class Settings(builtins.dict)
 |  Settings(user_settings)
 |  
 |  Context manager for temporarily changing global settings.
 |  
 |  Parameters
 |  ----------
 |  user_settings : dict
 |      a dictionary with only the fields from the original_settings that
 |      are to be changed
 |  
 |  Method resolution order:
 |      Settings
 |      builtins.dict
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __enter__(self)
 |  
 |  __exit__(self, *args, **kwargs)
 |  
 |  __init__(self, user_settings)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from builtins.dict:
 |  
 |  __contains__(self, key, /)
 |      True if the dictionary has the specified key, else False.
 |  
 |  __delitem__(self, key, /)
 |      Delete self[key].
 |  
 |  __eq__(self, value, /)
 |      Return self==value.
 |  
 |  __ge__(self, value, /)
 |      Return self>=value.
 |  
 |  __getattribute__(self, name, /)
 |      Return getattr(self, name).
 |  
 |  __getitem__(...)
 |      x.__getitem__(y) <==> x[y]
 |  
 |  __gt__(self, value, /)
 |      Return self>value.
 |  
 |  __ior__(self, value, /)
 |      Return self|=value.
 |  
 |  __iter__(self, /)
 |      Implement iter(self).
 |  
 |  __le__(self, value, /)
 |      Return self<=value.
 |  
 |  __len__(self, /)
 |      Return len(self).
 |  
 |  __lt__(self, value, /)
 |      Return self<value.
 |  
 |  __ne__(self, value, /)
 |      Return self!=value.
 |  
 |  __or__(self, value, /)
 |      Return self|value.
 |  
 |  __repr__(self, /)
 |      Return repr(self).
 |  
 |  __reversed__(self, /)
 |      Return a reverse iterator over the dict keys.
 |  
 |  __ror__(self, value, /)
 |      Return value|self.
 |  
 |  __setitem__(self, key, value, /)
 |      Set self[key] to value.
 |  
 |  __sizeof__(...)
 |      D.__sizeof__() -> size of D in memory, in bytes
 |  
 |  clear(...)
 |      D.clear() -> None.  Remove all items from D.
 |  
 |  copy(...)
 |      D.copy() -> a shallow copy of D
 |  
 |  get(self, key, default=None, /)
 |      Return the value for key if key is in the dictionary, else default.
 |  
 |  items(...)
 |      D.items() -> a set-like object providing a view on D's items
 |  
 |  keys(...)
 |      D.keys() -> a set-like object providing a view on D's keys
 |  
 |  pop(...)
 |      D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
 |      
 |      If key is not found, default is returned if given, otherwise KeyError is raised
 |  
 |  popitem(self, /)
 |      Remove and return a (key, value) pair as a 2-tuple.
 |      
 |      Pairs are returned in LIFO (last-in, first-out) order.
 |      Raises KeyError if the dict is empty.
 |  
 |  setdefault(self, key, default=None, /)
 |      Insert key with a value of default if key is not in the dictionary.
 |      
 |      Return the value for key if key is in the dictionary, else default.
 |  
 |  update(...)
 |      D.update([E, ]**F) -> None.  Update D from dict/iterable E and F.
 |      If E is present and has a .keys() method, then does:  for k in E: D[k] = E[k]
 |      If E is present and lacks a .keys() method, then does:  for k, v in E: D[k] = v
 |      In either case, this is followed by: for k in F:  D[k] = F[k]
 |  
 |  values(...)
 |      D.values() -> an object providing a view on D's values
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from builtins.dict:
 |  
 |  __class_getitem__(...) from builtins.type
 |      See PEP 585
 |  
 |  fromkeys(iterable, value=None, /) from builtins.type
 |      Create a new dictionary with keys from iterable and values set to value.
 |  
 |  ----------------------------------------------------------------------
 |  Static methods inherited from builtins.dict:
 |  
 |  __new__(*args, **kwargs) from builtins.type
 |      Create and return a new object.  See help(type) for accurate signature.
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes inherited from builtins.dict:
 |  
 |  __hash__ = None
```

### rnavigate.styles.set_defaults

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function set_defaults in module rnavigate.styles

set_defaults(context='paper', style='ticks', colors='default', dpi=140)
    Set or reset the major global style settings.
    
    Parameters
    ----------
    context : str, defaults to "paper"
        Passed to seaborn.set_context
        Defaults to "paper"
    style : str, defaults to "ticks"
        Passed to seaborn.set_style
    colors : str, defaults to "default"
        Passed to seaborn.set_palette
    dpi : int, defaults to 140
        Sets the dots-per-inch for inline and exported images
```

### rnavigate.styles.get_nt_color

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_nt_color in module rnavigate.styles

get_nt_color(nt, colors=None)
    Get the RNAvigate color for a given nucleotide.
    
    Invalid nucleotides are set to gray
    
    Parameters
    ----------
    nt : str
        a nucleotide letter
    colors "rnavigate" or "old", defaults to settings["sequence_colors"]
        "rnavigate" uses blue, light blue, red, light red for "AUCG"
        "old" uses traditional red, yellow, blue, green for "AUCG"
    
    Returns:
    color : str
        a hex color string
```

### rnavigate.styles.get_nt_cmap

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_nt_cmap in module rnavigate.styles

get_nt_cmap(colors=None)
    Get an rnavigate color map for nucleotides.
    
    Parameters
    ----------
    colors "rnavigate" or "old", defaults to settings["sequence_colors"]
        "rnavigate" uses blue, light blue, red, light red for "AUCG"
        "old" uses traditional red, yellow, blue, green for "AUCG"
    
    Returns
    -------
    cmap : rnavigate.data.colors.ScalarMappable (matplotlib.cm.ScalarMappable)
        a color map for nucleotides
```

### rnavigate.styles.apply_style

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function apply_style in module rnavigate.styles

apply_style(style_dict)
    Decorator for applying matplotlib style settings to a function.
```

## rnavigate.transcriptomics

### rnavigate.transcriptomics.BedFile

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class BedFile in module rnavigate.transcriptomics.bed

class BedFile(builtins.object)
 |  BedFile(bedfile)
 |  
 |  Methods defined here:
 |  
 |  __init__(self, bedfile)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_annotation(self, transcript, **kwargs)
 |  
 |  get_annotations(self, transcripts, **kwargs)
 |  
 |  get_density_profile(self, transcript, **kwargs)
 |  
 |  get_profile(self, transcript, **kwargs)
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.transcriptomics.NarrowPeak

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class NarrowPeak in module rnavigate.transcriptomics.bed

class NarrowPeak(BedFile)
 |  NarrowPeak(bedfile)
 |  
 |  Method resolution order:
 |      NarrowPeak
 |      BedFile
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, bedfile)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from BedFile:
 |  
 |  get_annotation(self, transcript, **kwargs)
 |  
 |  get_annotations(self, transcripts, **kwargs)
 |  
 |  get_density_profile(self, transcript, **kwargs)
 |  
 |  get_profile(self, transcript, **kwargs)
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from BedFile:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.transcriptomics.Transcriptome

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Transcriptome in module rnavigate.transcriptomics.transcriptome

class Transcriptome(builtins.object)
 |  Transcriptome(genome, annotation, path=WindowsPath('C:/Users/psirv/OneDrive - University of North Carolina at Chapel Hill/GitHub/RNAvigate/reference_data'), chr_ids=None)
 |  
 |  Methods defined here:
 |  
 |  __init__(self, genome, annotation, path=WindowsPath('C:/Users/psirv/OneDrive - University of North Carolina at Chapel Hill/GitHub/RNAvigate/reference_data'), chr_ids=None)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_sequence(self, chromosome, coordinates, strand)
 |  
 |  get_sequences(self, chromosomes, coordinates, strands)
 |  
 |  get_transcript(self, transcript_id)
 |  
 |  get_transcripts(self, transcript_ids)
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.transcriptomics.Transcript

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Transcript in module rnavigate.transcriptomics.transcriptome

class Transcript(rnavigate.data.data.Sequence)
 |  Transcript(parent, name, sequence, chromosome, strand, coordinates, tx_info, cds_coors=None, other_features=None)
 |  
 |  Method resolution order:
 |      Transcript
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, parent, name, sequence, chromosome, strand, coordinates, tx_info, cds_coors=None, other_features=None)
 |      Initialize the Sequence object.
 |  
 |  get_cds_annotation(self, **kwargs)
 |  
 |  get_cds_domains(self)
 |  
 |  get_coordinate_df(self)
 |  
 |  get_exon_annotation(self, exon_number, **kwargs)
 |  
 |  get_exon_domains(self)
 |  
 |  get_junctions_annotation(self, **kwargs)
 |  
 |  get_tx_coordinate(self, coordinate)
 |  
 |  get_tx_range(self, start, stop)
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return the name of the sequence.
 |  
 |  get_aligned_data(self, alignment)
 |      Get a copy of the sequence positionally aligned to another sequence.
 |      
 |      Parameters
 |      ----------
 |      alignment : rnavigate.data.Alignment
 |          the alignment to use
 |      
 |      Returns
 |      -------
 |      aligned_sequence : rnavigate.data.Sequence
 |          the aligned sequence
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get colors and colormap representing information about the sequence.
 |      
 |      Parameters
 |      ----------
 |      source : str, list, or matplotlib color-like
 |          the source of the color information
 |          if a string, must be one of:
 |              "sequence", "position", "profile", "structure", "annotations"
 |          if a list, must be a list of matplotlib color-like values, colormap
 |              will be None.
 |          if a matplotlib color-like value, all nucleotides will be colored
 |              that color, colormap will be None.
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors if source is "position"
 |      profile : rnavigate.data.Profile, optional
 |          the profile to use to get colors if source is "profile"
 |      structure : rnavigate.data.SecondaryStructure, optional
 |          the structure to use to get colors if source is "structure"
 |      annotations : list of rnavigate.data.Annotations, optional
 |          the annotations to use to get colors if source is "annotations"
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_annotations(self, annotations, default_color='gray')
 |      Get colors and colormap representing sequence annotations.
 |      
 |      Parameters
 |      ----------
 |      annotations : list of rnavigate.data.Annotations
 |          the annotations to use to get colors.
 |      default_color : matplotlib color-like, defaults to "gray"
 |          the color to use for nucleotides not in any annotation
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |      Get colors and colormap representing the nucleotide position.
 |      
 |      Parameters
 |      ----------
 |      pos_cmap : str, defaults to "rainbow"
 |          cmap used for position colors
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_profile(self, profile)
 |      Get colors and colormap representing per-nucleotide data.
 |      
 |      Parameters
 |      ----------
 |      profile : rnavigate.data.Profile
 |          the profile to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_sequence(self)
 |      Get a colors and colormap representing the nucleotide sequence.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_colors_from_structure(self, structure)
 |      Get colors and colormap representing base-pairing status.
 |      
 |      Parameters
 |      ----------
 |      structure : rnavigate.data.SecondaryStructure
 |          the structure to use to get colors.
 |      
 |      Returns
 |      -------
 |      colors : numpy array
 |          one matplotlib color-like value for each nucleotide in self.sequence
 |      colormap : rnavigate.data.ScalarMappable
 |          a colormap used for creating a colorbar
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Parameters
 |      ----------
 |      dataframe : pandas.DataFrame
 |          must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Parameters
 |      ----------
 |      t_or_u : "T", "U", or False, defaults to "U"
 |          "T" converts "U"s to "T"s
 |          "U" converts "T"s to "U"s
 |          False does nothing.
 |      uppercase : bool, defaults to True
 |          Whether to make sequence all uppercase
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence.
 |      
 |      Parameters
 |      ----------
 |      fasta : string
 |          path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns
 |      -------
 |      length : int
 |          the length of self.sequence
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Sequence:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.transcriptomics.eCLIPDatabase

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class eCLIPDatabase in module rnavigate.transcriptomics.eclip

class eCLIPDatabase(builtins.object)
 |  Methods defined here:
 |  
 |  __init__(self)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_annotation(self, transcript, cell_line, target)
 |  
 |  get_eclip_data(self)
 |  
 |  get_eclip_density(self, transcript, cell_line, targets=None)
 |  
 |  get_profile(self, transcript, cell_line, target)
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```

### rnavigate.transcriptomics.download_eclip_peaks

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function download_eclip_peaks in module rnavigate.transcriptomics.eclip

download_eclip_peaks(assembly='GRCh38', outpath=WindowsPath('C:/Users/psirv/OneDrive - University of North Carolina at Chapel Hill/GitHub/RNAvigate/reference_data/eCLIP_downloads'))
    download eCLIP bed files from ENCODE database
    
    Args:
        assembly (string, optional): reference genome ("h19" or "GRCh38")
            Defaults to "GRCh38"
        outpath (string, optional): output directory path
            Defaults to "eCLIP_downloads".
```

