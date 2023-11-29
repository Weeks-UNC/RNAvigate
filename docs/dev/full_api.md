
RNAvigate full API
==================

This is RNAvigate's (mostly) full API. I apologize that it is a bit ugly at the
moment, and that some docstrings are missing. I am currently working hard to
complete the following top priorities, in this order:

1. Fill in any missing docstrings.
2. Update current docstrings to Google-style.
3. Implement mkdocstrings for prettier automatic documentation.

Table of contents:

- rnavigate
    - [Sample](#rnavigatesample)
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
 |  Methods defined here:
 |  
 |  __init__(self, sample, inherit=None, keep_inherited_defaults=True, **data_keywords)
 |      Creates a Sample.
 |      
 |      Required arguments:
 |          sample (string)
 |              An arbitrary name. This will be used as a label in plot legends
 |              and titles to differentiate it from other samples
 |      
 |      Optional arguments:
 |          inherit (Sample or list of Samples)
 |              Data keywords and associated data from other samples become the
 |              data keywords and associated data from this sample. This does
 |              not make additional copies of the data: i.e. operations that
 |              make changes to inherited data change the original sample, and
 |              any other samples that inherited that data. This can be useful
 |              to save time and memory on operations and large data structures
 |              that are shared between samples.
 |          keep_inherited_defaults (True or False)
 |              whether to keep inherited default keywords
 |              defaults to True
 |      
 |      Data keywords:
 |          There are many built-in data keywords with different expectations
 |          and behaviors. For a full list with expected input formats and
 |          output behavior, visit:
 |      
 |          https://rnavigate.readthedocs.io/en/latest/loading-data/
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and filters interactions data.
 |      
 |      Args:
 |          interactions (rnavigate.data.Interactions | str):
 |              Interactions object to be filtered. If a string, value is
 |              replaced with self.get_data(interactions)
 |          metric (str, optional):
 |              column of interactions data to be used as metric for coloring
 |              interactions.
 |              "Distance" will compute 3D distance in "pdb", defaulting to
 |              2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |              those atoms to compute distance.
 |          cmap (str | list, optional):
 |              sets the interactions colormap, used to color interactions
 |              according to metric values.
 |          normalization (str, optional): 
 |              `'norm'`: extreme values in colormap are given to the extreme
 |                  values of interactions metric data
 |              `'bins'`: data are colored according to which bin they fall into
 |                  `values` defines bins (list, length = 2 less than cmap)
 |              `'min_max'`: extreme values in cmap are given to values beyond
 |                  minimum and maximum, defined by `values`
 |          values:
 |              behavior depends on normalization
 |              `'norm'`: values are not needed
 |              `'bins'`: list of floats containing the boundaries between bins
 |                  One fewer than the number of categories
 |              `'min_max'`: list of floats containing the minimum and maximum
 |          **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Required arguments:
 |          data_keyword (Data or data keyword or list/dict of such types)
 |              If None, returns None.
 |              If a data keyword, returns associated data from sample
 |              If Data, returns that data.
 |              If a list or dictionary, returns list or dictionary with
 |                  data keyword values replaced with associated Data
 |          data_class (RNAvigate Data class)
 |              If provided, ensures that returned data is of this type.
 |      
 |      Returns:
 |          Same type as data_keyword argument, but data keywords are replaced
 |              with associated data
 |      
 |      Raises:
 |          ValueError:
 |              if data is not found in sample
 |          ValueError:
 |              if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Args:
 |          inherit (RNAvigate Sample or list of Samples)
 |              Other samples from which to inherit data and data keywords
 |          keep_inherited_defaults (True or False)
 |              Use default values from inherited samples
 |          overwrite (True or False)
 |              whether to overwrite any existing keywords
 |      
 |      Raises:
 |          ValueError: if inherit is not a Sample or list of Samples
 |  
 |  print_data_keywords(self)
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              The data keyword to set as the default
 |      
 |      Optional arguments:
 |          overwrite (True or False)
 |              whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name',
 |              data_keyword=inputs)
 |      
 |          is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name')
 |          my_sample.add_data(
 |              'data_keyword', inputs)
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              a data keyword (arbitrary or standard) used to store and/or
 |              parse the inputs
 |          inputs (dictionary or RNAvigate Data)
 |              a dictionary used to create the data object
 |      
 |      Optional arguments:
 |          overwrite (bool)
 |              whether to overwrite a pre-existing data_keyword
 |              Defaults to False.
 |      
 |      Raises:
 |          ValueError:
 |              the data keyword already exists and overwrite is False
 |          ValueError:
 |              there was an issue parsing the data
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

## rnavigate.plot_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_alignment in module rnavigate.plotting_functions

plot_alignment(data1, data2, labels=None, plot_kwargs=None)
    Plots the sequence alignment used to compare two sequences
    
    Required arguments:
        data1 (tuple (rnavigate Sample, data keyword))
            a sample and data keyword to retrieve a sequence
        data2 (tuple (rnavigate Sample, data keyword))
            another sample and data keyword to retrieve a second sequence
    
    Optional display arguments:
        labels (list of 2 strings)
            Labels used for each sample
            Defaults to "sample.sample: data keyword" for each data input
        plot_kwargs (dict)
            passed to matplotlib.pyplot.subplots()
            Defaults to {}.
    
    Returns:
        rnavigate.plots.Alignment: the Alignment plot object
```

## rnavigate.plot_arcs

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_arcs in module rnavigate.plotting_functions

plot_arcs(samples, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, annotations=None, domains=None, labels=None, nt_ticks=(20, 5), profile_scale_factor=1, plot_error=False, annotation_mode='track', panels=None, seqbar=True, region='all', colorbars=True, title=True, plot_kwargs=None)
    Plots interactions and/or base-pairs as arcs.
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        sequence (data, data keyword, or raw sequence)
            All data are mapped to this sequence before plotting
            If a data keyword, data from the first sample will be used
    
    Optional data input arguments:
        structure (data or data keyword)
            secondary structure to plot as arcs
            Defaults to None
        structure2 (data or data keyword)
            another secondary structure to compare with the first structure
            arcs will be colored depending on which structure they are in
            Defaults to None
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot as arcs, no filtering performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
        interactions2 (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot as arcs, no filtering performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            Defaults to None
        profile (data or data keyword)
            Profile from which values will be plotted
            Defaults to None
        annotations (list of data or data keywords)
            Annotations used to highlight regions or sites of interest
            Defaults to [].
        domains (data or data keyword)
            domains to label along x-axis
            Defaults to None
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        nt_ticks (tuple of two integers)
            first integer is the gap between major tick marks
            second integer is the gap between minor tick marks
            defaults to (20, 5)
        profile_scale_factor (number)
            small profile values will be hard to see
            large profile values will overwhelm the plot
            e.g. use 1/10 to scale values down 10-fold, use 10 to scale up
            Defaults to 1
        plot_error (True or False)
            Whether to plot error bars, values are determined by profile.metric
            Defaults to False
        annotation_mode ('track' or 'bars')
            'track' will highlight annotations along the x-axis
            'bars' will use a vertical transparent bar over the plot
            Defaults to 'track'
        panels (dictionary)
            a dictionary of whether plot elements are displayed on the 'top'
            (above x-axis) or 'bottom' (below x-axis)
            Only the values you wish to change from the default are needed
            defaults to {'interactions': 'bottom',
                         'interactions2': 'bottom',
                         'structure': 'top',
                         'profile': 'top'}
        seqbar (True or False)
            whether to display the sequence along the x-axis
            Defaults to True
        region (list of 2 integers)
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, length of sequence]
    
    Optional plot display arguments:
        colorbars (True or False)
            Whether to plot colorbars for all plot elements
            Defaults to True
        title (True or False)
            Whether to display titles for each axis
            Defaults to True
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.AP: the ArcPlot object
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
    
    Required arguments:
        samples (list of 2 rnavigate Samples)
            samples used to retrieve data
            This plotting function can only compare two samples at a time
        sequence (data keyword)
            All data are mapped to this sequence taken from their respective
            sample before plotting
    
    Optional data input arguments:
        structure (data keyword)
            secondary structure to plot as arcs
            Defaults to None
        structure2 (data keyword)
            another secondary structure to compare with the first structure
            arcs will be colored depending on which structure they are in
            Defaults to None
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot as arcs, no filtering performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
        interactions2 (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot as arcs, no filtering performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            Defaults to None
        profile (data or data keyword)
            Profile from which values will be plotted
            Defaults to None
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        profile_scale_factor (number)
            small profile values will be hard to see
            large profile values will overwhelm the plot
            e.g. use 1/10 to scale values down 10-fold, use 10 to scale up
            Defaults to 1
        plot_error (True or False)
            Whether to plot error bars, values are determined by profile.metric
            Defaults to False
        region (list of 2 integers)
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, length of sequence]
    
    Optional plot display arguments:
        colorbars (True or False)
            Whether to plot color scales for all plot elements
            Defaults to True
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.AP plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
```

## rnavigate.plot_circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_circle in module rnavigate.plotting_functions

plot_circle(samples, sequence, structure=None, structure2=None, interactions=None, interactions2=None, annotations=None, profile=None, colors=None, nt_ticks=(20, 5), gap=30, labels=None, colorbars=True, plot_kwargs=None)
    Creates a figure containing a circle plot for each sample given.
    
    Data that can be plotted on circle plots includes annotations (highlights
    regions around the edge.Generates a multipanel secondary structure drawing with optional
    coloring by per-nucleotide data and display of inter-nucleotide data and/or
    sequence annotations. Each plot may display a unique sample and/or
    inter-nucleotide data filtering scheme.
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        sequence (data or data keyword)
            All data are mapped to this sequence before plotting
    
    Optional data input arguments:
        structure (data or data keyword)
            Structure used to plot base-pairs on circle plot
        structure2 (data or data keyword or list of either)
            Structures to compare with Structure. Each base-pair is colored by
            which structure contains it or how many structures contain it.
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot on cirle plot, no filtering
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
        interactions2 (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot on circle plot, no filtering
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            Defaults to None
        annotations (list of data or data keywords)
            Annotations used to highlight regions or sites of interest
            Defaults to [].
        profile (data or data keyword)
            Profile used for coloring if "profile" used in colors dictionary
            Defaults to None
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        colors (dictionary)
            a dictionary of element: value pairs that determines how colors
            will be applied to each plot element and if that element is plotted
            only the elements you wish to change need to be included
            Defaults to {'sequence': None, 'nucleotides': 'sequence'}
            value options and what the colors represent:
                None: don't plot this elelement
                'sequence': nucleotide identity
                'position': position in sequence
                'annotations': sequence annotations
                'profile': per-nucleotide data from profile
                    profile argument must be provided
                'structure': base-pairing status
                matplotlib color: all positions plotted in this color
                array of colors: a different color for each position
                    must be the same length as structure
            'sequence' may also use 'contrast' which automatically chooses
                white or black, which ever contrasts better with 'nucleotide'
                color
        nt_ticks (tuple of two integers)
            first integer is the gap between major tick marks
            second integer is the gap between minor tick marks
            defaults to (20, 5)
        gap (float)
            Width of gap between 5' and 3' end in degrees
            Defaults to 30
    
    Optional plot display arguments:
        colorbars (True or False)
            Whether to plot color scales for all plot elements
            Defaults to True
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.Circle plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
```

## rnavigate.plot_disthist

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_disthist in module rnavigate.plotting_functions

plot_disthist(samples, structure, interactions, bg_interactions=None, labels=None, same_axis=False, atom="O2'", rows=None, cols=None, plot_kwargs=None)
    Calculates 3D distance of nucleotides in inter-nucleotide data and plots
    the distribution of these distances. Compares this to a 'background'
    distribution consisting of either all pairwise distances in structure, or
    those defined by bg_interactions and bg_interactions_filter
    
    Required arguments:
        samples (list of rnavigate Samples)
            Samples from which to retreive data
            There will be one panel for each sample unless same_axis is True
        structure (str)
            secondary structure or 3D structure to calculate inter-nucleotide
            contact distance or 3D distance, respectively
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions used to calculate distance histogram, no filtering
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
    
    Optional data input arguments:
        bg_interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to calculate background distance histogram, no
                filtering is performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            if not provided, background distance histogram is calculated from
            all pairwise distances in structure
            Defaults to None
    
    Optional data display arguments:
        labels (str)
            Labels to be used as titles, must be same length as samples list
            Defaults to sample.sample for each sample
        atom (string or dictionary)
            from which atoms to calculate distances
            for DMS reactive atoms (N1 for A and G, N3 for U and C) use "DMS"
            use a dictionary to specify a different atom for each nucleotide
                e.g. "DMS" == {'A': 'N1', 'G': 'N1', 'U': 'N3', 'C': 'N3'}
            Defaults to "O2'"
    
    Optional plot display arguments:
        rows (integer)
            number of rows of plots
            Defaults to None (determined automatically)
        cols (integer)
            number of columns of plots
            Defaults to None (determined automatically)
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}
    
    Returns:
        rnavigate.plots.DistHist: object containing matplotlib figure and axes
            with additional plotting and file saving methods
```

## rnavigate.plot_heatmap

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_heatmap in module rnavigate.plotting_functions

plot_heatmap(samples, sequence, structure=None, interactions=None, regions=None, labels=None, levels=None, interpolation='nearest', atom="O2'", plot_type='heatmap', weights=None, rows=None, cols=None, plot_kwargs=None)
    Generates a multipanel plot displaying a heatmap of inter-nucleotide
    data (nucleotide resolution of 2D KDE) and/or contour map of pdb
    distances. Each plot may display a unique sample and/or filtering scheme.
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        sequence (data, data keyword or sequence string)
            All data are mapped to this sequence before plotting
    
    Optional data input arguments:
        structure (data or data keyword)
            secondary structure or 3D structure used to plot contour lines
            contour lines are drawn according to levels argument
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot as a heatmap, no filtering performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
        regions (list of lists of 4 integers)
            each inner list defines two regions of the RNA that are interacting
            a box will be drawn around this interaction on the heatmap
            e.g. [[10, 20, 50, 60], [35, 45, 70, 80]] draws 2 boxes
                the first box will connect nucleotides 10-20 and 50-60
                the second box will connect nucleotides 35-45 and 70-80
    
    Optional data display arguments:
        labels (string)
            Labels to be used as titles, must be same length as samples list
            Defaults to sample.sample for each sample
        levels (list of numbers)
            contours are drawn separating nucleotides above and below these
            distances
            if structure argument is a secondary structure
                distance refers to contact distance
                Defaults to [5]
            if structure argument is a 3D structure
                distance refers to spatial distance in angstroms
                Defaults to [20]
        interpolation (string)
            one of matplotlib's interpolations for heatmap (used with imshow)
            'nearest' works well for shorter RNAs (under 300 nt)
            'none' works well for longer RNAs (over 1200 nt)
            defaults to None (uses default)
        atom (string or dictionary)
            from which atoms to calculate distances
            for DMS reactive atoms (N1 for A and G, N3 for U and C) use "DMS"
            use a dictionary to specify a different atom for each nucleotide
                e.g. "DMS" == {'A': 'N1', 'G': 'N1', 'U': 'N3', 'C': 'N3'}
            Defaults to "O2'"
        plot_type ("heatmap" or "kde")
            how to plot interactions data
            "heatmap" will plot raw data, each interaction is a pixel in a grid
            "kde" will calculate a kernel density estimate and plot 5 levels
            Defaults to "heatmap"
        weights (string)
            weights to be used in kernel density estimation
            must be a column of interactions data
            Defaults to None
    
    Optional plot display arguments:
        rows (integer)
            number of rows of plots
            Defaults to None (determined automatically)
        cols (integer)
            number of columns of plots
            Defaults to None (determined automatically)
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}
    
    Returns:
        rnavigate.plots.Heatmap plot: object containing matplotlib figure and
            axes with additional plotting and file saving methods
```

## rnavigate.plot_linreg

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_linreg in module rnavigate.plotting_functions

plot_linreg(samples, profile, sequence=None, structure=None, annotations=None, labels=None, kde=False, scale='linear', regression='pearson', colors='sequence', column=None, region='all', colorbars=True, plot_kwargs=None)
    Performs linear regression analysis and generates scatter plots of all
    sample-to-sample profile vs. profile comparisons. Colors nucleotides by
    identity or base-pairing status.
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        profile (data or data keyword)
            per-nucleotide data to perform linear regression
            all data are mapped to the sequence of the profile data from the
            first sample before plotting, unless sequence is supplied
    
    Optional data input arguments:
        sequence (data or data keyword)
            a sequence from which to align all profiles
            if a data keyword, uses data from the first sample
            Defaults to None
        structure (data or data keyword)
            Structure used for coloring if colors argument is "structure"
            Defaults to None
        annotations (list of data or data keywords)
            Annotations used for coloring if colors argument is "annotations"
            Defaults to [].
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        kde (True or False)
            whether to plot kde (density) instead of a scatter plot
        scale ('linear' or 'log')
            'linear' performs regression on raw values, displays linear units
            'log' performs regression on log10(values), displays log10 units
            Defaults to 'linear'
        regression ('pearson' or 'spearman')
            'pearson' calculates Pearson R-squared (standard)
            'spearman' calculates Spearman R-squared (rank-order)
            Defaults to 'pearson'
        colors (string or list)
            value options and what the colors represent:
                'sequence': nucleotide identity
                'position': position in sequence
                'annotations': sequence annotations
                'profile': per-nucleotide data from profile
                    profile argument must be provided
                'structure': base-pairing status
                matplotlib color: all positions plotted in this color
                array of colors: a different color for each position
                    must be the same length as structure
            Defaults to 'sequence'
        column (string)
            column name of values from profile to use in regression
            Defaults to profile.metric
        region (list of 2 integers)
            start and end nucleotide positions to include. 1-indexed, inclusive
            Defaults to [1, length of sequence]
    
    Optional plot display arguments:
        colorbars (True or False)
            Whether to plot colorbars for scatter plot colors
            Defaults to True
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.LinReg: object containing matplotlib figure and axes
            with additional plotting and file saving methods
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
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        structure (data or data keyword)
            3D structure to view as interactive molecule
            All data are mapped to this sequence before plotting
    
    Optional data input arguments:
        profile (data or data keyword)
            Profile used to color nucleotides if colors="profile"
            Defaults to None
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot on molecule, no filtering performed
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot titles
            Defaults to sample.sample for each sample
        style ('cartoon', 'cross', 'line', 'sphere' or 'stick')
            sets the py3Dmol style for drawing the molecule
            Defaults to 'cartoon'
        hide_cylinders (True or False)
            whether to hide nucleotide cylinders (only shows backbone ribbon)
            Defaults to False
        colors (string or list)
            value options and what the colors represent:
                'sequence': nucleotide identity
                'position': position in sequence
                'annotations': sequence annotations
                'profile': per-nucleotide data from profile
                    profile argument must be provided
                'structure': base-pairing status
                matplotlib color: all positions plotted in this color
                array of colors: a different color for each position
                    must be the same length as structure
            Defaults to 'grey'
        atom (string or dictionary)
            which atoms to draw interactions between
            for DMS reactive atoms (N1 for A and G, N3 for U and C) use "DMS"
            use a dictionary to specify a different atom for each nucleotide
                e.g. "DMS" == {'A': 'N1', 'G': 'N1', 'U': 'N3', 'C': 'N3'}
            Defaults to "O2'"
        rotation (dictionary)
            axis-degrees pairs for setting the starting orientation of the
            molecule, only the axes to be rotated are needed
            e.g. {'x': 180} flips the molecule on the x-axis
            Defaults to None
        orientation (list of floats)
            set the precise starting orientation
            see get_orientation for more details
            Defaults to None
        get_orientation (True or False)
            allows getting the orientation for use with orientation argument
            all other arguments will be ignored and a larger, single panel view
            window is displayed with no title
                1. adjust the molecule to the desired orientation
                2. click on the molecule to display the orientation vector
                3. copy this orientation vector (manually)
                4. provide this list of values to the orientation argument
            Defaults to False
    
    Optional viewer display arguments:
        title (True or False)
            whether to display the title
            Defaults to True
        colorbars (True or False)
            Whether to plot color scales for all plot elements
            Defaults to True
        width (integer)
            width of view window in pixels
            Defaults to 400
        height (integer)
            height of view window in pixels
            Defaults to 400
        rows (integer)
            the number of rows in the view window
            Defaults to None (set automatically)
        cols (integer)
            the number of columns in the view window
            Defaults to None (set automatically)
        background_alpha (float)
            the opacity of the view window, must be between 0 and 1
            Defaults to 1 (completely opaque)
        show (True or False)
            whether to display the viewer object
            Defaults to True
    
    Returns:
        rnavigate.plots.Mol plot: object containing py3dmol viewer with
            additional plotting and file saving methods
```

## rnavigate.plot_ntdist

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_ntdist in module rnavigate.plotting_functions

plot_ntdist(samples, profile, labels=None, column=None, plot_kwargs=None)
    Plots the distributions of values at A, U, C, and G.
    
    Calculates the kernel density estimate (KDE) for each nucleobase and plots
    them on one axis per sample.
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        profile (data or data keyword)
            per-nucleotide data to plot per-nt-identity distributions
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        column (string)
            which column of data to use for KDE
            defaults to 'AUCG'
    
    Optional plot display arguments:
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}
    
    Returns:
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
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        profile (data or data keyword)
            Profile from which values will be plotted
    
    Optional data input arguments:
        sequence (data, data keyword, or raw sequence)
            All data are mapped to this sequence before plotting
            If a data keyword, data from the first sample will be used
            Defaults to the value of the profile argument
        annotations (list of data or data keywords)
            Annotations used to highlight regions or sites of interest
            Defaults to [].
        domains (data or data keyword)
            domains to label along x-axis
            Defaults to None
    
    Optional data display arguments:
        labels (list of strings)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        nt_ticks (tuple of two integers)
            first integer is the gap between major tick marks
            second integer is the gap between minor tick marks
            defaults to (20, 5)
        column (string)
            column name of values from profile to plot
            Defaults to profile.metric
        plot_error (True or False)
            Whether to plot error bars, values are determined by profile.metric
            Defaults to True
        annotations_mode ('track' or 'bars')
            'track' will highlight annotations along the x-axis
            'bars' will use a vertical transparent bar over the plot
            Defaults to 'track'
        seqbar (True or False)
            whether to display the sequence along the x-axis
            Defaults to True
        region (list of 2 integers)
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, length of sequence]
    
    Optional plot display arguments:
        colorbars (True or False)
            Whether to plot color scales for per-nucleotide data
            Defaults to True
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.Profile: the Profile plot object
```

## rnavigate.plot_qc

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_qc in module rnavigate.plotting_functions

plot_qc(samples, profile, labels=None)
    Creates a multipanel quality control plot displaying mutations per
    molecule, read length distribution, and mutation rate distributions for
    modified and unmodified samples.
    
    Required arguments:
        samples (list of rnavigate.Sample)
            samples to retrieve data from
        profile (data or data keyword)
            ShapeMaP or similar data for plotting reactivity distributions
            Must contain data from ShapeMapper log file
    
    Optional display arguments:
        labels (list of str)
            labels to be used on legends, must be same length as samples list
            Defaults to sample.sample for each sample.
    
    Returns:
        rnavigate.plots.QC: the quality control plot object
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
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        structure (data or data keyword)
            secondary structure to use as classifier (paired or unpaired)
            profile data for each sample is first aligned to this structure
        profile (data or data keyword)
            per-nucleotide data to perform ROC analysis
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        nts (string)
            which nucleotides to plot nucleotide-type ROC plots
            defaults to 'AUCG'
    
    Optional plot display arguments:
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}
    
    Returns:
        rnavigate.plots.ROC: object containing matplotlib figure and axes with
            additional plotting and file saving methods
```

## rnavigate.plot_shapemapper

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_shapemapper in module rnavigate.plotting_functions

plot_shapemapper(sample, profile, label=None, panels=None)
    Makes a standard ShapeMapper2 profile plot with 3 panels: Normalized
    Reactivities, mutation rates, and read depths.
    
    Required arguments:
        sample (rnavigate Sample)
            The sample from which data profile and label will be retreived
        profile (data or data keyword)
            ShapeMaP or similar data for plotting profiles
    
    Optional display arguments:
        label (string)
            A label to use as the title of the figure
        panels (list)
            Which of the three panels to include.
            Defaults to ["profile", "rates", "depth"]
    
    Returns:
        rnavigate.plots.SM: the ShapeMapper2 plot object
```

## rnavigate.plot_skyline

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_skyline in module rnavigate.plotting_functions

plot_skyline(samples, profile, sequence=None, annotations=None, domains=None, labels=None, nt_ticks=(20, 5), columns=None, errors=None, annotations_mode='track', seqbar=True, region='all', plot_kwargs=None)
    Plots multiple per-nucleotide datasets on a single axis.
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        profile (data or data keyword)
            Profile from which values will be plotted
    
    Optional data input arguments:
        sequence (data, data keyword, or raw sequence string)
            All data are mapped to this sequence before plotting
            If a data keyword, data from the first sample will be used
            Defaults to the value of the profile argument
        annotations (list of data or data keywords)
            Annotations used to highlight regions or sites of interest
            Defaults to [].
        domains (data or data keyword)
            domains to label along x-axis
            Defaults to None
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        nt_ticks (tuple of two integers)
            first integer is the gap between major tick marks
            second integer is the gap between minor tick marks
            defaults to (20, 5)
        columns (string or list of strings)
            columns names of values from profile to plot
            Defaults to profile.metric
        errors (string or list of strings)
            column names of error values for plotting error bars
            Defaults to None (no error bars)
        annotations_mode ('track' or 'bars')
            'track' will highlight annotations along the x-axis
            'bars' will use a vertical transparent bar over the plot
            Defaults to 'track'
        seqbar (True or False)
            whether to display the sequence along the x-axis
            Defaults to True
        region (list of 2 integers)
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, length of sequence]
    
    Optional plot display arguements:
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.Skyline: the skyline plot object
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
    
    Required arguments:
        samples (list of rnavigate Samples)
            samples used to retrieve data
        structure (data or data keyword)
            secondary structure to plot as arcs
            All data are mapped to this sequence before plotting
    
    Optional data input arguments:
        profile (data or data keyword)
            Profile used for coloring if "profile" used in colors dictionary
            Defaults to None
        annotations (list of data or data keywords)
            Annotations used to highlight regions or sites of interest
            Defaults to [].
        interactions (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot on secondary structure, no filtering
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            format 3 (list of format 2 dictionaries)
                This format allows multiple filtering schemes to be applied,
                each will be plotted on a seperate axis
            Defaults to None
        interactions2 (one of the formats below)
            format 1 (data or data keyword)
                Interactions to plot on secondary structure, no filtering
            format 2 (dictionary)
                e.g. {'interactions': format 1}
                additional filtering options can be added to the dictionary
            Defaults to None
    
    Optional data display arguments:
        labels (list of str)
            list containing Labels to be used in plot legends
            Defaults to sample.sample for each sample
        colors (dictionary)
            a dictionary of element: value pairs that determines how colors
            will be applied to each plot element and if that element is plotted
            only the elements you wish to change need to be included
            value options and what the colors represent:
                None: don't plot this elelement
                'sequence': nucleotide identity
                'position': position in sequence
                'annotations': sequence annotations
                'profile': per-nucleotide data from profile
                    profile argument must be provided
                'structure': base-pairing status
                matplotlib color: all positions plotted in this color
                array of colors: a different color for each position
                    must be the same length as structure
            'sequence' may also use 'contrast' which automatically chooses
                white or black, which ever contrasts better with 'nucleotide'
                color
            Defaults to {'sequence': None,
                         'nucleotides': 'sequence',
                         'structure': 'grey',
                         'basepairs': 'grey'}
        nt_ticks (integer)
            gap between major tick marks
            defaults to None (no position labels)
        bp_style ('dotted', 'line', or 'conventional')
            'dotted' plots basepairs as a dotted line
            'line' plots basepairs as a solid line
            'conventional' plots basepairs using Leontis-Westhof conventions
                for canonical and wobble pairs ('G-A' plotted as solid dot)
    
    Optional plot display arguments:
        colorbars (True or False)
            Whether to plot color scales for all plot elements
            Defaults to True
        plot_kwargs (dictionary)
            Keyword-arguments passed to matplotlib.pyplot.subplots
            Defaults to {}.
    
    Returns:
        rnavigate.plots.SS plot: object containing matplotlib figure and axes
            with additional plotting and file saving methods
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
 |  First, the log10(modified/untreated) rate is calculated. These values a
 |  scaled to minimize the median of the absolute value of the difference
 |  between samples. The standard error in these values is computed for each
 |  replicate. Z-scores between samples are calculated. The results are plotted
 |  in two panels: (1) the scaled log10(modified/untreated) rate for each
 |  sample with error bars, and (2) the difference between samples, colored by
 |  z-score.
 |  
 |  Methods:
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
 |  Methods:
 |      __init__: performs the analysis
 |      plot_lowss: displays the result and returns plot object
 |  
 |  Attributes:
 |      sample (str): the new label for this Sample's data on plots
 |      parent (rnavigate.Sample): the sample from which data is retrieved
 |      window (int): size of the windows, must be odd
 |      median_shape (float): global median SHAPE reactivity
 |      median_entropy (float): global median Shannon entropy
 |      data (dictionary): dictionary of data keyword: Data objects, keys are:
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
 |      Perform Low SHAPE, low Shannon entropy analysis given a sample with
 |      1) reactivities 2) MFE structure 3) pairing probabilities.
 |      
 |      Required arguments:
 |          sample (rnav.Sample)
 |              sample with shapemap, pairing probabilities and MFE structure
 |      
 |      Optional arguments:
 |          window (integer)
 |              Window size for calculating median SHAPE and Shannon entropy
 |              Defaults to 55
 |          shapemap (data keyword)
 |              data keyword that points to SHAPE-MaP data
 |              defaults to "shapemap"
 |          pairprob (data keyword)
 |              data keyword that points to pairing probabilities data
 |              defaults to "pairprob"
 |          structure (data keyword)
 |              data keyword that points to MFE structure data
 |              defaults to "ss"
 |  
 |  plot_lowss(self, region=None, colorbars=True)
 |      Visualize LowSS analysis over the given region.
 |      
 |      Optional arguments:
 |          region (integer or list of 2 integers)
 |              lowSS region number or start and end positions to plot.
 |              For region numbers, +/- 150 nts are shown.
 |              Defaults to entire sequence.
 |          colorbars (True or False)
 |              whether to plot colorbars for pairing probability
 |      
 |      Returns:
 |          rnavigate.plots.AP: LowSS visualization
 |  
 |  reset_lowss(self, maximum_shape=None, maximum_entropy=0.08)
 |      Generates an annotation of lowSS regions. Stored as self.lowSS
 |      
 |      Args:
 |          maximum_shape (float, optional): maximum normalized SHAPE
 |              reactivity to be called lowSS. Defaults to median reactivity.
 |          maximum_entropy (float, optional): maximum shannon entropy to be
 |              called lowSS. Defaults to 0.8.
 |  
 |  reset_window(self, window=None)
 |      Resets the window size and recalculates windowed SHAPE reactivities
 |      and shannon entropies and lowSS region annotations.
 |      
 |      Optional arguments:
 |          window (integer)
 |              window size, must be an odd integer
 |              Defaults to value of window attribute
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.rnavigate.Sample:
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and filters interactions data.
 |      
 |      Args:
 |          interactions (rnavigate.data.Interactions | str):
 |              Interactions object to be filtered. If a string, value is
 |              replaced with self.get_data(interactions)
 |          metric (str, optional):
 |              column of interactions data to be used as metric for coloring
 |              interactions.
 |              "Distance" will compute 3D distance in "pdb", defaulting to
 |              2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |              those atoms to compute distance.
 |          cmap (str | list, optional):
 |              sets the interactions colormap, used to color interactions
 |              according to metric values.
 |          normalization (str, optional): 
 |              `'norm'`: extreme values in colormap are given to the extreme
 |                  values of interactions metric data
 |              `'bins'`: data are colored according to which bin they fall into
 |                  `values` defines bins (list, length = 2 less than cmap)
 |              `'min_max'`: extreme values in cmap are given to values beyond
 |                  minimum and maximum, defined by `values`
 |          values:
 |              behavior depends on normalization
 |              `'norm'`: values are not needed
 |              `'bins'`: list of floats containing the boundaries between bins
 |                  One fewer than the number of categories
 |              `'min_max'`: list of floats containing the minimum and maximum
 |          **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Required arguments:
 |          data_keyword (Data or data keyword or list/dict of such types)
 |              If None, returns None.
 |              If a data keyword, returns associated data from sample
 |              If Data, returns that data.
 |              If a list or dictionary, returns list or dictionary with
 |                  data keyword values replaced with associated Data
 |          data_class (RNAvigate Data class)
 |              If provided, ensures that returned data is of this type.
 |      
 |      Returns:
 |          Same type as data_keyword argument, but data keywords are replaced
 |              with associated data
 |      
 |      Raises:
 |          ValueError:
 |              if data is not found in sample
 |          ValueError:
 |              if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Args:
 |          inherit (RNAvigate Sample or list of Samples)
 |              Other samples from which to inherit data and data keywords
 |          keep_inherited_defaults (True or False)
 |              Use default values from inherited samples
 |          overwrite (True or False)
 |              whether to overwrite any existing keywords
 |      
 |      Raises:
 |          ValueError: if inherit is not a Sample or list of Samples
 |  
 |  print_data_keywords(self)
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              The data keyword to set as the default
 |      
 |      Optional arguments:
 |          overwrite (True or False)
 |              whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name',
 |              data_keyword=inputs)
 |      
 |          is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name')
 |          my_sample.add_data(
 |              'data_keyword', inputs)
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              a data keyword (arbitrary or standard) used to store and/or
 |              parse the inputs
 |          inputs (dictionary or RNAvigate Data)
 |              a dictionary used to create the data object
 |      
 |      Optional arguments:
 |          overwrite (bool)
 |              whether to overwrite a pre-existing data_keyword
 |              Defaults to False.
 |      
 |      Raises:
 |          ValueError:
 |              the data keyword already exists and overwrite is False
 |          ValueError:
 |              there was an issue parsing the data
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
 |  Citation: (doi:10.1021/acs.biochem.5b00977)
 |  
 |  Algorithm:
 |      1. Extract SHAPE-MaP sequence, normalized profile, and normalized
 |         standard error from given samples
 |      2. Calculated smoothed profiles (mean) and propagate standard errors
 |         over rolling windows
 |      3. Subtract raw and smoothed normalized profiles and propogate errors
 |      4. Calculate Z-factors for smoothed data. This is the magnitude of the
 |         difference relative to the standard error
 |      5. Calculate Z-scores for smoothed data. This is the magnitude of the
 |         difference in standard deviations from the mean difference
 |      6. Call sites. Called sites must have # nucleotides that pass Z-factor
 |         and Z-score thresholds per window.
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
 |      Required Arguments:
 |          sample1 (rnavigate.Sample)
 |              First sample to compare
 |          sample2 (rnavigate.Sample)
 |              Second sample to compare
 |      
 |      Optional Arguments:
 |          profile (string)
 |              Data keyword pointing to SHAPE-MaP data in samples 1 and 2
 |              Defaults to 'shapemap'
 |          smoothing_window (integer)
 |              Size of windows for data smoothing
 |              Defaults to 3
 |          zf_coeff (float)
 |              Sites must have a difference more than zf_coeff standard errors
 |              Defaults to 1.96 (95% confidence interval)
 |          ss_thresh (float)
 |              Sites must have a difference that is ss_thresh standard
 |              deviations from the mean difference
 |              Defaults to 1
 |          site_window (integer)
 |              Number of nucleotides to include when calling sites
 |              Defaults to 3
 |          site_nts (integer)
 |              Number of nts within site_window that must pass thresholds
 |              Defaults to 2
 |  
 |  calculate_deltashape(self, smoothing_window=3, zf_coeff=1.96, ss_thresh=1, site_window=2, site_nts=3)
 |      Calculate or recalculate deltaSHAPE profile and called sites
 |      
 |      Optional Arguments:
 |          smoothing_window (integer)
 |              Size of windows for data smoothing
 |              Defaults to 3
 |          zf_coeff (float)
 |              Sites must have a difference more than zf_coeff standard errors
 |              Defaults to 1.96 (95% confidence interval)
 |          ss_thresh (float)
 |              Sites must have a difference that is ss_thresh standard
 |              deviations from the mean difference
 |              Defaults to 1
 |          site_window (integer)
 |              Number of nucleotides to include when calling sites
 |              Defaults to 3
 |          site_nts (integer)
 |              Number of nts within site_window that must pass thresholds
 |              Defaults to 2
 |  
 |  plot(self, region='all')
 |      Plot the deltaSHAPE result
 |      
 |      Optional arguments:
 |          region (list of 2 integers)
 |              start and end positions to plot
 |              Defaults to 'all'.
 |      
 |      Returns:
 |          rnav.plots.Profile: The plot object
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.rnavigate.Sample:
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and filters interactions data.
 |      
 |      Args:
 |          interactions (rnavigate.data.Interactions | str):
 |              Interactions object to be filtered. If a string, value is
 |              replaced with self.get_data(interactions)
 |          metric (str, optional):
 |              column of interactions data to be used as metric for coloring
 |              interactions.
 |              "Distance" will compute 3D distance in "pdb", defaulting to
 |              2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |              those atoms to compute distance.
 |          cmap (str | list, optional):
 |              sets the interactions colormap, used to color interactions
 |              according to metric values.
 |          normalization (str, optional): 
 |              `'norm'`: extreme values in colormap are given to the extreme
 |                  values of interactions metric data
 |              `'bins'`: data are colored according to which bin they fall into
 |                  `values` defines bins (list, length = 2 less than cmap)
 |              `'min_max'`: extreme values in cmap are given to values beyond
 |                  minimum and maximum, defined by `values`
 |          values:
 |              behavior depends on normalization
 |              `'norm'`: values are not needed
 |              `'bins'`: list of floats containing the boundaries between bins
 |                  One fewer than the number of categories
 |              `'min_max'`: list of floats containing the minimum and maximum
 |          **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Required arguments:
 |          data_keyword (Data or data keyword or list/dict of such types)
 |              If None, returns None.
 |              If a data keyword, returns associated data from sample
 |              If Data, returns that data.
 |              If a list or dictionary, returns list or dictionary with
 |                  data keyword values replaced with associated Data
 |          data_class (RNAvigate Data class)
 |              If provided, ensures that returned data is of this type.
 |      
 |      Returns:
 |          Same type as data_keyword argument, but data keywords are replaced
 |              with associated data
 |      
 |      Raises:
 |          ValueError:
 |              if data is not found in sample
 |          ValueError:
 |              if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Args:
 |          inherit (RNAvigate Sample or list of Samples)
 |              Other samples from which to inherit data and data keywords
 |          keep_inherited_defaults (True or False)
 |              Use default values from inherited samples
 |          overwrite (True or False)
 |              whether to overwrite any existing keywords
 |      
 |      Raises:
 |          ValueError: if inherit is not a Sample or list of Samples
 |  
 |  print_data_keywords(self)
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              The data keyword to set as the default
 |      
 |      Optional arguments:
 |          overwrite (True or False)
 |              whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name',
 |              data_keyword=inputs)
 |      
 |          is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name')
 |          my_sample.add_data(
 |              'data_keyword', inputs)
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              a data keyword (arbitrary or standard) used to store and/or
 |              parse the inputs
 |          inputs (dictionary or RNAvigate Data)
 |              a dictionary used to create the data object
 |      
 |      Optional arguments:
 |          overwrite (bool)
 |              whether to overwrite a pre-existing data_keyword
 |              Defaults to False.
 |      
 |      Raises:
 |          ValueError:
 |              the data keyword already exists and overwrite is False
 |          ValueError:
 |              there was an issue parsing the data
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
 |      Args:
 |          input_data (tuple of RNAvigate Profiles or Pandas Dataframe)
 |              if tuple of Profiles, the unified Dataframe will be created
 |  
 |  calculate_deltashape(self, smoothing_window=3, zf_coeff=1.96, ss_thresh=1, site_window=3, site_nts=2)
 |      Calculate the deltaSHAPE profile metrics
 |      
 |      Args:
 |          smoothing_window (int, optional): Defaults to 3.
 |          zf_coeff (float, optional): Defaults to 1.96.
 |          ss_thresh (int, optional): Defaults to 1.
 |          site_window (int, optional): Defaults to 3.
 |          site_nts (int, optional): Defaults to 2.
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
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from rnavigate.data.profile.Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.profile.Profile:
 |  
 |  recreation_kwargs
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Class for computing and displaying windowed AUROC analysis. This
 |  analysis computes the ROC curve over a sliding window for the performance
 |  of per-nucleotide data (usually SHAPE-MaP or DMS-MaP Normalized reactivity)
 |  in predicting the base-pairing status of each nucleotide. The area under
 |  this curve (AUROC) is displayed compared to the median across the RNA.
 |  Below, an arc plot displays the secondary structure and per-nucleotide
 |  profile.
 |  
 |  Citation:
 |  Lan, T.C.T., Allan, M.F., Malsick, L.E. et al. Secondary structural
 |      ensembles of the SARS-CoV-2 RNA genome in infected cells. Nat Commun
 |      13, 1128 (2022). https://doi.org/10.1038/s41467-022-28603-2
 |  
 |  Methods:
 |      __init__: Computes the AUROC array and AUROC median.
 |      plot_auroc: Displays the AUROC analysis over the given region.
 |          Returns Plot object
 |  
 |  Attributes:
 |      sample: an rnavigate.Sample to retrieve profile and secondary structure
 |      structure: sample.data[structure]
 |      profile: sample.data[profile]
 |      sequence: the sequence string of sample.data[structure]
 |      window: the size of the windows
 |      nt_length: the length of sequence string
 |      auroc: the auroc numpy array, length = nt_length, padded with np.nan
 |      median_auroc: the median of the auroc array
 |  
 |  Methods defined here:
 |  
 |  __init__(self, sample, window=81, profile='default_profile', structure='default_structure')
 |      Compute the AUROC for all windows. AUROC is a measure of how well a
 |      reactivity profile predicts paired vs. unpaired nucleotide status.
 |      
 |      Args:
 |          sample (rnav.Sample): Your rnavigate sample
 |          window (int, optional): number of nucleotides to include in window
 |              Defaults to 81.
 |          profile (str, optional): data keyword of provided sample pointing
 |              to a profile.
 |              Defaults to "default_profile"
 |          structure (str, optional): data keyword of provided sample pointing
 |              to a secondary structure.
 |              Defaults to "default_structure"
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
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
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
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.profile.Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from rnavigate.data.profile.Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Required arguments:
 |          sample (string)
 |              An arbitrary name. This will be used as a label in plot legends
 |              and titles to differentiate it from other samples
 |      
 |      Optional arguments:
 |          inherit (Sample or list of Samples)
 |              Data keywords and associated data from other samples become the
 |              data keywords and associated data from this sample. This does
 |              not make additional copies of the data: i.e. operations that
 |              make changes to inherited data change the original sample, and
 |              any other samples that inherited that data. This can be useful
 |              to save time and memory on operations and large data structures
 |              that are shared between samples.
 |          keep_inherited_defaults (True or False)
 |              whether to keep inherited default keywords
 |              defaults to True
 |      
 |      Data keywords:
 |          There are many built-in data keywords with different expectations
 |          and behaviors. For a full list with expected input formats and
 |          output behavior, visit:
 |      
 |          https://rnavigate.readthedocs.io/en/latest/loading-data/
 |  
 |  plot_scatter(self, column='Modified_rate')
 |      Generates scatter plots useful for fragmapper quality control.
 |      
 |      Args:
 |          column (str, optional): Dataframe column containing data to plot 
 |                                  (must be avalible for the sample and control).
 |                                  Defaults to 'Modified_rate'.
 |      
 |      Returns:
 |          scatter_plot: Scatter plot with control values on the x-axis,
 |                        sample values on the y-axis, and each point
 |                        representing a nucleotide not filtered out in the
 |                        fragmapper pipeline.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.rnavigate.Sample:
 |  
 |  filter_interactions(self, interactions, metric=None, cmap=None, normalization=None, values=None, **kwargs)
 |      sets coloring properties and filters interactions data.
 |      
 |      Args:
 |          interactions (rnavigate.data.Interactions | str):
 |              Interactions object to be filtered. If a string, value is
 |              replaced with self.get_data(interactions)
 |          metric (str, optional):
 |              column of interactions data to be used as metric for coloring
 |              interactions.
 |              "Distance" will compute 3D distance in "pdb", defaulting to
 |              2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
 |              those atoms to compute distance.
 |          cmap (str | list, optional):
 |              sets the interactions colormap, used to color interactions
 |              according to metric values.
 |          normalization (str, optional): 
 |              `'norm'`: extreme values in colormap are given to the extreme
 |                  values of interactions metric data
 |              `'bins'`: data are colored according to which bin they fall into
 |                  `values` defines bins (list, length = 2 less than cmap)
 |              `'min_max'`: extreme values in cmap are given to values beyond
 |                  minimum and maximum, defined by `values`
 |          values:
 |              behavior depends on normalization
 |              `'norm'`: values are not needed
 |              `'bins'`: list of floats containing the boundaries between bins
 |                  One fewer than the number of categories
 |              `'min_max'`: list of floats containing the minimum and maximum
 |          **kwargs: Other arguments are passed to interactions.filter()
 |  
 |  get_data(self, data_keyword, data_class=None)
 |      Replaces data keyword with data object, even if nested.
 |      
 |      Required arguments:
 |          data_keyword (Data or data keyword or list/dict of such types)
 |              If None, returns None.
 |              If a data keyword, returns associated data from sample
 |              If Data, returns that data.
 |              If a list or dictionary, returns list or dictionary with
 |                  data keyword values replaced with associated Data
 |          data_class (RNAvigate Data class)
 |              If provided, ensures that returned data is of this type.
 |      
 |      Returns:
 |          Same type as data_keyword argument, but data keywords are replaced
 |              with associated data
 |      
 |      Raises:
 |          ValueError:
 |              if data is not found in sample
 |          ValueError:
 |              if the data retrieved is not of the specified data_class
 |  
 |  inherit_data(self, inherit, keep_inherited_defaults, overwrite)
 |      retrieves and stores data and data keywords from other samples
 |      
 |      Args:
 |          inherit (RNAvigate Sample or list of Samples)
 |              Other samples from which to inherit data and data keywords
 |          keep_inherited_defaults (True or False)
 |              Use default values from inherited samples
 |          overwrite (True or False)
 |              whether to overwrite any existing keywords
 |      
 |      Raises:
 |          ValueError: if inherit is not a Sample or list of Samples
 |  
 |  print_data_keywords(self)
 |  
 |  set_as_default(self, data_keyword, overwrite=True)
 |      Set the given data keyword as the default for its data class
 |      
 |      It's data class is determined automatically. Only one default exists
 |      per data class and per Sample object.
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              The data keyword to set as the default
 |      
 |      Optional arguments:
 |          overwrite (True or False)
 |              whether to overwrite a pre-existing default data keyword
 |  
 |  set_data(self, data_keyword, inputs, overwrite=False)
 |      Add data to Sample using the given data keyword and inputs
 |      
 |      This methods works similarly to the data keywords arguments used
 |      during Sample initialization:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name',
 |              data_keyword=inputs)
 |      
 |          is equivalent to:
 |      
 |          my_sample = rnavigate.Sample(
 |              sample='name')
 |          my_sample.add_data(
 |              'data_keyword', inputs)
 |      
 |      Required arguments:
 |          data_keyword (string)
 |              a data keyword (arbitrary or standard) used to store and/or
 |              parse the inputs
 |          inputs (dictionary or RNAvigate Data)
 |              a dictionary used to create the data object
 |      
 |      Optional arguments:
 |          overwrite (bool)
 |              whether to overwrite a pre-existing data_keyword
 |              Defaults to False.
 |      
 |      Raises:
 |          ValueError:
 |              the data keyword already exists and overwrite is False
 |          ValueError:
 |              there was an issue parsing the data
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
 |  Attributes:
 |      samples: the list of samples
 |      sequences: a list of all unique sequence strings stored in the list of
 |          samples. These are converted to an all uppercase RNA alphabet.
 |      keywords: a list of unique data keywords stored in the list of samples.
 |      which_sequences: a dataframe of samples and keywords and which
 |          sequences each contains.
 |  
 |  Methods defined here:
 |  
 |  __init__(self, samples)
 |      Creates an instance of SequenceChecker given a list of samples
 |      
 |      Arguments:
 |          samples (list of rnav.Sample)
 |              samples for which to compare data keywords and sequences.
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
 |      Optional arguments:
 |          print_format (string)
 |              What format to print the alignments in:
 |              "cigar" prints the cigar string
 |              "short" prints the numbers of mismatches and indels
 |              "long" prints the location and nucleotide identity of all
 |                  mismatches, insertions and deletions.
 |              Defaults to "long".
 |          which (pair of integers)
 |              two sequence IDs to compare.
 |              Defaults to every pairwise comparison.
 |  
 |  print_mulitple_sequence_alignment(self, base_sequence)
 |      Print the multiple sequence alignment with nice formatting.
 |      
 |      Required arguments:
 |          base_sequence (string)
 |              a sequence string that represents the longest common sequence.
 |              Usually, this is the return value from:
 |                  rnav.data.set_multiple_sequence_alignment()
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
 |      Required arguments:
 |          filename (string)
 |              path to a new file to which fasta entries are written
 |      
 |      Optional arguments:
 |          which (list of integers)
 |              Sequence IDs to write to file.
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
 |  Methods defined here:
 |  
 |  __init__(self, input_data, name=None)
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Method resolution order:
 |      Data
 |      Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric, metric_defaults, read_table_kw=None, name=None)
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Parent classes:
 |      rnav.data.Data(rnav.data.Sequence)
 |  
 |  Attributes:
 |      data (str | pandas.DataFrame): dataframe storing base-pairs
 |          One row for every nucleotide position
 |          Required columns:
 |              "Nucleotide"   the nucleotide position
 |              "Sequence"     the nucleotide letter
 |              "Pair"         the nucleotide position of the base-pair
 |                             0 indicates a single stranded nucleotide
 |          Optional columns:
 |              "X_coordinate" the x position of the nucleotide in the drawing
 |              "Y_coordinate" the y position of the nucleotide in the drawing
 |      sequence (str): sequence string
 |      nts (numpy.array): "Nucleotide" column of data
 |      pair_nts (numpy.array): "Pair" column of data
 |      header (str): header information from CT file
 |      xcoordinates (numpy array): "X_coordinate" column of data
 |      ycoordinates (numpy array): "X_coordinate" column of data
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
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              A dataframe or filepath containing a secondary structure
 |              DataFrame should contain these columns:
 |                  ["Nucleotide", "Sequence", "Pair"]
 |              "Pair" column must be redundant.
 |              Filepath parsing is determined by file extension:
 |                  varna, xrna, nsd, cte, ct, dbn, bracket, json (R2DT), forna
 |  
 |  __str__(self)
 |      print the filepath and length of the RNA
 |  
 |  add_pairs(self, pairs, break_conflicting_pairs=False)
 |      Add base pairs to current secondary structure
 |      
 |      Args:
 |          pairs (list): 1-indexed list of paired residues.
 |              e.g. [(1, 20), (2, 19)]
 |  
 |  as_interactions(self, structure2=None)
 |      returns list of i, j basepairs as rnavigate.Interactions data.
 |      
 |      Args:
 |          structure2 (rnavigate.data.SecondaryStructure, optional):
 |              If provided, basepairs from both structures are included and
 |              colored by structure (left, right, or both)
 |              defaults to None.
 |  
 |  break_noncanonical_pairs(self)
 |      Removes non-canonical basepairs from the secondary structure.
 |      
 |      WARNING: this deletes information.
 |  
 |  break_pairs_nts(self, nt_positions)
 |      break base pairs at the given list of positions.
 |      
 |      Args:
 |          nt_positions (list of int): 1-indexed positions
 |  
 |  break_pairs_region(self, start, end, break_crossing=True, inverse=False)
 |      Removes pairs from the specified region (1-indexed, inclusive)
 |      
 |      WARNING: this deletes information
 |      
 |      Args:
 |          start (int): start position
 |          end (int): end position
 |          break_crossing (bool, optional): whether to keep pairs that cross
 |              over the specified region. Defaults to True.
 |          inverse (bool, optional): invert the behavior. Defaults to False.
 |  
 |  break_singleton_pairs(self)
 |      Removes singleton basepairs from the secondary structure.
 |      
 |      WARNING: This deletes information.
 |  
 |  compute_ppv_sens(self, structure2, exact=True)
 |      Compute the PPV and sensitivity between self and another
 |      SecondaryStructure object.
 |      
 |      Args:
 |          structure2 (SecondaryStructure): The SecondaryStructure to compare to.
 |          exact (bool, optional): True requires BPs to be exactly correct.
 |                                  False allows +/-1 bp slippage.
 |                                  Defaults to True.
 |      
 |      Returns:
 |          float, float, tuple: sensitivity, PPV, (TP, TP+FP, TP+FN)
 |  
 |  contact_distance(self, i, j)
 |      Returns the contact distance between positions i and j
 |  
 |  copy(self)
 |  
 |  extractPK(self, fill_mismatches=True)
 |      Returns the pk1 and pk2 pairs from the secondary structure.
 |      
 |      Ignores single base pairs. PK1 is defined as the helix crossing the
 |      most other bps. if there is a tie, the most 5' helix is called pk1
 |      returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
 |  
 |  fill_mismatches(self, mismatch=1)
 |      Adds base pairs to fill 1,1 and optionally 2,2 mismatches.
 |      
 |      mismatch = 1 will fill only 1,1 mismatches
 |      mismatch = 2 will fill 1,1 and 2,2
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new SecondaryStructure matching the alignment target.
 |      
 |      Args:
 |          alignment (data.Alignment): an alignment object used to map values
 |  
 |  get_distance_matrix(self, recalculate=False)
 |      Based on Tom's contact_distance function, but instead returns
 |      the all pairs shortest paths matrix, and stores it as an attribute. If
 |      the attribute has already been set, it returns the attribute. This is
 |      faster than calling contact_distance pairwise to fill the matrix.
 |      
 |      Args:
 |          recalculate (bool, optional): Set to true to recalculate the matrix
 |              even if the attribute is set. In case changes to the structure
 |              have been made.
 |  
 |  get_dotbracket(self)
 |      Returns a dotbracket notation representation of the secondary
 |      structure.
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
 |  get_helices(self, fill_mismatches=True, split_bulge=True, keep_singles=False)
 |      Returns a dictionary of helices from the secondary structure.
 |      Keys are equivalent to list indices. Values are lists of paired
 |      nucleotides (1-indexed) in that helix. e.g. {0:[(1,50),(2,49),(3,48)}
 |      
 |      Args:
 |          fill_mismatches (bool, optional):
 |              Whether 1-1 and 2-2 bulges are replaced with base pairs.
 |              Defaults to True.
 |          split_bulge (bool, optional):
 |              Whether to split helices on bulges.
 |              Defaults to True.
 |          keep_singles (bool, optional):
 |              Whether to return helices that contain only 1 base-pair.
 |              Defaults to False.
 |  
 |  get_human_dotbracket(self)
 |      Returns dotbracket notation string representing SecondaryStructure
 |      object. This is an experimental format designed to be more human
 |      readable, i.e. no counting of brackets required.
 |      
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
 |      
 |      Pseudoknot levels:
 |          1: Aa, Bb, Cc, etc.
 |          2: [], 3: {}, 4: <>
 |  
 |  get_interactions_df(self)
 |  
 |  get_junction_nts(self)
 |      Returns a list of residues at junctions (paired, but adjacent to
 |      an unpaired residue or the end of a chain)
 |  
 |  get_nonredundant_ct(self)
 |      Returns the ct attribute in non-redundant form - only pairs i<j
 |  
 |  get_paired_nts(self)
 |      Returns a list of residues that are paired.
 |  
 |  get_pairs(self)
 |      Returns a non-redundant list of base pairs i < j as a array of tuples.
 |      e.g., [(19,50),(20,49)....]
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
 |      Returns a list of residues that are paired.
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_ct(self, structure_number=0)
 |      Loads secondary structure information from a given ct file.
 |      
 |      Requires a properly formatted header!
 |      
 |      Args:
 |          structure_number (int, optional): If ct file contains multiple structures,
 |              uses the given structure. Defaults to 0.
 |  
 |  read_cte(self)
 |      Generates SecondaryStructure object data from a CTE file, including
 |      nucleotide x and y coordinates.
 |  
 |  read_dotbracket(self)
 |      Generates SecondaryStructure object data from a dot-bracket notation
 |      file, including nucleotide x and y coordinates.
 |  
 |  read_forna(self)
 |      Generates SecondaryStructure object data from a FORNA JSON file,
 |      including nucleotide x and y coordinates.
 |  
 |  read_nsd(self)
 |      Generates SecondaryStructure object data from an NSD file
 |      (format for RNAStructure StructureEditor), including nucleotide x and
 |      y coordinates.
 |  
 |  read_r2dt(self)
 |      Generates SecondaryStructure object data from an R2DT JSON file,
 |      including nucleotide x and y coordinates.
 |  
 |  read_varna(self)
 |      Generates SecondaryStructure object data from a VARNA file,
 |      including nucleotide x and y coordinates.
 |  
 |  read_xrna(self)
 |      Generates SecondaryStructure object data from an XRNA file,
 |      including nucleotide x and y coordinates.
 |  
 |  transform_coordinates(self, flip=None, scale=None, center=None, rotate_degrees=None)
 |      Perform transformation of structure coordinates
 |      
 |      Args:
 |          flip (str, optional):
 |              "horizontal" or "vertical".
 |              Defaults to None.
 |          scale (float, optional):
 |              new median distance of basepairs.
 |              Defaults to None.
 |          center (tuple of floats, optional):
 |              new center x and y coordinate.
 |              Defaults to None.
 |          rotate_degrees (float, optional):
 |              number of degrees to rotate structure.
 |              Defaults to None.
 |  
 |  write_ct(self, out_file)
 |      Writes a ct file from the SecondaryStructure object.
 |  
 |  write_cte(self, outputPath)
 |      writes the current structure out to CTE format for Structure Editor.
 |      
 |      Args:
 |          outputPath (string): path to output cte file to be created
 |  
 |  write_dbn(self, filename, rna_name)
 |      Writes dot-bracket notation of current structure to a new file.
 |      
 |      Args:
 |          filename (str): name of new file
 |          rna_name (str): name of this RNA for the dbn header
 |  
 |  write_sto(self, outfile, name='seq')
 |      "write structure file out into Stockholm (STO) file format
 |      for use for infernal searches
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_pairs_list(input_data, sequence) from builtins.type
 |      Creates a SecondaryStructure from a list of pairs and a sequence.
 |      
 |      Args:
 |          pairs (list): 1-indexed list of base pairs. e.g. [(1, 20), (2, 19)]
 |          sequence (str): sequence string. e.g., "AUCGUGUCAUGCUA"
 |  
 |  from_sequence(input_data) from builtins.type
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
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Methods defined here:
 |  
 |  __init__(self, x, y, pairs=None)
 |      initialize structure with arrays holding x and y coordinates. A list
 |      of base-paired positions is required only if scaling coordinates.
 |      
 |      Args:
 |          x (numpy.array): x coordinates
 |          y (numpy.array): y coordinates
 |          pairs (list of pairs, optional): list of base-paired positions.
 |  
 |  center(self, x=0, y=0)
 |      Center structure on given x, y coordinate
 |      
 |      Args:
 |          x (int, optional): x coordinate of structure center. Defaults to 0.
 |          y (int, optional): y coordinate of structure center. Defaults to 0.
 |  
 |  flip(self, horizontal=True)
 |      Flip structure vertically or horizontally
 |      
 |      Args:
 |          horizontal (bool, optional):
 |              whether to flip structure horizontally, otherwise vertically
 |              Defaults to True.
 |  
 |  get_center_point(self)
 |      return the x, y coordinates for the center of structure
 |      
 |      Returns:
 |          tuple: x and y coordinates of structure center as floats
 |  
 |  rotate(self, degrees)
 |      Rotate structure on current center point
 |      
 |      Args:
 |          degrees (float): number of degrees to rotate structure
 |  
 |  scale(self, median_bp_distance=1)
 |      Scale structure such that median base-pair distance is constant.
 |      
 |      Args:
 |          median_bp_distance (int, optional):
 |              New median distance between all base-paired nucleotides.
 |              Defaults to 1.
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
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              A dataframe or filepath containing a secondary structure
 |              DataFrame should contain these columns:
 |                  ["Nucleotide", "Sequence", "Pair"]
 |              "Pair" column must be redundant.
 |              Filepath parsing is determined by file extension:
 |                  varna, xrna, nsd, cte, ct, dbn, bracket, json (R2DT), forna
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from SecondaryStructure:
 |  
 |  __str__(self)
 |      print the filepath and length of the RNA
 |  
 |  add_pairs(self, pairs, break_conflicting_pairs=False)
 |      Add base pairs to current secondary structure
 |      
 |      Args:
 |          pairs (list): 1-indexed list of paired residues.
 |              e.g. [(1, 20), (2, 19)]
 |  
 |  as_interactions(self, structure2=None)
 |      returns list of i, j basepairs as rnavigate.Interactions data.
 |      
 |      Args:
 |          structure2 (rnavigate.data.SecondaryStructure, optional):
 |              If provided, basepairs from both structures are included and
 |              colored by structure (left, right, or both)
 |              defaults to None.
 |  
 |  break_noncanonical_pairs(self)
 |      Removes non-canonical basepairs from the secondary structure.
 |      
 |      WARNING: this deletes information.
 |  
 |  break_pairs_nts(self, nt_positions)
 |      break base pairs at the given list of positions.
 |      
 |      Args:
 |          nt_positions (list of int): 1-indexed positions
 |  
 |  break_pairs_region(self, start, end, break_crossing=True, inverse=False)
 |      Removes pairs from the specified region (1-indexed, inclusive)
 |      
 |      WARNING: this deletes information
 |      
 |      Args:
 |          start (int): start position
 |          end (int): end position
 |          break_crossing (bool, optional): whether to keep pairs that cross
 |              over the specified region. Defaults to True.
 |          inverse (bool, optional): invert the behavior. Defaults to False.
 |  
 |  break_singleton_pairs(self)
 |      Removes singleton basepairs from the secondary structure.
 |      
 |      WARNING: This deletes information.
 |  
 |  compute_ppv_sens(self, structure2, exact=True)
 |      Compute the PPV and sensitivity between self and another
 |      SecondaryStructure object.
 |      
 |      Args:
 |          structure2 (SecondaryStructure): The SecondaryStructure to compare to.
 |          exact (bool, optional): True requires BPs to be exactly correct.
 |                                  False allows +/-1 bp slippage.
 |                                  Defaults to True.
 |      
 |      Returns:
 |          float, float, tuple: sensitivity, PPV, (TP, TP+FP, TP+FN)
 |  
 |  contact_distance(self, i, j)
 |      Returns the contact distance between positions i and j
 |  
 |  copy(self)
 |  
 |  extractPK(self, fill_mismatches=True)
 |      Returns the pk1 and pk2 pairs from the secondary structure.
 |      
 |      Ignores single base pairs. PK1 is defined as the helix crossing the
 |      most other bps. if there is a tie, the most 5' helix is called pk1
 |      returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
 |  
 |  fill_mismatches(self, mismatch=1)
 |      Adds base pairs to fill 1,1 and optionally 2,2 mismatches.
 |      
 |      mismatch = 1 will fill only 1,1 mismatches
 |      mismatch = 2 will fill 1,1 and 2,2
 |  
 |  get_aligned_data(self, alignment)
 |      Returns a new SecondaryStructure matching the alignment target.
 |      
 |      Args:
 |          alignment (data.Alignment): an alignment object used to map values
 |  
 |  get_distance_matrix(self, recalculate=False)
 |      Based on Tom's contact_distance function, but instead returns
 |      the all pairs shortest paths matrix, and stores it as an attribute. If
 |      the attribute has already been set, it returns the attribute. This is
 |      faster than calling contact_distance pairwise to fill the matrix.
 |      
 |      Args:
 |          recalculate (bool, optional): Set to true to recalculate the matrix
 |              even if the attribute is set. In case changes to the structure
 |              have been made.
 |  
 |  get_dotbracket(self)
 |      Returns a dotbracket notation representation of the secondary
 |      structure.
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
 |  get_helices(self, fill_mismatches=True, split_bulge=True, keep_singles=False)
 |      Returns a dictionary of helices from the secondary structure.
 |      Keys are equivalent to list indices. Values are lists of paired
 |      nucleotides (1-indexed) in that helix. e.g. {0:[(1,50),(2,49),(3,48)}
 |      
 |      Args:
 |          fill_mismatches (bool, optional):
 |              Whether 1-1 and 2-2 bulges are replaced with base pairs.
 |              Defaults to True.
 |          split_bulge (bool, optional):
 |              Whether to split helices on bulges.
 |              Defaults to True.
 |          keep_singles (bool, optional):
 |              Whether to return helices that contain only 1 base-pair.
 |              Defaults to False.
 |  
 |  get_human_dotbracket(self)
 |      Returns dotbracket notation string representing SecondaryStructure
 |      object. This is an experimental format designed to be more human
 |      readable, i.e. no counting of brackets required.
 |      
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
 |      
 |      Pseudoknot levels:
 |          1: Aa, Bb, Cc, etc.
 |          2: [], 3: {}, 4: <>
 |  
 |  get_interactions_df(self)
 |  
 |  get_junction_nts(self)
 |      Returns a list of residues at junctions (paired, but adjacent to
 |      an unpaired residue or the end of a chain)
 |  
 |  get_nonredundant_ct(self)
 |      Returns the ct attribute in non-redundant form - only pairs i<j
 |  
 |  get_paired_nts(self)
 |      Returns a list of residues that are paired.
 |  
 |  get_pairs(self)
 |      Returns a non-redundant list of base pairs i < j as a array of tuples.
 |      e.g., [(19,50),(20,49)....]
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
 |      Returns a list of residues that are paired.
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_ct(self, structure_number=0)
 |      Loads secondary structure information from a given ct file.
 |      
 |      Requires a properly formatted header!
 |      
 |      Args:
 |          structure_number (int, optional): If ct file contains multiple structures,
 |              uses the given structure. Defaults to 0.
 |  
 |  read_cte(self)
 |      Generates SecondaryStructure object data from a CTE file, including
 |      nucleotide x and y coordinates.
 |  
 |  read_dotbracket(self)
 |      Generates SecondaryStructure object data from a dot-bracket notation
 |      file, including nucleotide x and y coordinates.
 |  
 |  read_forna(self)
 |      Generates SecondaryStructure object data from a FORNA JSON file,
 |      including nucleotide x and y coordinates.
 |  
 |  read_nsd(self)
 |      Generates SecondaryStructure object data from an NSD file
 |      (format for RNAStructure StructureEditor), including nucleotide x and
 |      y coordinates.
 |  
 |  read_r2dt(self)
 |      Generates SecondaryStructure object data from an R2DT JSON file,
 |      including nucleotide x and y coordinates.
 |  
 |  read_varna(self)
 |      Generates SecondaryStructure object data from a VARNA file,
 |      including nucleotide x and y coordinates.
 |  
 |  read_xrna(self)
 |      Generates SecondaryStructure object data from an XRNA file,
 |      including nucleotide x and y coordinates.
 |  
 |  transform_coordinates(self, flip=None, scale=None, center=None, rotate_degrees=None)
 |      Perform transformation of structure coordinates
 |      
 |      Args:
 |          flip (str, optional):
 |              "horizontal" or "vertical".
 |              Defaults to None.
 |          scale (float, optional):
 |              new median distance of basepairs.
 |              Defaults to None.
 |          center (tuple of floats, optional):
 |              new center x and y coordinate.
 |              Defaults to None.
 |          rotate_degrees (float, optional):
 |              number of degrees to rotate structure.
 |              Defaults to None.
 |  
 |  write_ct(self, out_file)
 |      Writes a ct file from the SecondaryStructure object.
 |  
 |  write_cte(self, outputPath)
 |      writes the current structure out to CTE format for Structure Editor.
 |      
 |      Args:
 |          outputPath (string): path to output cte file to be created
 |  
 |  write_dbn(self, filename, rna_name)
 |      Writes dot-bracket notation of current structure to a new file.
 |      
 |      Args:
 |          filename (str): name of new file
 |          rna_name (str): name of this RNA for the dbn header
 |  
 |  write_sto(self, outfile, name='seq')
 |      "write structure file out into Stockholm (STO) file format
 |      for use for infernal searches
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from SecondaryStructure:
 |  
 |  from_pairs_list(input_data, sequence) from builtins.type
 |      Creates a SecondaryStructure from a list of pairs and a sequence.
 |      
 |      Args:
 |          pairs (list): 1-indexed list of base pairs. e.g. [(1, 20), (2, 19)]
 |          sequence (str): sequence string. e.g., "AUCGUGUCAUGCUA"
 |  
 |  from_sequence(input_data) from builtins.type
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
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
    
    Required arguments:
        sequence1 (string)
            the first sequence
        sequence2 (string)
            the second sequence
        alignment1 (string)
            first sequence, plus dashes "-" indicating indels
        alignment2 (string)
            second sequence, plus dashes "-" indicating indels
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
    
    Required arguments:
        fasta (string)
            location of Pearson fasta file
    
    Optional arguments:
        set_pairwise (True or False)
            whether to set every pairwise alignment as well as the multiple
            sequence alignment.
            Defaults to False
```

### rnavigate.data.lookup_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function lookup_alignment in module rnavigate.data.alignments

lookup_alignment(sequence1, sequence2, t_or_u='U')
    look up a previously set alignment in the _alignments_cache
    
    Required arguments:
        sequence1 (string)
            The first sequence to align
        sequence2 (string)
            The second sequence to be aligned to
    
    Optional arguments:
        t_or_u ("T", "U", or False)
            "T" converts "U"s to "T"s
            "U" converts "U"s to "T"s
            False does nothing
            defaults to "U"
    
    Returns:
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
 |  Attributes:
 |      sequence1 (str): the sequence to be aligned
 |      sequence2 (str): the sequence to align to
 |      alignment1 (str): the alignment string matching sequence1 to sequence2
 |      alignment2 (str): the alignment string matching sequence2 to sequence1
 |      starting_sequence (str): sequence1
 |      target_sequence(str): sequence2 if full is False, else alignment2
 |      mapping (numpy.array): the alignment map array.
 |          index of starting_sequence is mapping[index] of target_sequence
 |  
 |  Methods:
 |      All map_functions map from starting sequence to target sequence.
 |      map_values: maps per-nucleotide values
 |      map_indices: maps a list of indices
 |      map_positions: maps a list of positions
 |      map_dataframe: maps a dataframe with multiple position columns
 |          (rows that cannot be mapped are dropped)
 |      map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
 |          (rows that cannot be mapped are dropped)
 |          (missing rows filled with NaN)
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
 |      Args:
 |          sequence1 (str): the starting sequence
 |          sequence2 (str): the target sequence
 |          full (bool, optional): whether to keep unmapped starting sequence
 |              positions. Defaults to False.
 |  
 |  __repr__(self)
 |      a nice text only representation of an alignment
 |  
 |  get_alignment(self)
 |      Gets an alignment that has either been user-defined or previously
 |      calculated or produces a new pairwise alignment between two sequences.
 |      
 |      Returns:
 |          (tuple of 2 str): alignment1 and alignment2
 |  
 |  get_inverse_alignment(self)
 |      Alignments require a method to get the inverted alignment
 |  
 |  get_mapping(self)
 |      Calculates a mapping from starting sequence to target sequence.
 |      
 |      Returns:
 |          numpy.array: an array of length of starting sequence that maps to
 |              an index of target sequence. Stored as self.mapping
 |              starting_sequence[idx] == target_sequence[self.mapping[idx]]
 |  
 |  print(self, print_format='full')
 |      Print the alignment in a human-readable format.
 |      
 |      Arguments:
 |          print_format (string)
 |              how to format the alignment.
 |              "full": the full length alignment with changes labeled "X"
 |              "cigar": the CIGAR string
 |              "long": locations and sequences of each change
 |              "short": total number of matches, mismatches, and indels
 |              Defaults to "full".
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
 |      Unmapped positions are dropped.
 |      
 |      Args:
 |          dataframe (pandas.DataFrame): a dataframe with position columns
 |          position_columns (list of str): a list of columns containing
 |              positions to map
 |      
 |      Returns:
 |          pandas.DataFrame: a new dataframe (copy) with position columns
 |              mapped or dropped
 |  
 |  map_indices(self, indices, keep_minus_one=True)
 |      Takes a list of indices (0-index) and maps them to target sequence
 |      
 |      Args:
 |          indices (int | list): a single or list of integer indices
 |          keep_minus_one (bool, optional): whether to keep unmapped starting
 |              sequence indices (-1) in the returned array. Defaults to True.
 |      
 |      Returns:
 |          numpy.array: the equivalent indices in target sequence
 |  
 |  map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide', sequence_column='Sequence')
 |      Takes a per-nt dataframe and map it to the target sequence.
 |      
 |      Dataframe must have 1 row per nucleotide in starting sequence,
 |      with a position column and a sequence column. Dataframe is
 |      mapped to have the same format, but for target sequence nucleotides and
 |      positions.
 |      
 |      Required arguments:
 |          dataframe (pandas.DataFrame)
 |              a per-nucleotide dataframe
 |      
 |      Optional arguments
 |          position_column (string)
 |              name of the position column.
 |              Defaults to "Nucleotide".
 |          sequence_column (string)
 |              name of the sequence column.
 |              Defaults to "Sequence".
 |      
 |      Returns:
 |          pandas.DataFrame: a new dataframe (copy) mapped to target sequence.
 |              Unmapped starting sequence positions are dropped and unmapped
 |              target sequence positions are filled.
 |  
 |  map_positions(self, positions, keep_zero=True)
 |      Takes a list of positions (1-index) and maps them to target sequence
 |      
 |      Args:
 |          positions (int | list): a single or list of integer positions
 |          keep_zero (bool, optional): whether to keep unmapped starting
 |              sequence positions (0) in the returned array. Defaults to True.
 |      
 |      Returns:
 |          numpy.array: the equivalent positions in target sequence
 |  
 |  map_values(self, values, fill=nan)
 |      Takes an array of length equal to starting sequence and maps them to
 |      target sequence, unmapped positions in starting sequence are dropped
 |      and unmapped positions in target sequence are filled with fill value.
 |      
 |      Args:
 |          values (iterable): values to map to target sequence
 |          fill (any, optional): a value for unmapped positions in target.
 |          Defaults to numpy.nan.
 |      
 |      Returns:
 |          numpy.array: an array of values equal in length to target sequence
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
 |  Attributes:
 |      alignments (list): the constituent alignments
 |      starting_sequence (str): starting sequence of alignments[0]
 |      target_sequence (str): target sequence of alignments[-1]
 |      mapping (numpy.array): a vector that maps from starting to target
 |          index of starting_sequence is mapping[index] of target sequence
 |  
 |  Methods:
 |      All map_functions map from starting sequence to target sequence.
 |      map_values: maps per-nucleotide values
 |      map_indices: maps a list of indices
 |      map_positions: maps a list of positions
 |      map_dataframe: maps a dataframe with multiple position columns
 |          (rows that cannot be mapped are dropped)
 |      map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
 |          (rows that cannot be mapped are dropped)
 |          (missing rows filled with NaN)
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
 |      Creates a single alignment from multiple alignments
 |      
 |      Raises:
 |          ValueError: if the target sequence of one alignment doesn't match
 |              the starting sequence of the next.
 |  
 |  get_inverse_alignment(self)
 |      Alignments require a method to get the inverted alignment
 |  
 |  get_mapping(self)
 |      combines mappings from each alignment.
 |      
 |      Returns:
 |          numpy.array: a mapping from initial starting sequence to final
 |              target sequence
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
 |      Unmapped positions are dropped.
 |      
 |      Args:
 |          dataframe (pandas.DataFrame): a dataframe with position columns
 |          position_columns (list of str): a list of columns containing
 |              positions to map
 |      
 |      Returns:
 |          pandas.DataFrame: a new dataframe (copy) with position columns
 |              mapped or dropped
 |  
 |  map_indices(self, indices, keep_minus_one=True)
 |      Takes a list of indices (0-index) and maps them to target sequence
 |      
 |      Args:
 |          indices (int | list): a single or list of integer indices
 |          keep_minus_one (bool, optional): whether to keep unmapped starting
 |              sequence indices (-1) in the returned array. Defaults to True.
 |      
 |      Returns:
 |          numpy.array: the equivalent indices in target sequence
 |  
 |  map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide', sequence_column='Sequence')
 |      Takes a per-nt dataframe and map it to the target sequence.
 |      
 |      Dataframe must have 1 row per nucleotide in starting sequence,
 |      with a position column and a sequence column. Dataframe is
 |      mapped to have the same format, but for target sequence nucleotides and
 |      positions.
 |      
 |      Required arguments:
 |          dataframe (pandas.DataFrame)
 |              a per-nucleotide dataframe
 |      
 |      Optional arguments
 |          position_column (string)
 |              name of the position column.
 |              Defaults to "Nucleotide".
 |          sequence_column (string)
 |              name of the sequence column.
 |              Defaults to "Sequence".
 |      
 |      Returns:
 |          pandas.DataFrame: a new dataframe (copy) mapped to target sequence.
 |              Unmapped starting sequence positions are dropped and unmapped
 |              target sequence positions are filled.
 |  
 |  map_positions(self, positions, keep_zero=True)
 |      Takes a list of positions (1-index) and maps them to target sequence
 |      
 |      Args:
 |          positions (int | list): a single or list of integer positions
 |          keep_zero (bool, optional): whether to keep unmapped starting
 |              sequence positions (0) in the returned array. Defaults to True.
 |      
 |      Returns:
 |          numpy.array: the equivalent positions in target sequence
 |  
 |  map_values(self, values, fill=nan)
 |      Takes an array of length equal to starting sequence and maps them to
 |      target sequence, unmapped positions in starting sequence are dropped
 |      and unmapped positions in target sequence are filled with fill value.
 |      
 |      Args:
 |          values (iterable): values to map to target sequence
 |          fill (any, optional): a value for unmapped positions in target.
 |          Defaults to numpy.nan.
 |      
 |      Returns:
 |          numpy.array: an array of values equal in length to target sequence
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
 |  Attributes:
 |      sequence1 (str): the sequence to be aligned
 |      sequence2 (str): the sequence to align to
 |      structure1
 |      structure2
 |      aa_sequence1
 |      aa_sequence2
 |      alignment1 (str): the alignment string matching sequence1 to sequence2
 |      alignment2 (str): the alignment string matching sequence2 to sequence1
 |      starting_sequence (str): sequence1
 |      target_sequence(str): sequence2 if full is False, else alignment2
 |      mapping (numpy.array): the alignment map array.
 |          index of starting_sequence is mapping[index] of target_sequence
 |  
 |  Methods:
 |      All map_functions map from starting sequence to target sequence.
 |      map_values: maps per-nucleotide values
 |      map_indices: maps a list of indices
 |      map_positions: maps a list of positions
 |      map_dataframe: maps a dataframe with multiple position columns
 |          (rows that cannot be mapped are dropped)
 |      map_nucleotide_dataframe: maps a dataframe with 1 row per nucleotide
 |          (rows that cannot be mapped are dropped)
 |          (missing rows filled with NaN)
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
 |      Args:
 |          structure1 (str): the starting structure
 |          structure2 (str): the target structure
 |          full (bool, optional): whether to keep unmapped starting sequence
 |              positions. Defaults to False.
 |  
 |  get_alignment(self)
 |      Aligns pseudo-amino-acid sequences according to RNAlign2D rules.
 |      
 |      Returns:
 |          (tuple of 2 str): alignment1 and alignment2
 |  
 |  get_inverse_alignment(self)
 |      Alignments require a method to get the inverted alignment
 |  
 |  get_mapping(self)
 |      Calculates a mapping from starting sequence to target sequence.
 |      
 |      Returns:
 |          numpy.array: an array of length of starting sequence that maps to
 |              an index of target sequence. Stored as self.mapping
 |              starting_sequence[idx] == target_sequence[self.mapping[idx]]
 |  
 |  set_as_default_alignment(self)
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
 |      Unmapped positions are dropped.
 |      
 |      Args:
 |          dataframe (pandas.DataFrame): a dataframe with position columns
 |          position_columns (list of str): a list of columns containing
 |              positions to map
 |      
 |      Returns:
 |          pandas.DataFrame: a new dataframe (copy) with position columns
 |              mapped or dropped
 |  
 |  map_indices(self, indices, keep_minus_one=True)
 |      Takes a list of indices (0-index) and maps them to target sequence
 |      
 |      Args:
 |          indices (int | list): a single or list of integer indices
 |          keep_minus_one (bool, optional): whether to keep unmapped starting
 |              sequence indices (-1) in the returned array. Defaults to True.
 |      
 |      Returns:
 |          numpy.array: the equivalent indices in target sequence
 |  
 |  map_nucleotide_dataframe(self, dataframe, position_column='Nucleotide', sequence_column='Sequence')
 |      Takes a per-nt dataframe and map it to the target sequence.
 |      
 |      Dataframe must have 1 row per nucleotide in starting sequence,
 |      with a position column and a sequence column. Dataframe is
 |      mapped to have the same format, but for target sequence nucleotides and
 |      positions.
 |      
 |      Required arguments:
 |          dataframe (pandas.DataFrame)
 |              a per-nucleotide dataframe
 |      
 |      Optional arguments
 |          position_column (string)
 |              name of the position column.
 |              Defaults to "Nucleotide".
 |          sequence_column (string)
 |              name of the sequence column.
 |              Defaults to "Sequence".
 |      
 |      Returns:
 |          pandas.DataFrame: a new dataframe (copy) mapped to target sequence.
 |              Unmapped starting sequence positions are dropped and unmapped
 |              target sequence positions are filled.
 |  
 |  map_positions(self, positions, keep_zero=True)
 |      Takes a list of positions (1-index) and maps them to target sequence
 |      
 |      Args:
 |          positions (int | list): a single or list of integer positions
 |          keep_zero (bool, optional): whether to keep unmapped starting
 |              sequence positions (0) in the returned array. Defaults to True.
 |      
 |      Returns:
 |          numpy.array: the equivalent positions in target sequence
 |  
 |  map_values(self, values, fill=nan)
 |      Takes an array of length equal to starting sequence and maps them to
 |      target sequence, unmapped positions in starting sequence are dropped
 |      and unmapped positions in target sequence are filled with fill value.
 |      
 |      Args:
 |          values (iterable): values to map to target sequence
 |          fill (any, optional): a value for unmapped positions in target.
 |          Defaults to numpy.nan.
 |      
 |      Returns:
 |          numpy.array: an array of values equal in length to target sequence
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
 |  Method resolution order:
 |      ScalarMappable
 |      matplotlib.cm.ScalarMappable
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, cmap, normalization, values, title='', tick_labels=None, **cbar_args)
 |      Parameters
 |      ----------
 |      norm : `matplotlib.colors.Normalize` (or subclass thereof)
 |          The normalizing object which scales data, typically into the
 |          interval ``[0, 1]``.
 |          If *None*, *norm* defaults to a *colors.Normalize* object which
 |          initializes its scaling based on the first data processed.
 |      cmap : str or `~matplotlib.colors.Colormap`
 |          The colormap used to map normalized data values to RGBA colors.
 |  
 |  get_cmap(self, cmap)
 |      Given a matplotlib color, list of colors, or colormap name, return
 |      a colormap object
 |      
 |      Args:
 |          cmap (str | tuple | float | list): A valid mpl color, list of valid
 |              colors or a valid colormap name
 |      
 |      Returns:
 |          matplotlib colormap: listed colormap matching the input
 |  
 |  get_norm(self, normalization, values, cmap)
 |  
 |  is_equivalent_to(self, cmap2)
 |  
 |  values_to_hexcolors(self, values, alpha=1.0)
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
 |  Method resolution order:
 |      Interactions
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, metric, metric_defaults, read_table_kw=None, window=1, name=None)
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
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
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Parses a deletions.txt file and stores data as a dataframe at
 |      self.data, sets self.window=1, and calculates a "Percentile" column.
 |      
 |      Args:
 |          input_data (str): path to deletions.txt file
 |          read_table_kw (dict, optional): kwargs passed to pandas.read_table().
 |              Defaults to None.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
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
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  data_specific_filter(self, positive_only=False, negative_only=False, **kwargs)
 |      Adds filters for "Sign" column to parent filter() function
 |      
 |      Args:
 |          positive_only (bool, optional): whether to require that sign is 1.
 |              Defaults to False.
 |          negative_only (bool, optional): whether to require that sign is -1.
 |              Defaults to False.
 |      
 |      Returns:
 |          dict: any additional keyword-argument pairs are returned
 |  
 |  get_sorted_data(self)
 |      Overwrites parent get_normalized_ij_data. Uses these values instead:
 |          self.data[self.metric] * self.data["Sign"]
 |      Except when self.metric is "Distance".
 |      
 |      Args:
 |          min_max (list of float): min and max values for normalization
 |      
 |      Returns:
 |          numpy array: normalized values for color mapping
 |  
 |  read_file(self, filepath, read_table_kw=None)
 |      Parses a correlations file and stores data as a dataframe at
 |      self.data, sets self.window=1, and renames "+/-" column to "Sign".
 |      
 |      Args:
 |          filepath (str):
 |              path to a RingMapper correlations file
 |          read_table_kw (dict, optional):
 |              kwargs passed to pandas.read_table().
 |              Defaults to None.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  data_specific_filter(self, all_pairs=False, **kwargs)
 |      Used by Interactions.filter(). By default, non-primary and
 |      -secondary pairs are removed. all_pairs=True changes this behavior.
 |      
 |      Args:
 |          all_pairs (bool, optional): whether to include all PAIRs.
 |              Defaults to False.
 |      
 |      Returns:
 |          dict: remaining kwargs are passed back to Interactions.filter()
 |  
 |  get_sorted_data(self)
 |      Same as parent function, unless metric is set to "Class", in which
 |      case ij pairs are returned in a different order.
 |      
 |      Args:
 |          min_max (list of int, length 2): minimum and maximum bounds for
 |              colormapping
 |      
 |      Returns:
 |          pandas DataFrame: Dataframe providing i, j, and normalized data
 |              values for plotting
 |  
 |  read_file(self, filepath, read_table_kw=None)
 |      Parses a pairmap.txt file and stores data as a dataframe at
 |      self.data, sets self.window (usually 3, from header).
 |      
 |      Args:
 |          filepath (str): path to a PAIR-MaP pairmap.txt file
 |          read_table_kw (dict, optional): kwargs passed to pandas.read_table().
 |              Defaults to None.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Constructs Interactions data from a pairing probability text file
 |      containing i, j, and -log10(P) values. Can be obtained using partition
 |      and ProbabilityPlot functions from RNAStructure (Matthews Lab).
 |      
 |      Args:
 |          input_data (str): path to pairing probability text file
 |          sequence (str, optional): Sequence string. Defaults to None.
 |  
 |  data_specific_filter(self, **kwargs)
 |      Used by parent filter function. By default, filters out pairs with
 |      probability less that 3%
 |      
 |      Returns:
 |          dict: keyword arguments are passed back to Interactions.filter()
 |  
 |  get_entropy_profile(self, print_out=False, save_file=None)
 |      Calculates per-nucleotide Shannon entropy and stores as self.entropy
 |      
 |      Args:
 |          print_out (bool, optional): whether to print the result.
 |              Defaults to False.
 |          save_file (str, optional): file to write result to.
 |              Defaults to None.
 |  
 |  read_file(self, filepath, read_table_kw=None)
 |      Parses a pairing probability text file to create a DataFrame
 |      containing i, j, -log10(P) and Probability (0-1).
 |      
 |      Args:
 |          filepath (str): path to pairing probability text file
 |          read_table_kw (None, optional): ignored. Defaults to None.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
 |  
 |  filter(self, prefiltered=False, reset_filter=True, structure=None, min_cd=None, max_cd=None, paired_only=False, ss_only=False, ds_only=False, profile=None, min_profile=None, max_profile=None, compliments_only=False, nts=None, max_distance=None, min_distance=None, exclude_nts=None, isolate_nts=None, resolve_conflicts=None, **kwargs)
 |      Convenience function that applies the above filters simultaneously.
 |      
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
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
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
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
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
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
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Given a dataframe or a data file, construct the interactions object
 |      
 |      Args:
 |          input_data (str | pandas.DataFrame):
 |              path to a file or dataframe containing interactions data.
 |              Must have at least "i" and "j" columns indicating the 5' and 3'
 |              ends of the interactions.
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |          default_metric (str, optional): column name to use as the default
 |              metric. Defaults to None.
 |          read_table_kw (dict, optional): other options for read_table.
 |              Defaults to {}.
 |          window (int, optional): 5' and 3' interactions windows.
 |              Defaults to 1.
 |          fasta (str, optional): path to fasta file. Defaults to None.
 |          fill (dict, optional): dictionary specifying a fill value (values)
 |              to use with a metric (keys). Defaults to {}.
 |          cmaps (dict, optional): specifies cmaps (values) to use with
 |              metrics (keys). Defaults to {}.
 |          mins_maxes (dict, optional): specifies minimum and maximum values
 |              as a list of floats (values) to use with given metric (keys).
 |              Defaults to {}.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Interactions:
 |  
 |  copy(self, apply_filter=False)
 |      Returns a deep copy of the Interactions.
 |      
 |      Optional arguments:
 |          apply_filter (True or False)
 |              whether to remove masked rows
 |              Defaults to False
 |      
 |      Returns:
 |          rnavigate.data.Interactions
 |              The same subclass as the original object
 |  
 |  count_filter(self, **kwargs)
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
 |      Args:
 |          prefiltered (bool, optional):
 |              passed to self.set_mask_update().
 |              Defaults to False.
 |          structure (SecondaryStructure, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          min_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          max_cd (int, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to None.
 |          paired_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ss_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          ds_only (bool, optional):
 |              passed to self.mask_on_structure().
 |              Defaults to False.
 |          profile (Profile or subclass, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          min_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          max_profile (float, optional):
 |              passed to self.mask_on_profile().
 |              Defaults to None.
 |          compliments_only (bool, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to False.
 |          nts (str, optional):
 |              passed to self.mask_on_sequence().
 |              Defaults to None.
 |          max_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          min_distance (int, optional):
 |              passed to self.mask_on_distance().
 |              Defaults to None.
 |          exclude_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          isolate_nts (list of int, optional):
 |              passed to self.mask_on_position().
 |              Defaults to None.
 |          resolve_conflicts (str, optional):
 |              passed to self.resolve_conflicts().
 |              Defaults to None.
 |          **kwargs:
 |              additional arguments are first passed to
 |              self.data_specific_filter(), remaining kwargs are passed to
 |              self.mask_on_values()
 |  
 |  get_aligned_data(self, alignment, apply_filter=True)
 |      Get a new copy of the data with i and j mapped to new positions
 |      using an alignment. Interactions in which i or j does not map are
 |      dropped.
 |  
 |  get_ij_colors(self)
 |      Gets i, j, and colors lists for plotting interactions. i and j are
 |      the 5' and 3' ends of each interaction, and colors is the color to use
 |      for each interaction. Values of self.data[self.metric] are normalized
 |      to 0 to 1, which correspond to self.min_max values. These are then
 |      mapped to a color using self.cmap.
 |      
 |      Returns:
 |          list, list, list: 5' and 3' ends of each pair, color for each pair
 |  
 |  get_sorted_data(self)
 |      Returns a sorted copy of interactions data.
 |      
 |      Returns:
 |          pandas.DataFrame: i, j, and metric values, sorted
 |  
 |  mask_on_distance(self, max_dist=None, min_dist=None)
 |      Mask interactions based on their primary sequence distance (j-i).
 |      
 |      Args:
 |          max_dist (int): maximum allowable distance
 |          min_dist (int): minimum allowable distance
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_position(self, exclude=None, isolate=None)
 |      Masks interactions based on position in sequence
 |      
 |      Args:
 |          exclude (list of int, optional): a list of nucleotide positions to
 |              exclude if i or j is in list
 |          isolate (list of int, optional): a list of nucleotide positions to
 |              isolate if i and j are not in list
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_profile(self, profile, min_profile=None, max_profile=None)
 |      Masks interactions based on per-nucleotide information. Positions
 |      that are not mapped to profile or are np.nan values in profile are not
 |      masked.
 |      
 |      Args:
 |          profile (Profile or subclass): a data object containing
 |              per-nucleotide information
 |          min_profile (float, optional): minimum allowable per-nucleotide
 |              value. Defaults to None.
 |          max_profile (float, optional): maximum allowable per-nucleotide
 |              value. Defaults to None.
 |      
 |      Returns:
 |          numpy array: the mask array
 |  
 |  mask_on_sequence(self, compliment_only=None, nts=None)
 |      Mask interactions based on sequence content
 |      
 |      Args:
 |          compliment_only (bool): require that i and j windows are reverse
 |              complimentary
 |          nts (str): require that all nucleotides in i and j windows are in
 |              nts
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_structure(self, structure, min_cd=None, max_cd=None, ss_only=False, ds_only=False, paired_only=False)
 |      Mask interactions based on secondary structure
 |      
 |      Args:
 |          structure (SecondaryStructure):
 |              A SecondaryStructure object or list of these objects
 |              Filters are applied based on all structures.
 |          min_cd (int, optional): minimum allowable contact distance.
 |              Defaults to None.
 |          max_cd (int, optional): maximum allowable contact distance.
 |              Defaults to None.
 |          ss_only (bool, optional): whether to require i and j to be single-
 |              stranded. Defaults to False.
 |          ds_only (bool, optional): whether to require i and j to be double-
 |              stranded. Defaults to False.
 |          paired_only (bool, optional): whether to require that i and j are
 |              base paired. Defaults to False.
 |      
 |      Returns:
 |          numpy array: the masking boolean array
 |  
 |  mask_on_values(self, **kwargs)
 |      Mask interactions on values in self.data. Each keyword should have
 |      the format "column_operator" where column is a valid column name of
 |      the dataframe and operator is one of:
 |          "ge": greater than or equal to
 |          "le": less than or equal to
 |          "gt": greater than
 |          "lt": less than
 |          "eq": equal to
 |          "ne": not equal to
 |      The values given to these keywords are then used in the comparison and
 |      False comparisons are filtered out. e.g.:
 |          self.mask_on_values(Statistic_ge=23) evaluates to:
 |          self.update_mask(self.data["Statistic"] >= 23)
 |  
 |  print_new_file(self, outfile=None)
 |      Prints a new file containing repositioned and filtered interactions
 |      in the original format
 |      
 |      Args:
 |          outfile (str, optional): path to an output file. If None, file
 |              string is printed to console. Defaults to None.
 |  
 |  reset_mask(self)
 |      Resets the mask to all True (removes previous filters)
 |  
 |  resolve_conflicts(self, metric=None)
 |      Resolves conflicting windows using the Maximal Weighted Independent
 |      Set. The weights are taken from the metric value. The graph is first
 |      broken into components to speed up the identification of the MWIS. Then
 |      the mask is updated to only include the MWIS.
 |  
 |  set_3d_distances(self, pdb, atom)
 |      Creates or overwrites values in self.data["Distance"] by calculating
 |      the distance between atoms in i and j in the PDB structure.
 |      
 |      Args:
 |          pdb (PDB): a data object containing atomic coordinates
 |          atom (str): an atom id
 |  
 |  update_mask(self, mask)
 |      Given a new masking array, the mask is updated
 |      
 |      Args:
 |          mask (numpy array of bool): the new mask to update on
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Stores a tertiary structure with atomic coordinates
 |  
 |  PDB data can come from a PDB or CIF file.
 |  Note: this is the only object that cannot be aligned to another sequence.
 |  
 |  Attributes:
 |      data (Bio.PDB object)
 |          the structure from the PDB or CIF file
 |      sequence (string)
 |          the sequence of the structure, not always in the same coordinates
 |      offset (integer)
 |          the difference between sequence and residue positions
 |      chain (string)
 |          the string chain identifier for the RNA
 |      path (string)
 |          the path to the PDB or CIF file
 |      distance_matrix (dictionary of 2D Numpy array)
 |          keys are atom identifiers, values are the pairwise atom distance
 |          matrices between residues. These are only stored if computed.
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
 |      Required arguments:
 |          input_data (string)
 |              path to a PDB or CIF file
 |          chain (string)
 |              chain identifier of RNA of interest
 |      
 |      Optional arguments:
 |          sequence (rnavigate.Sequence or string)
 |              A sequence to use as the reference sequence.
 |              This is required if the sequence cannot be found in the header
 |              Defaults to None.
 |  
 |  get_distance(self, i, j, atom="O2'")
 |      Get the atomic distance between nucleotides i and j (1-indexed).
 |  
 |  get_distance_matrix(self, atom="O2'")
 |      Get the pairwise atomic distance matrix for all residues.
 |  
 |  get_pdb_idx(self, seq_idx)
 |      Return the PDB index given the sequence index.
 |  
 |  get_seq_idx(self, pdb_idx)
 |      Return the sequence index given the PDB index.
 |  
 |  get_sequence(self, pdb)
 |      Find the sequence in the provided CIF or PDB file.
 |  
 |  get_sequence_from_seqres(self, seqres)
 |      Used by get_sequence to parse the SEQRES entries.
 |  
 |  get_xyz_coord(self, nt, atom)
 |      Return the x, y, and z coordinates for a given residue and atom.
 |  
 |  is_valid_idx(self, pdb_idx=None, seq_idx=None)
 |      Determines if a PDB or sequence index is in the PDB structure.
 |  
 |  read_pdb(self, pdb)
 |      Read a PDB or CIF file into the data structure.
 |  
 |  set_indices(self)
 |      Uses self.data and self.sequence to set self.offset
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Method resolution order:
 |      Profile
 |      rnavigate.data.data.Data
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, metric='default', metric_defaults=None, read_table_kw=None, sequence=None, name=None)
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  recreation_kwargs
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  read_log(self, log)
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_rnaframework(input_data, normalize=None) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Profile:
 |  
 |  recreation_kwargs
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  read_file(self, input_data, read_table_kw={})
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties defined here:
 |  
 |  recreation_kwargs
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from SHAPEMaP:
 |  
 |  read_log(self, log)
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from SHAPEMaP:
 |  
 |  from_rnaframework(input_data, normalize=None) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Profile:
 |  
 |  recreation_kwargs
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from Profile:
 |  
 |  calculate_gini_index(self, values)
 |      Calculate the Gini index of an array of values.
 |  
 |  calculate_windows(self, column, window, method='median', new_name=None, minimum_points=None, mask_na=True)
 |      calculates a windowed operation over a column of self.data and
 |      stores the result as a new column. Value of each window is assigned to
 |      the center position of the window.
 |      
 |      Args:
 |          column (str): name of column to perform operation on
 |          window (int): window size, must be an odd number
 |          method (str, optional): operation to perform over windows, must be
 |              one of 'median', 'mean', 'minimum', 'maximum'
 |              Defaults to 'median'.
 |          new_name (str, optional): name of new column for stored result.
 |              Defaults to f"{method}_{window}_nt", e.g. "median_55_nt".
 |          minimum_points (int, optional): minimum number of points within
 |              each window.
 |              Defaults to the size of the window.
 |  
 |  copy(self)
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_plotting_dataframe(self)
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
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_eDMS(self, values)
 |      Returns normalization factors for normalize values following eDMS
 |      pernt scheme in ShapeMapper 2.2
 |      
 |      Args:
 |          values (1D numpy array): values to scale
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  norm_percentiles(self, values, lower_bound=90, upper_bound=99, median_or_mean='mean')
 |      Calculates profile scaling factors and error propagation by scaling
 |      the median between upper and lower bound percentiles to 1.
 |      
 |      Args:
 |          values (1D numpy.array): values to scale
 |          lower_bound (int or float, optional): percentile of lower bound
 |              Defaults to 90
 |          upper_bound (int or float, optional): percentile of upper bound
 |              Defaults to 99
 |      
 |      Returns:
 |          (float, float): scaling factor and error propagation factor
 |  
 |  normalize(self, profile_column=None, new_profile=None, error_column=None, new_error=None, norm_method=None, nt_groups=None, profile_factors=None, **norm_kwargs)
 |      Normalize values in a column, and store in a new column.
 |      
 |      By default, performs ShapeMapper2 boxplot normalization on self.metric
 |      and stores the result as "Norm_profile".
 |      
 |      Optional arguments:
 |          profile_column (string)
 |              column name of values to normalize
 |              Defaults to self.metric
 |          new_profile (string)
 |              column name of new normalized values
 |              Defaults to "Norm_profile"
 |          error_column (string)
 |              column name of error values to propagate
 |              Defaults to self.error_column
 |          new_error (string)
 |              column name of new propagated error values
 |              Defaults to "Norm_error"
 |          norm_method (string)
 |              normalization method to use.
 |              "DMS" uses self.norm_percentile and nt_groups=['AC', 'UG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  As and Cs are normalized seperately from Us and Gs
 |              "eDMS" uses self.norm_eDMS and  nt_groups=['A', 'U', 'C', 'G']
 |                  Applies the new eDMS-MaP normalization.
 |                  Each nucleotide is normalized seperately.
 |              "boxplot" uses self.norm_boxplot and nt_groups=['AUCG']
 |                  removes outliers (> 1.5 iqr) and scales median to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              "percentile" uses self.norm_percentile and nt_groups=['AUCG']
 |                  scales the median of 90th to 95th percentiles to 1
 |                  scales nucleotides together unless specified with nt_groups
 |              Defaults to "boxplot": the default normalization of ShapeMapper
 |          nt_groups (list of strings)
 |              A list of nucleotides to group
 |              e.g. ['AUCG'] groups all nts together
 |                   ['AC', 'UG'] groups As with Cs and Us with Gs
 |                   ['A', 'C', 'U', 'G'] scales each nt seperately
 |              Default depends on norm_method
 |          profile_factors (dictionary)
 |              a scaling factor (float) for each nucleotide. keys must be:
 |                  'A', 'C', 'U', 'G'
 |              Note: using this argument overrides any calculation of scaling
 |              Defaults to None
 |          **norm_kwargs: these are passed to the norm_method function
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_external(self, profiles, **kwargs)
 |      normalize reactivities using other profiles to normfactors.
 |      
 |      Args:
 |          profiles (list of rnavigate.data.Profile): a list of other profiles
 |              used to compute scaling factors
 |      
 |      Returns:
 |          dict: the new profile scaling factors dictionary
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  winsorize(self, column, lower_bound=None, upper_bound=None)
 |      Winsorize the data between bounds.
 |      
 |      If either bound is set to None, one-sided Winsorization is performed.
 |      
 |      Required arguments:
 |          column (string)
 |              the column of data to be winsorized
 |          lower_bound (Number or None)
 |              Data below this value is set to this value.
 |              If None, no lower bound is applied.
 |          upper_bound (Number or None)
 |              Data above this value is set to this value.
 |              If None, no upper bound is applied.
 |  
 |  ----------------------------------------------------------------------
 |  Class methods inherited from Profile:
 |  
 |  from_array(input_data, sequence, **kwargs) from builtins.type
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from Profile:
 |  
 |  recreation_kwargs
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Data:
 |  
 |  add_metric_defaults(self, metric_defaults)
 |  
 |  read_file(self, filepath, read_table_kw)
 |      Convert data file to pandas dataframe and store as self.data
 |      
 |      Args:
 |          filepath (str): path to data file containing interactions
 |          read_table_kw (dict): kwargs dictionary passed to pd.read_table
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Data:
 |  
 |  cmap
 |  
 |  color_column
 |  
 |  colors
 |  
 |  error_column
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from rnavigate.data.data.Data:
 |  
 |  metric
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Attributes:
 |      data (Pandas DataFrame): stores the list of sites or regions
 |      name (string): the label for this annotation for use on plots
 |      color (valid matplotlib color): color to represent annotation on plots
 |      sequence (string): the reference sequence string
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
 |      Create an annotation from a list of sites or regions.
 |      
 |      Args:
 |          input_data (list):
 |              List will be treated according to annotation_type argument.
 |              Expected behaviors for each value of annotation_type:
 |              "sites" or "group": 1-indexed location of sites of interest
 |              "spans": 1-indexed, inclusive locations of spans of interest
 |                  e.g. [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30
 |              "primers": Similar to spans, but 5'/3' direction is preserved.
 |                  e.g. [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20
 |          annotation_type (str):
 |              "group", "sites", "spans", or "primers".
 |          sequence (str | pandas.DataFrame):
 |              Nucleotide sequence, path to fasta file, or dataframe
 |              containing a "Sequence" column.
 |          name (str, optional): Name of annotation.
 |              Defaults to None.
 |          color (matplotlib color-like, optional): Color to be used for
 |              displaying this annotation on plots.
 |              Defaults to "blue".
 |  
 |  __iter__(self)
 |  
 |  __len__(self)
 |  
 |  from_sites(self, sites)
 |  
 |  from_spans(self, spans)
 |      Create the self.data dataframe from a list of spans.
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_sites(self)
 |      Returns a list of nucleotide positions and colors based on these
 |      sequence annotations.
 |      
 |      Returns:
 |          tuple: a list of nucleotide positions and a list of colors
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
 |      Required arguments:
 |          values (list of True or False)
 |              the boolean array
 |          sequence (string or rnav.data.Sequence)
 |              the sequence of the Annotation
 |          annotation_type ("spans", "sites", "primers", or "group")
 |              the type of the new annotation
 |              If "spans" or "primers", adjacent True values, or values within
 |              window are collapse to a region.
 |          name (string): a name for labelling the annotation.
 |      
 |      Optional arguments:
 |          color (string)
 |              a color for plotting the annotation
 |              Defaults to "blue"
 |          window (integer)
 |              a window around True values to include in the annotation
 |              Defaults to 1
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Method resolution order:
 |      Motif
 |      Annotation
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, sequence, name=None, color='blue')
 |      Creates a Motif annotation, which acts like a span Annotation, for
 |      highlighting a sequence motif of interest, given with conventional
 |      nucleotide codes. e.g. "DRACH"
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence to be searched.
 |              Defaults to None.
 |          motif (str, optional): sequence motif to be searched for.
 |              Defaults to None.
 |          name (str, optional): name of this annotation.
 |              Defaults to None.
 |          color (str, optional): color used to display these motif locations.
 |              Defaults to "blue".
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_spans_from_motif(self, sequence, motif)
 |      Returns a list of spans [[start, end], [start, end]] for each
 |      location of motif found within sequence, using conventional nucleotide
 |      codes.
 |      
 |      Args:
 |          sequence (str): sequence to be searched
 |          motif (str): sequence motif to be searched for.
 |      
 |      Returns:
 |          _type_: _description_
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
 |  
 |  from_spans(self, spans)
 |      Create the self.data dataframe from a list of spans.
 |  
 |  get_sites(self)
 |      Returns a list of nucleotide positions and colors based on these
 |      sequence annotations.
 |      
 |      Returns:
 |          tuple: a list of nucleotide positions and a list of colors
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
 |      Required arguments:
 |          values (list of True or False)
 |              the boolean array
 |          sequence (string or rnav.data.Sequence)
 |              the sequence of the Annotation
 |          annotation_type ("spans", "sites", "primers", or "group")
 |              the type of the new annotation
 |              If "spans" or "primers", adjacent True values, or values within
 |              window are collapse to a region.
 |          name (string): a name for labelling the annotation.
 |      
 |      Optional arguments:
 |          color (string)
 |              a color for plotting the annotation
 |              Defaults to "blue"
 |          window (integer)
 |              a window around True values to include in the annotation
 |              Defaults to 1
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
 |  Method resolution order:
 |      ORFs
 |      Annotation
 |      rnavigate.data.data.Sequence
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, input_data, name=None, sequence=None, color='blue')
 |      Create an annotation from a list of sites or regions.
 |      
 |      Args:
 |          input_data (list):
 |              List will be treated according to annotation_type argument.
 |              Expected behaviors for each value of annotation_type:
 |              "sites" or "group": 1-indexed location of sites of interest
 |              "spans": 1-indexed, inclusive locations of spans of interest
 |                  e.g. [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30
 |              "primers": Similar to spans, but 5'/3' direction is preserved.
 |                  e.g. [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20
 |          annotation_type (str):
 |              "group", "sites", "spans", or "primers".
 |          sequence (str | pandas.DataFrame):
 |              Nucleotide sequence, path to fasta file, or dataframe
 |              containing a "Sequence" column.
 |          name (str, optional): Name of annotation.
 |              Defaults to None.
 |          color (matplotlib color-like, optional): Color to be used for
 |              displaying this annotation on plots.
 |              Defaults to "blue".
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_spans_from_orf(self, sequence, which='all')
 |      Given a sequence string, returns spans for specified ORFs
 |      
 |      Args:
 |          sequence (str): RNA nucleotide sequence
 |          which (str): "all" returns all spans, "longest" returns the longest
 |              defaults to "all"
 |      
 |      Returns:
 |          list of tuples: (start, end) position of each ORF
 |              1-indexed, inclusive
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
 |  
 |  from_spans(self, spans)
 |      Create the self.data dataframe from a list of spans.
 |  
 |  get_sites(self)
 |      Returns a list of nucleotide positions and colors based on these
 |      sequence annotations.
 |      
 |      Returns:
 |          tuple: a list of nucleotide positions and a list of colors
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
 |      Required arguments:
 |          values (list of True or False)
 |              the boolean array
 |          sequence (string or rnav.data.Sequence)
 |              the sequence of the Annotation
 |          annotation_type ("spans", "sites", "primers", or "group")
 |              the type of the new annotation
 |              If "spans" or "primers", adjacent True values, or values within
 |              window are collapse to a region.
 |          name (string): a name for labelling the annotation.
 |      
 |      Optional arguments:
 |          color (string)
 |              a color for plotting the annotation
 |              Defaults to "blue"
 |          window (integer)
 |              a window around True values to include in the annotation
 |              Defaults to 1
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from rnavigate.data.data.Sequence:
 |  
 |  __str__(self)
 |      Return str(self).
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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
```

## rnavigate.plots

### rnavigate.plots.Alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Alignment in module rnavigate.plots.alignment

class Alignment(rnavigate.plots.plots.Plot)
 |  Alignment(num_samples, rows=None, cols=1, **kwargs)
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
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_data(self, alignment, label, ax=None)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=0.03, width_ax_rel=0.03, width_ax_in=None, height_ax_in=None, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      AP
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, nt_length, region='all', track_labels=True, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_data(self, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, annotations=None, domains=None, label='', ax=None, seqbar=True, title=True, panels=None, annotation_mode='track', track_height=None, profile_scale_factor=1, plot_error=False, nt_ticks=(20, 5))
 |  
 |  set_axis(self, ax, sequence, track_height=0, nt_ticks=(20, 5), max_height=300, yticks=None, ylabels=None)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=0.03, width_ax_rel=0.03, width_ax_in=None, height_ax_in=None, height_gap_in=0.5, width_gap_in=0.5, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      Circle
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_data(self, sequence, structure=None, structure2=None, interactions=None, interactions2=None, profile=None, annotations=None, label=None, colors=None, gap=30, nt_ticks=(20, 5))
 |  
 |  set_axis(self, ax, label, seq_circle, gap, nt_ticks)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=0.035, width_ax_rel=0.035, width_ax_in=None, height_ax_in=None, height_gap_in=1, width_gap_in=1, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      DistHist
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **plot_kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_data(self, structure, interactions, bg_interactions, label, atom="O2'", ax=None)
 |  
 |  plot_experimental_distances(self, ax, structure, interactions, atom, histtype='bar')
 |  
 |  plot_structure_distances(self, ax, structure, atom)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=1, width_gap_in=0.4, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      Heatmap
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, structure, **plot_kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_contour_distances(self, ax, levels, atom)
 |  
 |  plot_contour_regions(self, ax, interactions, regions)
 |  
 |  plot_data(self, interactions, label, levels=None, regions=None, interpolation=None, atom="O2'", plot_type='heatmap', weights=None)
 |  
 |  plot_heatmap_data(self, ax, interactions, interpolation)
 |  
 |  plot_kde_data(self, ax, interactions, weights=None, **kwargs)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      LinReg
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, scale='linear', regression='pearson', kde=False, region='all')
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  finalize(self)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_data(self, structure, profile, annotations, label, column=None, colors='sequence')
 |  
 |  plot_regression(self, i, j)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=0.3, width_gap_in=0.3, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      NucleotideDistribution
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **plot_kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_data(self, profile, label, column=None, normalize=None, ax=None)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=0.2, width_gap_in=0.4, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  add_lines(self, i, j, color, viewer, atom)
 |  
 |  get_orientation(self)
 |  
 |  get_viewer(self, i=None)
 |  
 |  hide_cylinders(self)
 |  
 |  plot_data(self, interactions=None, profile=None, label=None, colors='grey', atom="O2'", title=True, get_orientation=False, viewer=None)
 |  
 |  plot_interactions(self, viewer, interactions, atom)
 |  
 |  save(self)
 |      output png image of viewer, which must already be instantiated
 |  
 |  set_colors(self, viewer, profile, colors)
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=None, height_ax_in=None, height_gap_in=None, width_gap_in=None, top_in=None, bottom_in=None, left_in=None, right_in=None)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  Method resolution order:
 |      Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, rows=None, cols=None, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  add_colorbar_args(self, cmap)
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  plot_data(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
 |              A file path to write to. The file format is provided by this
 |              file extension (svg, pdf, or png).
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=None, height_ax_in=None, height_gap_in=None, width_gap_in=None, top_in=None, bottom_in=None, left_in=None, right_in=None)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  Method resolution order:
 |      ColorBar
 |      Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  plot_data(self, colorbar)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=3, height_ax_in=0.1, height_gap_in=0.75, width_gap_in=0.5, top_in=None, bottom_in=None, left_in=None, right_in=None)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  add_colorbar_args(self, cmap)
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      Profile
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, nt_length, region='all', **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_data(self, profile, annotations, domains, label, plot_error=True, column=None, seqbar=True, annotations_mode='track', nt_ticks=(20, 5))
 |  
 |  set_axis(self, ax, sequence, nt_ticks=(20, 5))
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=0.03, width_ax_in=None, height_ax_in=2, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
 |  
 |  set_labels(self, ax, axis_title='Reactivity Profile', xlabel='Nucleotide Position', ylabel='Reactivity')
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
 |  
 |  get_ax(self, i=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      QC
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  make_boxplot(self, labels)
 |  
 |  plot_data(self, profile, label)
 |  
 |  plot_mutations_per_molecule(self, profile, label, upper_limit=12)
 |  
 |  plot_read_lengths(self, profile, label, upper_limit=12)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=2, height_ax_in=2, height_gap_in=1, width_gap_in=1, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      ROC
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_data(self, structure, profile, label, nts='AUCG')
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=None, width_ax_in=1.5, height_ax_in=1.5, height_gap_in=0.3, width_gap_in=0.2, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      Skyline
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, nt_length, region='all', **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_data(self, profile, annotations=None, domains=None, label=None, columns=None, seqbar=True, errors=None, annotations_mode='track', nt_ticks=(20, 5))
 |  
 |  set_axis(self, ax, sequence, nt_ticks)
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=0.03, width_ax_in=None, height_ax_in=2, height_gap_in=1, width_gap_in=0.5, top_in=1, bottom_in=1, left_in=1, right_in=1)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
 |  
 |  set_labels(self, ax, axis_title='Raw Reactivity Profile', legend_title='Samples', xlabel='Nucleotide Position', ylabel='Profile')
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
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      SM
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, nt_length, region=None, panels=['profile', 'rates', 'depth'])
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  metric_abbreviate(self, num)
 |      takes a large number and applies an appropriate abbreviation
 |      
 |      Args:
 |          num (int): number to be abbreviated
 |      
 |      Returns:
 |          str: abbreviated number
 |  
 |  plot_data(self, profile, label)
 |      Creates a figure with the three classic Shapemapper plots.
 |  
 |  plot_sm_depth(self, ax, profile)
 |      Plots classic ShapeMapper read depth on the given ax
 |      
 |      Args:
 |          ax (pyplot ax): ax on which to add plot
 |  
 |  plot_sm_profile(self, ax, profile)
 |      Plots classic ShapeMapper normalized reactivity on the given ax
 |      
 |      Args:
 |          ax (pyplot ax): ax on which to add plot
 |  
 |  plot_sm_rates(self, ax, profile)
 |      Plots classic ShapeMapper mutation rates on the given ax
 |      
 |      Args:
 |          ax (pyplot ax): ax on which to add plot
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=None, width_ax_rel=0.03, width_ax_in=None, height_ax_in=2, height_gap_in=0.5, width_gap_in=0.5, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
 |  Method resolution order:
 |      SS
 |      rnavigate.plots.plots.Plot
 |      abc.ABC
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __init__(self, num_samples, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  plot_data(self, structure, interactions=None, interactions2=None, profile=None, annotations=None, label='', colors=None, nt_ticks=None, bp_style='dotted')
 |  
 |  set_figure_size(self, fig=None, ax=None, rows=None, cols=None, height_ax_rel=0.2, width_ax_rel=0.2, width_ax_in=None, height_ax_in=None, height_gap_in=0.5, width_gap_in=0.2, top_in=1, bottom_in=0.5, left_in=0.5, right_in=0.5)
 |      Sets figure size so that axes sizes are always consistent.
 |      
 |      Args:
 |          height_ax_rel (float, optional): ax unit to inches ratio for the
 |              y-ax.
 |          width_ax_rel (float, optional): ax unit to inches ration for the
 |              x-ax.
 |          width_ax_in (float, optional): fixed width of each ax in inches
 |          height_ax_in (float, optional): fixed height of each ax in inches
 |          width_gap_in (float, optional): fixed width of gaps between each
 |              ax in inches
 |          height_gap_in (float, optional): fixed height of gaps between each
 |              ax in inches
 |          top_in (float, optional): fixed height of top margin in inches
 |          bottom_in (float, optional): fixed height of bottom margin in inches
 |          left_in (float, optional): fixed width of left margin in inches
 |          right_in (float, optional): fixed width of right margin in inches
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
 |  
 |  get_ax(self, i=None)
 |  
 |  get_rows_columns(self, rows=None, cols=None)
 |  
 |  plot_colorbars(self)
 |  
 |  save(self, filename)
 |      Saves the figure to a file
 |      
 |      Args:
 |          filename (str):
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
```

### rnavigate.plots.adjust_spines

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function adjust_spines in module rnavigate.plots.functions.functions

adjust_spines(ax, spines)
```

### rnavigate.plots.clip_spines

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function clip_spines in module rnavigate.plots.functions.functions

clip_spines(ax, spines)
```

### rnavigate.plots.get_nt_ticks

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_nt_ticks in module rnavigate.plots.functions.functions

get_nt_ticks(sequence, region, gap)
```

### rnavigate.plots.set_nt_ticks

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function set_nt_ticks in module rnavigate.plots.functions.functions

set_nt_ticks(ax, sequence, region, major, minor)
```

### rnavigate.plots.box_xtick_labels

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function box_xtick_labels in module rnavigate.plots.functions.functions

box_xtick_labels(ax)
```

### rnavigate.plots.plot_interactions_arcs

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_interactions_arcs in module rnavigate.plots.functions.functions

plot_interactions_arcs(ax, interactions, panel, yvalue=0, region='all')
```

### rnavigate.plots.plot_profile_bars

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_profile_bars in module rnavigate.plots.functions.functions

plot_profile_bars(ax, profile, scale_factor=1, plot_error=True, bottom=0, region='all')
```

### rnavigate.plots.plot_profile_skyline

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_profile_skyline in module rnavigate.plots.functions.functions

plot_profile_skyline(ax, profile, label, columns, errors)
```

### rnavigate.plots.plot_sequence_alignment

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_sequence_alignment in module rnavigate.plots.functions.functions

plot_sequence_alignment(ax, alignment, labels, top=5, bottom=-5, ytrans='data')
```

### rnavigate.plots.plot_annotation_track

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_annotation_track in module rnavigate.plots.functions.tracks

plot_annotation_track(ax, annotation, yvalue, height, mode, region='all', ytrans='data')
```

### rnavigate.plots.plot_domain_track

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_domain_track in module rnavigate.plots.functions.tracks

plot_domain_track(ax, spans, yvalue, height, region='all', ytrans='data')
```

### rnavigate.plots.plot_sequence_track

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_sequence_track in module rnavigate.plots.functions.tracks

plot_sequence_track(ax, sequence, yvalue=-0.05, height=0.05, ytrans='data', verticalalignment='bottom', region='all')
    # 1-dimensional x-axis tracks
```

### rnavigate.plots.plot_annotation_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_annotation_ss in module rnavigate.plots.functions.ss

plot_annotation_ss(ax, structure, annotation)
```

### rnavigate.plots.plot_basepairs_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_basepairs_ss in module rnavigate.plots.functions.ss

plot_basepairs_ss(ax, structure, bp_style)
```

### rnavigate.plots.plot_interactions_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_interactions_ss in module rnavigate.plots.functions.ss

plot_interactions_ss(ax, structure, interactions)
```

### rnavigate.plots.plot_nucleotides_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_nucleotides_ss in module rnavigate.plots.functions.ss

plot_nucleotides_ss(ax, structure, colors)
```

### rnavigate.plots.plot_positions_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_positions_ss in module rnavigate.plots.functions.ss

plot_positions_ss(ax, structure, xticks=20)
```

### rnavigate.plots.plot_sequence_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_sequence_ss in module rnavigate.plots.functions.ss

plot_sequence_ss(ax, structure, colors)
```

### rnavigate.plots.plot_structure_ss

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_structure_ss in module rnavigate.plots.functions.ss

plot_structure_ss(ax, structure, colors)
    # Secondary structure diagram ploting functions
```

### rnavigate.plots.plot_annotation_circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_annotation_circle in module rnavigate.plots.functions.circle

plot_annotation_circle(ax, seq_circle, annotation, offset=1)
```

### rnavigate.plots.plot_interactions_circle

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function plot_interactions_circle in module rnavigate.plots.functions.circle

plot_interactions_circle(ax, seq_circle, interactions)
```

## rnavigate.styles

### rnavigate.styles.update_copy

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function update_copy in module rnavigate.styles

update_copy(original_settings, user_settings)
```

### rnavigate.styles.Settings

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: class Settings in module rnavigate.styles

class Settings(builtins.dict)
 |  Settings(user_settings)
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
```

### rnavigate.styles.get_nt_color

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_nt_color in module rnavigate.styles

get_nt_color(nt, colors=None)
```

### rnavigate.styles.get_nt_cmap

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function get_nt_cmap in module rnavigate.styles

get_nt_cmap()
```

### rnavigate.styles.apply_style

[back to top](#rnavigate-full-api)

```text
Python Library Documentation: function apply_style in module rnavigate.styles

apply_style(style_dict)
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
 |      Constructs a Data object given a sequence string, fasta file, or
 |      dataframe containing a "Sequence" column.
 |      
 |      Args:
 |          sequence (str | pandas.DataFrame):
 |              sequence string, fasta file, or a pandas dataframe containing
 |              a "Sequence" column.
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
 |      Return str(self).
 |  
 |  get_aligned_data(self, alignment)
 |  
 |  get_colors(self, source, pos_cmap='rainbow', profile=None, structure=None, annotations=None)
 |      Get a numpy array of colors that fits the current sequence.
 |      
 |      Args:
 |          source (str | array of color-like): One of the following:
 |              "position": colors represent position in sequence
 |              "sequence": colors represent nucleotide identity
 |              "annotations": colors represent sequence annotations
 |              "profile": colors represent per-nucleotide data
 |              "structure": colors represent base-pairing status
 |              matplotlib color-like: all colors are this color
 |              array of color like: must match length of sequence
 |          pos_cmap (str, optional): cmap used if source="position".
 |              Defaults to 'rainbow'.
 |          profile (Profile or subclass, optional): Data object containing
 |              per-nucleotide information. Defaults to None.
 |          structure (SecondaryStructure or subclass, optional): Data object
 |              containing secondary structure information.
 |              Defaults to None.
 |          annotations (list of Annotations or subclass, optional): list of
 |              Data objects containing annotations. Defaults to None.
 |      
 |      Returns:
 |          numpy array: one matplotlib color-like value for each nucleotide in
 |              self.sequence
 |  
 |  get_colors_from_annotations(self, annotations)
 |  
 |  get_colors_from_positions(self, pos_cmap='rainbow')
 |  
 |  get_colors_from_profile(self, profile)
 |  
 |  get_colors_from_sequence(self)
 |  
 |  get_colors_from_structure(self, structure)
 |  
 |  get_seq_from_dataframe(self, dataframe)
 |      Parse a dataframe for the sequence string, store as self.sequence.
 |      
 |      Args:
 |          dataframe (pandas DataFrame): must contain a "Sequence" column
 |  
 |  normalize_sequence(self, t_or_u='U', uppercase=True)
 |      Converts sequence to all uppercase nucleotides and corrects T or U.
 |      
 |      Optional arguments:
 |          t_or_u ("T", "U", or False)
 |              "T" converts "U"s to "T"s
 |              "U" converts "T"s to "U"s
 |              False does nothing.
 |              Defaults to "U"
 |          uppercase (True or False)
 |              Whether to make sequence all uppercase
 |              Defaults to True
 |  
 |  read_fasta(self, fasta)
 |      Parse a fasta file for the first sequence. Store the sequence name
 |      as self.gene and the sequence string as self.sequence.
 |      
 |      Args:
 |          fasta (str): path to fasta file
 |  
 |  ----------------------------------------------------------------------
 |  Readonly properties inherited from rnavigate.data.data.Sequence:
 |  
 |  length
 |      Get the length of the sequence
 |      
 |      Returns:
 |          int: the length of self.sequence
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

