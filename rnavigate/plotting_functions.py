"""Contains all rnavigate convenience plotting functions"""

from rnavigate import plots
from rnavigate import data
from rnavigate.helper_functions import PlottingArgumentParser, _parse_plot_kwargs, fit_data
from rnavigate.data_loading import get_sequence


__all__ = [
    "plot_qc",
    "plot_shapemapper",
    "plot_skyline",
    "plot_profile",
    "plot_alignment",
    "plot_arcs",
    "plot_arcs_compare",
    "plot_ss",
    "plot_mol",
    "plot_heatmap",
    "plot_circle",
    "plot_linreg",
    "plot_roc",
    "plot_disthist",
    "plot_ntdist",
    ]

def plot_qc(
        # required
        samples,
        profile,
        # optional display
        labels=None
        ):
    """Creates a multipanel quality control plot displaying mutations per
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
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        profile=profile,
    )
    plot = plots.QC(num_samples=len(samples))
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict)
    plot.set_figure_size()
    return plot


def plot_shapemapper(
        # required
        sample,
        profile,
        # optional display
        label=None,
        panels=None):
    """Makes a standard ShapeMapper2 profile plot with 3 panels: Normalized
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
    """
    if panels is None:
        panels = ["profile", "rates", "depth"]
    if label is None:
        label = sample.sample
    profile = sample.get_data(profile)
    plot = plots.SM(profile.length, panels=panels)
    plot.plot_data(profile=profile, label=label)
    plot.set_figure_size()
    return plot


def plot_skyline(
        # required
        samples,
        profile,
        # optional data inputs
        sequence=None,
        annotations=None,
        domains=None,
        # optional data display
        labels=None,
        nt_ticks=(20, 5),
        columns=None,
        errors=None,
        annotations_mode='track',
        seqbar=True,
        region="all",
        # optional plot display
        plot_kwargs=None,
        ):
    """Plots multiple per-nucleotide datasets on a single axis.

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
    """
    sequence = get_sequence(sequence, samples[0], profile)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence.null_alignment,
        profile=profile,
        annotations=annotations,
        domains=domains,
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Skyline")
    plot = plots.Skyline(
        num_samples=len(samples), nt_length=sequence.length, region=region,
        **plot_kwargs
        )
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(
            **data_dict, columns=columns, errors=errors, nt_ticks=nt_ticks,
            annotations_mode=annotations_mode, seqbar=seqbar,
            )
    plot.set_figure_size()
    plot.plot_colorbars()
    return plot


def plot_profile(
        # required
        samples,
        profile,
        # optional data inputs
        sequence=None,
        annotations=None,
        domains=None,
        # optional data display
        labels=None,
        nt_ticks=(20, 5),
        column=None,
        plot_error=True,
        annotations_mode='track',
        seqbar=True,
        region="all",
        # optional plot display
        colorbars=True,
        plot_kwargs=None,
        ):
    """Aligns reactivity profiles by sequence and plots them on seperate axes.

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
    """
    sequence = get_sequence(sequence, samples[0], profile)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence.null_alignment,
        annotations=annotations,
        domains=domains,
        profile=profile)
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Profile")
    plot = plots.Profile(
        num_samples=len(samples), nt_length=sequence.length, region=region,
        **plot_kwargs
        )
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(
            **data_dict, column=column, plot_error=plot_error, seqbar=seqbar,
            annotations_mode=annotations_mode, nt_ticks=nt_ticks,
            )
    plot.set_figure_size()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_alignment(
        # required
        data1,
        data2,
        # optional display
        labels=None,
        plot_kwargs=None
        ):
    """Plots the sequence alignment used to compare two sequences

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
    """
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Alignment")
    if labels is None:
        labels = [f"{s.sample}: {seq}" for s, seq in [data1, data2]]
    plot = plots.Alignment(num_samples=1, **plot_kwargs)
    alignment = data.SequenceAlignment(
        data1[0].data[data1[1]], data2[0].data[data2[1]])
    plot.plot_data(alignment=alignment, label=labels)
    plot.set_figure_size()
    return plot


def plot_arcs(
        # required
        samples,
        sequence,
        # optional data inputs
        structure=None,
        structure2=None,
        interactions=None,
        interactions2=None,
        profile=None,
        annotations=None,
        domains=None,
        # optional data display
        labels=None,
        nt_ticks=(20, 5),
        profile_scale_factor=1,
        plot_error=False,
        annotation_mode='track',
        panels=None,
        seqbar=True,
        region="all",
        # optional plot display
        colorbars=True,
        title=True,
        plot_kwargs=None,
        ):
    """Plots interactions and/or base-pairs as arcs.

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
    """
    sequence = get_sequence(sequence, samples[0], structure)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence.null_alignment,
        annotations=annotations,
        domains=domains,
        structure=structure,
        structure2=structure2,
        profile=profile,
        interactions2=interactions2,
        interactions=interactions
        )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.AP")
    parsed_args.update_rows_cols(plot_kwargs)
    # initialize plot
    plot = plots.AP( # TODO: plots should take alignments instead of region + length
        num_samples=parsed_args.num_samples,
        nt_length=sequence.length,
        region=region, **plot_kwargs)
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(
            **data_dict, sequence=sequence, title=title, panels=panels,
            plot_error=plot_error, annotation_mode=annotation_mode,
            seqbar=seqbar, profile_scale_factor=profile_scale_factor,
            nt_ticks=nt_ticks,
            )
    plot.set_figure_size()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_arcs_compare(
        # required
        samples,
        sequence,
        # optional data inputs
        structure=None,
        structure2=None,
        interactions=None,
        interactions2=None,
        profile=None,
        # optional data display
        labels=None,
        profile_scale_factor=1,
        plot_error=False,
        region="all",
        # optional plot display
        colorbars=True,
        plot_kwargs=None,
        ):
    """Generates a single arc plot displaying combinations of secondary
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
    """
    if len(samples) != 2:
        raise ValueError("Only 2 samples can be compared.")
    seq1 = get_sequence(sequence, samples[0], structure)
    seq2 = get_sequence(sequence, samples[1], structure)
    alignment = data.SequenceAlignment(seq1, seq2, full=True)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        structure=structure,
        structure2=structure2,
        profile=profile,
        interactions2=interactions2,
        interactions=interactions,
        )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.AP")
    # initialize plot
    plot = plots.AP(
        num_samples=1, nt_length=len(alignment.target_sequence), region=region,
        **plot_kwargs
        )
    # loop through samples and filters, adding each as a new axis
    labels = []
    for data_dict, seq, panel in zip(parsed_args.data_dicts, [1, 2], ["top", "bottom"]):
        if seq == 1:
            seq, other_seq = seq1, seq2
        else:
            seq, other_seq = seq2, seq1
        alignment = data.SequenceAlignment(seq, other_seq, full=True)
        panels = {
            'structure': panel, 
            'interactions': panel,
            'interactions2': panel,
            'profile': panel}
        labels.append(data_dict.pop('label'))
        data_dict = fit_data(data_dict, alignment)
        plot.plot_data(
            ax=0, sequence=seq, track_height=6, label='', seqbar=False,
            panels=panels, **data_dict, plot_error=plot_error,
            profile_scale_factor=profile_scale_factor,
            )
    plots.plot_sequence_alignment(
        ax=plot.axes[0, 0], alignment=alignment, labels=labels, top=6,
        bottom=0,
        )
    plot.set_figure_size()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_ss(
        # required
        samples,
        structure,
        # optional data inputs
        profile=None,
        annotations=None,
        interactions=None,
        interactions2=None,
        # optional data display
        labels=None,
        colors=None,
        nt_ticks=None,
        bp_style='dotted',
        # optional plot display
        colorbars=True,
        plot_kwargs=None,
        ):
    """Generates a multipanel secondary structure drawing with optional
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
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        structure=structure,
        annotations=annotations,
        profile=profile,
        interactions=interactions,
        interactions2=interactions2,
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.SS")
    parsed_args.update_rows_cols(plot_kwargs)
    # initialize plot using all structure drawings
    plot = plots.SS(num_samples=parsed_args.num_samples, **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        data_dict = fit_data(data_dict, data_dict['structure'].null_alignment)
        plot.plot_data(
            **data_dict, colors=colors, nt_ticks=nt_ticks, bp_style=bp_style
            )
    plot.set_figure_size()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_mol(
        # required
        samples,
        structure,
        # optional data inputs
        profile=None,
        interactions=None,
        # optional data display
        labels=None,
        style='cartoon',
        hide_cylinders=False,
        colors='grey',
        atom="O2'",
        rotation=None,
        orientation=None,
        get_orientation=False,
        # optional viewer display
        title=True,
        colorbars=True,
        width=400,
        height=400,
        rows=None,
        cols=None,
        background_alpha=1,
        show=True,
        ):
    """Generates a multipanel interactive 3D molecular rendering of a PDB
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
    """
    sequence = get_sequence(structure, samples[0])
    if get_orientation:
        plot = plots.Mol(num_samples=1, pdb=sequence, width=800, height=800)
        plot.plot_data(title=False, get_orientation=True)
        plot.view.show()
        return plot
    parsed_args = PlottingArgumentParser(
        samples=samples,
        alignment=sequence.null_alignment,
        labels=labels,
        profile=profile,
        interactions=interactions,
    )
    # initialize plot using 1st 3D structure (applies to all samples)
    structure = samples[0].get_data(structure)
    plot = plots.Mol(
        num_samples=parsed_args.num_samples, pdb=structure, width=width,
        height=height, background_alpha=background_alpha, rotation=rotation,
        orientation=orientation, style=style, rows=rows, cols=cols
        )
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(
            **data_dict, colors=colors, atom=atom, title=title
            )
    # hide nucleotide cylinders in all viewers
    if hide_cylinders:
        plot.hide_cylinders()
    # show viewer grid
    if show:
        plot.view.show()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_heatmap(
        # required
        samples,
        sequence,
        # optional data input
        structure=None,
        interactions=None,
        regions=None,
        # optional data display
        labels=None,
        levels=None,
        interpolation='nearest',
        atom="O2'",
        plot_type="heatmap",
        weights=None,
        # optional plot display
        rows=None,
        cols=None,
        plot_kwargs=None,
        ):
    """Generates a multipanel plot displaying a heatmap of inter-nucleotide
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
    """
    # TODO: this ones a little tricky to get working with arg parse function
    sequence = get_sequence(sequence, samples[0], structure)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence.null_alignment,
        interactions=interactions
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Heatmap")
    # initialize plot using 1st 3D pdb (applies to all samples)
    structure = samples[0].data[structure]
    plot = plots.Heatmap(
        parsed_args.num_samples, structure, rows=rows, cols=cols, **plot_kwargs
        )
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(
            **data_dict, regions=regions, levels=levels, weights=weights,
            interpolation=interpolation, atom=atom, plot_type=plot_type,
            )
    plot.set_figure_size()
    return plot

# TODO: refactor to match profile, arcs, skyline, ss, then update docstring
def plot_circle(
        # required
        samples,
        sequence,
        # optional data inputs
        structure=None,
        structure2=None,
        interactions=None,
        interactions2=None,
        annotations=None,
        profile=None,
        # optional data display
        colors=None,
        nt_ticks=(20, 5),
        gap=30,
        labels=None,
        # optional plot display
        colorbars=True,
        plot_kwargs=None):
    """Creates a figure containing a circle plot for each sample given.

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
    """

    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        sequence=sequence,
        structure=structure,
        structure2=structure2,
        interactions2=interactions2,
        profile=profile,
        annotations=annotations,
        interactions=interactions
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Circle")
    parsed_args.update_rows_cols(plot_kwargs)
    plot = plots.Circle(num_samples=parsed_args.num_samples, **plot_kwargs)
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        data_dict = fit_data(data_dict, data_dict['sequence'].null_alignment)
        plot.plot_data(
            **data_dict, colors=colors, nt_ticks=nt_ticks,
            gap=gap
            )
    plot.set_figure_size()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_linreg(
        # required
        samples,
        profile,
        # optional data inputs
        sequence=None,
        structure=None,
        annotations=None,
        # optional data display
        labels=None,
        kde=False,
        scale='linear',
        regression='pearson',
        colors='sequence',
        column=None,
        region='all',
        # optional plot display
        colorbars=True,
        plot_kwargs=None,
        ):
    """Performs linear regression analysis and generates scatter plots of all
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
    """
    sequence = get_sequence(sequence, samples[0], profile)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence.null_alignment,
        structure=structure,
        profile=profile,
        annotations=annotations
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.LinReg")
    plot = plots.LinReg(
        num_samples=parsed_args.num_samples, scale=scale, regression=regression,
        kde=kde, region=region, **plot_kwargs
        )
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, colors=colors, column=column)
    plot.set_figure_size()
    if colorbars:
        plot.plot_colorbars()
    return plot


def plot_roc(
        # required
        samples,
        structure,
        profile,
        # optional data display
        labels=None,
        nts='AUCG',
        # optional plot display
        plot_kwargs=None
        ):
    """Performs receiver operator characteristic analysis (ROC), calculates
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
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        structure=structure,
        profile=profile,
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.ROC")
    plot = plots.ROC(parsed_args.num_samples, **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        data_dict = fit_data(data_dict, data_dict['structure'].null_alignment)
        plot.plot_data(**data_dict, nts=nts)
    plot.set_figure_size()
    return plot


def plot_disthist(
        # required
        samples,
        structure,
        interactions,
        # optional data inputs
        bg_interactions=None,
        # optional data display
        labels=None,
        same_axis=False,
        atom="O2'",
        # optional plot display
        rows=None,
        cols=None,
        plot_kwargs=None
        ):
    """Calculates 3D distance of nucleotides in inter-nucleotide data and plots
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
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        structure=structure,
        bg_interactions=bg_interactions,
        interactions=interactions,
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.DistHist")
    # initialize plot
    if same_axis:
        plot = plots.DistHist(
            parsed_args.num_samples, rows=1, cols=1, **plot_kwargs
            )
        ax = plot.axes[0, 0]
    else:
        plot = plots.DistHist(
            parsed_args.num_samples, rows=rows, cols=cols, **plot_kwargs
            )
        ax = None
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, ax=ax, atom=atom)
    plot.set_figure_size()
    return plot


def plot_ntdist(
        # required
        samples,
        profile,
        # optional data display
        labels=None,
        column=None,
        # optional plot display
        plot_kwargs=None
        ):
    """Plots the distributions of values at A, U, C, and G.

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
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        profile=profile,
    )
    plot_kwargs = _parse_plot_kwargs(
        plot_kwargs, "rnavigate.plots.NucleotideDistribution"
        )
    plot = plots.NucleotideDistribution(parsed_args.num_samples, **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, column=column)
    plot.set_figure_size()
    return plot
