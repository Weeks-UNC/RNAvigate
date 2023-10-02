"""Contains all rnavigate convenience plotting functions"""

from rnavigate import plots
from rnavigate import data
from rnavigate.helper_functions import PlottingArgumentParser, _parse_plot_kwargs, fit_data
from rnavigate.data_loading import get_sequence


def plot_qc(
        samples,
        labels=None,
        plot_kwargs=None,
        **kwargs):
    """Creates a multipanel quality control plot displaying mutations per
    molecule, read length distribution, and mutation rate distributions for
    modified and unmodified samples.

    Args:
        samples (list of rnavigate.Sample): samples holding ShapeMapper logs.
        labels (list of str, optional): same length as samples list. labels to
            to be used on plot legends.
            Defaults to sample.sample for each sample.
        plot_kwargs (dict, optional): passed to rnavigate.plots.QC().
            Defaults to {}.
        **kwargs: passed to QC.plot_data

    Returns:
        rnavigate.plots.QC: the quality control plot object
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        log="log",
        profile="default_profile"
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.QC")
    plot = plots.QC(num_samples=len(samples), **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    return plot


def plot_shapemapper(
        sample,
        label=None,
        profile="default_profile",
        panels=None):
    """Makes a standard ShapeMapper2 profile plot with 3 panels: Normalized
    Reactivities, mutation rates, and read depths.

    Args:
        profile (str, optional):
            a data keyword or rnavigate.data.ShapeMaP object
        panels (list, optional):
            Which of the three panels to include.
            Defaults to ["profile", "rates", "depth"].

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
        samples,
        sequence=None,
        profile="default_profile",
        labels=None,
        annotations=None,
        region="all",
        plot_kwargs=None,
        **kwargs):
    """Plots multiple per-nucleotide datasets on a single axis.

    Args:
        samples (list of rnavigate.Sample):
            samples to plot.
        sequence (str | rnavigate.data.Sequence, optional):
            a data keyword, a sequence string, or a Sequence object.
            All data are mapped to this sequence before plotting.
            Defaults to the value of the profile argument below.
        profile (str, optional):
            a data keyword to retreive a Profile object.
            Defaults to "default_profile".
        labels (list of str, optional):
            list containing Labels to be used in plot legends for each sample.
            Defaults to sample.sample for each sample
        annotations (list of str, optional):
            a list of data keywords to retreive a list of Annotations objects.
            Defaults to [].
        region (list of int: length 2, optional):
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, sequence length].
        plot_kwargs (dict, optional): passed to rnavigate.plots.Skyline().
            Defaults to {}.
        **kwargs: passed to Skyline.plot_data.
            see rnavigate.plots.Skyline.plot_data for more detail.

    Returns:
        rnavigate.plots.Skyline: the skyline plot object
    """
    sequence = get_sequence(sequence, samples[0], profile)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence._alignment,
        profile=profile,
        annotations=annotations
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Skyline")
    plot = plots.Skyline(num_samples=len(samples),
                         nt_length=sequence.length,
                         region=region,
                         **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    return plot


def plot_profile(
        samples,
        sequence=None,
        profile="default_profile",
        labels=None,
        annotations=None,
        region="all",
        plot_kwargs=None,
        **kwargs):
    """Aligns reactivity profiles by sequence and plots them on seperate axes.

    Args:
        samples (list of rnavigate.Sample):
            samples to plot.
        sequence (str | rnavigate.data.Sequence, optional):
            a data keyword, a sequence string, or a Sequence object.
            All data are mapped to this sequence before plotting.
            Defaults to the value of the profile argument below.
        profile (str | rnavigate.data.Profile, optional):
            a data keyword or a Profile object.
            Defaults to "default_profile".
        labels (list of str, optional):
            list containing Labels to be used in plot legends for each sample.
            Defaults to sample.sample for each sample
        annotations (list of str, optional):
            a list of data keywords to retreive a list of Annotations objects.
            Defaults to [].
        region (list of int: length 2, optional):
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, sequence length].
        plot_kwargs (dict, optional): passed to rnavigate.plots.Profile().
            Defaults to {}.
        **kwargs: passed to Profile.plot_data.
            see rnavigate.plots.skyline.Profile.plot_data for more detail.

    Returns:
        rnavigate.plots.Profile: the Profile plot object
    """
    sequence = get_sequence(sequence, samples[0], profile)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence._alignment,
        annotations=annotations,
        profile=profile)
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Profile")
    plot = plots.Profile(num_samples=len(samples),
                                 nt_length=sequence.length,
                                 region=region,
                                 **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    return plot


def plot_alignment(
        data1,
        data2,
        labels=None,
        plot_kwargs=None,
        **kwargs):
    """Plots the sequence alignment used to compare two data objects.

    Args:
        data1 (tuple (rnavigate.Sample, str)): a sample and data keyword
        data2 (tuple (rnavigate.Sample, str)): a sample and data keyword
        labels (list of str: length 2, optional): Labels used for each sample.
            Defaults to sample.sample + data keyword for each sample
        plot_kwargs (dict, optional): passed to rnavigate.plots.Alignment().
            Defaults to {}.
        **kwargs: passed to Alignment.plot_data.
            help(rnavigate.plots.Alignment.plot_data) for more detail.

    Returns:
        rnavigate.plots.Alignment: the Alignment plot object
    """
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Alignment")
    if labels is None:
        labels = [f"{s.sample}: {seq}" for s, seq in [data1, data2]]
    plot = plots.Alignment(num_samples=1, **plot_kwargs)
    alignment = data.SequenceAlignment(
        data1[0].data[data1[1]], data2[0].data[data2[1]])
    plot.plot_data(alignment=alignment, label=labels, **kwargs)
    plot.set_figure_size()
    return plot


def plot_arcs(
        samples,
        sequence=None,
        structure="default_structure",
        structure2=None,
        interactions=None,
        interactions2=None,
        profile="default_profile",
        annotations=None,
        labels=None,
        region="all",
        plot_kwargs=None,
        colorbar=True,
        **kwargs):
    """Plots interactions and/or base-pairs as arcs.

    Args:
        samples (list of rnavigate.Sample):
            samples to plot.
        sequence (str | rnavigate.data.Sequence, optional):
            a data keyword, a sequence string, or a Sequence object.
            All data are mapped to this sequence before plotting.
            Defaults to the value of the profile argument below.
        structure (str| rnavigate.data.SecondaryStructure, optional):
            a data keyword or a SecondaryStructure object
            Defaults to "default_structure".
        structure2 (str | rnavigate.data.SecondaryStructure, optional):
            same as structure. base-pairs from structure and structure2 are
            overlayed and color-coded (shared, structure only, structure2 only)
            Defaults to None
        interactions (str | rnavigate.data.Interactions | dict | list, optional):
            parsed with parse_interactions
            defaults to [].
        interactions2 (str, optional):
            parsed with parse_interactions
            Defaults to {}.
        profile (str | rnavigate.data.Profile, optional):
            a data keyword or a Profile object.
            Defaults to "default_profile".
        labels (list of str, optional):
            list containing Labels to be used in plot legends for each sample.
            Defaults to sample.sample for each sample
        annotations (list of str, optional):
            a list of data keywords to retreive a list of Annotations objects.
            Defaults to [].
        region (list of int: length 2, optional):
            start and end positions to plot. 1-indexed, inclusive.
            Defaults to [1, sequence length].
        plot_kwargs (dict, optional): passed to rnavigate.plots.AP().
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to AP.plot_data.
            see rnavigate.plots.AP.plot_data for more detail.

    Returns:
        rnavigate.plots.AP: the ArcPlot object
    """
    sequence = get_sequence(sequence, samples[0], structure)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence._alignment,
        annotations=annotations,
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
        plot.plot_data(**data_dict, seq=sequence, **kwargs)
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_arcs_compare(
        samples,
        sequence=None,
        structure="default_structure",
        structure2=None,
        interactions=None,
        interactions2=None,
        profile=None,
        region="all",
        labels=None,
        plot_kwargs=None,
        colorbar=True,
        **kwargs):
    """Generates a single arc plot displaying combinations of secondary
    structures, per-nucleotide data, inter-nucleotide data, and sequence
    annotations. The first sample will be on top, the second on the bottom.
    Center shows how these sequences are being aligned.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
        sequence (str | rnavigate.Sequence, optional):
            a data keyword from sample.data, a sequence string, or an
            rnavigate.Sequence object. All data will be aligned and mapped to
            positions in this sequence.
            Defaults to the value of the structure argument below.
        structure (str | rnavigate.SecondaryStructure, optional):
            an rnavigate.SecondaryStructure object or a data keyword pointing
            to an rnavigate.SecondaryStructure.
            Defaults to "default_structure".
        structure2 (str | rnavigate.SecondaryStructure, optional):
            Same as structure.
            basepairs from structure and structure2 will be
            plotted on the top half of each panel. Basepairs are colored by
            which structure contains them (shared, structure only, structure2 only).
            Defaults to None.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to sequence coordinates,
            filtered using interactions_filter arguments, and displayed on the
            bottom half of each panel.
            Defaults to None.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are displayed in center of each panel.
            Defaults to "default_profile".
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        region (list of int: length 2, optional): start and end position of
            sequence to be plotted. 1-indexed, inclusive.
            Defaults to [0, sequence length].
        plot_kwargs (dict, optional): kwargs passed to AP(). See
            rnavigate.plots.AP for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to AP.plot_data.
            see rnavigate.plots.AP.plot_data for more detail.

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
    plot = plots.AP(num_samples=1, nt_length=len(alignment.target_sequence),
                    region=region, **plot_kwargs)
    # loop through samples and filters, adding each as a new axis
    labels = []
    for data_dict, seq, panel in zip(parsed_args.data_dicts, [1, 2], ["top", "bottom"]):
        if seq == 1:
            seq, other_seq = seq1, seq2
        else:
            seq, other_seq = seq2, seq1
        alignment = data.SequenceAlignment(seq, other_seq, full=True)
        panels = {
            'ct_panel': panel, 
            'interactions_panel': panel,
            'interactions2_panel': panel,
            'profile_panel': panel}
        labels.append(data_dict.pop('label'))
        data_dict = fit_data(data_dict, alignment)
        plot.plot_data(ax=0, seq=None, annotation_gap=10, label='',
                       seqbar=False, **panels, **data_dict, **kwargs)
    plots.alignment.plot_alignment(
        plot=plot, ax=plot.axes[0, 0], alignment=alignment, label=labels,
        center=-5, offset=4, spines_positions={"top": 0, "bottom": -10})
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_ss(
        samples,
        structure="default_structure",
        profile="default_profile",
        annotations=None,
        interactions=None,
        interactions2=None,
        labels=None,
        plot_kwargs=None,
        colorbar=True,
        **kwargs):
    """Generates a multipanel secondary structure drawing with optional
    coloring by per-nucleotide data and display of inter-nucleotide data and/or
    sequence annotations. Each plot may display a unique sample and/or
    inter-nucleotide data filtering scheme.

    Args:
        samples (list of rnavigate.Sample):
            A list containing rnavigate.Samples to visualize.
        structure (str, optional):
            A data keyword pointing to a secondary structure to visualize.
            Defaults to "default_structure"
        profile (str, optional):
            A data keyword pointing to profile data.
            Defaults to "default_profile".
        annotations (list, optional):
            A list of data keywords pointing to annotations data to visualize.
            Defaults to [].
        interactions (str, optional):
            A data keyword pointing to internucleotide data.
            These data are first filtered using interactions_filter arguments.
            Defaults to None.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to SS(). See
            rnavigate.plots.SS for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to SS.plot_data.
            see rnavigate.plots.SS.plot_data for more detail.

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
        data_dict = fit_data(data_dict, data_dict['structure']._alignment)
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_mol(
        samples,
        structure="default_pdb",
        interactions=None,
        profile="default_profile",
        labels=None,
        show=True,
        hide_cylinders=False,
        colorbar=True,
        custom_function=None,
        plot_kwargs=None,
        **kwargs):
    """Generates a multipanel interactive 3D molecular rendering of a PDB
    structure. Nucleotides may be colored by per-nucleotide data or custom
    color lists. Inter-nucleotide data may be displayed as cylinders connecting
    atoms or residues. Each plot may display a unique sample and/or filtering
    scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        structure (str, optional): a key from sample.data to retrieve PDB data
            with atomic coordinates.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to structre sequence
            coordinates, filtered using interactions_filter arguments, and
            displayed as cylinders connecting nucleotides in the 3D structure.
            Defaults to None.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data may be used to color nucleotides.
            Defaults to "default_profile".
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends.
            Defaults to default sample name.
        show (bool, optional): whether to display the interactive rendering.
            Defaults to True
        hide_cylinders (bool, optional): whether to display cylinders
            representing nucleoside orientation. Setting to false will display
            only the backbone as a ribbon.
            Defaults to False.
        plot_kwargs (dict, optional): kwargs passed to Mol(). See
            rnavigate.plots.Mol for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Mol.plot_data.
            see rnavigate.plots.Mol.plot_data for more detail.

    Returns:
        rnavigate.plots.Mol plot: object containing py3dmol viewer with
            additional plotting and file saving methods
    """
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        profile=profile,
        interactions=interactions,
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.mol")
    parsed_args.update_rows_cols(plot_kwargs)
    # initialize plot using 1st 3D structure (applies to all samples)
    structure = samples[0].get_data(structure)
    plot = plots.Mol(num_samples=parsed_args.num_samples,
                     pdb=structure, **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    # apply custom function
    if custom_function is not None:
        custom_function(plot)
    # hide nucleotide cylinders in all viewers
    if hide_cylinders:
        plot.hide_cylinders()
    # show viewer grid
    if show:
        plot.view.show()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_heatmap(
        samples,
        structure=None,
        interactions=None,
        labels=None,
        plot_kwargs=None,
        **kwargs):
    """Generates a multipanel plot displaying a heatmap of inter-nucleotide
    data (nucleotide resolution of 2D KDE) and/or contour map of pdb
    distances. Each plot may display a unique sample and/or filtering scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        structure (str, optional): a key from sample.data to retrieve either
            PDB data with atomic coordinates (contours outline 3D distance) or
            secondary pdb data (contours outline contact distance).
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to sequence coordinates,
            filtered using interactions_filter arguments, and displayed as
            either nucleotide-resolution heatmaps, or 2D kernel density
            estimate.
            Defaults to None.
        labels (str, optional): Same length as samples list. Labels to
            be used in plot legends. Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to Heatmap(). See
            rnavigate.plots.Heatmap for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Heatmap.plot_data.
            see rnavigate.plots.Heatmap.plot_data for more detail.

    Returns:
        rnavigate.plots.Heatmap plot: object containing matplotlib figure and
            axes with additional plotting and file saving methods
    """
    # TODO: this ones a little tricky to get working with arg parse function
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        interactions=interactions
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Heatmap")
    parsed_args.update_rows_cols(plot_kwargs)
    # initialize plot using 1st 3D pdb (applies to all samples)
    structure = samples[0].data[structure]
    plot = plots.Heatmap(parsed_args.num_samples, structure, **plot_kwargs)
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    return plot


def plot_circle(
        samples,
        sequence=None,
        structure=None,
        structure2=None,
        interactions=None,
        interactions2=None,
        annotations=None,
        profile="default_profile",
        labels=None,
        colorbar=True,
        plot_kwargs=None,
        **kwargs):
    """Generates a multipanel circle plot displaying combinations of secondary
    structures, per-nucleotide data, inter-nucleotide data, and sequence
    annotations. Each plot may display a unique sample and/or filtering scheme.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list, unless filters
            argument below is also used.
        sequence (str or data object, optional): a key from sample.data, a
            sequence string, or a Data object. All data will be mapped to this
            sequence using a user-defined or pairwise sequence alignment.
            Defaults to the value of the structure argument below.
        structure (str, optional): a key from sample.data to retreive a secondary
            structure. Basepairs are plotted as grey arcs within the circle.
            Defaults to "structure".
        structure2 (str, optional): same as structure. basepairs from structure and structure2 will be
            plotted. Basepairs are colored by which structure contains them
            (shared, structure only, structure2 only).
            Defaults to None.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to sequence coordinates,
            filtered using interactions_filter arguments, and displayed as arcs
            within the circle.
            Defaults to None.
        interactions2 (str, optional): same as interactions above.
            Defaults to None.
        annotations (list, optional): a list of keys from sample.data used to
            retrieve sequence annotations. These annotations are displayed next
            to the sequence, outside of the circle.
            Defaults to [].
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data may be used to color nucleotides.
            Defaults to "default_profile".
        labels (str, optional): Same length as samples list. Labels to
            be used as titles.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to Circle(). See
            rnavigate.plots.Circle for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to Circle.plot_data.
            see rnavigate.plots.Circle.plot_data for more detail.

    Returns:
        rnavigate.plots.Circle: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    sequence = get_sequence(sequence, samples[0], structure)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence._alignment,
        structure=structure,
        structure2=structure2,
        interactions2=interactions2,
        profile=profile,
        annotations=annotations,
        interactions=interactions
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.Circle")
    parsed_args.update_rows_cols(plot_kwargs)
    plot = plots.Circle(num_samples=parsed_args.num_samples,
                        sequence=sequence, **plot_kwargs)
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    if colorbar:
        plot.plot_colorbars()
    return plot


def plot_linreg(
        samples,
        sequence=None,
        structure=None,
        profile="default_profile",
        annotations=None,
        labels=None,
        scale='linear',
        regression='pearson',
        plot_kwargs=None,
        **kwargs):
    """Performs linear regression analysis and generates scatter plots of all
    sample-to-sample profile vs. profile comparisons. Colors nucleotides by
    identity or base-pairing status.

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list squared.
        structure (str, optional): a key from sample.data to retreive a secondary
            structure. Scatter plot points may be colored by base-pairing
            status in this structure.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are used for the linear regression.
        labels (str, optional): Same length as samples list. Labels to
            be used in titles.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to LinReg(). See
            rnavigate.plots.LinReg for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to LinReg.plot_data.
            see rnavigate.plots.LinReg.plot_data for more detail.

    Returns:
        rnavigate.plots.LinReg: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    # TODO: maybe LinReg should align each profile individually?
    sequence = get_sequence(sequence, samples[0], profile)
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        alignment=sequence._alignment,
        structure=structure,
        profile=profile,
        annotations=annotations
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.LinReg")
    plot = plots.LinReg(num_samples=parsed_args.num_samples, scale=scale,
                        regression=regression, **plot_kwargs)
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    return plot


def plot_roc(
        samples,
        structure="default_structure",
        profile="default_profile",
        labels=None,
        plot_kwargs=None,
        **kwargs):
    """Performs receiver operator characteristic analysis (ROC), calculates
    area under ROC curve (AUC), and generates ROC plots to assess how well
    per-nucleotide data predicts base-paired status. Does this for all
    positions as well as positions categorized by nucleotide
    5 plots: All, A, U, C, G

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            All samples are plotted on the same set of axes.
        structure (str, optional): a key from sample.data to retreive a secondary
            structure. Base-pairing status retreived from this data.
        profile (str, optional): a key from sample.data used to retrieve per-
            nucleotide data. These data are used for the ROC/AUC analysis.
        labels (str, optional): Same length as samples list. Labels to
            be used in legends.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to ROC(). See
            rnavigate.plots.ROC for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to ROC.plot_data.
            see rnavigate.plots.ROC.plot_data for more detail.

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
        data_dict = fit_data(data_dict, data_dict['structure']._alignment)
        plot.plot_data(**data_dict, **kwargs)
    plot.set_figure_size()
    return plot


def plot_disthist(
        samples,
        structure=None,
        interactions=None,
        bg_interactions=None,
        labels=None,
        same_axis=False,
        plot_kwargs=None,
        **kwargs):
    """Calculates 3D distance of nucleotides in inter-nucleotide data and plots
    the distribution of these distances. Compares this to a 'background'
    distribution consisting of either all pairwise distances in structure, or
    those defined by bg_interactions and bg_interactions_filter

    Args:
        samples (list of rnavigate.Sample): Samples to retreive data from.
            number of panels will equal the length of this list unless
            same_axis is set to True.
        structure (str, optional): a key from sample.data to retreive a PDB
            structure with atomic coordinates.
        interactions (str, optional): a key from sample.data to retrieve inter-
            nucleotide data. These data are mapped to PDB sequence coordinates,
            filtered using interactions_filter arguments, and used to calculate
            distance distribution histograms.
            Defaults to None.
        bg_interactions (str, optional): same as interactions above. used to
            calulate background distance distribution histograms.
            Defaults to None.
        labels (str, optional): Same length as samples list. Labels to
            be used as titles.
            Defaults to default sample name.
        plot_kwargs (dict, optional): kwargs passed to DistHist(). See
            rnavigate.plots.DistHist for more detail.
            Defaults to {}.
        **kwargs: additional keyword arguments are passed to DistHist.plot_data
            see rnavigate.plots.DistHist.plot_data for more detail.

    Returns:
        rnavigate.plots.DistHist: object containing matplotlib figure and axes
            with additional plotting and file saving methods
    """
    # TODO: seperate args for structures and pdbs
    parsed_args = PlottingArgumentParser(
        samples=samples,
        labels=labels,
        structure=structure,
        bg_interactions=bg_interactions,
        interactions=interactions,
    )
    plot_kwargs = _parse_plot_kwargs(plot_kwargs, "rnavigate.plots.DistHist")
    parsed_args.update_rows_cols(plot_kwargs)
    # initialize plot
    if same_axis:
        plot_kwargs |= {'rows': 1, 'cols': 1}
        plot = plots.DistHist(parsed_args.num_samples, **plot_kwargs)
        ax = plot.axes[0, 0]
    else:
        parsed_args.update_rows_cols(plot_kwargs)
        plot = plots.DistHist(parsed_args.num_samples, **plot_kwargs)
        ax = None
    # loop through samples and interactions, adding each as a new axis
    for data_dict in parsed_args.data_dicts:
        plot.plot_data(**data_dict, ax=ax, **kwargs)
    plot.set_figure_size()
    return plot
