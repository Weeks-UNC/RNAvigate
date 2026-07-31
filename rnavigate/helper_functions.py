"""Contains PlottingArgumentParser for plotting_functions.py, fit_data for
retreiving aligned data objects, and resolve_data for resolving flexible
data-keyword arguments into fully-resolved Data objects."""

from rnavigate import Sample, data
from rnavigate.data_loading import get_sequence

__all__ = ["fit_data", "resolve_data", "PlottingArgumentParser"]


def _parse_plot_kwargs(plot_kwargs, plot):
    if plot_kwargs is None:
        plot_kwargs = {}
    elif not isinstance(plot_kwargs, dict):
        raise ValueError(
            "plot_kwargs must be a dictionary of keyword-arguments "
            f"to be passed to {plot}."
        )
    return plot_kwargs


def fit_data(data_object, alignment):
    """Given a sample and list of sample.data keys, Data objects are mapped to
    sequence

    Parameters
    ----------
    data_object : Data object or list or dict of Data objects
        Data object to be mapped to sequence
    alignment : SequenceAlignment
        SequenceAlignment object to be used for mapping

    Returns
    -------
    Data object or list or dict of Data objects
        A new copy of each Data object mapped to a new sequence
    """
    if alignment is None:
        return data_object
    if isinstance(data_object, dict):
        return {k: fit_data(v, alignment) for k, v in data_object.items()}
    if isinstance(data_object, list):
        return [fit_data(v, alignment) for v in data_object]
    if isinstance(data_object, data.Sequence):
        if isinstance(data_object, data.PDB):
            return data_object
        alignment = data.AlignmentChain(
            data.SequenceAlignment(data_object, alignment.starting_sequence), alignment
        )
        return data_object.get_aligned_data(alignment)
    else:
        return data_object


def resolve_data(value, data_class=data.Sequence, sample=None, alignment=None):
    """Resolves a flexible data argument into a fully-resolved Data object.

    This is the single entry point for RNAvigate's flexible data-keyword
    input: a bare Data object, a data keyword string (or list/dict of these,
    resolved via ``sample.get_data``), or -- for `rnavigate.data.Profile` and
    `rnavigate.data.Interactions` -- a dictionary requesting coloring and/or
    filtering.

    Parameters
    ----------
    value : rnavigate.data.Data, str, dict, list, or None
        The value to resolve.
        If None, None is returned.
        If a Data object, it is validated against `data_class` and returned
        (aligned, see `alignment` below).
        If a string (or list/dict of these), it is looked up via
        ``sample.get_data(value, data_class)``.
        If a dict (only supported when `data_class` is
        `rnavigate.data.Profile` or `rnavigate.data.Interactions`, or a
        subclass of either): the dict must contain a "profile" or
        "interactions" key (matching `data_class`) whose value is any of the
        above (a Data object, string, list, or None). The remaining dict keys
        configure coloring and filtering:
            "metric", "cmap", "normalization", "values" : set the resolved
                object's metric/coloring, equivalent to setting `.metric`
                with a dict.
            "metric_kwargs" : dict, additional metric dictionary keys (e.g.
                "color_column", "title", "ticks", "extend"), merged in after
                "metric"/"cmap"/"normalization"/"values".
            "structure" (and "profile", for Interactions) : sibling Data used
                by structure/profile-based filters (`ss_only`, `min_profile`,
                etc). If omitted, defaults to the sample's
                "default_structure"/"default_profile", silently falling back
                to None if the sample has no such default.
            "sequence" : the coordinate frame for `exclude_nts`/`isolate_nts`,
                resolved the same way as the top-level `sequence` plotting
                argument (a data keyword string, a literal sequence string, or
                a Data object) via `get_sequence`, then passed to `.filter()`.
            All other keys are passed to the resolved object's `.filter()`
            method (e.g. `exclude_nts`, `isolate_nts`, `nts`, column
            comparison filters like `Statistic_ge`).
        For Interactions, "metric" also accepts "Distance" or
        "Distance_<atom>", which computes 3D distance using the sample's
        "default_pdb".
    data_class : rnavigate.data.Data class or subclass, or tuple of these, optional
        The expected type of the resolved data. Used to validate bare Data
        objects, to look up data keywords via `sample.get_data`, and to
        determine whether dict input is supported. Defaults to
        `rnavigate.data.Sequence` (i.e. any RNAvigate data object).
    sample : rnavigate.Sample, optional
        The sample used to resolve data keywords (strings) and sibling
        defaults ("default_structure", "default_profile", "default_pdb").
        Required unless `value` is already a Data object or None.
    alignment : rnavigate.data.SequenceAlignment, optional
        If provided, the resolved data is aligned to this alignment's target
        sequence via `fit_data`. If None, no alignment is performed.

    Returns
    -------
    rnavigate.data.Data, list, dict, or None
        The resolved (and aligned) data, matching the structure of `value`.
    """
    if value is None:
        return None
    elif isinstance(value, list):
        return [
            resolve_data(v, data_class, sample=sample, alignment=alignment)
            for v in value
        ]
    elif isinstance(value, data.Sequence):
        if not isinstance(value, data_class):
            raise ValueError(f"{value} is not {data_class}")
    elif isinstance(value, dict):
        if "sequence" in value:
            value["sequence"] = get_sequence(value["sequence"], sample)
        return data_class.resolve_from_dict(value, sample)
    elif isinstance(value, str):
        if sample is None:
            raise ValueError(f"Cannot resolve data keyword {value} without a sample.")
        if not isinstance(sample, Sample):
            raise ValueError(f"sample must be a rnavigate.Sample, not {type(sample)}")
        value = sample.get_data(value, data_class)
    return fit_data(value, alignment)


class PlottingArgumentParser:
    """Parse arguments for high-level plotting functions.

    Given samples list and data keywords, returns a dictionary of data objects
    for each sample, optionally aligned to a target sequence.
    """

    def __init__(self, samples, labels, alignment=None, **data_dict):
        self.samples = self._parse_samples(samples)
        self.labels = self._parse_labels(labels)
        self.rows = len(samples)
        self.cols = 1
        # parse special formats
        for key, value in data_dict.items():
            if key == "interactions":
                data_dict[key] = self._parse_interactions(value)
                self.cols = len(data_dict[key])
            elif "interactions" in key:
                data_dict[key] = self._parse_interactions(value, False)
            elif key == "annotations":
                data_dict[key] = self._parse_annotations(value)
        self.data_dicts = []
        classes = {
            # in order of approximate usage
            "profile": data.Profile,
            "annotations": data.Annotation,
            "structure": (data.SecondaryStructure, data.PDB),
            "interactions": data.Interactions,
            "domains": list,
        }
        for sample, label in zip(self.samples, self.labels):
            this_data_dict = {"label": label}
            for key, value in data_dict.items():
                if key == "sequence":
                    new_value = get_sequence(value, sample)
                    new_value = fit_data(new_value, alignment)
                    this_data_dict[key] = new_value
                    continue
                if key == "interactions":
                    # skip interactions in this loop in case it is a list.
                    continue
                for name, data_class in classes.items():
                    if name in key:
                        this_data_dict[key] = resolve_data(
                            value, data_class, sample=sample, alignment=alignment
                        )
                        break
                else:
                    this_data_dict[key] = value
            if "interactions" not in data_dict:
                self.data_dicts.append(this_data_dict)
            else:
                for each_interaction in data_dict["interactions"]:
                    new_value = resolve_data(
                        each_interaction,
                        data.Interactions,
                        sample=sample,
                        alignment=alignment,
                    )
                    self.data_dicts.append(this_data_dict | {"interactions": new_value})

    def update_rows_cols(self, plot_kwargs):
        if self.rows > 1 and self.cols > 1:
            plot_kwargs["rows"] = self.rows
            plot_kwargs["cols"] = self.cols

    @property
    def num_samples(self):
        return self.rows * self.cols

    def _parse_samples(self, samples):
        error = ValueError("`samples` must be a list of rnavigate.Sample")
        if isinstance(samples, Sample):
            samples = list(samples)
        elif not isinstance(samples, list):
            raise error
        return samples

    def _parse_labels(self, labels):
        """Ensures the value of labels is a list, same length as samples."""
        error = ValueError(
            "labels must be a list of strings of length equal to length of sample list."
        )
        if labels is None:
            labels = [sample.sample for sample in self.samples]
        elif not isinstance(labels, list):
            raise error
        elif len(labels) != len(self.samples):
            raise error
        return labels

    def _parse_annotations(self, annotations):
        """Ensures the value of annotations is a list.

        Parameters
        ----------
        annotations : None, list of str or rnav.data.Annotations
            None is converted to []
            a single string or data object is enclosed in a list

        Returns
        -------
        annotations list
            a list of annotations data or data keywords (or an empty list)
        """
        error = ValueError(
            "annotations must be a list containing data keywords or objects"
        )
        if annotations is None:
            return []
        if isinstance(annotations, list):
            for annotation in annotations:
                if not isinstance(annotation, (data.Annotation, str)):
                    raise error
            return annotations
        if isinstance(annotations, (data.Annotation, str)):
            return [annotations]
        raise error

    def _parse_interactions(self, interactions, return_list=True):
        """Ensure interactions is the expected format.

        interactions must follow one of the following formats:
        format 1: a data keyword or data object
        format 2: dictionary containing format 1:
                    {'interactions': format 1, 'filter': True}
        format 3 (if return_list is True): list of format 2 dictionaries:
                    [format 2, format 2]
        1 is converted to 2 if return_list is False.
        1 and 2 are converted to 3 if return_list is True.

        Parameters
        ----------
        interactions : format 1, 2, or 3
        return_list : bool, defaults to True
            Whether to return format 3, otherwise returns format 2

        Returns
        -------
        interactions dict
            format 2 if return_list is False
        interactions list
            format 3 if return_list is True
        """
        error = ValueError(
            "interactions must follow one of the following formats:\n"
            "format 1: a data keyword or data object\n"
            "format 2: dictionary containing format 1:\n"
            "            {'interactions': format 1, 'filter': True}\n"
            "format 3: list of format 2 dictionaries:\n"
            "            [format 2, format 2]\n"
        )
        if interactions is None:
            if return_list:
                return [{"interactions": None}]
            else:
                return {"interactions": None}
        if isinstance(interactions, (data.Interactions, str)):
            interactions = {"interactions": interactions}
        if isinstance(interactions, dict):
            if return_list:
                return [interactions]
            else:
                return interactions
        if isinstance(interactions, list) and return_list:
            return interactions
        raise error
