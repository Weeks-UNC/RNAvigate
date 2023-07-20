from rnavigate import data, Sample

def _parse_plot_kwargs(plot_kwargs, plot):
    error = ValueError("plot_kwargs must be a dictionary of keyword-arguments "
                    f"to be passed to {plot}.")
    if plot_kwargs is None:
        plot_kwargs = {}
    elif not isinstance(plot_kwargs, dict):
        raise error
    return plot_kwargs


def get_sequence(sequence, sample=None, default=None):
    """Flexible function that returns a Data object containing a sequence.

    Args:
        sequence (any): a sequence string, a data keyword, or a data object
        sample (Sample, optional): sample holding sequence data.
            Defaults to None.
        default (any, optional): same format as sequence above. uses this
            as sequence value if sequence is None.
            Defaults to None.

    Raises:
        ValueError: If sequence and default are both None
        ValueError: If Data object could not be retreived based on inputs

    Returns:
        rnavigate.Sequence: containing the desired sequence
    """
    if sequence is None and default is None:
        raise ValueError("A sequence must be provided.")
    elif sequence is None:
        sequence = default
    if isinstance(sequence, data.Sequence):
        return sequence
    elif isinstance(sequence, str) and os.path.isfile(sequence):
        return data.Sequence(sequence)
    elif isinstance(sequence, str) and sequence in sample.data.keys():
        return sample.data[sequence]
    elif isinstance(sequence, str) and all([nt.upper() in "AUCGT." for nt in sequence]):
        return data.Sequence(sequence=sequence)
    else:
        raise ValueError(f"Cannot find sequence from {sequence}")


def fit_data(data_object, alignment):
    """Given a sample and list of sample.data keys, Data objects are mapped to
    sequence

    Args:
        sample (rnavigate.Sample): sample to retrieve data from
        data_list (list): list of sample.data keys or None
        sequence (rnavigate.Sequence): Data object with a sequence to fit to
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
            data.SequenceAlignment(
                data_object,
                alignment.starting_sequence),
            alignment)
        return data_object.get_aligned_data(alignment)
    else:
        return data_object


class PlottingArgumentParser():
    def __init__(self, samples, labels, alignment=None, **data_dict):
        self.samples = self._parse_samples(samples)
        self.labels = self._parse_labels(labels)
        self.rows = len(samples)
        self.cols = 1
        # parse special formats
        for key, value in data_dict.items():
            if "interactions" in key:
                return_list = key == "interactions"
                data_dict[key] = self._parse_interactions(value, return_list)
                if return_list:
                    self.cols = len(value)
            elif key == "annotations":
                data_dict[key] = self._parse_annotations(value)
        self.data_dicts = []
        classes = { # in order of approximate usage
            'profile': data.Profile,
            'annotations': data.Annotation,
            'structure': data.SecondaryStructure,
            # 'interactions': data.Interactions,
            'pdb': data.PDB,
            'log': data.Log,
        }
        for sample, label in zip(self.samples, self.labels):
            this_data_dict = {"label": label}
            for key, value in data_dict.items():
                for name, data_class in classes.items():
                    if name in key:
                        new_value = sample.get_data(value, data_class)
                        new_value = fit_data(new_value, alignment)
                        this_data_dict[key] = new_value
                        break
                    elif ("interactions" in key) and key != "interactions":
                        sample.filter_interactions(**value)
                        new_value = sample.get_data(
                            value["interactions"], data.Interactions)
                        new_value = fit_data(new_value, alignment)
                        this_data_dict[key] = new_value
                        break
                else:
                    this_data_dict[key] = value
            if "interactions" not in data_dict:
                self.data_dicts.append()
            else:
                for each_interaction in data_dict["interactions"]:
                    sample.filter_interactions(**each_interaction)
                    new_value = each_interaction["interaction"]
                    new_value = sample.get_data(new_value, data.Interactions)
                    new_dict = {"interactions": fit_data(new_value, alignment)}
                    self.data_dicts.append(this_data_dict | new_dict)

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
        error = ValueError("labels must be a list of strings of length equal "
                           "to length of sample list.")
        if labels is None:
            labels = [sample.sample for sample in self.samples]
        elif not isinstance(labels, list):
            raise error
        elif len(labels) != len(self.samples):
            raise error
        return labels

    def _parse_annotations(self, annotations):
        error = ValueError(
            "annotations must be a list of data keywords or Annotation objects.")
        if isinstance(annotations, (data.Annotation, str)):
            annotations = list(annotations)
        elif isinstance(annotations, list):
            for annotation in annotations:
                if not isinstance(annotation, (data.Annotation, str)):
                    raise error
        return annotations

    def _parse_interactions(self, interactions, return_list=True):
        error = ValueError("""
            interactions must follow one of the following formats:
            format 1: a data keyword or data object
            format 2: dictionary containing format 1:
                        {'interactions': format 1, 'filter': True}
            format 3: list of format 2 dictionaries:
                        [format 2, format 2]
            """)
        if interactions is None:
            if return_list:
                return []
            else:
                return None
        if isinstance(interactions, (data.Interactions, str)):
            interactions = {"interactions": interactions}
        if isinstance(interactions, dict):
            if return_list:
                return list(interactions)
            else:
                return interactions
        if isinstance(interactions, list):
            if return_list:
                return interactions
            else:
                raise error
        else:
            raise error
