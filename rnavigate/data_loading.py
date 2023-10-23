import os
from rnavigate import data

_required = object()

data_keyword_defaults = {
    "sequence": {"data_class": data.Sequence},
    "shapemap": {"data_class": data.SHAPEMaP},
    "dancemap": {"data_class": data.DanceMaP, "component": _required},
    "rnpmap": {"data_class": data.RNPMaP},
    "ringmap": {"data_class": data.RINGMaP, "sequence": "default_profile"},
    "pairmap": {"data_class": data.PAIRMaP, "sequence": "default_profile"},
    "shapejump": {"data_class": data.SHAPEJuMP, "sequence": _required},
    "pairprob": {"data_class": data.PairingProbability,
                    "sequence": "default_profile"},
    "ss": {"data_class": data.SecondaryStructure},
    "pdb": {"data_class": data.PDB, "chain": _required},
    "allpossible": {"data_class": data.AllPossible, "sequence": _required},
    "motif": {"data_class": data.Motif, "name": _required,
                "sequence": _required, "color": _required},
    "orfs": {"data_class": data.ORFs, "name": _required,
                "sequence": _required, "color": _required},
    "spans": {"data_class": data.Annotation, "name": _required,
                "sequence": _required, "color": _required},
    "sites": {"data_class": data.Annotation, "name": _required,
                "sequence": _required, "color": _required},
    "group": {"data_class": data.Annotation, "name": _required,
                "sequence": _required, "color": _required},
    "primers": {"data_class": data.Annotation, "name": _required,
                "sequence": _required, "color": _required},
    "domains": {"data_class": data.domains, "sequence": _required,
                "colors": _required, "names": _required},
    }

def create_data(data_keyword, inputs, sample=None):
    """Convenience function for creating rnavigate.data objects. This function
    is used to parse **data_keywords passed to rnavigate.Sample, but can also
    be used on it's own using the same syntax.

    Required arguments:
        data_keyword (string)
            a standard data keyword
            if inputs is a dictionary, this can be an arbitrary data keyword
        inputs (any)
            If a dictionary, the first key must be the standard data keyword
            Otherwise, the expected type depends on the data_keyword

    Optional arguments:
        sample (rnavigate.Sample)
            Requried only if inputs is a dictionary containing a 'sequence' key
            with a value that is a standard data keyword
                {'sequence': 'data_keyword'}
                    is replaced with:
                {'sequence': sample.get_data('data_keyword').sequence}
            Defaults to None.

    Returns: new RNAvigate data object

    Usage:
        In this example we create a data object from a fasta file. There are
        two valid syntaxes:

            create_data(
                data_keyword="sequence",
                inputs="my_sequence.fa"
                )

        OR

            create_data(
                data_keyword="arbitrary_string",
                inputs={'sequence': 'my_sequence.fa'}
                )

        Returns: rnavigate.data.Sequence(input_data="my_sequence.fa")
    """
    # if given existing rnavigate.data object return
    if isinstance(inputs, data.Sequence):
        return inputs

    # data_keyword='A', inputs='B'
    # OR
    # data_keyword='arbitrary', inputs={'A': 'B', etc.}
    # becomes
    # data_keyword='A', inputs={'input_data': 'B', etc.}
    if isinstance(inputs, dict):
        first_key = list(inputs.keys())[0]
        data_keyword=first_key
        inputs["input_data"] = inputs.pop(first_key)
    else:
        inputs = {'data_keyword': data_keyword, 'input_data': inputs}

    # handle special cases
    if data_keyword == 'allpossible':
        inputs['sequence'] = inputs.pop('input_data')
    elif data_keyword in ['sites', 'spans', 'primers', 'group']:
        inputs['annotation_type'] = data_keyword

    # retrieve defaults for given data_class
    try:
        inputs = data_keyword_defaults[data_keyword] | inputs
        data_constructor = inputs.pop('data_class')
    except KeyError as exception:
        raise KeyError(
            f"{data_keyword} is not a valid data keyword.") from exception

    # get sequence from another object if appropriate
    if 'sequence' in inputs:
        inputs["sequence"] = get_sequence(inputs["sequence"], sample)

    # check if any required arguments are not provided
    required_arguments = []
    for key, value in inputs.items():
        if value is _required:
            required_arguments.append(key)
    if len(required_arguments) > 0:
        raise ValueError(
            f"{data_keyword} missing required arguments:\n"
            ", ".join(required_arguments))

    # instantiate and return the data class
    try:
        return data_constructor(**inputs)
    except BaseException as exception:
        raise ValueError(
            f"data_keyword={data_keyword}\ninputs={inputs}"
            ) from exception


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
    try:
        sequence = sample.get_data(sequence)
    except:
        pass
    if isinstance(sequence, data.Sequence):
        return sequence
    if isinstance(sequence, str) and os.path.isfile(sequence):
        return data.Sequence(sequence)
    if isinstance(sequence, str) and all([nt.upper() in "AUCGT." for nt in sequence]):
        return data.Sequence(sequence)
    else:
        raise ValueError(f'Cannot find sequence from "{sequence}"')
