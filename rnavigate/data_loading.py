import os
from rnavigate import data

_required = object()
data_keyword_defaults = {
    "fasta": {
        "data_class": data.Sequence},
    "log": {
        "data_class": data.Log},
    "shapemap": {
        "data_class": data.SHAPEMaP},
    "dmsmap": {
        "data_class": data.SHAPEMaP,
        "dms": True},
    "dancemap": {
        "data_class": data.DanceMaP,
        "component": _required},
    "rnpmap": {
        "data_class": data.RNPMaP},
    "ringmap": {
        "data_class": data.RINGMaP,
        "sequence": "default_profile"},
    "pairmap": {
        "data_class": data.PAIRMaP,
        "sequence": "default_profile"},
    "allcorrs": {
        "data_class": data.RINGMaP,
        "sequence": "default_profile"},
    "shapejump": {
        "data_class": data.SHAPEJuMP,
        "sequence": _required},
    "pairprob": {
        "data_class": data.PairingProbability,
        "sequence": "default_profile"},
    "ss": {
        "data_class": data.SecondaryStructure},
    "pdb": {
        "data_class": data.PDB,
        "chain": _required},
    "allpossible": {
        "data_class": data.AllPossible,
        "sequence": _required},
    "motif": {
        "data_class": data.Motif,
        "sequence": _required,
        "color": _required},
    "orfs": {
        "data_class": data.ORFs,
        "sequence": _required,
        "color": _required},
    "spans": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "sites": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "group": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "primers": {
        "data_class": data.Annotation,
        "sequence": _required,
        "color": _required},
    "domains": {
        "data_class": data.domains,
        "sequence": _required,
        "colors": _required,
        "names": _required},
}

def create_data(sample=None, **data_keyword):
    """Convenience function for creating rnavigate.data objects. This function
    is used to parse **data_keywords passed to rnavigate.Sample, but can also
    be used on it's own using the same syntax.

    Args:
        sample (rnavigate.Sample, optional):
            Requried only if 'sequence' is provided as a data keyword.
                {'sequence': 'data_keyword'}
                    is replaced with:
                {'sequence': sample.data['data_keyword'].sequence}
            Defaults to None.
        **data_keyword:
            There must only one additional argument provided.
            This argument is used to create the rnavigate.data object, it's
            syntax is flexible, see the example below.

    Usage:
        In this example we create a data object from a fasta file using
        three different syntaxes, starting with the most explicit:

            get_data(arbitrary_name={
                        'data_class':'fasta',
                        'input_data': 'my_sequence.fa'})

        This can be simplified by using the value of 'data_class' to replace
        the key 'input_data':

            get_data(arbitrary_name={
                        'fasta': 'my_sequence.fa'})

        Above, an arbitrary_name is used. This allows Data to be assigned to
        arbitrary data keywords of the given Sample. This name can be replaced
        by our data class:

            get_data(fasta='my_sequence.fa')

        In all three cases, the result returned is equivalent to:

            rnavigate.data.Sequence(input_data="my_sequence.fa")
    """
    # Parse data_keyword for name and inputs
    if len(data_keyword) > 1:
        raise ValueError('Only one data_keyword can be provided')

    name = list(data_keyword.keys())[0]
    inputs = data_keyword[name]

    # if given existing rnavigate.data object return
    if isinstance(inputs, data.Sequence):
        return inputs

    # convert case 1 or 2 into case 3:
    # 1) name='A', inputs='B'
    # 2) inputs={'A': 'B'}
    # 3) inputs={'data_keyword': 'A', 'input_data': 'B'}
    if not isinstance(inputs, dict):
        if name in data_keyword_defaults:
            inputs = {
                'data_keyword': name,
                'input_data': inputs}
    elif all([name in data_keyword_defaults,
              name not in inputs,
              'data_keyword' not in inputs]):
        inputs['data_keyword'] = name
    else:
        for key in list(inputs.keys()):
            if key in data_keyword_defaults:
                inputs['data_keyword'] = key
                inputs['input_data'] = inputs.pop(key)

    # handle special cases
    if inputs['data_keyword'] == 'allpossible':
        inputs['sequence'] = inputs.pop('input_data')
    elif inputs['data_keyword'] in ['sites', 'spans', 'primers', 'group']:
        inputs['annotation_type'] = inputs['data_keyword']

    # 'filepath' or 'dataframe' may be used in place of 'input_data'
    for alt_input_data in ['filepath', 'dataframe']:
        if alt_input_data in inputs:
            inputs['input_data'] = inputs.pop(alt_input_data)
    if 'fasta' in inputs:
        inputs['sequence'] = inputs.pop('fasta')
    elif 'seq_source' in inputs:
        inputs['sequence'] = inputs.pop('seq_source')

    # retrieve defaults for given data_class
    data_class = inputs.pop('data_keyword')
    try:
        inputs = data_keyword_defaults[data_class] | inputs
        data_class = inputs.pop('data_class')
    except KeyError as exception:
        raise KeyError(
            f"{data_class} is not a valid data keyword.") from exception

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
            f"Required arguments for {data_class} were not provided:\n"
            ", ".join(required_arguments))

    # instantiate and return the data class
    try:
        return data_class(**inputs)
    except BaseException as exception:
        raise ValueError(
            f"data_class={data_class}\ninputs={inputs}") from exception


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
        raise ValueError(f"Cannot find sequence from {sequence}")
