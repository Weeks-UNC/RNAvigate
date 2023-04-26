import re
from .data import Data


def is_list_of_ints(my_list, length=None):
    """Checks if the list is a 1-D list of integers. If length is provided,
    also checks that the list has this length. For confirming that annotations
    lists have the correct format.

    Args:
        my_list (list): a list to check for formatting
        length (int, optional): if provided, checks that len(list)==length.
            Defaults to None.

    Returns:
        bool: whether my_list matches the format (or is an empty list)
    """
    if isinstance(my_list, list):
        if length is not None and (len(my_list) != length):
            return False
        for elem in my_list:
            if not isinstance(elem, int):
                return False
        return True
    else:
        return False


def is_list_of_lists_of_ints(my_list, length=None):
    """Checks if inputted list is a 2-D list of lists of integers. If length is
    provided, checks if inner-most lists have the given length. For confirming
    that annotations lists have the correct format.

    Args:
        my_list (list): a list to check for formatting
        length (int, optional): if provided, checks that each inner list has
            this length. Defaults to None.

    Returns:
        bool: whether my_list matches the format (or is an empty list)
    """
    if isinstance(my_list, list):
        return all(is_list_of_ints(elem, length) for elem in my_list)
    else:
        return False


class Annotation(Data):
    def __init__(self,
                 name=None,
                 datatype="annotation",
                 filepath="",
                 fasta=None,
                 sequence=None,
                 annotation_type=None,
                 sites=None,
                 spans=None,
                 groups=None,
                 primers=None,
                 annotations=None,
                 color="blue"):
        """Base annotation class to store 1D features of an RNA. This can
        include groups of separted nucleotides (e.g. binding pocket), spans of
        nucleotides (e.g. coding sequence, Alu elements), list of sites (e.g.
        m6A locations) or primer binding sites.

        Args:
            name (str, optional): Name of annotation.
                Defaults to None.
            datatype (str, optional): annotation datatype.
                Defaults to "annotation".
            filepath (str, optional): path to file containing annotations.
                Defaults to "".
            fasta (str, optional): path to fasta file containing one sequence.
                Defaults to None.
            sequence (str, optional): Nucleotide sequence.
                Defaults to None.
            sites (list of int, optional): 1-indexed location of sites of
                interest within sequence or fasta.
                Defaults to None.
            spans (list of pairs, optional): 1-indexed locations of spans
                of interest within sequence or fasta.
                e.g. [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30.
                Defaults to None.
            groups (list of lists, optional): 1-indexed locations of groups of
                sites of interest.
                Defaults to None.
            primers (list of pairs, optional): Similar to spans above,
                but reverse primers are reverse ordered.
                e.g. [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20.
                Defaults to None.
            annotations (list, optional): catchall for the above list formats.
                List will be treated according to annotation_type argument.
                Defaults to None.
            annotation_type (str, optional): "groups", "sites", "spans", or
                "primers". Must match the type used above. Required for
                annotations argument to be correctly parsed.
                Defaults to None.
            color (matplotlib color-like, optional): Color to be used for
                displaying this annotation on plots.
                Defaults to "blue".
        """
        self.name = name
        self.color = color
        self.datatype = datatype
        super().__init__(filepath=fasta, sequence=sequence)

        # make sure whichever list is stored matches expected format
        if (annotations is not None) and (annotation_type is not None):
            self.annotation_type = annotation_type
            if self.annotation_type == "sites":
                sites = annotations
            elif self.annotation_type == "spans":
                spans = annotations
            elif self.annotation_type == "primers":
                primers = annotations
            elif self.annotation_type == "groups":
                groups = annotations
        elif (annotations is not None) and (annotation_type is None):
            raise ValueError("annotation requires annotation_type")

        if sites is not None:
            if not is_list_of_ints(sites):
                raise ValueError(f"sites must be a list of integers:\n{sites}")
            self.annotation_type = "sites"
            self._list = sites
        if spans is not None:
            if not is_list_of_lists_of_ints(spans, length=2):
                raise ValueError("spans format is must be a list of lists "
                                 f"containing 2 integers:\n{spans}")
            self.annotation_type = "spans"
            self._list = spans
        if primers is not None:
            if not is_list_of_lists_of_ints(primers, length=2):
                raise ValueError("primers format is must be a list of lists "
                                 f"containing 2 integers:\n{primers}")
            self.annotation_type = "primers"
            self._list = primers
        if groups is not None:
            for group in groups:
                if not set(group.keys()) == set(["sites", "colors"]):
                    raise ValueError("Groups must have 'sites' and 'colors'")
                if not is_list_of_ints(group["sites"]):
                    raise ValueError("'sites' must be a list of lists")
            self.annotation_type = "groups"
            self._list = groups

    def fit_to(self, fit_to, store=True):
        """Creates a new Annotation, stored as self.fitted, which maps the
        indices to a new sequence. For sites and groups, deleted nucleotides
        will be dropped. For spans and primers, a deletion/insertion within a
        region will shrink/expand the annotated region, respectively. If either
        end of a region is deleted, that region is dropped. For different
        behavior, create a new annotation for the target sequence.

        Args:
            fit_to (rnavigate.data.Data): A data object containing a sequence.
            store (bool, optional): whether to store result as self.fitted, if
                False, returns the value instead.
        """
        am = self.get_alignment_map(fit_to=fit_to)

        def recursive_fit_to(indices):
            # given nested list of indices, returns same shape with indices
            # mapped to fit_to sequence
            new_list = []
            for idx in indices:
                if isinstance(idx, list):
                    new_list.append(recursive_fit_to(idx))
                elif isinstance(idx, int):
                    new_idx = int(am[idx-1]+1)
                    if new_idx != 0:
                        new_list.append(int(am[idx-1]+1))
            return new_list

        if self.annotation_type == "groups":
            new_list = [{"sites": recursive_fit_to(d["sites"]),
                        "color": d["color"]} for d in self._list]
        else:
            new_list = recursive_fit_to(self._list)
        if self.annotation_type in ["spans", "primers"]:
            new_list = [x for x in new_list if len(x) == 2]

        fitted = Annotation(
            name=self.name,
            color=self.color,
            sequence=fit_to,
            annotation_type=self.annotation_type,
            annotations=new_list)
        if store:
            self.fitted = fitted
        else:
            return fitted

    def get_sites_colors(self, fit_to=None):
        """Returns a list of nucleotide positions and colors based on these
        sequence annotations. If a fit_to data object is provided, annotations
        are first fitted using self.fit_to()

        Returns:
            tuple: a list of nucleotide positions and a list of colors
        """
        if fit_to is not None:
            return self.fit_to(fit_to, store=False).get_sites_colors()
        sites = []
        colors = []
        if self.annotation_type in ["spans", "primers"]:
            for start, stop in self._list:
                if start > stop:
                    start, stop = stop, start
                for nt in range(start, stop + 1):
                    sites.append(nt)
                    colors.append(self.color)
        elif self.annotation_type == "sites":
            for site in self._list:
                sites.append(site)
                colors.append(self.color)
        elif self.annotation_type == "groups":
            for group in self._list:
                colors.extend([group["color"]] * len(group["sites"]))
                sites.extend[group["sites"]]
        return sites, colors

    def __getitem__(self, i):
        return self._list[i]

    def __len__(self):
        return len(self._list)


class Motif(Annotation):
    def __init__(self,
                 name=None,
                 filepath="",
                 fasta=None,
                 sequence=None,
                 motif=None,
                 color="blue"):
        """Creates a Motif annotation, which acts like a span Annotation, for
        highlighting a sequence motif of interest, given with conventional
        nucleotide codes. e.g. "DRACH"

        Args:
            name (str, optional): name of this annotation.
                Defaults to None.
            filepath(str, optional): Defaults to "".
            fasta (str, optional): path to fasta file containing sequence.
                Defaults to None.
            sequence (str, optional): sequence to be searched.
                Defaults to None.
            motif (str, optional): sequence motif to be searched for.
                Defaults to None.
            color (str, optional): color used to display these motif locations.
                Defaults to "blue".
        """
        self.motif = motif
        span_list = self.get_spans_from_motif(sequence, motif)
        super().__init__(name=name, fasta=fasta, sequence=sequence,
                         annotation_type='spans',
                         spans=span_list, color=color)

    def get_spans_from_motif(self, sequence, motif):
        """Returns a list of spans [[start, end], [start, end]] for each
        location of motif found within sequence, using conventional nucleotide
        codes.

        Args:
            sequence (str): sequence to be searched
            motif (str): sequence motif to be searched for.

        Returns:
            _type_: _description_
        """
        nuc_codes = {"A": "A", "T": "T", "U": "U", "G": "G", "C": "C",
                     "B": "[CGTU]", "D": "[ATUG]", "H": "[ATUC]", "V": "[ACG]",
                     "W": "[ATU]", "S": "[CG]",  # strong and weak
                     "M": "[AC]", "K": "[GTU]",  # amino and ketone
                     "R": "[AG]", "Y": "[CTU]",  # purine and pyrimidine
                     "N": "[ATUGC]"}  # any nuc
        re_pattern = ''.join([nuc_codes[n] for n in motif])
        spans = []
        for match in re.finditer(re_pattern, sequence):
            start, end = match.span()
            spans.append([start+1, end])
        return spans

    def fit_to(self, fit_to, store=True):
        """Creates a new annotation, stored as self.fitted. Fit_to sequence is
        searched for motif matches

        Args:
            fit_to (Data object): data containing new sequence to be fitted
            store (bool, optional): If True, new annotations are stored as
                self.fitted. Else, new annotations are returned.
                Defaults to True."""
        fitted = Motif(
            name=self.name,
            color=self.color,
            sequence=fit_to.sequence,
            motif=self.motif,
        )
        if store:
            self.fitted = fitted
        else:
            return fitted


class ORFs(Annotation):
    def __init__(self,
                 name=None,
                 filepath=None,
                 fasta=None,
                 sequence=None,
                 color="blue"):
        span_list = self.get_spans_from_orf(sequence)
        super().__init__(name=name, fasta=fasta, sequence=sequence,
                         annotation_type='spans',
                         spans=span_list, color=color)

    def get_spans_from_orf(self, sequence):
        spans = []
        stop_codons = "UAA|UAG|UGA"
        stop_sites = []
        for match in re.finditer(stop_codons, sequence):
            stop_sites.append(match.span()[1])
        start_codon = "AUG"
        start_sites = []
        for match in re.finditer(start_codon, sequence):
            start_sites.append(match.span()[0]+1)
        for start in start_sites:
            for stop in stop_sites:
                if ((stop-start) % 3 == 2) and (start < stop):
                    spans.append([start, stop])
        return spans
