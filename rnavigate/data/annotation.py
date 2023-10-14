"""annotations.py contains Annotations and subclasses."""
import re
import pandas as pd
from rnavigate import data


class Annotation(data.Sequence):
    def __init__(self, input_data, annotation_type,
                 sequence, name=None, color="blue"):
        """Base annotation class to store 1D features of an RNA. This can
        include a group of separted nucleotides (e.g. binding pocket), spans of
        nucleotides (e.g. coding sequence, Alu elements), list of sites (e.g.
        m6A locations) or primer binding sites.

        Args:
            input_data (list):
                List will be treated according to annotation_type argument.
                Expected behaviors for each value of annotation_type:
                "sites" or "group": 1-indexed location of sites of interest
                "spans": 1-indexed, inclusive locations of spans of interest
                    e.g. [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30
                "primers": Similar to spans, but 5'/3' direction is preserved.
                    e.g. [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20
            annotation_type (str):
                "group", "sites", "spans", or "primers".
            sequence (str | pandas.DataFrame):
                Nucleotide sequence, path to fasta file, or dataframe
                containing a "Sequence" column.
            name (str, optional): Name of annotation.
                Defaults to None.
            color (matplotlib color-like, optional): Color to be used for
                displaying this annotation on plots.
                Defaults to "blue".
        """
        self.name = name
        self.color = color
        super().__init__(sequence)

        # make sure input data matches expected format
        valid_types = ["sites", "spans", "group", "primers"]
        self.annotation_type = annotation_type
        if isinstance(input_data, pd.DataFrame):
            self.data = input_data
        elif annotation_type in ['spans', 'primers']:
            self.data = self.from_spans(input_data)
        elif annotation_type in ['sites', 'group']:
            self.data = self.from_sites(input_data)
        else:
            raise ValueError(f"annotation_type not in {valid_types}")

    @classmethod
    def from_boolean_array(self, values, window, sequence, annotation_type,
                           name=None, color='blue'):
        annotations = []
        current_annotation = None
        pad = window // 2
        for i, value in enumerate(values):
            position = i + 1
            if value and annotation_type in ['sites', 'group']:
                annotations.append(position)
            elif value and annotation_type == 'spans':
                start = max(1, position - pad)
                stop = min(len(sequence), position + pad)
                if current_annotation is None:
                    current_annotation = [start, stop]
                elif start <= current_annotation[1]+1:
                    current_annotation[1] = stop
                elif start > current_annotation[1]:
                    annotations.append(current_annotation)
                    current_annotation = [start, stop]
        if annotation_type == 'spans':
            annotations.append(current_annotation)
        return Annotation(
            input_data=annotations, annotation_type=annotation_type,
            color=color, sequence=sequence, name=name)

    def from_spans(self, spans):
        data_dict = {'start':[], 'end':[]}
        for span in spans:
            if ((len(span) != 2)
                    and any(not isinstance(pos, int) for pos in span)):
                raise ValueError(
                    f'{self.annotation_type} must be a list of pairs of '
                    f'integers:\n{spans}'
                )
            start, end = span
            data_dict['start'].append(int(start))
            data_dict['end'].append(int(end))
        return pd.DataFrame(data_dict)

    def from_sites(self, sites):
        if any(not isinstance(site, int) for site in sites):
            raise ValueError(
                f'{self.annotation_type} must be a list of integers:\n{sites}')
        return pd.DataFrame({'site': sites})

    def get_aligned_data(self, alignment):
        new_input_data = alignment.map_dataframe(self.data, self.data.columns)
        return Annotation(
            name=self.name,
            color=self.color,
            sequence=alignment.target_sequence,
            annotation_type=self.annotation_type,
            input_data=new_input_data,
        )

    def get_sites(self):
        """Returns a list of nucleotide positions and colors based on these
        sequence annotations.

        Returns:
            tuple: a list of nucleotide positions and a list of colors
        """
        if self.annotation_type in ["spans", "primers"]:
            sites = []
            for _, row in self.data.iterrows():
                sites.extend(list(range(row['start'], row['end']+1)))
        elif self.annotation_type in ["sites", "group"]:
            sites = list(self.data['site'].values)
        return sites

    def __getitem__(self, idx):
        return self.data.loc[idx]

    def __iter__(self):
        for _, row in self.data.iterrows():
            yield row

    def __len__(self):
        return len(self._list)


class Motif(Annotation):
    def __init__(self, input_data, sequence,
                 name=None, color="blue"):
        """Creates a Motif annotation, which acts like a span Annotation, for
        highlighting a sequence motif of interest, given with conventional
        nucleotide codes. e.g. "DRACH"

        Args:
            sequence (str | pandas.DataFrame):
                sequence to be searched.
                Defaults to None.
            motif (str, optional): sequence motif to be searched for.
                Defaults to None.
            name (str, optional): name of this annotation.
                Defaults to None.
            color (str, optional): color used to display these motif locations.
                Defaults to "blue".
        """
        self.motif = input_data
        span_list = self.get_spans_from_motif(sequence, input_data)
        super().__init__(
            name=name,
            sequence=sequence,
            annotation_type="spans",
            input_data=span_list,
            color=color,
        )

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
        nuc_codes = {
            "A": "A",
            "T": "T",
            "U": "U",
            "G": "G",
            "C": "C",
            "B": "[CGTU]",
            "D": "[ATUG]",
            "H": "[ATUC]",
            "V": "[ACG]",
            "W": "[ATU]",
            "S": "[CG]",  # strong and weak
            "M": "[AC]",
            "K": "[GTU]",  # amino and ketone
            "R": "[AG]",
            "Y": "[CTU]",  # purine and pyrimidine
            "N": "[ATUGC]",
        }  # any nuc
        re_pattern = "".join([nuc_codes[n] for n in motif])
        spans = []
        for match in re.finditer(re_pattern, sequence):
            start, end = match.span()
            spans.append([start + 1, end])
        return spans

    def get_aligned_data(self, alignment):
        return Motif(
            input_data=self.motif,
            name=self.name,
            color=self.color,
            sequence=alignment.target_sequence,
            )


class ORFs(Annotation):
    def __init__(self, name=None, sequence=None, color="blue"):
        span_list = self.get_spans_from_orf(sequence)
        super().__init__(
            name=name, sequence=sequence, annotation_type="spans",
            input_data=span_list, color=color,
            )

    def get_spans_from_orf(self, sequence):
        """Given a sequence string, returns all possible ORFs

        Args:
            sequence (str): RNA nucleotide sequence

        Returns:
            list of tuples: (start, end) position of each ORF
                1-indexed, inclusive
        """
        spans = []
        stop_codons = "UAA|UAG|UGA"
        stop_sites = []
        for match in re.finditer(stop_codons, sequence):
            stop_sites.append(match.span()[1])
        start_codon = "AUG"
        start_sites = []
        for match in re.finditer(start_codon, sequence):
            start_sites.append(match.span()[0] + 1)
        for start in start_sites:
            for stop in stop_sites:
                if ((stop - start) % 3 == 2) and (start < stop):
                    spans.append([start, stop])
        return spans

    def get_aligned_data(self, alignment):
        return ORFs(
            name=self.name,
            color=self.color,
            sequence=alignment.target_sequence)

def domains(input_data, names, colors, sequence):
    return [
        Annotation(
            input_data=[span], annotation_type='spans', name=name,
            color=color, sequence=sequence,
            ) for span, name, color in zip(input_data, names, colors)
        ]
