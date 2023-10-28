from os.path import isfile
import Bio.SeqIO
import numpy as np
import pandas as pd
import matplotlib.colors as mpc
from rnavigate import styles
from rnavigate import data



class Sequence():
    def __init__(self, input_data):
        """Constructs a Data object given a sequence string, fasta file, or
        dataframe containing a "Sequence" column.

        Args:
            sequence (str | pandas.DataFrame):
                sequence string, fasta file, or a pandas dataframe containing
                a "Sequence" column.
        """
        self.other_info = dict()
        if isinstance(input_data, str):
            if isfile(input_data):
                self.sequence = self.read_fasta(input_data)
            else:
                self.sequence = input_data
        elif isinstance(input_data, pd.DataFrame):
            self.sequence = self.get_seq_from_dataframe(input_data)
        elif isinstance(input_data, Sequence):
            self.sequence = input_data.sequence
        self.null_alignment = data.SequenceAlignment(self, self)

    def read_fasta(self, fasta):
        """Parse a fasta file for the first sequence. Store the sequence name
        as self.gene and the sequence string as self.sequence.

        Args:
            fasta (str): path to fasta file
        """
        with open(fasta, 'r') as file:
            fasta = list(Bio.SeqIO.parse(file, 'fasta'))
        return str(fasta[0].seq).replace("T", "U")

    def get_seq_from_dataframe(self, dataframe):
        """Parse a dataframe for the sequence string, store as self.sequence.

        Args:
            dataframe (pandas DataFrame): must contain a "Sequence" column
        """
        sequence = ''.join(dataframe["Sequence"].values)
        return sequence.replace("T", "U").replace("t", "u")

    @property
    def length(self):
        """Get the length of the sequence

        Returns:
            int: the length of self.sequence
        """
        return len(self.sequence)

    def get_aligned_data(self, alignment):
        return Sequence(alignment.target_sequence)

    def get_colors_from_sequence(self):
        colors = np.array([styles.get_nt_color(nt) for nt in self.sequence])
        colormap = data.ScalarMappable(
            cmap=[styles.get_nt_color(nt) for nt in 'AUGC'],
            normalization='none', values=None, extend='neither',
            title='Nucleotide identity', alpha=1, ticks=[0, 1, 2, 3],
            tick_labels=['A', 'U', 'G', 'C']
            )
        return colors, colormap

    def get_colors_from_positions(self, pos_cmap='rainbow'):
        colormap = data.ScalarMappable(
            cmap=pos_cmap, normalization='min_max', values=[1, self.length],
            extend='neither', title='Nucleotide position', alpha=1
            )
        colors = colormap.values_to_hexcolors(np.arange(self.length))
        return colors, colormap

    def get_colors_from_profile(self, profile):
        alignment = data.SequenceAlignment(profile, self)
        colors = alignment.map_values(profile.colors, fill="#808080")
        colormap = profile.cmap
        return colors, colormap

    def get_colors_from_annotations(self, annotations):
        colors = np.full(self.length, 'gray', dtype='<U16')
        cmap = ['gray']
        tick_labels = ['other']
        for annotation in annotations:
            cmap.append(annotation.color)
            tick_labels.append(annotation.name)
            annotation = annotation.get_aligned_data(
                data.SequenceAlignment(annotation, self))
            for site in annotation.get_sites():
                colors[site-1] = annotation.color
        colormap = data.ScalarMappable(
            cmap=cmap, normalization='none', values=None, extend='neither',
            title='Annotations', alpha=1, ticks=list(range(len(annotations))),
            tick_labels=tick_labels,
            )
        return colors, colormap

    def get_colors_from_structure(self, structure):
        cmap = ['darkOrange', 'darkOrchid', 'gray']
        ct_colors = [cmap[int(pair == 0)] for pair in structure.pair_nts]
        alignment = data.SequenceAlignment(structure, self)
        colors = alignment.map_values(ct_colors, fill='gray')
        colormap = data.ScalarMappable(
            cmap=cmap, normalization='none', values=None, extend='neither',
            title='Base-pairing status', alpha=1, ticks=[0, 1, 2],
            tick_labels=['unpaired', 'paired', 'unaligned'],
            )
        return colors, colormap

    def get_colors(self, source, pos_cmap='rainbow', profile=None,
                   structure=None, annotations=None):
        """Get a numpy array of colors that fits the current sequence.

        Args:
            source (str | array of color-like): One of the following:
                "position": colors represent position in sequence
                "sequence": colors represent nucleotide identity
                "annotations": colors represent sequence annotations
                "profile": colors represent per-nucleotide data
                "structure": colors represent base-pairing status
                matplotlib color-like: all colors are this color
                array of color like: must match length of sequence
            pos_cmap (str, optional): cmap used if source="position".
                Defaults to 'rainbow'.
            profile (Profile or subclass, optional): Data object containing
                per-nucleotide information. Defaults to None.
            structure (CT of subclass, optional): Data object containing secondary
                structure information. Defaults to None.
            annotations (list of Annotations or subclass, optional): list of
                Data objects containing annotations. Defaults to None.

        Returns:
            numpy array: one matplotlib color-like value for each nucleotide in
                self.sequence
        """
        if isinstance(source, str) and (source == "sequence"):
            return self.get_colors_from_sequence()
        elif isinstance(source, str) and (source == "position"):
            return self.get_colors_from_positions(pos_cmap=pos_cmap)
        elif isinstance(source, str) and (source == "profile"):
            return self.get_colors_from_profile(profile=profile)
        elif isinstance(source, str) and (source == "annotations"):
            return self.get_colors_from_annotations(annotations=annotations)
        elif isinstance(source, str) and (source == "structure"):
            return self.get_colors_from_structure(structure)
        elif mpc.is_color_like(source):
            return np.full(self.length, source, dtype="<U16"), None
        elif ((len(source) == self.length)
              and all(mpc.is_color_like(c) for c in source)):
            return np.array(list(source)), None
        else:
            print('Invalid colors:\n\tchoices are "profile", "sequence", '
                  '"position", "structure", "annotations", a list of mpl '
                  'colors, or a single mpl color.\nDefaulting to sequence.')
            return self.get_colors_from_sequence()


class Data(Sequence):
    def __init__(self, input_data, sequence, metric, metric_defaults,
                 read_table_kw=None):
        if read_table_kw is None:
            read_table_kw = {}
        # assign data
        if isinstance(input_data, pd.DataFrame):
            self.data = input_data
            self.filepath = 'dataframe'
        elif isfile(input_data):
            self.data = self.read_file(input_data, read_table_kw)
            self.filepath = input_data
        else:
            print(f"{self} initialized without data.")
        # assign sequence
        if sequence is None:
            sequence = self.data
        super().__init__(sequence)
        # assign metrics
        self.metric_defaults = {
            'default': {
                'metric_column': 'Profile',
                'error_column': None,
                'color_column': None,
                'cmap': 'viridis',
                'normalization': '0_1',
                'values': None,
                'extend': 'neither',
                'title': 'Type: metric',
                'alpha': 0.7},
            "Distance": {
                "metric_column": "Distance",
                'error_column': None,
                'color_column': None,
                "cmap": "cool",
                "normalization": "min_max",
                "values": [5, 50],
                'extend':'both',
                'title':'3D distance',
                'alpha': 0.7},
            }
        self.add_metric_defaults(metric_defaults)
        self.default_metric = metric
        self._metric = None
        self.metric = metric

    def add_metric_defaults(self, metric_defaults):
        default_defaults = self.metric_defaults['default']
        for metric, defaults in metric_defaults.items():
            self.metric_defaults[metric] = default_defaults | defaults

    @property
    def metric(self):
        return self._metric['metric_column']

    @metric.setter
    def metric(self, value):
        if isinstance(value, str):
            try:
                self._metric = self.metric_defaults[value]
                return
            except KeyError as exception:
                if value not in self.data.columns:
                    print(f"metric ({value}) not found in data:\n"
                          + str(list(self.data.columns)))
                    raise exception
                self._metric = self.metric_defaults['default']
                self._metric['metric_column'] = value
                return
        try:
            defaults = self.metric_defaults[value["metric_column"]]
        except KeyError:
            defaults = self.metric_defaults['default']
        for key in value:
            if key not in defaults:
                raise ValueError(f'{key} is not an expected value for metric '
                                 f'setting:\n{list(defaults.keys())}')
        self._metric = defaults | value

    @property
    def error_column(self):
        return self._metric['error_column']

    @property
    def color_column(self):
        if self._metric['color_column'] is None:
            return self._metric['metric_column']
        else:
            return self._metric['color_column']

    @property
    def cmap(self):
        cmap_kwargs = {}
        for k, v in self._metric.items():
            if not k.endswith("column"):
                cmap_kwargs[k] = v
        return data.ScalarMappable(**cmap_kwargs)

    @property
    def colors(self):
        return self.cmap.values_to_hexcolors(
            np.ma.masked_invalid(self.data[self.color_column])
            )

    def read_file(self, filepath, read_table_kw):
        """Convert data file to pandas dataframe and store as self.data

        Args:
            filepath (str): path to data file containing interactions
            read_table_kw (dict): kwargs dictionary passed to pd.read_table
        """
        return pd.read_table(filepath, **read_table_kw)
