import pandas as pd
import numpy as np
from Bio import SeqIO
from rnavigate import data
from pathlib import Path



data_path = Path(__file__).parent.parent.parent / 'reference_data'


class Transcript(data.Sequence):
    def __init__(self, sequence, chromosome, strand, coordinates, tx_info,
                 cds_coors=None, other_features=None):
        super().__init__(sequence)
        self.chromosome = chromosome
        self.strand = strand
        self.coordinates = coordinates
        self.tx_info = tx_info
        self.cds_coors=cds_coors
        self.other_features = other_features
        self.coordinate_df = self.get_coordinate_df()

    def get_coordinate_df(self):
        coordinates = [np.arange(x, y+1) for x, y in zip(*self.coordinates)]
        coordinates = np.concatenate(coordinates)
        if self.strand == '-':
            coordinates = coordinates[::-1]
        df = pd.DataFrame({
            'Sequence': list(self.sequence),
            'Nucleotide': np.arange(self.length)+1,
            'Coordinate': coordinates,
        })
        return df

    def get_tx_coordinate(self, coordinate):
        df = self.coordinate_df
        coordinate = df['Coordinate'] == coordinate
        if sum(coordinate) == 0:
            return None
        return df.loc[coordinate, 'Nucleotide'].values[0]

    def get_tx_range(self, start, stop):
        df = self.coordinate_df
        nts = df.loc[df['Coordinate'].between(start, stop), 'Nucleotide']
        return [nts.min(), nts.max()]

    def get_cds_annotation(self, **kwargs):
        cds = self.get_tx_range(*sorted(self.cds_coors))
        return data.Annotation(
            input_data=[cds],
            annotation_type='spans',
            sequence=self.sequence,
            name="CDS",
            **kwargs)

    def get_junctions_annotation(self, **kwargs):
        junctions = [self.coordinates[0][1:], self.coordinates[1][:-1]]
        junctions = [self.get_tx_range(*sorted(junc)) for junc in zip(*junctions)]
        return data.Annotation(
            input_data=junctions,
            annotation_type='spans',
            sequence=self.sequence,
            name="exon junctions",
            **kwargs)

    def get_exon_annotation(self, exon_number, **kwargs):
        if self.strand == '-':
            exon_number = len(self.coordinates[0]) + 1 - exon_number
        try:
            exon = [
                self.coordinates[0][exon_number-1],
                self.coordinates[1][exon_number-1]]
        except IndexError as e:
            raise ValueError(f"exon number {exon_number} is too high") from e
        exon = self.get_tx_range(*sorted(exon))
        return data.Annotation(
            input_data=[exon],
            annotation_type='spans',
            sequence=self.sequence,
            name=f"exon {exon_number}",
            **kwargs)


class Transcriptome():
    def __init__(self, genome, annotation, path=data_path):
        if path is not data_path:
            path = Path(path)
        self.genome = data_path / genome
        self.annotation = data_path / annotation
        self.start_coors = []
        self.end_coors = []

    def get_transcript(self, transcript_id):
        def parse_entry(entry):
            entry = entry.strip().split('\t')
            if entry[2] == "transcript":
                info = entry[8].split(';')
                info = [item.strip().split(' ') for item in info[:-1]]
                info = {k: v.strip('"') for k, v in info}
            else:
                info = {}
            return {
                'chromosome': entry[0],
                'db': entry[1],
                'feature': entry[2],
                'start': entry[3],
                'stop': entry[4],
                'strand': entry[6],
                'info': info
            }

        with open(self.annotation) as txome:
            start_coors = []
            stop_coors = []
            cds_coors = []
            other_features = []
            for entry in txome.readlines():
                if f'transcript_id "{transcript_id}"' not in entry:
                    continue
                tx = parse_entry(entry)
                if tx['feature'] == 'transcript':
                    tx_info = tx['info']
                    chromosome = tx['chromosome']
                    strand = tx['strand']
                elif tx['feature'] == 'exon':
                    start_coors.append(int(tx['start']))
                    stop_coors.append(int(tx['stop']))
                elif tx['feature'] in ['start_codon', 'stop_codon']:
                    cds_coors.append(int(tx['start']))
                    cds_coors.append(int(tx['stop']))
                elif tx['feature'] in ['UTR', 'CDS']:
                    continue
                else:
                    other_features.append(
                        {key: tx[key] for key in ['feature', 'start', 'stop']})
        cds_coors = [min(cds_coors), max(cds_coors)]
        coordinates = [sorted(start_coors), sorted(stop_coors)]
        sequence = self.get_sequence(chromosome, coordinates, strand)
        return Transcript(
            sequence=sequence,
            chromosome=chromosome,
            strand=strand,
            coordinates=coordinates,
            tx_info=tx_info,
            cds_coors=cds_coors,
            other_features=other_features)

    # create fasta file from genomic coordinate using pybedtools get fasta function
    def get_sequence(self, chrom, coordinates, strand):
        seq = ''
        chrom = {
            'chr1':  'NC_000001',
            'chr2':  'NC_000002',
            'chr3':  'NC_000003',
            'chr4':  'NC_000004',
            'chr5':  'NC_000005',
            'chr6':  'NC_000006',
            'chr7':  'NC_000007',
            'chr8':  'NC_000008',
            'chr9':  'NC_000009',
            'chr10': 'NC_000010',
            'chr11': 'NC_000011',
            'chr12': 'NC_000012',
            'chr13': 'NC_000013',
            'chr14': 'NC_000014',
            'chr15': 'NC_000015',
            'chr16': 'NC_000016',
            'chr17': 'NC_000017',
            'chr18': 'NC_000018',
            'chr19': 'NC_000019',
            'chr20': 'NC_000020',
            'chr21': 'NC_000021',
            'chr22': 'NC_000022',
            'chrX':  'NC_000023',
            'chrY':  'NC_000024',
            'chrMT': 'NC_012920'}[chrom]
        for record in SeqIO.parse(self.genome, 'fasta'):
            if record.id.startswith(chrom):
                for start, end in zip(*coordinates):
                    seq += str(record.seq[start-1:end])
        seq = seq.upper().replace('T', 'U')
        if strand == '-':
            seq = ''.join(['AUCG'['UAGC'.index(nt)] for nt in seq[::-1]])
        return seq
