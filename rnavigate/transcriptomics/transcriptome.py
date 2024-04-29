"""This submodule defines the Transcriptome and Transcript classes.

Transcriptome objects require genome fasta and annotation gtf files. Then, if
provided with a transcript ID, it will return a Transcript object.

Transcript objects contain a transcript sequence and genome coordinates.
It can return annotations for CDS, UTR, and exon junctions.
It can be used with Bed objects to convert genome coordinate data to transcript
coordinate RNAvigate data classes.
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
from rnavigate import data
from pathlib import Path


human_chromosome_ids = {
    "chr1": "NC_000001",
    "chr2": "NC_000002",
    "chr3": "NC_000003",
    "chr4": "NC_000004",
    "chr5": "NC_000005",
    "chr6": "NC_000006",
    "chr7": "NC_000007",
    "chr8": "NC_000008",
    "chr9": "NC_000009",
    "chr10": "NC_000010",
    "chr11": "NC_000011",
    "chr12": "NC_000012",
    "chr13": "NC_000013",
    "chr14": "NC_000014",
    "chr15": "NC_000015",
    "chr16": "NC_000016",
    "chr17": "NC_000017",
    "chr18": "NC_000018",
    "chr19": "NC_000019",
    "chr20": "NC_000020",
    "chr21": "NC_000021",
    "chr22": "NC_000022",
    "chrX": "NC_000023",
    "chrY": "NC_000024",
    "chrMT": "NC_012920",
}


class Transcript(data.Sequence):
    """Transcript object for a single transcript.

    Parameters
    ----------
    parent : Transcriptome
        Parent Transcriptome object
    name : str
        Transcript ID
    sequence : str
        Transcript sequence
    chromosome : str
        Chromosome ID
    strand : str
        Strand of the transcript
    coordinates : tuple
        Tuple of two lists of genome coordinates for the transcript, e.g.:
        [(start1, start2, ...), (stop1, stop2, ...)]
    tx_info : dict
        Dictionary of transcript information from the GTF file
    cds_coors : list
        List of genome coordinates for the CDS
    other_features : list
        List of dictionaries of other features from the GTF file
    """

    def __init__(
        self,
        parent,
        name,
        sequence,
        chromosome,
        strand,
        coordinates,
        tx_info,
        cds_coors=None,
        other_features=None,
    ):
        """Initialize the Transcript object."""
        super().__init__(sequence)
        self.name = name
        self.parent = parent
        self.chromosome = chromosome
        self.strand = strand
        self.coordinates = coordinates
        self.tx_info = tx_info
        self.cds_coors = cds_coors
        self.other_features = other_features
        self.coordinate_df = self.get_coordinate_df()

    def get_coordinate_df(self):
        """Return a DataFrame of transcript coordinates."""
        coordinates = [np.arange(x, y + 1) for x, y in zip(*self.coordinates)]
        coordinates = np.concatenate(coordinates)
        if self.strand == "-":
            coordinates = coordinates[::-1]
        df = pd.DataFrame(
            {
                "Sequence": list(self.sequence),
                "Nucleotide": np.arange(self.length) + 1,
                "Coordinate": coordinates,
            }
        )
        return df

    def get_tx_coordinate(self, coordinate):
        """Return the transcript coordinate for a genome coordinate."""
        df = self.coordinate_df
        coordinate = df["Coordinate"] == coordinate
        if sum(coordinate) == 0:
            return None
        return df.loc[coordinate, "Nucleotide"].values[0]

    def get_tx_range(self, start, stop):
        """Return the transcript coordinates for a genome range."""
        df = self.coordinate_df
        nts = df.loc[df["Coordinate"].between(start, stop), "Nucleotide"]
        return [nts.min(), nts.max()]

    def get_cds_annotation(self, **kwargs):
        """Return an Annotation object for the CDS.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the Annotation object
        """
        cds = self.get_tx_range(*sorted(self.cds_coors))
        return data.Annotation(
            input_data=[cds],
            annotation_type="spans",
            sequence=self.sequence,
            name="CDS",
            **kwargs,
        )

    def get_cds_domains(self):
        """Return a Domains object for the 5' UTR, CDS and 3' UTR."""
        cds = self.get_tx_range(*sorted(self.cds_coors))
        spans = [[1, cds[0] - 1], cds, [cds[1] + 1, self.length]]
        data.domains(
            input_data=spans,
            names=["5' UTR", "CDS", "3' UTR"],
            colors=["skyblue", "palegreen", "orchid"],
            sequence=self.sequence,
        )

    def get_junctions_annotation(self, **kwargs):
        """Return an Annotation object for the exon junctions.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the Annotation object
        """
        junctions = [self.coordinates[0][1:], self.coordinates[1][:-1]]
        junctions = [
            self.get_tx_range(*sorted(junction)) for junction in zip(*junctions)
        ]
        return data.Annotation(
            input_data=junctions,
            annotation_type="spans",
            sequence=self.sequence,
            name="exon junctions",
            **kwargs,
        )

    def get_exon_annotation(self, exon_number, **kwargs):
        """Return an Annotation object for a single exon.

        Parameters
        ----------
        exon_number : int
            Exon number
        **kwargs
            Additional keyword arguments for the Annotation object
        """
        if self.strand == "-":
            exon_number = len(self.coordinates[0]) + 1 - exon_number
        try:
            exon = [
                self.coordinates[0][exon_number - 1],
                self.coordinates[1][exon_number - 1],
            ]
        except IndexError as e:
            raise ValueError(f"exon number {exon_number} is too high") from e
        exon = self.get_tx_range(*sorted(exon))
        return data.Annotation(
            input_data=[exon],
            annotation_type="spans",
            sequence=self.sequence,
            name=f"exon {exon_number}",
            **kwargs,
        )

    def get_exon_domains(self):
        """Return a Domains object for the exons."""
        exons = [self.get_tx_range(*sorted(e)) for e in zip(*self.coordinates)]
        data.domains(
            input_data=exons,
            names=[f"exon {i}" for i in range(len(exons))],
            colors=[["lightsalmon", "gold"][i % 2] for i in range(len(exons))],
            sequence=self.sequence,
        )


class Transcriptome:
    """Transcriptome object for a genome and annotation file.

    Parameters
    ----------
    genome : str
        Path to the genome fasta file
    annotation : str
        Path to the annotation gtf file
    path : str or Path
        Path to the genome and annotation files
    chr_ids : dict
        Dictionary of chromosome IDs
    """

    def __init__(self, genome, annotation, path, chr_ids=None):
        """Initialize the Transcriptome object."""
        path = Path(path)
        if chr_ids is None:
            chr_ids = human_chromosome_ids
        self.genome = path / genome
        self.annotation = path / annotation
        self.chr_ids = chr_ids

    def get_transcript(self, transcript_id):
        """Return a Transcript object for a single transcript ID."""
        return self.get_transcripts([transcript_id])[transcript_id]

    def get_transcripts(self, transcript_ids):
        """Return a dictionary of Transcript objects for a list of transcript IDs."""

        def parse_entry(entry):
            entry = entry.strip().split("\t")
            if entry[2] == "transcript":
                info = entry[8].split(";")
                info = [item.strip().split(" ") for item in info[:-1]]
                info = {k: v.strip('"') for k, v in info}
            else:
                info = {}
            return {
                "chromosome": entry[0],
                "db": entry[1],
                "feature": entry[2],
                "start": entry[3],
                "stop": entry[4],
                "strand": entry[6],
                "info": info,
            }

        with open(self.annotation) as txome:
            tx_info, chromosome, strand = {}, {}, {}
            start_coors = {tx_id: [] for tx_id in transcript_ids}
            stop_coors = {tx_id: [] for tx_id in transcript_ids}
            cds_coors = {tx_id: [] for tx_id in transcript_ids}
            other_features = {tx_id: [] for tx_id in transcript_ids}
            for entry in txome.readlines():
                tx_id = None
                for this_id in transcript_ids:
                    if f'transcript_id "{this_id}"' in entry:
                        tx_id = this_id
                if tx_id is None:
                    continue
                entry = parse_entry(entry)
                if entry["feature"] == "transcript":
                    tx_info[tx_id] = entry["info"]
                    chromosome[tx_id] = entry["chromosome"]
                    strand[tx_id] = entry["strand"]
                elif entry["feature"] == "exon":
                    start_coors[tx_id].append(int(entry["start"]))
                    stop_coors[tx_id].append(int(entry["stop"]))
                elif entry["feature"] in ["start_codon", "stop_codon"]:
                    cds_coors[tx_id].append(int(entry["start"]))
                    cds_coors[tx_id].append(int(entry["stop"]))
                elif entry["feature"] in ["UTR", "CDS"]:
                    continue
                else:
                    other_features[tx_id].append(
                        {key: entry[key] for key in ["feature", "start", "stop"]}
                    )
        transcripts = {}
        coordinates = {
            id: (sorted(start_coors[id]), sorted(stop_coors[id]))
            for id in transcript_ids
        }
        sequences = self.get_sequences(chromosome, coordinates, strand)
        for tx_id in transcript_ids:
            cds_coors[tx_id] = [min(cds_coors[tx_id]), max(cds_coors[tx_id])]
            transcripts[tx_id] = Transcript(
                parent=self,
                name=tx_id,
                sequence=sequences[tx_id],
                chromosome=chromosome[tx_id],
                strand=strand[tx_id],
                coordinates=coordinates[tx_id],
                tx_info=tx_info[tx_id],
                cds_coors=cds_coors[tx_id],
                other_features=other_features[tx_id],
            )
        return transcripts

    def get_sequences(self, chromosomes, coordinates, strands):
        """Return a dictionary of transcript sequences."""
        sequences = {tx_id: "" for tx_id in chromosomes}
        for record in SeqIO.parse(self.genome, "fasta"):
            for tx_id, chromosome in chromosomes.items():
                if record.id.startswith(self.chr_ids[chromosome]):
                    for start, end in zip(*coordinates[tx_id]):
                        sequences[tx_id] += str(record.seq[start - 1 : end])
        for tx_id, seq in sequences.items():
            seq = seq.upper().replace("T", "U")
            # reverse compliment if on minus strand
            if strands[tx_id] == "-":
                seq = "".join(["AUCG"["UAGC".index(nt)] for nt in seq[::-1]])
            sequences[tx_id] = seq
        return sequences

    def get_sequence(self, chromosome, coordinates, strand):
        """Return a transcript sequence for a single transcript."""
        seq = ""
        chromosome = self.chr_ids[chromosome]
        for record in SeqIO.parse(self.genome, "fasta"):
            if record.id.startswith(chromosome):
                for start, end in zip(*coordinates):
                    seq += str(record.seq[start - 1 : end])
        seq = seq.upper().replace("T", "U")
        # reverse compliment if on minus strand
        if strand == "-":
            seq = "".join(["AUCG"["UAGC".index(nt)] for nt in seq[::-1]])
        return seq
