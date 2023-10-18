from rnavigate.rnavigate import *
from rnavigate.plotting_functions import *
from rnavigate import analysis
from rnavigate import data
from rnavigate import plots
from rnavigate import styles


standard_data_keywords = """
The following is a list of standard data keywords, what data type they
return and example input formats.

Some keywords require an external sequence, because the input data or file
format does not provide one.

ringmap, pairmap, allcorrs, shapejump, pairprob
    By default, the sequence for these data will be taken from the default
    profile, which is the first data keyword that returns per-nucleotide
    data. If this is not the desired behavior, a sequence can be specified.
motif, orfs, spans, sites, group, primers, domains
    These require a sequence to be specified.
pdb
    This keyword requires a sequence to be specified only if this
    information cannot be found in the header of the pdb or cif file.

To specify a sequence for any of these keywords, use the following format:
    keyword={"keyword": normal input,
                "sequence": same input format as fasta keyword}

Data keywords:

fasta
    data type: rnav.data.Sequence
    input explaination:
        Input should be a fasta file, a sequence string, or another data
        keyword. If another data keyword is provided, the sequence from
        that data is retrieved.
    example inputs:
        fasta file:
            fasta="path/to/my_sequence.fa"
            arbitrary={"fasta": "path/to/my_sequence.fa"}
        sequence string:
            fasta="AUCAGCGCUAUGACUGCGAUGACUGA"
            arbitrary={"fasta": "AUCAGCGCUAUGACUGCGAUGACUGA"}
        data keyword:
            fasta="another_data_keyword"
            arbitrary={"fasta": "another_data_keyword"}
log
    data type: rnav.data.Log
    input explaination:
        Input should be a ShapeMapper2 log.txt file. If ShapeMapper2 was
        run using the --per-read-histograms flag, this file contains
        information about fragment length distributions and per-read
        mutation counts distribution.
    example inputs:
        log="myshapemap_log.txt"
        arbitrary={"log": "myshapemap_log.txt"}
shapemap
    data type: rnav.data.SHAPEMaP
    input explaination:
        Input should be a ShapeMapper2 profile.txt file. This file contains
        the most complete per-nucleotide data from a ShapeMapper2 run.
    example inputs:
        shapemap="myshapemap_profile.txt"
        arbitrary={"shapemap": "myshapemap_profile.txt"}
dmsmap
    data type: rnav.data.SHAPEMaP
    input explaination:
        Input should be the same as for shapemap keyword. RNAvigate will
        automatically renormalize the profile according to the convention
        described here: https://doi.org/10.1073/pnas.1905491116
    example inputs:
        dmsmap="mydmsmap_profile.txt"
        arbitrary={"dmsmap": "mydmsmap_profile.txt"}
dancemap
    data type: rnav.data.DanceMaP
    input explaination:
        Input should be a dictionary containing:
            "dancemap": the DanceMapper reactivities.txt file
            "component": which component of the DANCE model to load
        Each component will need to be loaded seperately.
    example inputs:
        dancemap={"dancemap": "mydancemap_reactivities.txt",
                    "component": 0}
        arbitrary={"dancemap": "mydancemap_reactivities.txt",
                    "component": 0}
rnpmap
    data type: rnav.data.RNPMaP
    input explaination:
        Input should be the output csv file from RNPMapper
    example inputs:
        rnpmap="myrnpmap_output.csv"
        arbitrary={"rnpmap": "myrnpmap_output.csv"}
ringmap
    data type: rnav.data.RINGMaP
    input explaination:
        Input should be the correlations file output from RingMapper.
        See note above about specifying a sequence.
    example inputs:
        ringmap="myringmap_corrs.txt"
        arbitrary={"ringmap": "myringmap_corrs.txt"}
pairmap
    data type: rnav.data.PAIRMaP
    input explaination:
        Input should be the pairmap.txt output file from PairMapper.
        See note above about specifying a sequence.
    example inputs:
        pairmap="mydata_pairmap.txt"
        arbitrary={"pairmap": "mydata_pairmap.txt"}
allcorrs
    data type: rnav.data.RINGMaP
    input explaination:
        Input should be the allcorrs.txt output file from Pairmapper.
        See note above about specifying a sequence.
    example inputs:
        allcorrs="mydata_allcorrs.txt"
        arbitrary={"allcorrs": "mydata_allcorrs.txt"}
shapejump
    data type: rnav.data.SHAPEJuMP
    input explaination:
        Input should be the deletions.txt output file from ShapeJumper
        See note above about specifying a sequence.
    example inputs:
        shapejump="mydata_deletions.txt"
        arbitrary={"shapejump": "mydata_deletions.txt"}
pairprob
    data type: rnav.data.PairingProbability
    input explaination:
        Input should be the text format output file from running RNAstructure.
        partition followed by ProbabilityPlot -t, e.g.:
            partition my_sequence.fa pair_probabilities.dp
            ProbabilityPlot -t pair_probabilities.dp pair_probabilities.txt
        See note above about specifying a sequence.
    example inputs:
        pairprob={"pairprob": "pair_probabilities.txt",
                  "sequence": "my_sequence.fa"}
        arbitrary={"pairprob": "pair_probabilities.txt",
                  "sequence": "my_sequence.fa"}
ss
    data type: rnav.data.SecondaryStructure
    input explaination:
        Input should be one of the following secondary structure formats:
            connection table (.ct) or dotbracket notation (.dot, .dbn, etc.)
        OR one of the following secondary structure diagram formats:
            StructureEditor (.nsd or .cte), XRNA (.xrna), VARNA (.varna)
            FORNA (.forna or .json), R2DT (.r2dt or .json)
        The file format to parse will be determined by the file extension.
        Since FORNA and R2DT json specifications are different, the extension
        should be provided, e.g.:
            ss={"ss": "my_structure.json",
                "extension": "forna"}
    example inputs:
        ss="my_structure.ct"
        arbitrary={"ss": "my_structure.ct"}
pdb
    data type: rnav.data.PDB
    input explaination:
        Input should be a standard PDB file (.pdb or .cif).
        A chain ID must also be provided. If the file header is missing, a
        sequence must also be provided, see note about specifying a sequence.
    example inputs:
        pdb={"pdb": "my_structure.pdb",
             "chain": "X"}
        arbitrary={"pdb": "my_structure.pdb",
                   "chain": "X"}
allpossible
    data type: rnav.data.AllPossible
    input explaination:
        This data type has the same expected inputs as fasta keyword.
        The data returned represents all possible inter-nucleotide pairings.
        For very long RNAs this can be a large data structure (length squared)
    example inputs:
        allpossible=same expected input as fasta keyword
        arbitrary={"allpossible": same expected input as fasta keyword}

Annotations keywords:

Annotations are a little bit different from the above because they do not come
in standard file formats. Annotations inputs should all be dictionaries that
include the standard data keyword, "colors", "sequence", and "name":
    The input for each data keyword is different, and described below
    "sequence" has the same expected format as the fasta keyword
    "color" input should be a valid matplotlib color or RGB hexcode
        e.g. "grey", "green", "darkOrchid", "#f801ec"
        Web search for "Matplotlib specifying colors"
    "name" input can be any string. This will be used as a label on plots.

motif
    data type: rnav.data.Motif
    input explaination:
        input is a dictionary containing "motif", "sequence", "color", "name"
        "motif" input should be a string that uses the nucleotide alphabet
            A, U, C, G, T matches that nucleotide
            B = not A, D = not C, H = not G, V = not U or T
            W = weak (A, U, T), S = strong (C, G)
            M = amino (A, C), K = ketone (G, U, T)
            R = purine (A, G), Y = pyrimidine (C, U, T)
            N = any (A, U, C, G, T)
            e.g.: "DRACH" for potential m6A modification sites
        "sequence", "color" and "name" are described above
    example inputs:
        motif={"motif": "DRACH",
               "sequence": same expected input as fasta keyword,
               "color": "blue",
               "name": "m6A motif"}
        arbitrary={"motif": "DRACH",
                   "sequence": same expected input as fasta keyword,
                   "color": "blue",
                   "name": "m6A motif"}
orfs
    data type: rnav.data.ORFs
    input explaination:
        input is a dictionary containing "orfs", "sequence", "color", "name"
        "orfs" should be either of:
            "all" annotate all open-reading frames
            "longest" annotate only the longest open reading frame
        "sequence", "color" and "name" are described above
    example inputs:
        orfs={"orfs": "longest",
              "sequence": same expected input as fasta keyword,
              "color": "green",
              "name": "Longest ORF"}
        arbitrary={"orfs": "longest",
                   "sequence": same expected input as fasta keyword,
                   "color": "green",
                   "name": "Longest ORF"}
spans
    data type: rnav.data.Annotation
    input explaination:
        input is a dictionary containing "spans", "sequence", "color", "name"
        "spans" should be list of lists. Each inner list should be a pair of
        integers specifying the start and end of a span (1-indexed, inclusive)
        "sequence", "color" and "name" are described above
    example inputs:
        spans={"spans": [[10, 13], [65, 72]],
              "sequence": same expected input as fasta keyword,
              "color": "purple",
              "name": "interesting sites"}
        arbitrary={"spans": [[10, 13], [65, 72]],
                   "sequence": same expected input as fasta keyword,
                   "color": "purple",
                   "name": "interesting regions"}
sites
    data type: rnav.data.Annotation
    input explaination:
        input is a dictionary containing "sites", "sequence", "color", "name"
        "sites" should be list positions (1-indexed) of sites of interest
        "sequence", "color" and "name" are described above
    example inputs:
        sites={"sites": [10, 13, 65, 72],
               "sequence": same expected input as fasta keyword,
               "color": "red",
               "name": "interesting sites"}
        arbitrary={"sites": [[10, 13], [65, 72]],
                   "sequence": same expected input as fasta keyword,
                   "color": "red",
                   "name": "interesting sites"}
group
    data type: rnav.data.Annotation
    input explaination:
        input is a dictionary containing "group", "sequence", "color", "name"
        "group" should be list positions (1-indexed) of related nucleotides
        "sequence", "color" and "name" are described above
    example inputs:
        group={"group": [10, 13, 65, 72],
              "sequence": same expected input as fasta keyword,
              "color": "orange",
              "name": "ligand-binding pocket"}
        arbitrary={"group": [10, 13, 65, 72],
                   "sequence": same expected input as fasta keyword,
                   "color": "orange",
                   "name": "ligand-binding pocket"}
primers
    data type: rnav.data.Annotation
    input explaination:
        input is a dictionary containing "primers", "sequence", "color", "name"
        "primers" should be list of lists. Each inner list should be a pair of
        integers specifying the start and end of a primer (1-indexed, inclusive)
        Unlike spans keyword, primers are directional, i.e. reverse primer
        start sites are higher than end sites.
        "sequence", "color" and "name" are described above
    example inputs:
        primers={"primers": [[1, 23], [305, 292]],
                 "sequence": same expected input as fasta keyword,
                 "color": "orange",
                 "name": "PCR primers"}
        arbitrary={"primers": [[1, 23], [305, 292]],
                   "sequence": same expected input as fasta keyword,
                   "color": "orange",
                   "name": "PCR primers"}
domains
    data type: rnav.data.Annotation
    input explaination:
        input is a dictionary containing "domains", "sequence", "colors", "names"
        "domains" should be list of lists. Each inner list should be a pair of
        integers specifying the start and end of a domain (1-indexed, inclusive)
        "sequence" is described above
        "colors" and "names" are similar to the other annotations, but should
        be lists specifying a color and name for each domain
    example inputs:
        domains={"domains": [[1, 50], [51, 200] [201, 305]],
                 "sequence": same expected input as fasta keyword,
                 "colors": ["orange", "blue", "green"],
                 "names": ["5' UTR", "CDS", "3' UTR"]}
        arbitrary={"domains": [[1, 50], [51, 200] [201, 305]],
                 "sequence": same expected input as fasta keyword,
                 "colors": ["orange", "blue", "green"],
                 "names": ["5' UTR", "CDS", "3' UTR"]}
"""
