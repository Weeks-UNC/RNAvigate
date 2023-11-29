Data types analyzed by RNAvigate
================================

RNA structure data can come from databases, experimental measurements, and
computational analyses or predictions. RNAvigate will accept text files or
python objects (usually a pandas.DataFrame) as inputs.

Click [here](dev/data-sources.md) for methods, databases, and software that
supply RNA structure information.

RNA structure data can be analyzed with RNAvigate if it falls into one of these
five categories. More information on these categories is provided below.

- [Annotations](#annotations)
- [Per-nucleotide measurements](#profiles)
- [Secondary structures](#secondary-structures)
- [Inter-nucleotide measurements](#interactions)
- [3D Structures](#3d-structures)

Annotations
-----------

RNAvigate data class: `rnav.data.Annotations`

Annotations define a set of related RNA features. These can be
individual nucleotides, regions, or discontinuous groups.

Types of features:

- individual sites:
    - modified nucleotides
        - m6A, Pseudouridine, m7G, and others
    - nucleotides with measurements above a threshold value
    - SNPs and riboSNitches
- regions:
    - protein binding sites, such as eCLIP peaks
    - UTRs, ORFs and codons
    - introns and exons
    - sequence motifs
    - primer binding sites
    - structural features
        - pseudoknots, kissing loops, G-quadruplexes and others
- discontinuous groups:
    - nucleotides proximal to a protein, ligand, or pocket

Profiles
--------

RNAvigate data class: `rnav.data.Profile`

Profiles define per-nucleotide measurements.

Types of measurements:
- chemical reactivity, such as mutational profiles (MaP) and RT stop signals
- read counts or enrichment scores, such as CLIP or RIP based methods
- sequence conservation
- shannon entropy and pairing probability
- structural proximity to some other feature or molecule

Secondary structures
--------------------

RNAvigate data class: `rnav.data.SecondaryStructure`

Secondary structures define a pattern of base-pairing. Additionally, these may
contain secondary structure diagram layout coordinates for each nucleotide.

Types of secondary structures:
- Experimentally determined from CryoEM or crystal structures
- computationally modeled de novo or informed by chemical probing data
- Secondary structure drawing layouts, such as from VARNA, XRNA, R2DT, etc.

Interactions
------------

RNAvigate data class: `rnav.data.Interactions`

Interactions define inter-nucleotide measurements. These can be between
individual nucleotides or uniform windows of nucleotides.

Types of interactions:
- Single molecule correlated events
- Interactions data from proximity ligation, SHARC, SHAPE-JuMP, etc.
- Base-pairing probabilities
- Sequence covariation

3D structures
-------------

RNAvigate data class: `rnav.data.PDB`

3D structures define the atomic coordinates of residues in an RNA.

Types of 3D structures:
- Experimentally determined from CryoEM or crystal structures
- computationally modeled de novo or informed by chemical probing data
