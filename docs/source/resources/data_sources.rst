Sources of data
===============

.. note::

   This document is a work in progress.


These methods, databases, and software provide RNA structure data that can be
used with RNAvigate.

To add a database, method, or software package for obtaining RNA structure data, submit an issue on [Github issues](https://github.com/Weeks-UNC/RNAvigate/issues). Include the name of the database, method or software, a citation, and a brief explanation.

Categories of RNA structure data:

- post-transcriptional modifications
- base-pairing
- through-space interactions
- solvent accessibility
- ligand binding
- ribosome binding
- protein binding
- 3D structure
- conservation and covariation
- folding landscapes
- multi-state ensembles

small molecule binding
----------------------

chem-CLIP-Map-Seq

Frag-MaP

ribosomal profiling
-------------------

Ribosome profiling https://cshperspectives.cshlp.org/content/early/2018/07/23/cshperspect.a032698

Intra-RNA and RNA-RNA interactions
----------------------------------

Note: RNAvigate is currently not well-suited for display or analysis of inter-molecular data

Proximity detection through chimeric reads

transcriptome-wide RNA-RNA inter- and intra-molecular base pairing

RPL, PARIS, SPLASH, LIGR-seq, MARIO, COMRADES (uses psoralen)
SHARC, SHAPE-JuMP (uses bifunctional acylation reagents)

single-molecule correlated chemical probing for intra-molecular structural communication

RING-MaP, PAIR-MaP (correlated adduct-formation)

secondary structure
-------------------

chemical probing of secondary structure, widely used to infer secondary structure

DMS, TMO, glyoxal, EDC (base-pairing status)
SHAPE reagents (backbone flexibility)
hydroxyl radicals, LASER (NAz) (solvent accessibility)

computational methods for secondary structure determination

RNAStructure, RNAfold (thermodynamics based)
Dynalign II, R-scape (covariation based)
RNAalifold, TurboFold II (thermodynamics and covariation based)
Pfold, ContraFold, ContextFold (machine learning based)
SPOT-RNA, MXfold2, CDPfold, DMfold, E2Efold, Ufold (deep learning based)

Protein binding
---------------

chemical probing for protein binding

footprinting SHAPE (fSHAPE), LASER/SHAPE (differential probing signal from protein binding)

RNP-MaP (per-nucleotide protein proximity)

enrichment quantification for protein binding

https://www.sciencedirect.com/science/article/pii/S0021925822003647

CLIP-seq/HITS-CLIP, PAR-CLIP, iCLIP, eCLIP, irCLIP, CLASH, vPAR-CL (CLIP-derivative methods)
RIP-seq	(IP without crosslinking)
APEX-seq/Proximity-CLIP (proximity biotinylation)
PIP-seq (Formaldehyde crosslink and RNase footprinting)

structure ensembles
-------------------

experimental ensemble deconvolution

SHAPE-seq, SPET-seq (co-transcriptional folding)
Mutate and map (M^2)

thermodynamic ensemble deconvolution

RSample, SLEQ, R2D2, REEFFIT, IRIS (based on fitting thermodynamic models to probing data)

smCCP-based ensemble deconvolution

DANCE-MaP, DRACO, DREEM (based on clustering of DMS-MaP sequencing reads)

Tertiary structure
------------------

https://www.nature.com/articles/s41592-022-01623-y

in-vitro experimental methods for determining 3D structure

X-ray crystalography, NMR, cryo-EM, SAXS

computational methods for tertiary structure determination

iFold, SimRNA (ab initio folding methods)
FARNA, MC-Sym, FARFAR2 (fragment assembly methods)
ARES (deep-learning method)

RNA post-transcriptional modifications
--------------------------------------

Ψ: Ψ-seq, Pseudo-seq, PSI-seq, CeU-seq, RBS-seq

m5C: Bisulfite sequencing, m5C-RIP, AZA-IP, m5C-miCLIP

m1A: m1A-meRIP, m1A-ID-seq, m1A-seq, m1A-MAP, m1a-Quant-seq

ac4C: acRIP, ac4C-seq

Nm: Nm-seq, RibOxi-seq, RiboMeth-seq, MeTH-seq

m7G: m7G-MAP-seq, m7G-(me)RIP seq, TRAC-seq, AlkAniline-seq, BoRed-seq, m7G-seq


Online databases
----------------

codon usage database https://www.kazusa.or.jp/codon/

The RNA families database https://rfam.org/

RNAcentral https://rnacentral.org/ (comprehensive database of non-coding RNA)

MODOMICS (RNA modifications, including their known locations in RNA sequences)

GENCODE high quality reference gene annotation and validation for human and mouse genomes

Ensembl is a genome browser for vertebrate genomes and model organisms

Online computational Tools
--------------------------

RBPmap http://rbpmap.technion.ac.il/ find RNA binding protein sequence motif matches in an RNA of interest
ViennaRNA Web Services https://rna.tbi.univie.ac.at/ (thermodynamic folding tools)

FORNA http://rna.tbi.univie.ac.at/forna/ (quickly produce an RNA structure diagram)

R2DT https://rnacentral.org/r2dt (create a secondary structure diagram in standardized orientations)

Clustal Omega http://www.ebi.ac.uk/Tools/msa/clustalo/ (multiple sequence alignments)
