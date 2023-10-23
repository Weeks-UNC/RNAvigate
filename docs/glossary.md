Glossary
========

Documentation terms
-------------------

* Sample: a python object created with `rnavigate.Sample`
* Data keyword: a string that points to data within a sample
* Data object: an rnavigate object representation of data
    - `my_sample.get_data("data_keyword")` returns a data object
* contact distance: the shortest path distance between 2 nucleotides in the
  secondary structure graph

Python arguments
----------------

* sequence: this variable can always be any of:
    1. sequence string: `"AUGCGUCGGUGUCUACUGA"`
    2. fasta filepath (1 sequence): `"./data/my_sequence.fasta"`
    3. any rnavigate data object
    4. a pandas.DataFrame which is:
        - 1 row per position
        - sorted by position
        - contains a "Sequence" column of single characters
* interactions(2) - rnavigate.Interactions data object
* structure(2) - rnavigate.SecondaryStructure data object
* profile - rnavigate.Profile data object

Python variables
----------------

* index, idx, or indices: python array indices (0-indexed)
* position, pos, or positions: nucleotide positions (1-indexed, inclusive)
