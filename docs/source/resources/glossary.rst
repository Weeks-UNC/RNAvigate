Glossary
========

.. glossary::

   Sample
      a python object created with `rnavigate.Sample`

   Data keyword
      a string that points to data within a sample

   Data object
      an rnavigate object representation of data

   contact distance
      the shortest path distance between 2 nucleotides in the secondary structure graph

   sequence
      This refers to an RNA sequence. As an argument, the value can be a sequence
      string (e.g.: ``"AUGCGUCGGUGUCUACUGA"``), a path to a fasta file containing only
      1 sequence (e.g.: ``"./data/my_sequence.fasta"``) any rnavigate data object, or a
      pandas.DataFrame which contains a "Sequence" column with the nucleotide identity
      sorted by position

   interactions
   interactions2
      An rnavigate.data.Interactions object

   structure
   structure2
      An rnavigate.data.SecondaryStructure object

   profile
      An rnavigate.data.Profile object

   index
   indices
   idx
      python array indices (0-indexed)

   position
   positions
   pos
      nucleotide positions (1-indexed, inclusive)
