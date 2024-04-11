``sequence``
~~~~~~~~~~~~

an RNA sequence

example uses:

- aligning data between sequences
- all data in RNAvigate is associated with a sequence and can be aligned to
  other data, or vice versa.

input explaination:

- Input should be a fasta file, a sequence string, or another data keyword. If
  another data keyword is provided, the sequence from that data is retrieved.

example inputs:

.. code-block:: python

   # fasta file
   my_sample = rnav.Sample(
      name="example",
      sequence="path/to/my_sequence.fa",
   )

   # sequence string
   my_sample = rnav.Sample(
      name="example",
      sequence="AUCAGCGCUAUGACUGCGAUGACUGA",
   )

   # data keyword
   my_sample = rnav.Sample(
      name="example",
      data_keyword="some_data_with_a_sequence"
      sequence="data_keyword",
   )

back to `Standard data keywords`_
