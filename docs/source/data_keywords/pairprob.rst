``pairprob``
~~~~~~~~~~~~

inter-nucleotide predicted pairing probabilities

`RNAStructure software`_

.. _RNAStructure software: https://rna.urmc.rochester.edu/RNAstructure.html

example uses:

- same as `interactions`_
- plus: calculating shannon entropy

input explaination:

- Input should be a dotplot plain text file from running RNAstructure
   ``partition`` followed by ``ProbabilityPlot`` with ``-t`` option

For example:

.. code-block:: bash

   partition my_sequence.fa pair_probabilities.dp
   ProbabilityPlot -t pair_probabilities.dp pair_probabilities.txt


example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      pairprob={"pairprob": "pair_probabilities.txt",
               "sequence": "my_sequence.fa"}
      )

other optional inputs:

- ``"sequence"`` is used to specify the sequence, and accepts the same inputs
  as the `sequence`_ keyword
- ``"metric"``, ``"metric_defaults"``, ``"read_table_kw"``, and ``"window"`` are
  explained in :doc:`/guides/custom_interactions`, but are generally not recommended
  for "" files

typical optional argument example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      pairprob={
         "pairprob": "pair_probabilities.txt",
         "sequence": "my_rna.fa",
      }
   )

back to `Standard data keywords`_
