``allpossible``
~~~~~~~~~~~~~~~

All possible inter-nucleotide pairings for a given sequence

example uses:

- same as `interactions`_
- plus: calculating the expected distance distribution of a filtering scheme

input explanation:

- This keyword has the same expected inputs as the `sequence`_ keyword.
- Note: the size of the data increases with the sequence length squared.

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      allpossible="my_rna.fa"
   )

other optional inputs:

- ``"window"`` is used to specify the window size of the interacting regions.

   - Default: ``"window": 1`` means nucleotide *i* to nucleotide *j*
   - ``"window": 3`` means nucleotides i:i+3 to nucleotides j:j+3

- ``"metric"``, ``"metric_defaults"``, and ``"read_table_kw"`` are explained
   in :doc:`/guides/custom_interactions`.

typical optional argument example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      allpossible={
         "allpossible": "my_rna.fa",
         "window": 3,
      }
   )

back to `Standard data keywords`_
