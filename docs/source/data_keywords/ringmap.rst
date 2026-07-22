``ringmap``
~~~~~~~~~~~

single-molecule correlations from a DMS- or eDMS-MaP experiment

`RingMapper software`_

.. _RingMapper software: https://github.com/Weeks-UNC/RingMapper

A similar data keyword, ``ringmap_N7G``, accepts a ``concat_rings.txt`` file produced
by the msDMS-MaP pipeline. Positions involving N7 of guanosine (N7G) are decoded back
to their true nucleotide indices by subtracting the RNA length from any i/j values that
are greater than the RNA length, and a ``Type`` column is added to classify each RING
as ``"N1N1"``, ``"N7N1"``, ``"N1N7"``, or ``"N7N7"``, where the order indicates which of
``i`` (first) and ``j`` (second) is the N7G position. Use ``Type_eq`` and ``Type_ne``
filter kwargs when plotting to separate backbone (N1/N1) and N7G-involving correlations.

example uses:

- same as `interactions`_

input explaination:

- Input should be the correlations file output from RingMapper.

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      ringmap="myringmap_corrs.txt",
      )

other optional inputs:

- ``"sequence"`` is used to specify the sequence, and accepts the same inputs
   as the `sequence`_ keyword. Defaults to ``"default_profile"``, the first profile
   added to the RNAvigate Sample.
- ``"metric"``, ``"metric_defaults"``, ``"read_table_kw"``, and ``"window"`` are
   explained in :doc:`/guides/custom_interactions`, but are generally not recommended
   for RingMapper files.

typical optional input example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      ringmap={
         "ringmap": "myringmap_corrs.txt"
         "sequence": "my_rna.fa"
         }
      )

back to `Standard data keywords`_
