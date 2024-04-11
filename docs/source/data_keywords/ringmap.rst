``ringmap``
~~~~~~~~~~~

single-molecule correlations from a DMS- or eDMS-MaP experiment

`RingMapper software`_

.. _RingMapper software: https://github.com/Weeks-UNC/RingMapper

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
   as the `sequence`_ keyword
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
