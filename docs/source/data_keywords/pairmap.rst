``pairmap``
~~~~~~~~~~~

single-molecule correlations from a DMS- or eDMS-MaP experiment reflective of
base pairing

`PairMapper software`_ (part of RingMapper)

.. _PairMapper software: https://github.com/Weeks-UNC/RingMapper

example uses:

- same as `interactions`_

input explaination:

- Input should be the pairmap.txt output file from PairMapper.

example inputs:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      pairmap="mydata_pairmap.txt",
      )

other optional inputs:

- ``"sequence"`` is used to specify the sequence, and accepts the same inputs
  as the `sequence`_ keyword
- ``"metric"``, ``"metric_defaults"``, ``"read_table_kw"``, and ``"window"`` are
  explained in :doc:`/guides/custom_interactions`, but are generally not recommended
  for PairMapper files

typtical optional input example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      pairmap={
         "pairmap": "mydata_pairmap.txt",
         "sequence": "my_rna.fa",
      }
   )

back to `Standard data keywords`_
