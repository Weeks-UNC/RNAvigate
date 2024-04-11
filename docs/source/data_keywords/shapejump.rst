``shapejump``
~~~~~~~~~~~~~

ShapeJump inter-nucleotide RT deletion events

`ShapeJumper software`_

.. _ShapeJumper software: https://github.com/Weeks-UNC/ShapeJumper_V1

example uses:

- same as `interactions`_

input explaination:

- Input should be the deletions.txt output file from ShapeJumper.

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      shapejump="mydata_deletions.txt",
      )

other optional inputs:

- ``"sequence"`` is used to specify the sequence, and accepts the same inputs
  as the `sequence`_ keyword
- ``"metric"``, ``"metric_defaults"``, ``"read_table_kw"``, and ``"window"`` are
  explained in :doc:`/guides/custom_interactions`, but are generally not recommended
  for ShapeJumper files.

typical optional argument example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      shapejump={
         "shapejump": "mydata_deletions.txt",
         "sequence": "my_rna.fa",
      }
   )

back to `Standard data keywords`_
