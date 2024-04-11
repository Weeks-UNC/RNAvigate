``rnpmap``
~~~~~~~~~~

RNP-MaP per-nucleotide reactivities.

`RNPMapper software`_

.. _RNPMapper software: https://github.com/Weeks-UNC/RNP-MaP

example uses:

- same as `profile`_

input explaination:

- Input should be the output csv file from RNPMapper

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      rnpmap="myrnpmap_output.csv",
      )

other optional inputs:

- ``"metric"``, ``"metric_defaults"``, ``"sequence"``, and ``"read_table_kw"`` are
  explained in :doc:`/guides/custom_interactions`, but are not recommended for standard
  RNPMapper files.

back to `standard data keywords`_
