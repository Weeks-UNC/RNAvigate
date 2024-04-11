``shapemap``
~~~~~~~~~~~~

SHAPE, DMS, or other reagent per-nucleotide reactivities

- `ShapeMapper2 software`_

.. _ShapeMapper2 software: https://github.com/Weeks-UNC/shapemapper2

Two similar data keywords:
- ``dmsmap`` applies DMS-MaP normalization to profile when loaded

   - `PAIR-MaP paper`_

- ``shapemap_rnaframework`` accepts an RNAframework xml file.

   - RNAframework files do not contain error estimates or raw data.
   - `RNAframework software`_

.. _RNAframework software: https://github.com/dincarnato/RNAFramework
.. _PAIR-MaP paper: https://doi.org/10.1073/pnas.1905491116

.. code-block:: python

   dmsmap="path/to/shapemap_profile.txt"

is equivalent to

.. code-block:: python

   dmsmap={'shapemap': "shapemap_profile.txt", "normalize": "DMS"}

example uses:

- same as `profile`_
- plus: visualizing quality control metrics if a log file is specified

input explaination:

- Input should be a ShapeMapper2 profile.txt file. This file contains the most
   complete per-nucleotide data from a ShapeMapper2 run.

example inputs:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      shapemap="shapemap_profile.txt"
   )

other optional inputs:

- ``"normalize"``:

   - Defaults to not performing any renormalization.
   - ``"DMS"`` will perform DMS-MaP renormalization
   - ``"eDMS"`` will perform eDMS-MaP renormalization
   - ``"boxplot"`` will perform ShapeMapper2 renormalization (with 1 improvement)
   - By default, renormalization is performed on the HQ_profile and HQ_stderr
      columns, and overwrites the Norm_profile and Norm_stderr columns
   - Normalization can also be done after ``rnav.Sample`` creation.

    - type: ``help(rnav.data.Profile.normalize)``

- ``"log"`` is used to specify a log file. If ShapeMapper2 was run with the
   ``--per-read-histograms`` flag, this file will contain read length distribution
   and mutations-per-read distribution. This data can then be visualized with
   ``rnav.plot_QC``.
- ``"metric"``, ``"metric_defaults"``, ``"sequence"``, and ``"read_table_kw"`` are
   explained in :doc:`/guides/custom_interactions`, but are not recommended for standard
   ShapeMapper2 files.

typical optional input example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      shapemap={
         "shapemap": "shapemap_profile.txt",
         "log": "shapemap_log.txt",
      }
   )

back to `Standard data keywords`_
