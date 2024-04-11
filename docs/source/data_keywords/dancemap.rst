``dancemap``
~~~~~~~~~~~~

Reactivity profile of a single component of a DanceMapper model

- `DanceMapper software`_
- `DANCE-MaP paper`_

.. _DanceMapper software: https://github.com/MustoeLab/DanceMapper
.. _DANCE-MaP paper: https://doi.org/10.1016/j.molcel.2022.02.009

example uses:

- same as `profile`_

input explaination:

- Input should be a dictionary containing:

  - ``"dancemap"``: the DanceMapper reactivities.txt file
  - ``"component"```: which component of the DANCE model to load

- This works best if each component is a seperate ``rnav.Sample``

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      dancemap={
         "dancemap": "mydancemap_reactivities.txt",
         "component": 0},
      )
   my_sample_2 = rnav.Sample(
      sample="example2",
      dancemap={
         "dancemap": "mydancemap_reactivities.txt",
         "component": 1},
      )

other optional inputs:

- ``"metric"``, ``"metric_defaults"``, ``"sequence"``, and ``"read_table_kw"`` are
  explained in :doc:`/guides/custom_profiles`, but are not recommended for standard
  DanceMapper files.

back to `Standard data keywords`_
