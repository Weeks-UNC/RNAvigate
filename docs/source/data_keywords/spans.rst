``spans``
~~~~~~~~~

Annotation of any regions of interest

example uses:

- same as `motif`_

input explaination:

- input is a dictionary containing

   - ``"spans"``: a list of lists of 2 integers. Each inner list specifies a
      start and end position of a span (1-indexed, inclusive)
   - ``"sequence"``: same as `sequence`_ keyword
   - ``"color"``: a valid color or hexcode, e.g. ``"blue"``, ``"grey"``, or ``"#fa4ce2"``
   - ``"name"``: an arbitrary name to use on plots

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      regions={
         "spans": [[10, 13], [65, 72]],
         "sequence": "my_sequence.fa",
         "color": "purple",
         "name": "interesting regions"
         }
      )

back to `Standard data keywords`_
