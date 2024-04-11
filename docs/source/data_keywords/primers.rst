``primers``
~~~~~~~~~~~

Annotation of primer binding sites

example uses:

- same as `motif`_

input explaination:

- input is a dictionary containing

   - ``"primers"``: a list of lists of 2 integers. Each inner list specifies a
      start and end position of a primer (1-indexed, inclusive). A reverse primer
      is specified by listing the 3' -> 5' start and end, e.g. [300, 278]
   - ``"sequence"``: same as `sequence`_ keyword
   - ``"color"``: a valid color or hexcode, e.g. ``"blue"``, ``"grey"``, or ``"#fa4ce2"``
   - ``"name"``: an arbitrary name to use on plots

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      pcr_primers={
         "primers": [[1, 22], [300, 278]],
         "sequence": "my_sequence.fa",
         "color": "purple",
         "name": "primer-binding sites"
         }
      )

back to `Standard data keywords`_
