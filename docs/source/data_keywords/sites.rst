``sites``
~~~~~~~~~

Annotation of any sites of interest

example uses:

- same as :doc:`/data_keywords/motif`

input explaination:

- input is a dictionary containing

   - ``"sites"``: a list of nucleotide positions (1-indexed, inclusive)
   - ``"sequence"``: same as :doc:`/data_keywords/sequence` keyword
   - ``"color"``: a valid color or hexcode, e.g. ``"blue"``, ``"grey"``, or ``"#fa4ce2"``
   - ``"name"``: an arbitrary name to use on plots

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      m6a_sites={
         "sites": [10, 13, 65, 72],
         "sequence": "my_sequence.fa",
         "color": "purple",
         "name": "m6A sites"
         }
      )

back to :ref:`Standard data keywords <standard-data-keywords>`
