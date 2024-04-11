``domains``
~~~~~~~~~~~

Annotation of RNA domains

example uses:

- same as `motif`_
- plus: labelling domains across the x-axis of skyline, profile, and arc plots

input explaination:

- input is a dictionary containing

   - ``"domains"``: a list of lists of 2 integers. Each inner list specifies a
      start and end position of a primer (1-indexed, inclusive).
   - ``"sequence"``: same as `sequence`_ keyword
   - ``"colors"``: a valid color or hexcode, e.g. ``"blue"``, ``"grey"``, or ``"#fa4ce2"``
   - ``"names"``: an arbitrary name to use on plots

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      mrna_domains={
         "domains": [[1, 62], [63, 205], [206,300]],
         "sequence": "my_rna.fa",
         "colors": ["purple", "green", "orange"],
         "names": ["5'UTR", "CDS", "3'UTR"],
         }
      )

back to `Standard data keywords`_
