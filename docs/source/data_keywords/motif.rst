``motif``
~~~~~~~~~

Annotation of occurances of a sequence motif

example uses:

- highlighting nucleotides in skyline, profile, arc, circle, or secondary
  structure diagram plots
- coloring nucleotides in circle plots, secondary structure diagrams, 3D
  molecule renderings, or linear regression scatter plots

input explaination:

- Input should be a dictionary containing:

   - ``"motif"``: a string that uses the nucleotide alphabet

      - e.g.: "DRACH" for potential m6A modification sites

+------------+------------+------------+
| alphabet   | meaning    | matches    |
+============+============+============+
| A, U, C, G | identity   | A, U, C, G |
+------------+------------+------------+
| B          | not A      | U/C/G      |
+------------+------------+------------+
| D          | not C      | A/U/G      |
+------------+------------+------------+
| H          | not G      | A/U/C      |
+------------+------------+------------+
| V          | not U      | A/C/G      |
+------------+------------+------------+
| W          | weak       | A/U        |
+------------+------------+------------+
| S          | strong     | C/G        |
+------------+------------+------------+
| M          | amino      | A/C        |
+------------+------------+------------+
| K          | ketone     | U/G        |
+------------+------------+------------+
| R          | purine     | A/G        |
+------------+------------+------------+
| Y          | pyrimidine | U/C        |
+------------+------------+------------+
| N          | any        | A/U/C/G    |
+------------+------------+------------+

   - ``"sequence"``: same as `sequence`_ keyword
   - ``"color"``: a valid color or hexcode, e.g. ``"blue"``, ``"grey"``, or ``"#fa4ce2"``
   - ``"name"``: an arbitrary name to use on plots

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      m6A={
         "motif": "DRACH",
         "sequence": "my_rna.fa",
         "color": "blue",
         "name": "m6A motif"
         }
      )

back to `Standard data keywords`_
