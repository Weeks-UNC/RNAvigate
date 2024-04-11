``orfs``
~~~~~~~~

Annotation of open-reading frames

example uses:

- same as `motif`_
- coming soon: displaying amino acid translation and codon usage scores

input explaination:

- Input should be a dictionary containing:

   - ``"orfs"``:

      - ``"all"`` annotates all open-reading frames
      - ``"longest"`` annotates only the longest open reading frame

   - ``"sequence"``: same as `sequence`_ keyword
   - ``"color"``: a valid color or hexcode, e.g. ``"blue"``, ``"grey"``, or ``"#fa4ce2"``
   - ``"name"``: an arbitrary name to use on plots

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      main_orf={
         "orfs": "longest",
         "sequence": "my_sequence.fa",
         "color": "green",
         "name": "Longest ORF"
         }
      )

back to `Standard data keywords`_
