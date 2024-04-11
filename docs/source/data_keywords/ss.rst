``ss``
~~~~~~

A secondary structure with optional diagram drawing

example uses:

- visualizing base pairs on arc plots, circle plots, and secondary structure
   diagrams
- calculating contact distances

   - the shortest path between nucleotides in a secondary structure graph

- determining how well per-nucleotide data predict base pairing status (AUROC)

input explaination:

- Input should be one of the following formats:

   - secondary structure files (no diagram)

      - connection table (.ct)
      - dotbracket notation (.dot, .dbn, etc.)

   - secondary structure diagram files

      - `StructureEditor`_ (.nsd or .cte)
      - `XRNA`_ (.xrna)
      - `VARNA`_ (.varna)
      - `FORNA`_ (.json)

         - click "add molecule" and paste in a dotbracket notation structure.
         - arrange it how you like
         - click the download button in the lower-right, then click "json"

      - `R2DT`_ (.json)

         - Type in an RNA sequence, R2DT creates the secondary structure
         - Click on the R2DT paper link to learn more about how it works
         - Once the structure is drawn, click "Edit in XRNA"
         - Arrange it how you like it
         - In the upper-left, type in a file name, choose "json", click "download"

- Note: The file format is determined by the file extension. Since FORNA and R2DT
  both produce json, the extension should be provided.

.. _StructureEditor: https://rna.urmc.rochester.edu/RNAstructure.html
.. _XRNA: http://rna.ucsc.edu/rnacenter/xrna/xrna.html
.. _VARNA: https://varna.lisn.upsaclay.fr/
.. _FORNA: http://rna.tbi.univie.ac.at/forna/
.. _R2DT: https://rnacentral.org/r2dt

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      ss="my_structure.ct"
      )

other optional inputs:

- ``"extension"`` is used to specify a file extension.
- ``"autoscale"`` is used to scale coordinates to look good in RNAvigate plots
- ``"structure_number"`` is used to specify which structure to load if the file
   contains multiple structures (0-indexed). Default is to load the first
   structure. This currently only works with "ct" files.

typical optional argument examples:

.. code-block:: python

   # specify r2dt vs forna json
   my_sample = rnav.Sample(
      sample="example",
      ss={
         "ss": "my_rna.json",
         "extension": "r2dt,
      }
   )

   #specify structure number
   my_sample = rnav.Sample(
      sample="example",
      ss={
         "ss": "my_rna.ct",
         "structure_number": 3,
      }
   )

back to `Standard data keywords`_
