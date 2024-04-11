``pdb``
~~~~~~~

A tertiary structure with atomic coordinates

example uses:

- rendering 3D molecules with data overlayed
- computing 3D distances between nucleotides

input explaination:

- Input should be a dictionary containing these keys:

   - "pdb": a standard PDB file (.pdb or .cif)
   - "chain": A chain ID
   - "sequence": same inputs as `sequence`_ keyword

      - This is not needed if a sequence is found in the file header.

example inputs:

.. code-block:: python

   my_sample_1 = rnav.Sample(
      sample="example1",
      pdb={
         "pdb": "my_structure.pdb",
         "chain": "X",
         "sequence": "my_rna.fa" # not needed if sequence in pdb header
         }
      )

back to `Standard data keywords`_
