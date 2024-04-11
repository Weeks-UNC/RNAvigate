``interactions``
~~~~~~~~~~~~~~~~

inter-nucleotide data that does not have a more specific data keyword:

- `ringmap`_ for RingMapper correlations
- `pairmap`_ for PairMapper correlations
- `shapejump`_ for ShapeJump deletion events
- `pairprob`_ for pairing probabilities
- `allpossible`_ for every possible nucleotide pairing from a sequence
- interactions: for everything else

example uses:

- visualizing interaction networks in arc and circle plots, secondary structure
   diagrams and 3D molecule renderings
- filtering interactions based on many different factors.

   - see :doc:'/guides/filters' guide.

- calculating a distance distribution histogram of a set of interactions

input explaination:

- These inputs allow a lot of customization in loading data.
- For a full explaination, see :doc:`/guides/custom_interactions`

back to `standard data keywords`_
