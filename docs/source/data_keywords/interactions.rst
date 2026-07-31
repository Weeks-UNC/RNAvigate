``interactions``
~~~~~~~~~~~~~~~~

inter-nucleotide data that does not have a more specific data keyword:

- :doc:`/data_keywords/ringmap` for RingMapper correlations
- :doc:`/data_keywords/pairmap` for PairMapper correlations
- :doc:`/data_keywords/shapejump` for ShapeJump deletion events
- :doc:`/data_keywords/pairprob` for pairing probabilities
- :doc:`/data_keywords/allpossible` for every possible nucleotide pairing from a sequence
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

back to :ref:`Standard data keywords <standard-data-keywords>`
