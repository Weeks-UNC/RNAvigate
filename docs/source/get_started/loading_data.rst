Loading data
============

RNAvigate is built around the `Sample`, which is a grouping of datasets that came from a single RNA studied under a single set of experimental conditions.
For example, a `Sample` could contain a sequence, primer location annotations, a ShapeMapper profile, and a predicted secondary structure for an *in-vitro* structure probing experiment.
A second `Sample` could contain the same data for an *in-vivo* experiment.

Creating a `Sample` and assigning data files to it using data keywords accomplishes 4 tasks:

1. The data are organized as a `Sample`.
2. The data are easy to access via the assigned data keywords.
3. The data keyword tells RNAvigate how to parse and represent the data as one of the data classes described below.
4. The data is then compatible with all of RNAvigate's visualization and analysis tools.

+---------------------+------------------------------------------------------+
| Data class          | Description                                          |
+=====================+======================================================+
| sequence            | an RNA sequence                                      |
+---------------------+------------------------------------------------------+
| annotation          | sites or regions of interest along an RNA sequence   |
+---------------------+------------------------------------------------------+
| secondary structure | the base-pairing pattern of an RNA sequence          |
+---------------------+------------------------------------------------------+
| tertiary structure  | the 3D atomic coordinates of an RNA sequence         |
+---------------------+------------------------------------------------------+
| profile             | per-nucleotide measurements along an RNA sequence    |
+---------------------+------------------------------------------------------+
| interactions        | inter-nucleotide measurements within an RNA sequence |
+---------------------+------------------------------------------------------+

Creating and using a Sample
---------------------------

Samples are created using ``rnav.Sample()``.

.. code-block:: python

   import rnavigate as rnav               # Load RNAvigate and give it the alias "rnav"

   my_sample = rnav.Sample(               # create a new sample
      sample="My sample name",            # provide a name for plot labels
      data_keyword="my_data.txt",         # load data file 1
      data_keyword2="my_other_data.txt",  # load data file 2
   )

Above, ``sample="My sample name"`` provides a label, to appear in plot titles and legends, for any data that came from this sample.
``"My sample name"`` should be replaced with any string that uniquely and succinctly identifies this sample.
A ``sample`` label is always required.

``data_keyword`` should be replaced with a data keyword appropriate for your specific data (see below).

Then, visualizing this data would look something like this:

.. code-block:: python

   plot = rnav.plot_arcs(         # represent my data as an arc plot
      samples=[my_sample],        # visualize my_sample
      sequence="data_keyword",    # positionally align all data to this sequence
      profile="data_keyword",     # display profile data
      structure="data_keyword2",  # display secondary structure
   )


``plot_arcs`` can be replaced with other plotting functions, which are introduced in the next guide: :doc:`visualizing_data`.

Before we get into data keywords, ``rnav.Sample`` accepts two other arguments: ``inherit`` and ``keep_inherited_defaults``.
These are used to share data between samples, e.g. a literature-accepted structure shared between experimental samples.
This sharing saves on memory and computation time.

Example usage:

.. code-block:: python

   shared_data = rnav.Sample(
      name='shared data',
      keyword1='big_structure.pdb')

   sample1 = rnav.Sample(
      name='knockout',
      inherit=shared_data,
      keyword2='sample1-data.txt')

   sample2 = rnav.Sample(
      name='control',
      inherit=shared_data,
      keyword2='sample2-data.txt')

``sample1`` and ``sample2`` now both have ``keyword1``, which is shared, and ``keyword2``, which is not.

At the moment, default keywords are only used to simplify data keyword inputs. For
example, the ``ringmap`` data keyword uses the sequence provided by ``default_profile``,
which is the first profile-type data provided to the Sample.

Data keywords
-------------

Data keywords can either be an arbitrary keyword or a standard keyword:

Arbitrary data keywords
```````````````````````

An arbitrary keyword is useful if you are loading 2 or more of the same data
type into a single sample. Arbitrary keywords must follow some simple rules:

1. Cannot conflict with a given sample's other data keywords.
2. Cannot be ``inherit`` or ``keep_inherited_defaults``
3. Cannot consist only of valid nucleotides: ``AUCGTaucgt``
4. Cannot start with a number: ``0123456789``
5. Must only contain numbers, letters and underscores.

If an arbitrary data keyword is used, a dictionary must be provided, specifying
the standard data keyword to use for parsing inputs.

Example:

.. code-block:: python

   my_sample = rnav.Sample(
      sample="example",
      standard_keyword="input_file_1.txt",
      arbitrary_keyword={"standard_keyword": "input_file_2.txt"}
   )

Standard data keywords
``````````````````````

.. contents::
   :local:

Sequence data
^^^^^^^^^^^^^

.. include:: ../data_keywords/sequence.rst

Annotation data
^^^^^^^^^^^^^^^

.. include:: ../data_keywords/motif.rst
.. include:: ../data_keywords/orfs.rst
.. include:: ../data_keywords/spans.rst
.. include:: ../data_keywords/sites.rst
.. include:: ../data_keywords/group.rst
.. include:: ../data_keywords/primers.rst
.. include:: ../data_keywords/domains.rst

Secondary structure data
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../data_keywords/ss.rst

Tertiary structure data
^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../data_keywords/pdb.rst

Profile data
^^^^^^^^^^^^

.. include:: ../data_keywords/profile.rst
.. include:: ../data_keywords/shapemap.rst
.. include:: ../data_keywords/dancemap.rst
.. include:: ../data_keywords/rnpmap.rst

Interactions data
^^^^^^^^^^^^^^^^^

.. include:: ../data_keywords/interactions.rst
.. include:: ../data_keywords/ringmap.rst
.. include:: ../data_keywords/pairmap.rst
.. include:: ../data_keywords/shapejump.rst
.. include:: ../data_keywords/pairprob.rst
.. include:: ../data_keywords/allpossible.rst
