Filtering
=========

Any visualization function that can display interactions (inter-nucleotide measurements)
can also reconfigure how that data is displayed. This includes how measurement values are mapped to colors,
and which interactions are shown. Thisese visualizations include arc plots, circle plots, secondary
structure diagrams, 3D molecules, distance distributions, and heatmaps.

The ``interactions`` argument of plotting functions can be tailored to automatically
filter one or more datasets on the fly. This argument can take one of three forms:

1. If provided with a data keyword, the full data set will be used without filtering, e.g.:

.. code-block:: python

   interactions="ringmap"

2. If provided with a dictionary, the key/value pairs in the dictionary define the filtering scheme and data representation, e.g.:

.. code-block:: python

   interactions={
      "interactions": "ringmap",
      "cmap": "Greens",
      "metric": "Zij",
      "ss_only": True,
   }

3. If provided a list of dictionaries, each will be interpretted as in option 2, and displayed on seperate axes, e.g.:

.. code-block:: python

   interactions=[
      {"interactions": "ringmap", "cmap": "Greens", "metric": "Zij", "ss_only": True},
      {"interactions": "shapejump", "metric": "Percentile", "Percentile_ge": 0.95},
   ]

Note that plotting functions typically produce a seperate axis for each sample as well.
Therefore, the total number of axes will be (number of samples * number of interactions).

Interactions configuration keys
-------------------------------

Below is an explanation of all of the available filters and their default behaviors.
I will use the shorthand *i* and *j* to describe the 5' and 3' ends of an interaction.
Some defaults depend on the interaction data object, which will be referred to as ``data``.

.. contents::
   :local:

Color display arguments
^^^^^^^^^^^^^^^^^^^^^^^

These 4 values determine how the data are to be represented as in color space. Breifly:

- metric: the column of the interactions dataframe to represent
- cmap: the colors to be used
- normalization: the method used to normalize data to fit the colormap
- values: the values to use with the normalization method

See below for more detail on each of these values.

The process of going from data values to colors is performed step-wise:

#. An array of values are taken from ``data.data[metric]``.
#. The values are normalized (integers for discrete colors, 0-1 for continuous colors).
#. The values are then mapped to colorspace using the colormap.

``"metric": None``
~~~~~~~~~~~~~~~~~~

- This can be set to any column name (string) from the interactions dataframe.
- The default is the value of ``data.default_metric``:

   - ``"Statistic"`` for RING-MaP data
   - ``"Class"`` for PAIR-MaP data
   - ``"Probability"`` for pairing probabilities
   - ``"Percentile"`` for SHAPE-JuMP data

- Each value in this column is mapped to a color.
- The column name can also be ``"Distance"``, which uses the 3D
  distance of the ``"O2'"`` atoms of *i* and *j* in the given PDB file.
- You can specifiy the atom that distances are calculated from by appending
  ``"_atom"`` , e.g.: ``"metric": "Distance_O3'"``
- ``"Distance_DMS"`` will use ``"N1"`` for A and C, and ``"N3"`` for U and G.

``"cmap": None``
~~~~~~~~~~~~~~~~

- This can be any of the following:

   - a single color, applied to all values
   - a list of colors, for a discrete colormap
   - the name of a matplotlib colormap, usually a continuous colormap. Search the web
     for "matplotlib colormaps" for colormap names.
   - For more control, users can provide a predefined colormap, usually the return
     value from a matplotlib or seaborn function. These can be discrete or continuous.

- Color and colormap names must be valid for matplotlib. e.g., ``"red", "#f0f0f0", "DodgerBlue", (0, 0.5, 0.5)``, etc.
- The default is the value of ``data.metric_defaults[data.default_metric]["cmap"]``:

  - ``"bwr"`` (blue to white to red, continuous) for RING-MaP data
  - ``matplotlib.colors.ListedColormap([[0.3, 0.3, 0.3, 0.2], [0.0, 0.0, 0.95, 0.6], [0.12, 0.76, 1.0, 0.6]])`` (grey, light blue, dark blue) for PAIR-MaP class.
  - ``seaborn.cubehelix_palette(10, 0.7, 0.9, 1.5, 2.5, 1, 0.4, False, True)`` for pairing probabilities.
  - ``"jet"`` (dark blue to light yellow to dark red) for 3D distances.
  - ``"YlGnBu"`` (yellow to green to blue) for SHAPE-JuMP data.

``"normalization": None`` and ``"values": None``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``"normalization"`` and ``"values"`` work together to fit the data values to the colormap.
- The behavior of ``"values"`` depends on the value of ``"normalization"``:

   - ``"normalization": "min_max", "values": [minumum, maximum]``

      - Data values are winsorized to be between ``minimum`` and ``maximum`` (floats).
      - Values are then scaled to between 0 and 1 to be used with a continuous colormap.

   - ``"normalization": "0_1", "values": None``

      - This behavior is similar to that above, but without winsorization.

   - ``"normalization": "bins", "values": [boundary_1, boundary_2, ..., boundary_n]``

      - Values are placed into integer bins based on the boundaries to use with a discrete colormap.
      - The number of boundaries must be 1 less than the number of discrete colors.
      - ``boundary_n`` is the value of the boundary between colors *n* and *n+1*.
      - The first and last bins extend to infinity.

   - ``"normalization": "none", "values": None``

      - No normalization is performed.
      - If the colormap is continuous, the data values must already be between 0 and 1.
      - If the colormap is discrete, the data values must be integers less than the length of the colormap.

- By default:

.. code-block:: python

   "normalization": data.metric_defaults[data.default_metric]["normalization"],
   "values": data.metric_defaults[data.default_metric]["values"],

Filtering arguments
^^^^^^^^^^^^^^^^^^^

``"prefiltered": False``
~~~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to retain previously applied filters.
- Use this only if you are manually filtering data outside of the plotting function.
- Do not use it to avoid re-typing a long list of filters.

Filtering on sequence
^^^^^^^^^^^^^^^^^^^^^

``"compliments_only": False``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to require that *i* and *j* are reverse complimentary

``"nts": None``
~~~~~~~~~~~~~~~

- A string containing valid nucleotides, e.g.: ``"nts": "AUC"`` excludes
  interactions in which *i* or *j* contains a G nucleotide

Filtering by position
^^^^^^^^^^^^^^^^^^^^^

``"max_distance": None`` and ``"min_distance": None``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A minimum and maximum primary sequence distance between *i* and *j*.

``"exclude_nts": None`` and  ``"isolate_nts": None``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A list of nucleotide positions to exclude or to isolate (meaning all other
  positions are excluded).

Filtering on per-nucleotide data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``"profile": None``
~~~~~~~~~~~~~~~~~~~

- A key of `sample.data`, *i* and *j* values will be taken from `profile.metric`.

``"min_profile": None`` and ``"max_profile": None``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The minimum and maximum allowable values of ``profile.metric`` at *i* and *j*.

Filtering on structure
^^^^^^^^^^^^^^^^^^^^^^

``"ct": None``
~~~~~~~~~~~~~~

- A key of `sample.data` pointing to a secondary structure.
- Contact distances and base-pairing status for the below filters are taken from this data.

``"min_cd: None"`` and ``"max_cd": None``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A minimum or maximum contact distance between *i* and *j*.
- Contact distance is the shortest path distance through the secondary structure graph.
- NOTE: Contact distance calculations can take a long time for large structures.

``"paired_only": False``
~~~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to require that *i* and *j* are base-paired to each other.

``"ss_only": False`` and ``"ds_only": False``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to require that *i* and *j* are both single-stranded or both double-stranded, respectively.

RING-MaP specific filters
^^^^^^^^^^^^^^^^^^^^^^^^^

``"positive_only": False`` and ``"negative_only": False``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to require that RING-MaP correlations are positive or negative.

PAIR-MaP specific filters
^^^^^^^^^^^^^^^^^^^^^^^^^

``"all_pairs": False``
~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to include PAIRs that are not classified as Primary or Secondary.

Filtering on data values
^^^^^^^^^^^^^^^^^^^^^^^^

kwargs: ``"ColumnName_operator": value``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``"ColumnName"`` is a valid column of interactions data
- ``"operator"`` is one of:

   - ``_ge`` (greater that or equal to)
   - ``_le`` (less than or equal to)
   - ``_gt`` (greater than)
   - ``_lt`` (less than)
   - ``_eq`` (equal to)
   - ``_ne`` (not equal to)

- Interactions data from this column is compared to the given value using the
  given operator. Interactions where that operation is ``True`` are kept.
- For example, ``"Statistic_ge": 23`` requires that the `"Statistic"` column is greater
  than 23.

Experimental filter
^^^^^^^^^^^^^^^^^^^

``"resolve_conflicts": None``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``True`` or ``False``, whether to reduce the number of interactions by finding
  the maximum weighted set of non-conflicting interactions (parallel
  correlations do not conflict)
