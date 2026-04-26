Getting Started
===============

Installation
------------

Install ``baltic`` from PyPI:

.. code-block:: bash

   pip install baltic

The package depends on ``numpy`` and ``matplotlib``. If you plan to load JSON trees from remote URLs, install ``requests`` as well.

Create a tree from a Newick string
----------------------------------

.. code-block:: python

   import baltic as bt

   tree_string = "((A:1.0,B:2.0):1.0,C:3.0);"
   ll = bt.make_tree(tree_string, treeType="divergence")
   ll.treeStats()

``treeType`` should be one of:

- ``"divergence"`` for branch lengths in relative units.
- ``"time"`` for time-calibrated trees.

Load trees from files
---------------------

The main loader functions are:

- ``bt.load_newick(...)``: for Newick files
- ``bt.load_nexus(...)``: for Nexus files (e.g. from BEAST analyses)

Examples:

.. code-block:: python

   import baltic as bt

   newick_tree = bt.io.load_newick("tree.nwk", treeType="divergence")
   nexus_tree = bt.io.load_nexus("tree.nex", treeType="time")

If sampling dates are encoded in tip labels, use ``tipRegex`` and ``dateFmt`` to extract them:

.. code-block:: python

   ll = bt.io.load_nexus(
       "example.tree",
       treeType="time",
       tipRegex=r"\|([0-9\-]+)$",
       dateFmt="%Y-%m-%d",
       absoluteTime=True,
   )

Extract features of a tree
--------------------------

Once loaded, a ``Tree`` exposes helpers for inspection and traversal:

.. code-block:: python

   ll.traverse_tree()

   tips = ll.get_external()
   internal_nodes = ll.get_internal()
   stats = ll.treeStatsDict()

   print(stats["treeHeight"])
   print(len(tips))

Plot a tree
-----------

``baltic`` integrates directly with ``matplotlib`` for plotting functionality.

``baltic`` plotting functions always require a ``matplotlib`` axis object, which can be created
with ``plt.subplots()`` or similar functions. They also allow users to pass in styling
arguments (e.g. for colors, line widths, marker sizes) and callback functions
(e.g. for filtering, sorting, coordinate adjustments) to customize the plot.

Styling functions can be passed as ``baltic`` keyword arguments (e.g. ``colorFxn``, ``sizeFxn``,
``lineWidthFxn``), or as standard ``matplotlib`` keyword arguments (e.g. ``color``, ``markersize``,
``linewidth``), which get passed onward as keyword arguments to the underlying ``matplotlib``
plotting functions.

.. code-block:: python

   import matplotlib.pyplot as plt
   import baltic as bt

   ll = bt.io.load_newick("tree.nwk", treeType="divergence")

   fig, ax = plt.subplots(figsize=(8, 10))
   # Plot the branches using `plot_tree`
   ll.plot_tree(ax)
   # Plot the tip labels using `plot_points`
   #   in this case, we use a lambda to filter for leaf nodes (tips) and plot their names
   ll.plot_text(ax, targetFxn=lambda k: k.is_leaf())
   # Plot the points for the tips using `plot_points`
   ll.plot_points(ax)

   ax.set_axis_off()
   fig.tight_layout()
   plt.show()

Useful plotting helpers include:

- ``plot_tree()`` for branch geometry
- ``plot_points()`` for node or tip markers
- ``plot_text()`` for labels
- ``plot_aligned_tip_labels()`` for aligned tip names
- ``plot_exploded_tree()`` for trait-partitioned layouts

The plotting ``treeType`` argument supports ``"rectangular"``, ``"circular"``, and ``"unrooted"`` layouts.

The ``baltic`` way
------------------

By convention, ``baltic`` is imported as ``bt`` and tree objects are often named ``ll`` (for linked-list) in examples.

.. code-block:: python

   import baltic as bt
   ll = bt.make_tree("((A:1.0,B:2.0):1.0,C:3.0);", treeType="divergence")

By convention, ``baltic`` also uses Python lambda functions in many places, especially for short filtering, sorting, coordinate, and styling callbacks. These can always be replaced with regular named functions, but ``baltic`` tends to use lambdas where possible because they are succinct and usually easy to read in context.

.. code-block:: python

   # Example of a lambda for filtering tips
   # This function returns a list of tips whose names start with "A"
   tips = ll.get_external(lambda k: k.name.startswith("A"))

   # The same function as a regular named function
   # these two code blocks do the same thing, but the lambda is more concise
   # and fits the baltic convention
   def filter_tips_starting_with_A(node):
       return node.name.startswith("A")
   tips = ll.get_external(filter_tips_starting_with_A)

API reference
-------------

After this guide, continue with the module reference:

- :doc:`baltic`
- :doc:`io`
- :doc:`tree`
- :doc:`bt_utils`
