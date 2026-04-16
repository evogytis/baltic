Getting Started
===============

Installation
------------

Install ``baltic`` from PyPI:

.. code-block:: bash

   pip install baltic

The package depends on ``numpy`` and ``matplotlib``. If you plan to load JSON trees from remote URLs, install ``requests`` as well.

Import convention
-----------------

By convention, ``baltic`` is imported as ``bt`` and tree objects are often named ``ll`` in examples.

.. code-block:: python

   import baltic as bt

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

- ``bt.load_newick(...)``
- ``bt.load_nexus(...)``
- ``bt.load_posterior_nexus(...)``
- ``bt.load_JSON(...)``

Examples:

.. code-block:: python

   import baltic as bt

   newick_tree = bt.io.load_newick("tree.nwk", treeType="divergence")
   nexus_tree = bt.io.load_nexus("tree.nex", treeType="time")
   json_tree, json_meta = bt.io.load_JSON("auspice.json", treeType="time")

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

``baltic`` integrates directly with ``matplotlib``.

.. code-block:: python

   import matplotlib.pyplot as plt
   import baltic as bt

   ll = bt.io.load_newick("tree.nwk", treeType="divergence")

   fig, ax = plt.subplots(figsize=(8, 10))
   ll.plot_tree(ax)
   ll.plot_text(ax, targetFxn=lambda k: k.is_leaf())

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

API reference
-------------

After this guide, continue with the module reference:

- :doc:`baltic`
- :doc:`io`
- :doc:`tree`
- :doc:`bt_utils`
