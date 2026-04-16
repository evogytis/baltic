``baltic``: the Backronymed Adaptable Lightweight Tree Import Code
==================================================================

``baltic`` is a Python library for parsing, manipulating, and visualizing phylogenetic trees. It is designed for exploratory phylogenetics work, especially with BEAST-style Nexus trees and Nextstrain/Auspice JSON exports, while staying lightweight enough for quick scripting.

Core capabilities
-----------------

1. Parse Newick, Nexus, and JSON tree formats into a shared in-memory tree model.
2. Traverse and manipulate trees programmatically.
3. Produce publication-ready plots with ``matplotlib``.

Start here
----------

New users should begin with :doc:`getting_started`, then use the API reference pages for module-level details.

.. toctree::
   :maxdepth: 2

   getting_started
   publications

.. toctree::
   :hidden:
   :maxdepth: 3

   baltic
   io
   branchLike
   node
   leaf
   clade
   reticulation
   tree
   bt_utils
   curonia
   samogitia
