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

Source code
-----------
The source code is available on GitHub: `dev-doc branch <https://github.com/evogytis/baltic/tree/dev-doc>`_.

Contributing
------------
``baltic`` is an open source/science tool and very open to contributions that fix and/or expand its utility. While phylogenetc data manipulation is one aspect of baltic it is understood that visualising such data will be the main application and as such the design philosophy of baltic follows closely that of ``matplotlib`` - functionality in ``baltic`` itself should be basic and generally applicable with more involved visualisations combining elements thereof. The most sought-after contributions are therefore bug fixes, code to manipulate trees and the occasional novel visualisation technique.

Note that baltic follows and enforces the `Contributor Covenant Code of Conduct <https://www.contributor-covenant.org/>`_.


If you'd like to contribute a bug fix feel free to start a pull request outlining the problem and its solution.
For adding code for a new tree manipulation function or for a novel visualisation please open an issue first and describe the goal and suggested approach first.
Even if suggested manipulations/visualisations don't make it into baltic's code base (*e.g.* if the application isn't broad) there will be room to share it under examples.


.. toctree::
   :hidden:
   :maxdepth: 3

   baltic
   io
   node
   leaf
   clade
   reticulation
   tree
   bt_utils
   curonia
   samogitia
   publications
