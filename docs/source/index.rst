.. BALTIC documentation master file, created by
   sphinx-quickstart on Thu Oct  2 13:42:26 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

baltic: the Backronymed Adaptable Lightweight Tree Import Code
===============================================================

**********
background
**********

BALTIC is a Python library for the parsing, manipulation, and visualization of phylogenetic trees. It has been developed to extract various statistics from molecular phylogenies generated using `BEAST X <https://beast.community>`_ in a customised way. This `influenza B virus reassortment paper <https://dx.doi.org/10.1093/molbev/msu287>`_ used early versions of baltic's code to look at how the human influenza B virus segment diversity is structured according to genomic background.

BALTIC's core functionality comprises two major components:

1) Phylogenetic tree parsing and manipulation
2) Manuscript-quality visualization leveraging `matplotlib <https://matplotlib.org/>`_

BALTIC has already been used in over forty :doc:`publications <publications>`!

************
Usage basics
************

By convention BALTIC is imported as ``bt``, and instances of BALTIC tree objects are given the name ``ll``. When called with a tree string the ``make_tree`` function return a baltic tree object:

.. code-block:: python

   import baltic as bt
   treeString='((A:1.0,B:2.0):1.0,C:3.0);'
   ll = bt.make_tree(treeString, treeType="divergence")

The resultant BALTIC tree can be examined via the ``treeStats`` method.

.. code-block:: python

   ll.treeStats()

.. code-block:: text

   Tree height: 3.000000
   Tree length: 7.000000
   strictly bifurcating tree

   Numbers of objects in tree: 5 (2 nodes and 3 leaves)

.. toctree::
   :hidden:
   :maxdepth: 2

   baltic
   io
   branchLike
   node
   leaf
   clade
   reticulation
   tree
   bt_utils
   publications
