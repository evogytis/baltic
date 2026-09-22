"""This module provides the ``baltic`` package's top-level import surface: the :class:`~baltic.tree.Tree` and :class:`~baltic.branchLike.BranchLike` subclasses, :func:`~baltic.baltic.make_tree` and :func:`~baltic.baltic.make_tree_JSON`, plus everything re-exported from :mod:`~baltic.io` and :mod:`~baltic.bt_utils`.

**Notes**

This version of ``baltic`` (v1.0 (Cedar)) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.
"""

from baltic.tree import Tree
from baltic.branchLike import BranchLike
from baltic.node import Node
from baltic.leaf import Leaf
from baltic.reticulation import Reticulation
from baltic.clade import Clade
from baltic.baltic import make_tree, make_tree_JSON
from baltic.io import *
from baltic.bt_utils import *
