"""This module provides the baltic ``Clade`` class.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Clade")


class Clade(BranchLike): ## clade class
    """``Clade`` class for terminal nodes that represent an entire monophyletic subtree.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str, default='leaf'
        DEPRECATED: should now be checked using ``is_leaf()``, or ``is_node()``.

    subtree : Tree
        Tree object containing all the branches that were collapses.

    leaves : list[:class:`.Leaf`]
        List of descendent leaves contained in the clade.

    name : str
        The pretend tip name for the clade.

    lastHeight : float
        The height of the highest (furthest in the past) tip in the collapsed tree.

    lastAbsoluteTime : float
        The absolute time of the highest (furthest in the past) tip in the collapsed clade.

    width : float
        Width value used in plotting of the collapsed clade.
    """

    def __init__(self,name, branchType='leaf', subtree=None, leaves=None, lastHeight=None, lastAbsoluteTime=None, width=1.0):
        self.name=name
        self.branchType=branchType
        self.subtree=subtree
        self.leaves=leaves
        self.lastHeight=lastHeight ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime=lastAbsoluteTime ## refers to the absolute time of the highest tip in the collapsed clade
        self.width=1.0

    def is_leaflike(self):
        """Returns True."""
        return True

    def is_leaf(self):
        """Returns False."""
        return False

    def is_node(self):
        """Returns False"""
        return False

    def __str__(self):
        return f"Clade({self.name})"
