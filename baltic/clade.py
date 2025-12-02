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
    subtree : bt.Tree, default=None
        Tree object containing all the branches that were collapses.
    leaves : list[bt.Leaf], default=None
        List of descendent leaves contained in the clade.
    name : str
        The pretend tip name for the clade.
    lastHeight : float, default=None
        The height of the highest (furthest in the past) tip in the collapsed tree.
    lastAbsoluteTime : float, default=None
        The absolute time of the highest (furthest in the past) tip in the collapsed clade
    width : float, default=1.0
        Width value used in plotting of the collapsed clade.
    """

    def __init__(self,givenName):
        self.branchType='leaf'
        self.subtree=None
        self.leaves=None
        self.name=givenName
        self.lastHeight=None ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime=None ## refers to the absolute time of the highest tip in the collapsed clade
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
