"""This module provides the baltic ``Clade`` class.

Notes
-----
This version of BALTIC (v0.1) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


Attributes
----------
logger : :external+python:py:class:`logging.Logger`
    Default logger which will be passed to other baltic functions.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Clade")


class Clade(BranchLike): ## clade class
    """
    ``Clade`` class for terminal nodes that represent an entire monophyletic subtree.

    Inherits from :class:`.BranchLike`.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str
        DEPRECATED: should now be checked using :meth:`~Clade.is_leaf()`, or :meth:`~Clade.is_node()`.

    subtree : Tree
        Tree object containing all the branches that were collapses.

    leaves : list[:class:`.Leaf`]
        List of descendent leaves contained in the clade.

    name : str
        The label for the clade.

    lastHeight : float
        The height of the highest (furthest in the past) tip in the collapsed tree.

    lastAbsoluteTime : float
        The absolute time of the highest (furthest in the past) tip in the collapsed clade.

    width : float
        Width value used in plotting of the collapsed clade.

        Default value is ``1.0``.
    """

    def __init__(
        self,
        name,
        subtree=None,
        leaves=None,
        lastHeight=None,
        lastAbsoluteTime=None,
        width=1.0,
        ):
        super().__init__() # Inherit traits from BranchLike
        self.name=name
        self.branchType='clade'
        self.subtree=subtree
        self.leaves=leaves
        self.lastHeight=lastHeight
        self.lastAbsoluteTime=lastAbsoluteTime
        self.width=width

    def is_leaflike(self):
        """Returns ``True``."""
        return True

    def is_leaf(self):
        """Returns ``False``."""
        return False

    def is_node(self):
        """Returns False"""
        return False

    def __str__(self):
        return f"Clade({self.name})"
