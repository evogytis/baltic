"""This module provides the baltic ``Node`` class.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Node")


class Node(BranchLike):
    """
    ``Node`` class for internal nodes that represent hypothetical common ancestors.

    Inherits from :class:`.BranchLike`.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str, default='node'
        DEPRECATED: should now be checked using :meth:`~Node.is_leaf()`, or :meth:`~Node.is_node()`.

    children : list[:class:`.BranchLike`]
        A list of the descendent branches of this node.

        By default this is an empty list, and should be populated during tree construction.

    childHeight : float
        The height of the youngest descendant of this node.

    leaves : set({:class:`.Leaf`})
        A set of all the tips that descend from this node (i.e. the full clade defined by this node).

        This set does not contain any internal nodes, only leaves.
    """

    def __init__(
        self,
        children=None,
        childHeight=None,
        leaves=None,
        ):
        super().__init__()
        self.branchType='node'
        self.children = children if children is not None else []
        self.childHeight = childHeight
        self.leaves = leaves if leaves is not None else set()


    def is_leaflike(self):
        """Returns ``False``."""
        return False

    def is_leaf(self):
        """Returns ``False``."""
        return False

    def is_node(self):
        """Returns ``True``."""
        return True

    def __str__(self):
        return f"Node({self.index})"
