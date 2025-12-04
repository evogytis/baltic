"""This module provides the baltic ``Node`` class.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Node")


class Node(BranchLike):
    """``Node`` class for internal nodes that represent hypothetical common ancestors.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str, default='node'
        DEPRECATED: should now be checked using ``is_leaf()``, or ``is_node()``.
    children : list[BranchLike], default=[]
        A list of the descendent branches of this node.
    childHeight : float, default=None
        The height of the youngest descendant of this node.
    leaves : set({bt.Leaf, bt.Clade}) : default={}
        A set of all the tips that descend from this node (i.e. the full clade defined by this node).
    """

    def __init__(self):
        super().__init__()
        self.branchType='node'
        self.children=[]
        self.childHeight=None
        self.leaves=set()


    def is_leaflike(self):
        """Returns False."""
        return False

    def is_leaf(self):
        """Returns False."""
        return False

    def is_node(self):
        """Returns True."""
        return True
