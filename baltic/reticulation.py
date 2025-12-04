"""This module provides the baltic ``Node`` class.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Reticulation")


class Reticulation(BranchLike):
    """``Reticulation`` class for external nodes that represent recombination, conversion, and reassortment events.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str, default='leaf'
        DEPRECATED: should now be checked using ``is_leaf()``, or ``is_node()``.
    name : str
        Name for the reticulation node
    width : float, default=0.5
        Plotting width for the reticulate branch.
    target : BranchLike : default=None
        The branch which the reticulation merges into.
    """

    def __init__(self,name):
        super().__init__() # Inherit traits from BranchLike
        self.branchType='leaf'
        self.name=name
        self.width=0.5
        self.target=None

    def is_leaflike(self):
        """Returns True."""
        return True

    def is_leaf(self):
        """Returns False."""
        return False

    def is_node(self):
        """Returns False."""
        return False
