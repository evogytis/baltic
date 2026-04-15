"""This module provides the baltic ``Leaf`` class.

Notes
-----
This version of ``baltic`` (v0.1) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


Attributes
----------
logger : ``logging.Logger``
    Default logger which will be passed to other baltic functions.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Leaf")


class Leaf(BranchLike):
    """
    ``Leaf`` class for terminal nodes that represent individual taxa.

    Inherits from :class:`.BranchLike`.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str
        DEPRECATED: should now be checked using :meth:`~Leaf.is_leaf()`, or :meth:`~Leaf.is_node()`.

    name : str
        The tip name for the leaf.
    """

    def __init__(self, name):
        """
        Initialize a terminal leaf.

        Parameters
        ----------
        name : str
            Tip label associated with the leaf.
        """
        super().__init__() # all of leaf's traits come from BranchLike
        self.branchType='leaf'
        self.name=name

    def is_leaflike(self):
        """Returns ``True``."""
        return True

    def is_leaf(self):
        """Returns ``True``."""
        return True

    def is_node(self):
        """Returns ``False``."""
        return False

    def __str__(self):
        """
        Return a compact string representation of the leaf.

        Returns
        -------
        str
            Human-readable representation containing the leaf name.
        """
        return f"Leaf({self.name})"
