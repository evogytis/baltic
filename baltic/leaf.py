"""This module provides the ``baltic`` :class:`~baltic.leaf.Leaf` class.

**Notes**

This version of ``baltic`` (v0.9 (Birch)) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Leaf")


class Leaf(BranchLike):
    """
    :class:`.Leaf` represents a terminal node corresponding to an individual taxon.

    **Note**

    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    **Attributes**

    branchType : str
        DEPRECATED: should now be checked using :meth:`~Leaf.is_leaf()`, or :meth:`~Leaf.is_node()`.

    name : str
        The tip name for the leaf.
    """

    def __init__(self, name):
        """
        Initialize a terminal leaf.

        This constructor specializes :class:`baltic.branchLike.BranchLike`
        for terminal taxa.

        **Parameters**

        name : str
            Tip label associated with the leaf.

        **Examples**

        >>> from baltic.leaf import Leaf
        >>> leaf = Leaf("A")
        >>> leaf.name
        'A'
        """
        super().__init__() # all of leaf's traits come from BranchLike
        self.branchType='leaf'
        self.name=name

    def is_leaflike(self):
        """Returns ``True`` for :class:`.Leaf` objects."""
        return True

    def is_leaf(self):
        """Returns ``True`` for :class:`.Leaf` objects."""
        return True

    def is_node(self):
        """Returns ``False`` because :class:`.Leaf` is not a :class:`baltic.node.Node`."""
        return False

    def __str__(self):
        """
        Return a compact string representation of the leaf.

        Matches the readable style used throughout :class:`.Leaf` handling.

        **Returns**

        str
            Human-readable representation containing the leaf name.
        """
        return f"Leaf({self.name})"
