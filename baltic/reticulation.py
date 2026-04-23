"""This module provides the ``baltic`` :class:`~baltic.reticulation.Reticulation` class.

**Notes**

This version of ``baltic`` (v0.1) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Reticulation")


class Reticulation(BranchLike):
    """
    :class:`.Reticulation` represents an external node for recombination, conversion, and reassortment events.

    **Note**

    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    **Attributes**

    branchType : str
        DEPRECATED: should now be checked using :meth:`~Reticulation.is_leaf()`, or :meth:`~Reticulation.is_node()`.

    name : str
        Name for the reticulate node.

    width : float
        Plotting width for the reticulate branch.

        Default value is ``0.5``.

    target : :class:`.BranchLike`, optional
        The branch which the reticulation merges into.
    """

    def __init__(
        self,
        name,
        width=0.5,
        target=None,
        ):
        """
        Initialize a reticulation branch.

        **Parameters**

        name : str
            Identifier for the reticulation event.

        width : float, optional
            Plotting width assigned to the reticulation. Defaults to ``0.5``.

        target : :class:`.BranchLike`, optional
            Destination branch into which the reticulation merges.

        **Examples**

        >>> from baltic.reticulation import Reticulation
        >>> ret = Reticulation("R1", width=1.0)
        >>> ret.name, ret.width
        ('R1', 1.0)
        """
        super().__init__() # Inherit traits from BranchLike
        self.name=name
        self.branchType='reticulation'
        self.width=width
        self.target=target

    def is_leaflike(self):
        """Returns ``True`` because :class:`.Reticulation` is handled as a terminal branch."""
        return True

    def is_leaf(self):
        """Returns ``False`` because :class:`.Reticulation` is not a true :class:`baltic.leaf.Leaf`."""
        return False

    def is_node(self):
        """Returns ``False`` because :class:`.Reticulation` is not an internal :class:`baltic.node.Node`."""
        return False

    def __str__(self):
        """
        Return a compact string representation of the reticulation.

        Matches the readable style used for :class:`.Reticulation` objects.

        **Returns**

        str
            Human-readable representation containing the reticulation name.
        """
        return f"Reticulation({self.name})"
