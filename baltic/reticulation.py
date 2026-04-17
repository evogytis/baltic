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
        """
        super().__init__() # Inherit traits from BranchLike
        self.name=name
        self.branchType='reticulation'
        self.width=width
        self.target=target

    def is_leaflike(self):
        """Returns ``True``."""
        return True

    def is_leaf(self):
        """Returns ``False``."""
        return False

    def is_node(self):
        """Returns ``False``."""
        return False

    def __str__(self):
        """
        Return a compact string representation of the reticulation.

        **Returns**

        str
            Human-readable representation containing the reticulation name.
        """
        return f"Reticulation({self.name})"
