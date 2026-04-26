"""This module provides the ``baltic`` :class:`~baltic.clade.Clade` class.

**Notes**

This version of ``baltic`` (v0.9 (Birch)) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Clade")


class Clade(BranchLike): ## clade class
    """
    :class:`.Clade` represents a terminal placeholder for an entire monophyletic subtree.

    **Note**

    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    **Attributes**

    branchType : str
        DEPRECATED: should now be checked using :meth:`~Clade.is_leaf()`, or :meth:`~Clade.is_node()`.

    subtree : :class:`baltic.tree.Tree`
        :class:`baltic.tree.Tree` object containing all the branches that were collapsed.

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
        """
        Initialize a collapsed clade placeholder.

        **Parameters**

        name : str
            Label for the collapsed clade.

        subtree : list[:class:`.BranchLike`], optional
            Branches contained in the collapsed subtree.

        leaves : list[:class:`.Leaf`], optional
            Descendant leaves represented by the clade.

        lastHeight : float, optional
            Height of the deepest descendant in the collapsed subtree.

        lastAbsoluteTime : float, optional
            Absolute time of the deepest descendant in the collapsed subtree.

        width : float, optional
            Plotting width assigned to the clade. Defaults to ``1.0``.

        **Examples**

        >>> from baltic.clade import Clade
        >>> clade = Clade("example", width=2.0)
        >>> clade.name, clade.width
        ('example', 2.0)
        """
        super().__init__() # Inherit traits from BranchLike
        self.name=name
        self.branchType='clade'
        self.subtree=subtree
        self.leaves=leaves
        self.lastHeight=lastHeight
        self.lastAbsoluteTime=lastAbsoluteTime
        self.width=width

    def is_leaflike(self):
        """Returns ``True`` because :class:`.Clade` behaves like a terminal placeholder."""
        return True

    def is_leaf(self):
        """Returns ``False`` because :class:`.Clade` is not a true :class:`baltic.leaf.Leaf`."""
        return False

    def is_node(self):
        """Returns ``False`` because :class:`.Clade` is not an internal :class:`baltic.node.Node`."""
        return False

    def __str__(self):
        """
        Return a compact string representation of the clade.

        Matches the readable style used for :class:`.Clade` placeholders.

        **Returns**

        str
            Human-readable representation containing the clade name.
        """
        return f"Clade({self.name})"
