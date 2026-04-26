"""This module provides the ``baltic`` :class:`~baltic.node.Node` class.

**Notes**

This version of ``baltic`` (v0.9 (Birch)) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""

import logging
from baltic.branchLike import BranchLike

logger = logging.getLogger("baltic.Node")


class Node(BranchLike):
    """
    :class:`.Node` represents an internal node corresponding to a hypothetical common ancestor.

    **Note**

    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    **Attributes**

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
        """
        Initialize an internal node.

        **Parameters**

        children : list[:class:`.BranchLike`], optional
            Descendant branches of the node.

        childHeight : float, optional
            Maximum descendant height observed below the node.

        leaves : set[str], optional
            Descendant leaf names associated with the node.

        **Examples**

        >>> from baltic.node import Node
        >>> node = Node(leaves={"A", "B"})
        >>> sorted(node.leaves)
        ['A', 'B']
        """
        super().__init__()
        self.branchType='node'
        self.children = children if children is not None else []
        self.childHeight = childHeight
        self.leaves = leaves if leaves is not None else set()


    def is_leaflike(self):
        """Returns ``False`` because :class:`.Node` is not terminal like :class:`baltic.leaf.Leaf`."""
        return False

    def is_leaf(self):
        """Returns ``False`` because :class:`.Node` is distinct from :class:`baltic.leaf.Leaf`."""
        return False

    def is_node(self):
        """Returns ``True`` for :class:`.Node` objects."""
        return True

    def __str__(self):
        """
        Return a compact string representation of the node.

        Matches the readable style used for :class:`.Node` debugging output.

        **Returns**

        str
            Human-readable representation containing the node index.
        """
        return f"Node({self.index})"
