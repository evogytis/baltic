"""This module provides the ``baltic`` :class:`~baltic.branchLike.BranchLike` superclass.

**Notes**

This version of ``baltic`` (v0.9 (Birch)) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""
import logging

logger = logging.getLogger("baltic.BranchLike")


class BranchLike:
    """
    :class:`.BranchLike` serves as the superclass inherited by
    :class:`.Clade`, :class:`.Leaf`, :class:`.Node`, and :class:`.Reticulation`.

    This class defines the recursive structures of which ``baltic`` :class:`.Tree` objects are built.

    **Note**

    Most attributes have null default values as they are either (1) set automatically
    during tree construction, or (2) defined by the specific inherited subclass.

    **Attributes**

    branchType : str
        DEPRECATED: should now be checked using :meth:`is_leaf()`, or :meth:`is_node()`.

    length : float, default=0.0
        The length of the incoming branch leading to the node/leaf.

    height : float
        Height, set by traversing the tree.

    absoluteTime : float
        Branch end point in absolute time.

    parent : :class:`baltic.branchLike.BranchLike`
        Reference to parent of ``self``.

    traits : dict, default={}
        Dictionary that will contain annotations from the tree string, e.g. ``{'posterior': 1.0}``.

    index : int
        Index of the character designating this object in the tree string.

        Uniquely identifies every object in the tree.

    x : float
        The x-coordinate of the node, according to the coordinate system of a plot.

    y : float
        The y-coordinate of the node, according to the coordinate system of a plot.
    """

    def __init__(
        self,
        length=0.0,
        height=None,
        absoluteTime=None,
        absoluteTimeRange=None,
        parent=None,
        traits=None,
        ):
        """
        Initialize a generic branch-like object.

        **Parameters**

        length : float, optional
            Length of the incoming branch. Defaults to ``0.0``.

        height : float, optional
            Height of the branch in tree coordinates.

        absoluteTime : float, optional
            Absolute time assigned to the branch endpoint.

        absoluteTimeRange : tuple[float, float], optional
            Uncertainty interval associated with *absoluteTime*.

        parent : :class:`.BranchLike`, optional
            Parent branch in the tree.

        traits : dict, optional
            Mapping of parsed annotations associated with the branch.

        **Examples**

        >>> from baltic.branchLike import BranchLike
        >>> branch = BranchLike(length=0.5, traits={"posterior": 1.0})
        >>> branch.length
        0.5
        """
        self.branchType=None
        self.length=length
        self.height=height
        self.absoluteTime=absoluteTime
        self.absoluteTimeRange=absoluteTimeRange ##todo: add this to docstring
        self.parent=parent
        self.traits=traits if traits is not None else {}

        # These three attributes are set later during tree construction / plotting
        self.index=None
        self.x=None
        self.y=None


    def get_path_to_root(self, path=None):
        """
        Recursively find the path from this node to the root, listed in reverse-chrnonological order (starting with current node, ending with the root).

        Operates by adding itself to the path, then recursively calling ``get_path_to_root`` on the current node's parent.

        **Parameters**

        path : list[:class:`.BranchLike`]
            The path that has been traversed so far.

        **Returns**

        path : list[:class:`.BranchLike`]

        **Examples**

        >>> import baltic as bt
        >>> ll = bt.make_tree("((A:1.0,B:1.0):1.0,C:1.0);", treeType="divergence")
        >>> [branch.index for branch in ll.get_leaf("A").get_path_to_root()]
        [2, 1]
        """
        if path is None: path = []

        #TODO: make this self.parent.parent or flag the superroot
        if self.parent.parent is None: # the root doesn't have a parent, so just return the path including the root // changed to .parent.parent because otherwise this function adds the "super root" - a fake parent to the root actually found in the tree string
            return path + [self]
        else: # if we aren't at the root, recurse
            newPath = path + [self]
            return self.parent.get_path_to_root(newPath)


    def get_siblings(self, include_self=False):
        """Get a list of the current node's siblings (i.e. the children of its parent).

        **Parameters**

        include_self : bool, default=False
            ``True`` if the list of siblings should include the current node.

        **Returns**

        list[:class:`.BranchLike`]

        **Examples**

        >>> import baltic as bt
        >>> ll = bt.make_tree("((A:1.0,B:1.0):1.0,C:1.0);", treeType="divergence")
        >>> siblings = ll.get_leaf("A").get_siblings()
        >>> [branch.name for branch in siblings]
        ['B']
        """
        if self.parent is None:
            logger.warning("Attempted to find siblings of root node. Returning empty list.")
            return []
        sibs = set(self.parent.children)  # using a set is easier because then we don't have to worry about order
        if include_self:
            return list(sibs)
        sibs.remove(self)
        return list(sibs)


    def is_leaflike(self) -> bool:
        """Returns ``True`` if the current node is a terminal node (i.e. :class:`.Leaf`, :class:`.Clade`, or :class:`.Reticulation`).

        **Returns**

        bool

        **Note**

        This is set to always return ``False``, however this behavior is overwritten by the appropriate subclasses.
        """
        return False


    def is_leaf(self) -> bool:
        """Returns ``True`` if the current node is a :class:`.Leaf` node.

        **Returns**

        bool

        **Note**

        This is set to always return ``False``, however this behavior is overwritten by the :class:`.Leaf` subclass.
        """
        return False


    def is_node(self) -> bool:
        """Returns ``True`` if the current node is an internal node.

        **Returns**

        bool

        **Note**

        This is set to always return ``False``, however this behavior is overwritten by the :class:`.Node` subclass.
        """
        return False


    def __str__(self):
        """
        Return a compact string representation of the branch.

        Mirrors the compact representations used by other
        :class:`.BranchLike` subclasses.

        **Returns**

        str
            Human-readable representation containing the branch index.
        """
        return f"BranchLike({self.index})"
