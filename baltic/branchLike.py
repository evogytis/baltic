"""This module provides the baltic ``BranchLike`` superclass.
"""
import logging

logger = logging.getLogger("baltic.BranchLike")


class BranchLike:
    """``BranchLike`` class serves as super class that gets inherited by ``Clade``, ``Leaf``, ``Node``, and ``Reticulation``.

    This class defines the recursive structures of which baltic ``Tree`` objects are built.

    Note
    ----
    Most attributes have null default values as they are either (1) set automatically during tree construction, or (2) defined by the specific inherited subclass.

    Attributes
    ----------
    branchType : str, default=None
        DEPRECATED: should now be checked using ``is_leaf()``, or ``is_node()``.
    length : float, default=0.0
        The length of the incoming branch leading to the node/leaf.
    height : float, default=None
        Height, set by traversing the tree.
    absoluteTime : float, default=None
        Branch end point in absolute time.
    parent : baltic.branchLike, default=None
        Reference to parent of ``self``.
    traits : dict, default={}
        Dictionary that will contain annotations from the tree string, e.g. ``{'posterior': 1.0}``.
    index : int, default=None
        Index of the character designating this object in the tree string. Uniquely identifies every object in the tree.
    x : float, default=None
        The x-coordinate of the node, according to the coordinate system of a plot.
    y : float, default=None
        The y-coordinate of the node, according to the coordinate system of a plot.
    """


    def __init__(self):
        self.branchType=None
        self.length=0.0
        self.height=None
        self.absoluteTime=None
        self.absoluteTimeRange=None ##todo: add this to docstring
        self.parent=None
        self.traits={}
        self.index=None
        self.x=None
        self.y=None


    def get_path_to_root(self, path=None):
        """Recursively find the path from this node to the root, listed in reverse-chrnonological order (starting with current node, ending with the root).

        Operates by adding itself to the path, then recursively calling ``get_path_to_root`` on the current node's parent.

        Parameters
        ----------
        path : str(baltic.branchLike), default=None
            The path that has been traversed so far.

        Returns
        -------
        path : str(baltic.branchLike)
            The full path from self to the root of the tree.
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

        Parameters
        ----------
        include_self : bool, default=False
            ``True`` if the list of siblings should include the current node.

        Returns
        -------
        list(baltic.branchLike)
        """
        if self.parent is None:
            logger.warning("Attempted to find siblings of root node. Returning empty list.")
            return []
        sibs = set(self.parent.children)  # using a set is easier because then we don't have to worry about order
        if include_self:
            return list(sibs)
        sibs.remove(self)
        return list(sibs)


    def is_leaflike(self):
        """Returns ``True`` if the current node is a terminal node (i.e. ``Leaf``, ``Clade``, or ``Reticulation``).

        Note
        ----
        This is set to always return ``False``, however this behavior is overwritten by the appropriate subclasses.
        """
        return False


    def is_leaf(self):
        """Returns ``True`` if the current node is a ``Leaf`` node.

        Note
        ----
        This is set to always return ``False``, however this behavior is overwritten by the ``Leaf`` subclass.
        """
        return False


    def is_node(self):
        """Returns ``True`` if the current node is an internal node.

        Note
        ----
        This is set to always return ``False``, however this behavior is overwritten by the ``Node`` subclass.
        """
        return False


    def __str__(self):
        return f"BranchLike({self.index})"
