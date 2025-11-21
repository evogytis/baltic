import copy
import json
import math
import re
import sys
from collections.abc import Callable
from functools import reduce
from operator import attrgetter
from pathlib import Path
from statistics import mean
from typing import Any, Literal

from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from typing_extensions import TypeIs

from baltic.utils import always_true, calendarDate, convertDate, decimalDate, initialized_property

__all__ = [
    "decimalDate",
    "convertDate",
    "calendarDate",
    "reticulation",  # make from baltic import * safe
    "clade",
    "leaf",
    "node",
    "tree",
    "make_tree",
    "make_treeJSON",
    "loadJSON",
    "loadNexus",
    "loadNewick",
    "untangle",
]

sys.setrecursionlimit(9001)


class Branch:
    """Parent class to tree components (nodes, tips, reticulation events and collapsed nodes)

    Attributes:
        branchType ("leaf" | "node"): Type of branch (defined in subclasses)
        x (float | None): x-coordinate for plotting, default is `None`
        y (float | None): y-coordinate for plotting, default is `None`
        length (float): Length of the branch, assigned in `make_tree()` or `collapseSubtree()`
        height (float): Height of the branch, assigned in `traverse_tree()` or `collapseSubtree()`
        parent (node): Parent node, assigned in `make_tree()` or `collapseSubtree()`
        traits (dict): Dictionary of traits associated with this object, assigned in `make_tree()` or `collapseSubtree()`
        index (int | str): The index of the character that defines this object in the tree string, or the parent node in `clade`
    """

    branchType: Literal["leaf", "node"]
    x: float | None
    y: float | None
    traits: dict
    absoluteTime: float | None

    def __init__(self, branchType: Literal["leaf", "node"]):
        self.branchType = branchType
        self.x = None
        self.y = None
        self.absoluteTime = None
        self.traits = {}
        self._height = None

    @initialized_property
    def length(self) -> float: ...

    @initialized_property
    def height(self) -> float: ...

    @initialized_property
    def parent(self) -> "node": ...

    @initialized_property
    def index(self) -> int | str: ...

    def is_node(self) -> bool:
        return isinstance(self, node)

    def is_leaf(self) -> bool:
        return isinstance(self, leaf)

    def is_leaflike(self) -> bool:
        return isinstance(self, (clade, leaf, reticulation))

    def __str__(self):
        return f"{self.__class__.__name__} {self.index!r}"


class reticulation(Branch):  ## reticulation class (recombination, conversion, reassortment)
    """
    Represents a reticulation event in a phylogenetic tree, such as recombination or reassortment.

    Attributes:
    branchType (str): The type of branch, default is 'leaf'.
    length (float): The length of the branch, assigned in `make_tree()`.
    height (float): The height of the branch, assigned in `traverse_tree()`.
    absoluteTime (float or None): The absolute time of the event, default is None, typically assigned in `setAbsoluteTime()`.
    parent (node): The parent node, assigned in `make_tree()`.
    traits (dict): Dictionary of traits associated with the reticulation event, default is empty.
    index (int): The index of the node in the tree string, assigned in `make_tree()`.
    name (str): The name of the reticulation event.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.
    width (float): The width of the node for plotting, default is 0.5.
    target (node or leaf): The target node where reticulation event lands, assigned in `make_tree()`.

    Docstring generated with ChatGPT 4o.
    """

    def __init__(self, name: str):
        super().__init__("leaf")

        self.name = name

        self.length = 0.0
        self.height = 0.0
        self.width = 0.5

    @initialized_property
    def target(self) -> "node | leaf": ...


class clade(Branch):  ## clade class
    """
    Represents a collapsed clade in a phylogenetic tree.

    Attributes:
    branchType (str): The type of branch, default is 'leaf'.
    subtree (list): The subtree containing all the branches that were collapsed, assigned in `collapseSubtree()`.
    leaves (set): Names of descendant tips in the collapsed clade, assigned in `collapseSubtree()`.
    length (float): The length of the branch, assigned in `collapseSubtree()`.
    height (float or None): The height of the branch, assigned in `collapseSubtree()`.
    absoluteTime (float or None): The absolute time of the event, assigned in `collapseSubtree()`.
    parent (node): The parent node, assigned in `collapseSubtree()`.
    traits (dict): Dictionary of traits associated with the clade, assigned in `collapseSubtree()`.
    index (int or None): The index of the node, assigned in `collapseSubtree()`.
    name (str): The name assigned to the clade when it was collapsed.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.
    lastHeight (float): The height of the highest tip in the collapsed clade, assigned in `collapseSubtree()`.
    lastAbsoluteTime (float or None): The absolute time of the highest tip in the collapsed clade, assigned in `collapseSubtree()`.
    width (float): The width of the node for plotting, default is 1.

    Docstring generated with ChatGPT 4o.
    """

    def __init__(self, givenName: str):
        super().__init__("leaf")
        self.name = givenName  ## the pretend tip name for the clade
        self.subtree: list[BranchType] = []  ## subtree will contain all the branches that were collapsed
        self.leaves: set[str] = set()
        self.lastHeight: float | None = None  ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime: float | None = (
            None  ## refers to the absolute time of the highest tip in the collapsed clade
        )
        self.width = 1


class node(Branch):  ## node class
    """
    Represents a node in a phylogenetic tree.

    Attributes:
    branchType (str): The type of branch, default is 'node'.
    length (float): The length of the branch, assigned in `make_tree()`.
    height (float): The height of the branch, assigned in `traverse_tree()`.
    absoluteTime (float or None): The branch endpoint in absolute time, assigned in `setAbsoluteTime()`.
    parent (node): The parent node, assigned in `make_tree()`.
    children (list): A list of descendant branches of this node, assigned in `make_tree()`.
    traits (dict): Dictionary containing annotations from the tree string, assigned in `make_tree()`.
    index (int): The index of the character designating this object in the tree string, a unique identifier.
    childHeight (float or None): The height of the youngest (last) descendant tip of this node, assigned in `traverse_tree()`.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.
    leaves (set): A set of tips that are descended from this node, assigned in `traverse_tree()`.

    Docstring generated with ChatGPT 4o.
    """

    yRange: tuple[float, float]
    """maximum extent of children's y coordinates"""

    def __init__(self):
        super().__init__("node")
        self.children: list[BranchType] = []  ## a list of descendent branches of this node
        self.childHeight: float | None = None  ## the youngest descendant tip of this node
        self.leaves: set[str] = set()  ## is a set of tips that are descended from it


class leaf(Branch):  ## leaf class
    """
    Represents a leaf in a phylogenetic tree.

    Attributes:
    branchType (str): The type of branch, default is 'leaf'.
    name (str): The name of the tip after translation, assigned in `make_tree()`.
    index (int): The index of the character that defines this object in the tree string, assigned in `make_tree()`.
    length (float): The length of the branch, assigned in `make_tree()`.
    absoluteTime (float or None): The position of the tip in absolute time, assigned in `setAbsoluteTime()`.
    height (float): The height of the tip, assigned in `traverse_tree()`.
    parent (node): The parent node, assigned in `make_tree()`.
    traits (dict): Dictionary containing traits associated with the leaf, assigned in `make_tree()`.
    x (float or None): The x-coordinate for plotting, default is None.
    y (float or None): The y-coordinate for plotting, default is None.

    Docstring generated with ChatGPT 4o.
    """

    def __init__(self):
        super().__init__("leaf")

    @initialized_property
    def name(self) -> str:
        """name of tip after translation, since BEAST trees will generally have numbers for taxa but will provide a map at the beginning of the file"""
        ...


BranchType = reticulation | clade | node | leaf
LeafLike = reticulation | clade | leaf


def is_node(obj: BranchType) -> TypeIs[node]:
    return obj.is_node()


def is_leaf(obj: BranchType) -> TypeIs[leaf]:
    return obj.is_leaf()


def is_leaflike(obj: BranchType) -> TypeIs[LeafLike]:
    return obj.is_leaflike()


class tree:  ## tree class
    """
    Represents a phylogenetic tree.

    Attributes:
    cur_node (node): The current node in the tree, initialized as a new instance of the node class when building the tree in `make_tree()`.
    root (node): The root of the tree.
    Objects (list): A flat list of all branches (nodes, leaves, reticulations) in the tree.
    tipMap (dict or None): A mapping of tip numbers to names used when importing trees in NEXUS format, assigned by `loadNexus()`.
    treeHeight (float): The height of the tree, defined as the distance between the root and the most recent tip.
    mostRecent (float or None): The age of the most recent node in the tree.
    ySpan (float): The vertical span of the tree for plotting.

    Docstring generated with ChatGPT 4o.
    """

    cur_node: node | leaf | reticulation | None
    Objects: list[BranchType]
    tipMap: dict[str, str] | None
    treeHeight: float
    mostRecent: float | None
    ySpan: float

    def __init__(self):
        """
        Initializes a new tree instance.

        Docstring generated with ChatGPT 4o.
        """
        self.cur_node = node()  ## current node is a new instance of a node class
        self.cur_node.index = "Root"  ## first object in the tree is the root to which the rest gets attached
        self.cur_node.length = 0.0  ## startind node branch length is 0
        self.cur_node.height = 0.0  ## starting node height is 0
        self.root = None  # self.cur_node ## root of the tree is current node
        self.Objects = []  ## tree objects have a flat list of all branches in them
        self.tipMap = None
        self.treeHeight = 0  ## tree height is the distance between the root and the most recent tip
        self.mostRecent = None
        self.ySpan = 0.0

    def add_reticulation(self, name: str):
        """
        Adds a reticulate branch to the tree.

        Parameters:
        name (str): The name of the reticulate branch.

        Docstring generated with ChatGPT 4o.
        """
        ret = reticulation(name)
        ret.index = name
        ret.parent = self.cur_node
        assert self.cur_node is not None and is_node(self.cur_node)
        self.cur_node.children.append(ret)
        self.Objects.append(ret)
        self.cur_node = ret

    def add_node(self, i: int):
        """
        Attaches a new node to the current node.

        Parameters:
        i (int): The index of the new node, representing its position along the tree string.

        Raises:
        AssertionError: If the current node is not a node (i.e., if it is a leaf or another non-node object).

        Docstring generated with ChatGPT 4o.
        """
        new_node = node()  ## new node instance
        new_node.index = i  ## new node's index is the position along the tree string
        if self.root is None:
            self.root = new_node
            self.root.length = 0.0

        new_node.parent = self.cur_node  ## new node's parent is current node
        assert self.cur_node is not None and is_node(self.cur_node), (
            "Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses or commas."
        )
        self.cur_node.children.append(new_node)  ## new node is a child of current node
        self.cur_node = new_node  ## current node is now new node
        self.Objects.append(self.cur_node)  ## add new node to list of objects in the tree

    def add_leaf(self, i: int, name: str):
        """
        Attaches a new leaf (tip) to the current node.

        Parameters:
        i (int): The index of the new leaf, representing its position along the tree string.
        name (str): The name of the new leaf.

        Docstring generated with ChatGPT 4o.
        """
        new_leaf = leaf()  ## new instance of leaf object
        new_leaf.index = i  ## index is position along tree string
        if self.root is None:
            self.root = new_leaf

        new_leaf.parent = self.cur_node  ## leaf's parent is current node
        assert self.cur_node is not None and is_node(self.cur_node), (
            "Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses."
        )
        self.cur_node.children.append(new_leaf)  ## assign leaf to parent's children
        new_leaf.name = name
        self.cur_node = new_leaf  ## current node is now new leaf
        self.Objects.append(self.cur_node)  ## add leaf to all objects in the tree

    def subtree(
        self,
        starting_node: node | None = None,
        traverse_condition: Callable[..., bool] | None = None,
        stem: bool = True,
    ):
        """
        Generate a subtree (as a baltic tree object) from a tree traversal starting from a provided node.

        Parameters:
        k (node or None): The starting branch for traversal. Default is None, which means the traversal starts from the root.
        traverse_condition (function or None): A function that determines whether a child branch should be visited.
                                                Default is None, which means all branches are visited.
        stem (bool): If True, includes the stem leading to the root of the subtree. Default is True.

        Returns:
        tree or None: A new baltic tree instance representing the subtree.
                      Returns None if the subtree is empty or contains no leaves.

        Notes:
        - Custom traversal functions can result in multitype trees.
          If this is undesired, call `singleType()` on the resulting subtree afterwards.

        Example:
        >>> subtree = original_tree.subtree(k=some_node, traverse_condition=lambda node: node.traits['some_state'] == some_node.traits['some_state']) ## will only traverse through children in the same trait state as the starting node

        Docstring generated with ChatGPT 4o.
        """

        if starting_node is None:
            assert self.root and is_node(self.root)
            starting_node = self.root
        if traverse_condition is None:
            traverse_condition = always_true

        node_: node = starting_node.parent if stem else starting_node  ## move up a node if we want the stem

        subtree_branches = self.traverse_tree(
            node_,
            include_condition=lambda k: True,
            traverse_condition=traverse_condition,
        )
        subtree_branches = [*copy.deepcopy(subtree_branches)]

        if stem:  ## using stem - need to prune subtrees from root now
            unwanted_branches = []
            for child in node_.children:  ## iterate over parent's children
                if child.index != starting_node.index:  ## not at focal branch (unwanted sibling)
                    unwanted_branches += self.traverse_tree(
                        child, include_condition=always_true
                    )  ## add all branches resulting from traversals of unwanted siblings

            remove = [
                k for k in subtree_branches if k.index in [ub.index for ub in unwanted_branches]
            ]  ## iterate over subtree branches, remember those that belong to unwanted subtrees
            for r in remove:  ## iterate over branches belong to unwanted subtrees
                subtree_branches.remove(r)  ## remove from list

        if subtree_branches is None or not any(
            node.is_leaf() for node in subtree_branches
        ):  ## nothing found or no leaf objects in traversal
            return None

        local_tree = tree()  ## create a new tree object where the subtree will be
        local_tree.Objects = subtree_branches  ## assign branches to new tree object

        local_tree.root = subtree_branches[0]  ## root is the beginning of the traversal
        local_tree.root.parent = None  ## tree begins strictly at node

        subtree_set = set(subtree_branches)  ## turn branches into set for quicker look up later

        if traverse_condition is not None:
            ## didn't use default traverse condition, might need to deal with hanging nodes and prune children
            for nd in local_tree.getInternal():  ## iterate over nodes
                nd.children = [
                    child for child in nd.children if child in subtree_set
                ]  ## only keep children seen in traversal
            local_tree.fixHangingNodes()

        if self.tipMap:  ## if original tree has a tipMap dictionary
            names = {w.name for w in local_tree.getExternal()}
            local_tree.tipMap = {
                k: v for k, v in self.tipMap.items() if v in names
            }  ## copy over the relevant tip translations

        return local_tree

    def singleType(self):
        """
        Removes any branches with a single child (multitype nodes) from the tree.

        This method simplifies the tree by removing nodes that have only one child, effectively
        merging these nodes with their single child to ensure each node has either no children (leaves)
        or multiple children.

        The process involves:
        - Identifying nodes with a single child.
        - Reassigning the child node to the grandparent of the original single child node.
        - Adjusting the branch lengths accordingly.
        - Removing the single child node from the tree.

        This process is repeated until no multitype nodes remain.

        Returns:
        None (modifies the baltic tree object in-place)

        Docstring generated with ChatGPT 4o.
        """
        while True:
            multitype_nodes = self.getInternal(lambda k: len(k.children) == 1)

            if not multitype_nodes:
                break

            for k in sorted(multitype_nodes, key=lambda x: -x.height):
                child = k.children[0]  ## fetch child
                grandparent = k.parent if k.parent.index else self.root  ## fetch grandparent
                assert isinstance(grandparent, node)

                child.parent = grandparent  ## child's parent is now grandparent
                grandparent.children.append(child)  ## add child to grandparent's children
                grandparent.children.remove(k)  ## remove old parent from grandparent's children
                grandparent.children = list(set(grandparent.children))

                child.length += k.length  ## adjust child length

                self.Objects.remove(k)  ## remove old parent from all objects
        self.sortBranches()

    def setAbsoluteTime(self, date: float):
        """
        Places all objects in absolute time by providing the date of the most recent tip.

        Parameters:
        date (float): The date of the most recent tip in the tree in decimal date format (see the `decimalDate()` function).

        This method calculates the absolute time for each object in the tree based on its height and
        the provided date of the most recent tip (most recent date - tree height + node/leaf height). The absolute time is stored in the `absoluteTime`
        attribute of each object.

        Example:
        >>> tree.setAbsoluteTime(2023.0)

        Docstring generated with ChatGPT 4o.
        """
        for k in self.Objects:  ## iterate over all objects
            k.absoluteTime = date - self.treeHeight + k.height  ## heights are in units of time from the root
        self.mostRecent = max(k.absoluteTime or 0 for k in self.Objects)

    def treeStats(self):
        """
        Provides information about the tree.

        This method traverses the tree to gather and print various statistics, including:
        - Tree height and length.
        - Whether the tree is strictly bifurcating, multitype, or a singleton (comprised of a single tip).
        - The presence of annotations (traits).
        - The number of nodes and leaves.

        The method prints the collected statistics to the console.

        Example:
        >>> tree.treeStats()

        Docstring generated with ChatGPT 4o.
        """
        self.traverse_tree()  ## traverse the tree
        obs = self.Objects  ## convenient list of all objects in the tree
        print(
            "\nTree height: %.6f\nTree length: %.6f" % (self.treeHeight, sum([x.length for x in obs]))
        )  ## report the height and length of tree

        nodes = self.getInternal()  ## get all nodes
        strictlyBifurcating = all(len(x.children) == 2 for x in nodes)  ## assume tree is not strictly bifurcating
        multiType = any(len(x.children) == 1 for x in nodes)
        singleton = len(nodes) == 0

        hasTraits = False  ## assume tree has no annotations
        max(len(k.traits) for k in obs)  ## check the largest number of annotations any branch has

        if strictlyBifurcating:
            print("strictly bifurcating tree")  ## report
        if multiType:
            print("multitype tree")  ## report
        if singleton:
            print("singleton tree")
        if hasTraits:
            print("annotations present")  ## report

        print(
            "\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n"
            % (len(obs), len(nodes), len(self.getExternal()))
        )  ## report numbers of different objects in the tree

    def traverse_tree(
        self,
        cur_node: BranchType | None = None,
        include_condition: Callable[[BranchType], bool] | None = None,
        traverse_condition: Callable[[BranchType], bool] | None = None,
        collect: list | None = None,
        verbose: bool = False,
    ) -> list[BranchType]:
        """
        Traverses the tree starting from a specified node and collects nodes based on conditions.

        Parameters:
        cur_node (node or None): The starting node for traversal. If None, starts from the root. Default is None.
        include_condition (function or None): A function that determines whether a node should be included in the `collect` list.
                                              Default includes leaf-like nodes.
        traverse_condition (function or None): A function that determines whether a child node should be traversed.
                                               Default traverses all nodes.
        collect (list): A list to collect nodes that meet the include condition. Default is None, which collects and returns only leaf and leaf-like (`clade` and `reticulation`) objects.
        verbose (bool): If True, prints verbose output during traversal. Default is False.

        Returns:
        list: A list of nodes that pass `include_condition` that were encountered during the traversal that meets the `traverse_condition`.

        Example:
        >>> descendant_tip_list = tree.traverse_tree(include_condition=lambda node: node.is_leaf()) ## will return a list of all leaf (but not leaf-like) objects encountered during the traversal

        Docstring generated with ChatGPT 4o.
        """
        if cur_node is None:  ## if no starting point defined - start from root
            if verbose:
                print("Initiated traversal from root")
            assert self.root is not None
            cur_node = self.root

            if traverse_condition is None and include_condition is None:  ## reset heights if traversing from scratch
                for k in self.Objects:  ## reset various parameters
                    if is_node(k):
                        k.leaves = set()
                        k.childHeight = None
                    k.height = None

        if include_condition is None:
            include_condition = is_leaflike

        if traverse_condition is None:
            traverse_condition = always_true

        ## initiate collect list if not initiated
        # assumes the side effect of passing collect along
        if collect is None:
            collect = []

        assert cur_node is not None
        if cur_node.parent and cur_node._height is None:
            ## cur_node has a parent - set height if it doesn't have it already
            cur_node.height = cur_node.length + cur_node.parent.height
        elif cur_node._height is None:
            ## cur_node does not have a parent (root), if height not set before it's zero
            cur_node.height = 0.0

        if verbose:
            print("at %s (%s)" % (cur_node.index, cur_node.branchType))

        if include_condition(cur_node):  ## test if interested in cur_node
            collect.append(cur_node)  ## add to collect list for reporting later

        if is_leaf(cur_node) and self.root != cur_node:  ## cur_node is a tip (and tree is not single tip)
            cur_node.parent.leaves.add(cur_node.name)  ## add to parent's list of tips

        elif is_node(cur_node):  ## cur_node is node
            for child in filter(
                traverse_condition, cur_node.children
            ):  ## only traverse through children we're interested
                if verbose:
                    print("visiting child %s" % (child.index))
                self.traverse_tree(
                    cur_node=child,
                    include_condition=include_condition,
                    traverse_condition=traverse_condition,
                    verbose=verbose,
                    collect=collect,
                )  ## recurse through children
                if verbose:
                    print("child %s done" % (child.index))
            assert len(cur_node.children) > 0, "Tried traversing through hanging node without children. Index: %s" % (
                cur_node.index
            )
            cur_node.childHeight = max(
                [child.childHeight or 0 if is_node(child) else child.height for child in cur_node.children]
            )

            if cur_node.parent:
                cur_node.parent.leaves = cur_node.parent.leaves.union(
                    cur_node.leaves
                )  ## pass tips seen during traversal to parent

            assert cur_node.childHeight is not None
            self.treeHeight = cur_node.childHeight  ## it's the highest child of the starting node
        return collect

    def renameTips(self, d: dict[str, str] | None = None):
        """
        Rename each tip using a dictionary.

        Parameters:
        d (dict or None): A dictionary mapping original tip names to new tip names.
                          If None, uses the tree's `tipMap` attribute if it exists (see `loadNexus()`). Default is to assume tipMap exists.

        This method iterates over all leaf objects in the tree and updates their names based on the provided dictionary.

        Example:
        >>> tree.renameTips({'tip1': 'new_name1', 'tip2': 'new_name2'})

        Docstring generated with ChatGPT 4o.
        """
        if d is None and self.tipMap is not None:
            d = self.tipMap

        if d is None:
            raise ValueError("No dictionary provided for renaming tips.")

        for k in self.getExternal():
            k.name = d.get(k.name, k.name)

    def sortBranches(
        self,
        descending: bool = True,
        sort_function: Callable[[BranchType], Any] | None = None,
        sortByHeight: bool = True,
    ):
        """
        Sort descendants of each node.

        Parameters:
        descending (bool): If True, sorts in descending order. Default is True.
        sort_function (function or None): A custom sorting function. Default sorts nodes by number of descendants and length.
        sortByHeight (bool): If True, sorts nodes by height and groups nodes and leaves together. Default is True.

        This method sorts the children of each internal node in the tree according to the specified sorting function
        and order. It then updates the x and y positions of each branch by calling `drawTree()`.

        Example:
        >>> tree.sortBranches(descending=False)

        Docstring generated with ChatGPT 4o.
        """

        mod = -1 if descending else 1
        if sort_function is None:

            def func(k):
                return (
                    (is_node(k), -len(k.leaves) * mod, k.length * mod) if is_node(k) else (is_node(k), k.length * mod)
                )

            sort_function = func

        if sortByHeight:  # Sort nodes by height and group nodes and leaves together
            """ Sort descendants of each node. """

            for k in self.getInternal():  ## iterate over nodes
                k.children = sorted(k.children, key=sort_function)
        else:  # Do not sort by height. Retain leaves at original positions. Only sort nodes
            for k in self.getInternal():
                leavesIdx = [
                    (i, ctr) for ctr, i in enumerate(k.children) if i.is_leaflike()
                ]  # Get original indices of leaves
                nodes = sorted(
                    [x for x in k.children if x.is_node()], key=sort_function
                )  # Sort nodes only by number of descendants
                children = nodes
                for i in leavesIdx:  # Insert leaves back into same positions
                    children.insert(i[1], i[0])
                k.children = children
        self.drawTree()  ## update x and y positions of each branch, since y positions will have changed because of sorting

    def drawTree(
        self,
        order: list[LeafLike] | None = None,
        width_function: Callable[..., float] | None = None,
        pad_nodes: dict | None = None,
        verbose: bool = False,
    ):
        """
        Assign x and y coordinates of each branch in the tree.

        Parameters:
        order (list or None): A list of tips recovered from a tree traversal to ensure they are plotted in the correct order along the vertical tree dimension.
                              If None, performs a pre-order traversal. Default is None.
        width_function (function or None): A function to determine the width of each branch.
                                           If None, uses default widths (1 unit for leaf objects, width + 1 for clades). Default is None.
        pad_nodes (dict or None): A dictionary specifying nodes to be padded with extra space around their descendants (`node` class : float or int). Default is None.
        verbose (bool): If True, prints verbose output during drawing. Default is False.

        This method assigns branch height as x and calculates the y coordinates (with adjustments, if any) for each branch in the tree for plotting purposes.

        Example:
        >>> tree.drawTree()

        Docstring generated with ChatGPT 4o.
        """
        if order is None:
            ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
            order = [k for k in self.traverse_tree() if is_leaflike(k)]
            if verbose:
                print("Drawing tree in pre-order")
        else:
            if verbose:
                print("Drawing tree with provided order")

        name_order = {x.name: i for i, x in enumerate(order)}
        assert len(name_order) == len(order), "Non-unique names present in tree"
        if width_function is None:
            if verbose:
                print("Drawing tree with default widths (1 unit for leaf objects, width+1 for clades)")
            skips = [1 if is_leaf(x) else x.width + 1 for x in order]
        else:
            skips = [width_function(k) for k in order]

        for k in self.Objects:  ## reset coordinates for all objects
            k.x = None
            k.y = None

        drawn = {}  ## drawn keeps track of what's been drawn
        for y_idx, k in enumerate(order):  ## iterate over tips
            x = k.height  ## x position is height
            y = sum(skips[y_idx:]) - skips[y_idx] / 2.0  ## sum across skips to find y position

            k.x = x  ## set x and y coordinates
            k.y = y
            drawn[k.index] = None  ## remember that this objects has been drawn

        if pad_nodes is not None:  ## will be padding nodes
            for n in pad_nodes:  ## iterate over nodes whose descendants will be padded
                idx = (
                    sorted([name_order[lf] for lf in n.leaves]) if n.is_node() else [order.index(n)]
                )  ## indices of all tips to be padded
                for i, k in enumerate(order):  ## iterate over all tips
                    if i < idx[0]:  ## tip below clade
                        k.y += pad_nodes[n]  ## pad

                    if (i - 1) < idx[-1]:  ## tip above clade
                        k.y += pad_nodes[n]  ## pad again

            all_ys = filter(None, self.getParameter("y"))  ## get all y positions in tree that aren't None
            minY = min(all_ys)  ## get min
            for k in self.getExternal():  ## reset y positions so tree starts at y=0.5
                k.y -= minY - 0.5

        assert len([k for k in self.Objects if is_leaflike(k)]) == len(order), (
            "Number of tips in tree does not match number of unique tips, check if two or more collapsed clades were assigned the same name."
        )
        storePlotted = 0

        while len(drawn) != len(self.Objects):  # keep drawing the tree until everything is drawn
            if verbose:
                print("Drawing iteration %d" % (len(drawn)))
            ## iterate through internal nodes that have not been drawn
            for k in self.getInternal(lambda w: w.index not in drawn):
                children_y_coords = [
                    q.y for q in k.children if q.y is not None
                ]  ## get all existing y coordinates of the node
                if len(children_y_coords) == len(k.children):  ## all y coordinates of children known
                    if verbose:
                        print("Setting node %s coordinates to" % (k.index))

                    k.x = k.height  ## x position is height
                    ## internal branch is in the middle of the vertical bar
                    k.y = mean(children_y_coords)
                    drawn[k.index] = None  ## remember that this objects has been drawn
                    if verbose:
                        print("%s (%s branches drawn)" % (k.y, len(drawn)))
                    minYrange = min(
                        [min(getattr(child, "yRange", [child.y])) for child in k.children]
                    )  ## get lowest y coordinate across children
                    maxYrange = max(
                        [max(getattr(child, "yRange", [child.y])) for child in k.children]
                    )  ## get highest y coordinate across children
                    k.yRange = (minYrange, maxYrange)  ## assign the maximum extent of children's y coordinates

            if len(self.Objects) > len(drawn):
                assert len(drawn) > storePlotted, (
                    "Got stuck trying to find y positions of objects (%d branches drawn this iteration, %d branches during previous iteration out of %d total)"
                    % (len(drawn), storePlotted, len(self.Objects))
                )
            storePlotted = len(drawn)  ## remember how many branches were drawn this iteration

        yvalues = [k.y or 0 for k in self.Objects]  ## all y values
        self.ySpan = max(yvalues) - min(yvalues) + min(yvalues) * 2  ## determine appropriate y axis span of tree

        assert self.root is not None
        if is_node(self.root):
            self.root.x = min(
                [q.x - q.length for q in self.root.children if q.x is not None]
            )  ## set root x and y coordinates
            children_y_coords = [q.y for q in self.root.children if q.y is not None]
            self.root.y = sum(children_y_coords) / float(len(children_y_coords))
        else:
            self.root.x = self.root.length

    def drawUnrooted(self, rotate: float = 0.0, n: BranchType | None = None, total: float | None = None):
        """
        Calculate x and y coordinates of each branch in an unrooted arrangement.

        This method arranges the branches of the tree in an unrooted, circular layout.
        The coordinates are calculated recursively for each node.

        Parameters:
        rotate (float): The initial rotation angle in radians. Default is 0.0.
        n (node or None): The current node being processed. If None, starts from the root. Default is None.
        total (int or None): The total number of tips or the sum of widths for clades. Default is None.

        Code translated from https://github.com/nextstrain/auspice/commit/fc50bbf5e1d09908be2209450c6c3264f298e98c, written by Richard Neher.

        Example:
        >>> tree.drawUnrooted(rotate=0.1)

        Docstring generated with ChatGPT 4o.
        """

        if n is None:
            total = sum([1 if is_leaf(x) else x.width + 1 for x in self.getExternal()])
            assert self.root is not None and is_node(self.root)
            n = self.root  # .children[0]
            for k in self.Objects:
                k.traits["tau"] = 2 * math.pi * rotate
                k.x = 0.0
                k.y = 0.0

        assert total is not None

        def width(n):
            return 2 * math.pi * (len(n.leaves) if isinstance(n, clade) else 1) / total

        if n.parent.x is None:
            n.parent.x = 0.0
            n.parent.y = 0.0

        n.x = n.parent.x + n.length * math.cos(n.traits["tau"] + width(n) * 0.5)
        n.y = n.parent.y + n.length * math.sin(n.traits["tau"] + width(n) * 0.5)
        eta = n.traits["tau"]

        if is_node(n):
            for ch in n.children:
                ch.traits["tau"] = eta
                eta += width(ch)
                self.drawUnrooted(rotate, ch, total)

    def commonAncestor(self, descendants):
        """
        Find the most recent node object that gave rise to a given list of descendant branches.

        Parameters:
        descendants (list): A list of descendant branches (as `node`, `leaf`, `clade` and/or `reticulation` classes) for which to find the most recent common ancestor.

        Returns:
        node: The most recent common ancestor node.

        Raises:
        AssertionError: If the number of descendants is less than 2.

        Example:
        >>> ancestor = tree.commonAncestor([descendant1, descendant2])

        Docstring generated with ChatGPT 4o.
        """
        assert len(descendants) > 1, "Not enough descendants to find common ancestor: %d" % (len(descendants))
        paths_to_root = {k.index: set() for k in descendants}  ## for every descendant create an empty set
        for k in descendants:  ## iterate through every descendant
            cur_node = k  ## start descent from descendant
            while cur_node:  ## while not at root
                paths_to_root[k.index].add(cur_node)  ## remember every node visited along the way
                cur_node = cur_node.parent  ## descend

        return sorted(reduce(set.intersection, paths_to_root.values()), key=lambda k: k.height)[
            -1
        ]  ## return the most recent branch that is shared across all paths to root

    def collapseSubtree(self, cl, givenName, verbose=False, widthFunction=lambda k: len(k.leaves)):
        """
        Collapse an entire subtree into a clade object.

        Parameters:
        cl (node): The node representing the root of the subtree to collapse.
        givenName (str): The name to assign to the new clade.
        verbose (bool): If True, prints verbose output during the process. Default is False.
        widthFunction (function): A function to determine the width of the clade when computing branch y coordinates in `drawTree()`. Default calculates width based on the number of leaves.

        Returns:
        clade: The newly created clade object representing the collapsed subtree.

        Raises:
        AssertionError: If the provided branch is not a `node` class or if attempting to collapse the entire tree.

        Example:
        >>> collapsed_clade = tree.collapseSubtree(node, "new_clade")

        Docstring generated with ChatGPT 4o.
        """
        assert is_node(cl), "Cannot collapse non-node class"
        collapsedClade = clade(givenName)
        collapsedClade.index = cl.index
        collapsedClade.leaves = cl.leaves
        collapsedClade.length = cl.length
        collapsedClade.height = cl.height
        collapsedClade.parent = cl.parent
        collapsedClade.absoluteTime = cl.absoluteTime
        collapsedClade.traits = cl.traits
        collapsedClade.width = widthFunction(cl)

        if verbose:
            print("Replacing node %s (parent %s) with a clade class" % (cl.index, cl.parent.index))
        parent = cl.parent

        remove_from_tree = self.traverse_tree(cl, include_condition=always_true)
        collapsedClade.subtree = remove_from_tree
        assert len(remove_from_tree) < len(self.Objects), "Attempted collapse of entire tree"
        collapsedClade.lastHeight = max([x.height for x in remove_from_tree])
        collapsedClade.lastAbsoluteTime = max(
            [x.absoluteTime for x in remove_from_tree if x.absoluteTime is not None], default=None
        )

        for k in remove_from_tree:
            self.Objects.remove(k)

        parent.children.remove(cl)
        parent.children.append(collapsedClade)
        self.Objects.append(collapsedClade)
        collapsedClade.parent = parent
        if self.tipMap is not None:
            self.tipMap[givenName] = givenName

        self.traverse_tree()
        self.sortBranches()
        return collapsedClade

    def uncollapseSubtree(self):
        """
        Uncollapse all collapsed subtrees in the tree.

        This method restores all previously collapsed clades back to their original subtree structures.
        It iterates through all objects in the tree, identifies clades, and replaces each clade with its
        corresponding subtree that was stored in the `clade` class.

        Example:
        >>> tree.uncollapseSubtree()

        Docstring generated with ChatGPT 4o.
        """
        while len([k for k in self.Objects if isinstance(k, clade)]) > 0:
            clades = [k for k in self.Objects if isinstance(k, clade)]
            for cl in clades:
                parent = cl.parent
                subtree = cl.subtree
                parent.children.remove(cl)
                parent.children.append(subtree[0])
                self.Objects += subtree
                self.Objects.remove(cl)
                if self.tipMap is not None:
                    self.tipMap.pop(cl.name, None)
        self.traverse_tree()

    def collapseBranches(
        self,
        collapseIf=lambda x: x.traits["posterior"] <= 0.5,
        designated_nodes: list[node] = [],
        verbose=False,
    ):
        """
        Collapse all branches that satisfy a function `collapseIf` (default is an anonymous function that returns true if posterior probability is <= 0.5).
        Alternatively, a list of nodes can be supplied to the script.
        A branch designated for deletion gets its descendants assigned to its parent with branch lengths adjusted accordingly before being pruned out of the tree.

        Parameters:
        collapseIf (function): A function that determines whether a branch should be collapsed. Default is a function that returns True if the posterior probability is <= 0.5.
        designated_nodes (list): A list of nodes to be collapsed. If empty, the collapseIf function is used to determine nodes to collapse. Default is an empty list.
        verbose (bool): If True, prints verbose output during the process. Default is False.

        Returns:
        tree: A deep copied version of the tree with the specified branches collapsed.

        Raises:
        AssertionError: If non-node classes are detected in the designated_nodes list or if the root node is designated for deletion.

        Example:
        >>> new_tree = tree.collapseBranches()

        Docstring generated with ChatGPT 4o.
        """
        newTree = copy.deepcopy(self)  ## work on a copy of the tree
        ## no nodes were designated for deletion - relying on anonymous function to collapse nodes
        if len(designated_nodes) == 0:
            ## fetch a list of all nodes who are not the root and who satisfy the condition
            nodes_to_delete = [n for n in newTree.Objects if is_node(n) and collapseIf(n) and n != newTree.root]
        else:
            assert len([w for w in designated_nodes if w.is_node()]) == len(designated_nodes), (
                "Non-node class detected in list of nodes designated for deletion"
            )
            assert len([w for w in designated_nodes if w != newTree.root]) == 0, "Root node was designated for deletion"

            ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
            nodes_to_delete = [w for w in newTree.Objects if w.index in [q.index for q in designated_nodes]]

        if verbose:
            print("%s nodes set for collapsing: %s" % (len(nodes_to_delete), [w.index for w in nodes_to_delete]))
        assert len(nodes_to_delete) < len(newTree.getInternal()) - 1, "Chosen cutoff would remove all branches"
        while len(nodes_to_delete) > 0:  ## as long as there are branches to be collapsed - keep reducing the tree
            if verbose:
                print("Continuing collapse cycle, %s nodes left" % (len(nodes_to_delete)))
            for k in sorted(nodes_to_delete, key=lambda x: -x.height):  ## start with branches near the tips
                assert is_node(k)
                zero_node = k.children  ## fetch the node's children
                k.parent.children += zero_node  ## add them to the zero node's parent
                old_parent = k  ## node to be deleted is the old parent
                new_parent = (
                    k.parent
                )  ## once node is deleted, the parent to all their children will be the parent of the deleted node
                if new_parent is None:
                    assert self.root is not None
                    new_parent = self.root
                if verbose:
                    print(
                        "Removing node %s, attaching children %s to node %s"
                        % (
                            old_parent.index,
                            [w.index for w in k.children],
                            new_parent.index,
                        )
                    )
                for w in (
                    newTree.Objects
                ):  ## assign the parent of deleted node as the parent to any children of deleted node
                    if w.parent == old_parent:
                        w.parent = new_parent
                        w.length += old_parent.length
                        if verbose:
                            print("Fixing branch length for node %s" % (w.index))
                k.parent.children.remove(
                    k
                )  ## remove traces of deleted node - it doesn't exist as a child, doesn't exist in the tree and doesn't exist in the nodes list
                newTree.Objects.remove(k)

                nodes_to_delete.remove(k)  ## in fact, the node never existed

                if len(designated_nodes) == 0:
                    nodes_to_delete = [n for n in newTree.Objects if is_node(n) and collapseIf(n) and n != newTree.root]
                else:
                    assert len([w for w in designated_nodes if w.is_node()]) == len(designated_nodes), (
                        "Non-node class detected in list of nodes designated for deletion"
                    )
                    assert len([w for w in designated_nodes if w != newTree.root]) == 0, (
                        "Root node was designated for deletion"
                    )
                    nodes_to_delete = [w for w in newTree.Objects if w.index in [q.index for q in designated_nodes]]

                if verbose:
                    print("Removing references to node %s" % (k.index))
        newTree.sortBranches()  ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree  ## return collapsed tree

    def toString(
        self,
        cur_node=None,
        traits=None,
        verbose=False,
        nexus=False,
        string_fragment=None,
        traverse_condition=None,
        rename=None,
        quotechar="'",
        json=False,
    ):
        """
        Output the topology of the tree with branch lengths and comments to a string.

        Parameters:
        cur_node (node or None): The starting point for traversal. Default is None, which starts at the root.
        traits (list or None): A list of keys to output entries in the traits dictionary of each branch. Default is all available traits.
        verbose (bool): If True, prints verbose output during the process. Default is False.
        nexus (bool): If True, outputs in NEXUS format. Default is False, which outputs in Newick format.
        string_fragment (list or None): A list of characters that comprise the tree string. Default is None.
        traverse_condition (function or None): A function that determines whether a child branch should be traversed. Default is None which traverses all children.
        rename (dict or None): A dictionary to rename tip names. Default is None.
        quotechar (str): The character to use for quoting tip names. Default is "'".
        json (bool): If True, outputs in auspice JSON format (somewhat experimental). Default is False.

        Returns:
        str: The tree string in the specified format.

        Example:
        >>> tree_string = tree.toString()

        Docstring generated with ChatGPT 4o.
        """
        if cur_node is None:
            assert self.root is not None
            cur_node = self.root  # .children[-1]
        if traits is None:
            traits = set(sum([list(k.traits.keys()) for k in self.Objects], []))  ## fetch all trait keys
        if string_fragment is None:
            string_fragment = []
            if nexus:
                assert not json, "Nexus format not a valid option for JSON output"
                if verbose:
                    print("Exporting to NEXUS format")
                string_fragment.append("#NEXUS\nBegin trees;\ntree TREE1 = [&R] ")
        if traverse_condition is None:
            traverse_condition = always_true

        comment = []  ## will hold comment
        if len(traits) > 0:  ## non-empty list of traits to output
            for tr in traits:  ## iterate through keys
                if tr in cur_node.traits:  ## if key is available
                    if verbose:
                        print(
                            "trait %s available for %s (%s) type: %s"
                            % (
                                tr,
                                cur_node.index,
                                cur_node.branchType,
                                type(cur_node.traits[tr]),
                            )
                        )
                    if isinstance(cur_node.traits[tr], str):  ## string value
                        comment.append('%s="%s"' % (tr, cur_node.traits[tr]))
                        if verbose:
                            print("adding string comment %s" % (comment[-1]))
                    elif isinstance(cur_node.traits[tr], float) or isinstance(
                        cur_node.traits[tr], int
                    ):  ## float or integer
                        comment.append("%s=%s" % (tr, cur_node.traits[tr]))
                        if verbose:
                            print("adding numeric comment %s" % (comment[-1]))
                    elif isinstance(cur_node.traits[tr], list):  ## lists
                        rangeComment = []
                        for val in cur_node.traits[tr]:
                            if isinstance(val, str):  ## string
                                rangeComment.append('"%s"' % (val))
                            elif isinstance(val, float) or isinstance(val, int):  ## float or integer
                                rangeComment.append("%s" % (val))
                            elif isinstance(val, list):  ## list of lists, example complete history annotated on tree
                                rangeComment.append("{{{}}}".format(",".join(val)))
                        comment.append("%s={%s}" % (tr, ",".join(rangeComment)))
                        if verbose:
                            print("adding range comment %s" % (comment[-1]))
                elif verbose:
                    print("trait %s unavailable for %s (%s)" % (tr, cur_node.index, cur_node.branchType))

        if is_node(cur_node):
            if verbose:
                print("node: %s" % (cur_node.index))
            string_fragment.append("(")
            traverseChildren = list(filter(traverse_condition, cur_node.children))
            assert len(traverseChildren) > 0, "Node %s does not have traversable children" % (cur_node.index)
            for c, child in enumerate(
                traverseChildren
            ):  ## iterate through children of node if they satisfy traverse condition
                if verbose:
                    print("moving to child %s of node %s" % (child.index, cur_node.index))
                self.toString(
                    cur_node=child,
                    traits=traits,
                    verbose=verbose,
                    nexus=nexus,
                    string_fragment=string_fragment,
                    traverse_condition=traverse_condition,
                    rename=rename,
                    quotechar=quotechar,
                )
                if (c + 1) < len(traverseChildren):  ## not done with children, add comma for next iteration
                    string_fragment.append(",")
            string_fragment.append(")")  ## last child, node terminates

        elif is_leaf(cur_node):
            if rename is None:
                treeName = cur_node.name  ## designated numName
            else:
                assert isinstance(rename, dict), 'Variable "rename" is not a dictionary'
                assert cur_node.name in rename, "Tip name %s not in rename dictionary" % (cur_node.name)
                treeName = rename[cur_node.name]

            if verbose:
                print("leaf: %s (%s)" % (cur_node.index, treeName))
            string_fragment.append("%s%s%s" % (quotechar, treeName, quotechar))

        if len(comment) > 0:
            if verbose:
                print("adding comment to %s" % (cur_node.index))
            comment = ",".join(comment)
            comment = "[&" + comment + "]"
            string_fragment.append("%s" % (comment))  ## end of node, add annotations

        if verbose:
            print("adding branch length to %s" % (cur_node.index))
        string_fragment.append(":%8f" % (cur_node.length))  ## end of node, add branch length

        if cur_node == self.root:  # .children[-1]:
            string_fragment.append(";")
            if nexus:
                string_fragment.append("\nEnd;")
            if verbose:
                print("finished")
            return "".join(string_fragment)

    def allTMRCAs(self):
        """
        Calculate the time to the most recent common ancestor (TMRCA) for all pairs of tips in the tree.

        This method creates a pairwise matrix of tips and iterates over all internal nodes to find the TMRCA (as `absoluteTime` attribute)
        for each pair of descendant tips. The matrix is symmetric, and the diagonal elements are set to 0.0
        as the TMRCA of a tip with itself is zero.

        Returns:
        dict: A dictionary where each key is a tip name and the corresponding value is another dictionary
              with tip names as keys and their TMRCA as values.

        Example:
        >>> tmrca_matrix = tree.allTMRCAs()

        Docstring generated with ChatGPT 4o.
        """
        tip_names: list[str] = [k.name for k in self.getExternal()]
        tmrcaMatrix = {
            x: {y: None if x != y else 0.0 for y in tip_names} for x in tip_names
        }  ## pairwise matrix of tips

        for k in self.getInternal():  ## iterate over nodes
            all_children = list(k.leaves)  ## fetch all descendant tips of node
            assert k.absoluteTime is not None

            for a, tipA in enumerate(all_children):
                for tipB in all_children[a + 1 :]:
                    v = tmrcaMatrix[tipA][tipB]
                    ## if node's time is more recent than previous entry - set new TMRCA value for pair of tips
                    if v is None or v <= k.absoluteTime:
                        tmrcaMatrix[tipA][tipB] = k.absoluteTime
                        tmrcaMatrix[tipB][tipA] = k.absoluteTime
        return tmrcaMatrix

    def reduceTree(self, keep: list[LeafLike], verbose: bool = False):
        """
        Reduce the tree to include only the branches tracking a specified set of tips to the root.

        Parameters:
        keep (list): A list of tip branches to retain in the reduced tree.
        verbose (bool): If True, prints verbose output during the process. Default is False.

        Returns:
        tree: A new tree object containing only the specified tips and the necessary branches to connect them to the root. Can result in a tree with multitype-like branches (nodes with a single child).

        Raises:
        AssertionError: If no tips are given to reduce the tree to, or if the list contains non-leaf-like branches.

        Example:
        >>> reduced_tree = tree.reduceTree([tip1, tip2, tip3])

        Docstring generated with ChatGPT 4o.
        """
        assert len(keep) > 0, "No tips given to reduce the tree to."
        assert len([k for k in keep if not k.is_leaflike()]) == 0, (
            "Embedding contains %d branches that are not leaf-like." % (len([k for k in keep if not k.is_leaflike()]))
        )
        if verbose:
            print("Preparing branch hash for keeping %d branches" % (len(keep)))
        branch_hash = {}
        for r in keep:
            branch_hash[r.index] = r
            # if isinstance(r, reticulation):
            #     # avoid dangling reticulations
            #     branch_hash[r.target.index] = r.target

        if verbose:
            print("Deep copying tree")
        reduced_tree = copy.deepcopy(self)  ## new tree object

        embedding = []
        for k in reduced_tree.Objects:  ## deep copy branches from current tree
            if k.index in branch_hash:  ## if branch is designated as one to keep
                cur_b = k
                if verbose:
                    print("Traversing to root from %s" % (cur_b.index))
                while cur_b != reduced_tree.root:  ## descend to root
                    if verbose:
                        print(f"at {cur_b.index} root: {cur_b == reduced_tree.root}")
                    embedding.append(cur_b)  ## keep track of the path to root
                    cur_b = cur_b.parent

        embedding.append(reduced_tree.root)  ## add root to embedding
        if verbose:
            print(
                "Finished extracting embedding with %s branches (%s tips, %s nodes)"
                % (
                    len(embedding),
                    len([w for w in embedding if is_leaf(w)]),
                    len([w for w in embedding if is_node(w)]),
                )
            )
        embedding = set(embedding)  ## prune down to only unique branches

        reduced_tree.Objects = sorted(
            list(embedding), key=lambda x: x.height
        )  ## assign branches that are kept to new tree's Objects

        if verbose:
            print("Pruning untraversed lineages")
        for k in reduced_tree.getInternal():  ## iterate through reduced tree
            k.children = [
                c for c in k.children if c in embedding
            ]  ## only keep children that are present in lineage traceback

        assert reduced_tree.root is not None and is_node(reduced_tree.root)
        reduced_tree.root.children = [c for c in reduced_tree.root.children if c in embedding]  ## do the same for root

        reduced_tree.fixHangingNodes()

        if verbose:
            print("Last traversal and branch sorting")
        reduced_tree.traverse_tree()  ## traverse
        reduced_tree.sortBranches()  ## sort

        return reduced_tree  ## return new tree

    def countLineages(self, t: float, attr: str = "absoluteTime", condition: Callable[..., bool] = always_true):
        """
        Count the number of lineages present at a specific time point.

        Parameters:
        t (float): The time point at which to count the lineages.
        attr (str): The attribute used to determine the time of the nodes. Default is `absoluteTime`.
        condition (function): A function that determines whether a lineage should be included in the count. Default is a function that always returns True.

        Returns:
        int: The number of lineages present at the specified time point (branches whose time is above and parent is below the time point provided).

        Example:
        >>> num_lineages = tree.countLineages(2020.5)

        Docstring generated with ChatGPT 4o.
        """
        return len(
            [
                k
                for k in self.Objects
                if getattr(k.parent, attr) is not None
                and getattr(k.parent, attr) < t <= getattr(k, attr)
                and condition(k)
            ]
        )

    def getExternal(self, secondFilter: Callable[[LeafLike], bool] | None = None):
        """
        Get all leaf-like branches (`leaf`, `clade`, and `reticulation` classes).

        Parameters:
        secondFilter (function or None): An optional function to further filter the leaf branches based on an additional property. Default is None.

        Returns:
        list: A list of leaf branches that optionally satisfy the secondFilter condition.

        Example:
        >>> leaves = tree.getExternal()
        >>> filtered_leaves = tree.getExternal(lambda x: x.absoluteTime >= 2023.0)

        Docstring generated with ChatGPT 4o.
        """
        if secondFilter is None:
            secondFilter = always_true
        externals = [k for k in self.Objects if is_leaflike(k) and secondFilter(k)]
        return externals

    def getInternal(self, secondFilter: Callable[[node], bool] | None = None):
        """
        Get all branches belonging to the `node` class.

        Parameters:
        secondFilter (function or None): An optional function to further filter the internal nodes based on an additional property. Default is None.

        Returns:
        list: A list of node branches that optionally satisfy the secondFilter condition.

        Example:
        >>> nodes = tree.getInternal()
        >>> filtered_nodes = tree.getInternal(lambda x: x.absoluteTime >= 2023.0)

        Docstring generated with ChatGPT 4o.
        """
        if secondFilter is None:
            secondFilter = always_true

        internals = [k for k in self.Objects if is_node(k) and secondFilter(k)]
        return internals

    def getBranches(self, attrs: Callable[[BranchType], bool] = always_true, warn: bool = True):
        """
        Get branches that satisfy a specified condition.

        Parameters:
        attrs (function): A function that determines whether a branch should be included. Default is a function that always returns True.
        warn (bool): If True, raises an exception if no branches satisfying the condition are found. Default is True.

        Returns:
        list or object: A list of branches that satisfy the condition (list is empty if warn is False). If only one branch satisfies the condition, returns that branch.

        Raises:
        Exception: If no branches satisfying the condition are found and warn is True.

        Example:
        >>> branches = tree.getBranches(lambda x: x.length > 0.5)
        >>> single_branch = tree.getBranches(lambda x: x.index == 'node1', warn=False)

        Docstring generated with ChatGPT 4o.
        """
        select = [k for k in self.Objects if attrs(k)]

        if len(select) == 0 and warn:
            raise Exception("No branches satisfying function were found amongst branches")
        elif len(select) == 0 and not warn:
            return []
        elif len(select) == 1:
            return select[-1]
        else:
            return select

    def getParameter(self, statistic: str, use_trait: bool = False, which: Callable[[BranchType], bool] = always_true):
        """
        Return a list of either branch trait or attribute states across branches.

        Parameters:
        statistic (str): The name of the trait or attribute to retrieve.
        use_trait (bool): If True, retrieves the trait from the branch's traits dictionary. If False, retrieves the attribute directly from branch attributes. Default is False (retrieves attributes).
        which (function): A function that determines which branches to include. Default includes all branches in the tree.

        Returns:
        list: A list of values for the specified trait or attribute across the selected branches.

        Note:
        - Branches that do not have the specified trait or attribute are skipped.

        Example:
        >>> branch_lengths = tree.getParameter('length')
        >>> posteriors = tree.getParameter('posterior', use_trait=True)
        >>> node_heights = tree.getParameter('height', which=lambda x: x.is_node())

        Docstring generated with ChatGPT 4o.
        """
        branches = [k for k in self.Objects if which(k)]

        if not use_trait:
            params = [getattr(k, statistic) for k in branches if hasattr(k, statistic)]
        elif use_trait:
            params = [k.traits[statistic] for k in branches if statistic in k.traits]
        return params

    def fixHangingNodes(self):
        """
        Remove internal nodes without any children. Used in `reduceTree()` and `subtree()` functions internally.

        This method iterates over all objects in the tree and removes nodes that have no children. It continues to check
        for and remove hanging nodes until none are left.

        Example:
        >>> tree.fixHangingNodes()

        Docstring generated with ChatGPT 4o.
        """
        reticulations = {r.target: r for r in self.Objects if isinstance(r, reticulation)}

        while True:
            hanging_nodes = [
                node for node in self.Objects if is_node(node) and not node.children
            ]  ## nodes without children (hanging nodes)
            if not hanging_nodes:
                break

            dangling_rets = [reticulations[n] for n in hanging_nodes if n in reticulations]
            if dangling_rets:
                print("Dangling reticulations:", dangling_rets)
            for node in hanging_nodes:
                node.parent.children.remove(node)
                self.Objects.remove(node)

    def addText(
        self,
        ax: Axes,
        target: Callable[[BranchType], bool] = is_leaf,
        x_attr: Callable[[BranchType], float] = attrgetter("x"),
        y_attr: Callable[[BranchType], float] = attrgetter("y"),
        text: Callable[[BranchType], str] = attrgetter("name"),
        zorder: int = 4,
        font_kwargs: Callable[[BranchType], dict[str, Any]] | None = None,
        **kwargs,
    ):
        """
        Add text annotations to the tree plot.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the text to.
        target (function): A function to select which branches to annotate. Default selects all `leaf` nodes.
        x_attr (function): A function to determine the x-coordinate for the text. Default uses the branch's x attribute.
        y_attr (function): A function to determine the y-coordinate for the text. Default uses the branch's y attribute.
        text (function): A function to determine the text content. Default uses the `leaf` name attribute.
        zorder (int): The z-order for the text. Default sets the z-order to 4.
        **kwargs: Additional keyword arguments to pass to the `ax.text` method.

        Returns:
        matplotlib.axes.Axes: The axes with the text annotations added.

        Example:
        >>> tree.addText(ax, target=lambda node: node.is_node(), x_attr=lambda node: node.x - 7/365, y_attr=lambda node: node.y - 0.25, text=lambda node: node.traits['posterior'], ha='right', va='top') ## adds posterior values to the left and below internal nodes

        Docstring generated with ChatGPT 4o.
        """
        local_kwargs = dict(kwargs)
        if "verticalalignment" not in local_kwargs:
            local_kwargs["verticalalignment"] = "center"

        if font_kwargs is None:

            def func(k):
                return {}

            font_kwargs = func

        for k in filter(target, self.Objects):
            x, y = x_attr(k), y_attr(k)
            ax.text(x, y, text(k), zorder=zorder, **font_kwargs(k), **local_kwargs)
        return ax

    def addTextUnrooted(
        self,
        ax: Axes,
        target: Callable[[BranchType], bool] = is_leaf,
        x_attr: Callable[[BranchType], float] = attrgetter("x"),
        y_attr: Callable[[BranchType], float] = attrgetter("y"),
        text: Callable[[BranchType], str] = attrgetter("name"),
        zorder: int = 4,
        **kwargs,
    ):
        """
        Add text annotations to an unrooted tree plot.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the text to.
        target (function): A function to select which branches to annotate. Default selects all `leaf` nodes.
        x_attr (function): A function to determine the x-coordinate for the text. Default uses the branch's x attribute.
        y_attr (function): A function to determine the y-coordinate for the text. Default uses the branch's y attribute.
        text (function): A function to determine the text content. Default uses the branch's name attribute.
        zorder (int): The z-order for the text. Default sets the z-order to 4.
        **kwargs: Additional keyword arguments to pass to the `ax.text` method.

        Returns:
        matplotlib.axes.Axes: The axes with the text annotations added.

        Example:
        >>> tree.addTextUnrooted(ax) ## adds tip names to the tree

        Docstring generated with ChatGPT 4o.
        """
        for k in filter(target, self.Objects):
            local_kwargs = dict(kwargs)

            x, y = x_attr(k), y_attr(k)
            z = zorder

            assert "tau" in k.traits, "Branch does not have angle tau computed by drawUnrooted()."

            rot = math.degrees(k.traits["tau"]) % 360

            if "horizontalalignment" not in local_kwargs:
                local_kwargs["horizontalalignment"] = "right" if 90 < rot < 270 else "left"
            if "verticalalignment" not in local_kwargs:
                local_kwargs["verticalalignment"] = "center"

            rot = rot + 180 if 90 < rot < 270 else rot

            ax.text(
                x,
                y,
                text(k),
                rotation=rot,
                rotation_mode="anchor",
                zorder=z,
                **local_kwargs,
            )

        return ax

    def addTextCircular(
        self,
        ax: Axes,
        target: Callable[[BranchType], bool] = is_leaf,
        x_attr: Callable[[BranchType], float] = attrgetter("x"),
        y_attr: Callable[[BranchType], float] = attrgetter("y"),
        text: Callable[[BranchType], str] = attrgetter("name"),
        circStart: float = 0.0,
        circFrac: float = 1.0,
        inwardSpace: float = 0.0,
        normaliseHeight: Callable[[float], float] | None = None,
        zorder: int = 4,
        **kwargs,
    ):
        """
        Add text annotations to a circular tree plot.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the text to.
        target (function): A function to select which branches to annotate. Default selects all `leaf` nodes.
        text (function): A function to determine the text content. Default uses the `leaf` name attribute.
        x_attr (function): A function to determine the x-coordinate for the text. Default uses the branch's x attribute.
        y_attr (function): A function to determine the y-coordinate for the text. Default uses the branch's y attribute.
        circStart (float): The starting angle (in fractions of 2*pi, i.e. radians) for the circular layout. Default is 0.0.
        circFrac (float): The fraction of the full circle to use for the layout. Default is 1.0.
        inwardSpace (float): Amount of space to leave in the middle of the tree (can be negative for inward-facing trees). Default is 0.0.
        normaliseHeight (function or None): A function to normalize the x-coordinates. Default is None, creates a normalisation that returns 0.0 at root and 1.0 at the most diverged tip.
        zorder (int): The z-order for the text. Default sets the z-order to 4.
        **kwargs: Additional keyword arguments to pass to the `ax.text` method.

        Returns:
        matplotlib.axes.Axes: The axes with the text annotations added.

        Example:
        >>> tree.addTextCircular(ax) ## adds `leaf` names to a circular tree plot

        Docstring generated with ChatGPT 4o.
        """
        circ_s = circStart * math.pi * 2
        circ = circFrac * math.pi * 2

        allXs = [x_attr(k) for k in self.Objects]
        if normaliseHeight is None:

            def func(value):
                return (value - min(allXs)) / (max(allXs) - min(allXs))

            normaliseHeight = func

        for k in filter(target, self.Objects):  ## iterate over branches
            local_kwargs = dict(kwargs)  ## copy global kwargs into a local version

            x = normaliseHeight(x_attr(k) + inwardSpace)  ## get branch x position
            y = y_attr(k)  ## get y position

            y = circ_s + circ * y / self.ySpan
            X = math.sin(y)
            Y = math.cos(y)

            rot = math.degrees(y) % 360

            if "horizontalalignment" not in local_kwargs:
                local_kwargs["horizontalalignment"] = (
                    "right" if 180 < rot < 360 else "left"
                )  ## rotate labels to aid readability
            if "verticalalignment" not in local_kwargs:
                local_kwargs["verticalalignment"] = "center"
            rot = 360 - rot - 90 if 180 < rot < 360 else 360 - rot + 90

            ax.text(
                X * x,
                Y * x,
                text(k),
                rotation=rot,
                rotation_mode="anchor",
                zorder=zorder,
                **local_kwargs,
            )

        return ax

    def plotPoints(
        self,
        ax: Axes,
        x_attr: Callable[[BranchType], float] = attrgetter("x"),
        y_attr: Callable[[BranchType], float] = attrgetter("y"),
        target: Callable[[BranchType], bool] = is_leaf,
        size: int | Callable[[BranchType], float] = 40,
        colour: str | Callable[[BranchType], Any] = "k",
        zorder: int = 3,
        outline: bool = True,
        outline_size: int | Callable[[BranchType], float] | None = None,
        outline_colour: str | Callable[[BranchType], Any] = "k",
        **kwargs,
    ):
        """
        Plot points on the tree plot.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the points to.
        x_attr (function): A function to determine the x-coordinate for the points. Default uses the branch's x attribute.
        y_attr (function): A function to determine the y-coordinate for the points. Default uses the branch's y attribute.
        target (function): A function to select which branches to annotate. Default selects all `leaf` objects.
        size (int or function): The size of the points. Default sets the size to 40.
        colour (str or function): The color of the points. Default sets the color to 'k' (black).
        zorder (int): The z-order for the points. Default sets the z-order to 3.
        outline (bool): If True, adds an outline to the points. Default sets the outline to True.
        outline_size (int or function or None): The size of the outline. Default is None, which sets the outline size to twice the size of the points.
        outline_colour (str or function): The color of the outline. Default sets the outline color to 'k' (black).
        **kwargs: Additional keyword arguments to pass to the `ax.scatter` method.

        Returns:
        matplotlib.axes.Axes: The axes with the points added.

        Example:
        >>> tree.plotPoints(ax, target=lambda node: node.traits['posterior'] >= 0.95) ## will add circles at nodes with greater than 0.95 posterior support

        Docstring generated with ChatGPT 4o.
        """
        if outline_size is None:

            def func(k):
                return size(k) * 2 if callable(size) else size * 2

            outline_size = func

        xs = []
        ys = []
        colours = []
        sizes = []

        outline_xs = []
        outline_ys = []
        outline_colours = []
        outline_sizes = []
        for k in filter(target, self.Objects):
            xs.append(x_attr(k))
            ys.append(y_attr(k))
            colours.append(colour(k)) if callable(colour) else colours.append(colour)
            sizes.append(size(k)) if callable(size) else sizes.append(size)

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outline_colour(k)) if callable(outline_colour) else outline_colours.append(
                    outline_colour
                )
                outline_sizes.append(outline_size(k)) if callable(outline_size) else outline_sizes.append(outline_size)

        ax.scatter(
            xs,
            ys,
            s=sizes,
            facecolor=colours,
            edgecolor="none",
            zorder=zorder,
            **kwargs,
        )  ## put a circle at each tip
        if outline:
            ax.scatter(
                outline_xs,
                outline_ys,
                s=outline_sizes,
                facecolor=outline_colours,
                edgecolor="none",
                zorder=zorder - 1,
                **kwargs,
            )  ## put a circle at each tip

        return ax

    def plotTree(
        self,
        ax: Axes,
        connection_type: Literal["baltic", "direct", "elbow"] = "baltic",
        target: Callable[[BranchType], bool] = always_true,
        x_attr: Callable[[BranchType], float] | str = "x",
        y_attr: Callable[[BranchType], float] | str = "y",
        width: int | Callable[[BranchType], float] = 2,
        colour: str | Callable[[BranchType], Any] = "k",
        autoscale: bool = True,
        **kwargs,
    ):
        """
        Plot the tree on a given matplotlib axes.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to plot the tree on.
        connection_type (str or None): The type of connection between nodes. Options are:
            - 'baltic' (parental branches are plotted as two straight lines - one horizontal, one vertical)
            - 'direct' (diagonal line that directly connects parent and child branches)
            - 'elbow' (each child has its own angled branch connecting it to the parent). Default is 'baltic'.
        target (function): A function to select which branches to plot. Default selects all branches.
        x_attr (function or str): A function or attribute name to determine the x-coordinate for the nodes. Default uses the branch's x attribute.
        y_attr (function or str): A function or attribute name to determine the y-coordinate for the nodes. Default uses the branch's y attribute.
        width (int or function): The width of the lines. Default sets the width to 2.
        colour (str or function): The color of the lines. Default sets the color to 'k' (black).
        **kwargs: Additional keyword arguments to pass to the LineCollection.

        Returns:
        matplotlib.axes.Axes: The axes with the tree plot added.

        Example:
        >>> tree.plotTree(ax)

        Docstring generated with ChatGPT 4o.
        """
        assert connection_type in ["baltic", "direct", "elbow"], f"Unrecognised drawing type {connection_type!r}"

        if isinstance(x_attr, str):
            x_attr = attrgetter(x_attr)

        if isinstance(y_attr, str):
            y_attr = attrgetter(y_attr)

        branches = []
        colours = []
        linewidths = []
        for k in filter(target, self.Objects):  ## iterate over branches
            x = x_attr(k)  ## get branch x position
            xp = x_attr(k.parent) if k.parent else x  ## get parent x position
            y = y_attr(k)  ## get y position

            colours.append(colour(k) if callable(colour) else colour)
            linewidths.append(width(k) if callable(width) else width)

            if connection_type == "baltic":
                ## each node has a single vertical line to which descendant branches are connected
                branches.append(((xp, y), (x, y)))
                if is_node(k):
                    children = [y_attr(ch) for ch in k.children]
                    yl, yr = min(children), max(children)  ## y positions of first and last child
                    branches.append(((x, yl), (x, yr)))
                    linewidths.append(linewidths[-1])
                    colours.append(colours[-1])
            elif connection_type == "elbow":
                ## more standard connection where each branch connects to its parent via a right-angled line
                yp = y_attr(k.parent) if k.parent else y  ## get parent x position
                branches.append(((xp, yp), (xp, y), (x, y)))
            elif connection_type == "direct":
                ## this gives triangular looking trees where descendants connect directly to their parents
                yp = y_attr(k.parent)  ## get y position
                branches.append(((xp, yp), (x, y)))
            else:
                pass  ## for now

        if "capstyle" not in kwargs:
            kwargs["capstyle"] = "projecting"
        line_segments = LineCollection(branches, lw=linewidths, color=colours, **kwargs)

        if autoscale:
            ax.add_collection(line_segments, autolim=True)
            ax.autoscale_view()
        else:
            ax.add_collection(line_segments)
        return ax

    def plotCircularTree(
        self,
        ax: Axes,
        target: Callable[[BranchType], bool] = always_true,
        x_attr: Callable[[BranchType], float] = attrgetter("x"),
        y_attr: Callable[[BranchType], float] = attrgetter("y"),
        width: float | Callable[[BranchType], float] = 2,
        colour: Callable[[BranchType], Any] | str = "k",
        circStart: float = 0.0,
        circFrac: float = 1.0,
        inwardSpace: float = 0.0,
        normaliseHeight: Callable[[float], float] | None = None,
        precision: int = 15,
        autoscale: bool = True,
        **kwargs,
    ):
        """
        Plot the tree in a circular layout on a given matplotlib axes.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to plot the tree on.
        target (function): A function to select which branches to plot. Default selects all branches.
        x_attr (function): A function to determine the x-coordinate for the nodes. Default uses the branch's x attribute.
        y_attr (function): A function to determine the y-coordinate for the nodes. Default uses the branch's y attribute.
        width (int or function): The width of the lines. Default is None, which sets the width to 2.
        colour (str or function): The color of the lines. Default is None, which sets the color to 'k' (black).
        circStart (float): The starting angle (in fractions of 2*pi, i.e. radians) for the circular layout. Default is 0.0.
        circFrac (float): The fraction of the full circle to use for the layout. Default is 1.0.
        inwardSpace (float): Amount of space to leave in the middle of the tree (can be negative for inward-facing trees). Default is 0.0.
        normaliseHeight (function or None): A function to normalize the x-coordinates. Default is None, creates a normalisation that returns 0.0 at root and 1.0 at the most diverged tip.
        precision (int): The number of points used to plot curved segments. Default is 15.
        **kwargs: Additional keyword arguments to pass to the LineCollection.

        Returns:
        matplotlib.axes.Axes: The axes with the circular tree plot added.

        Example:
        >>> tree.plotCircularTree(ax)

        Docstring generated with ChatGPT 4o.
        """
        if inwardSpace < 0:
            inwardSpace -= self.treeHeight

        branches = []
        colours = []
        linewidths = []

        circ_s = circStart * math.pi * 2
        circ = circFrac * math.pi * 2

        allXs = [x_attr(k) for k in self.Objects]
        if normaliseHeight is None:

            def func(value):
                return (value - min(allXs)) / (max(allXs) - min(allXs))

            normaliseHeight = func

        def linspace(start: float, stop: float, n: int):
            return [start + ((stop - start) / (n - 1)) * i for i in range(n)] if n > 1 else stop

        for k in filter(target, self.Objects):  ## iterate over branches
            x = normaliseHeight(x_attr(k) + inwardSpace)  ## get branch x position
            xp = normaliseHeight(x_attr(k.parent) + inwardSpace) if k.parent.parent else x  ## get parent x position
            y = y_attr(k)  ## get y position

            colours.append(colour(k) if callable(colour) else colour)
            linewidths.append(width(k) if callable(width) else width)

            y = circ_s + circ * y / self.ySpan
            X = math.sin(y)
            Y = math.cos(y)
            branches.append(((X * xp, Y * xp), (X * x, Y * x)))

            if is_node(k):
                yl, yr = (
                    y_attr(k.children[0]),
                    y_attr(k.children[-1]),
                )  ## get leftmost and rightmost children's y coordinates
                yl = circ_s + circ * yl / self.ySpan  ## transform y into a fraction of total y
                yr = circ_s + circ * yr / self.ySpan
                ybar = linspace(yl, yr, precision)  ## what used to be vertical node bar is now a curved line
                assert isinstance(ybar, list)

                xs = [yx * x for yx in map(math.sin, ybar)]  ## convert to polar coordinates
                ys = [yy * x for yy in map(math.cos, ybar)]

                branches += tuple(zip(zip(xs, ys), zip(xs[1:], ys[1:])))  ## add curved segment

                linewidths += [linewidths[-1] for q in zip(ys, ys[1:])]  ## repeat linewidths
                colours += [colours[-1] for q in zip(ys, ys[1:])]  ## repeat colours

        line_segments = LineCollection(
            branches, lw=linewidths, ls="-", color=colours, capstyle="projecting", zorder=1, **kwargs
        )  ## create line segments
        if autoscale:
            ax.add_collection(line_segments, autolim=True)
            ax.autoscale_view()
        else:
            ax.add_collection(line_segments)
        return ax

    def plotCircularPoints(
        self,
        ax,
        x_attr: Callable[[BranchType], float] = attrgetter("x"),
        y_attr: Callable[[BranchType], float] = attrgetter("y"),
        target: Callable[[BranchType], bool] = is_leaf,
        size: float | Callable[[BranchType], float] = 40,
        colour: Callable[[BranchType], Any] | str = "k",
        circStart: float = 0.0,
        circFrac: float = 1.0,
        inwardSpace: float = 0.0,
        normaliseHeight: Callable[[float], float] | None = None,
        zorder: int = 3,
        outline: bool = True,
        outline_size: float | Callable[[BranchType], float] | None = None,
        outline_colour: Callable[[BranchType], Any] | str = "k",
        **kwargs,
    ):
        """
        Plot points on a circular tree plot.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the points to.
        x_attr (function or None): A function to determine the x-coordinate for the points. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the points. Default is None, which uses the branch's y attribute.
        target (function or None): A function to select which branches to annotate. Default is None, which selects all `leaf` nodes.
        size (int or function or None): The size of the points. Default is None, which sets the size to 40.
        colour (str or function or None): The color of the points. Default is None, which sets the color to 'k' (black).
        circStart (float): The starting angle (in fractions of 2*pi, i.e. radians) for the circular layout. Default is 0.0.
        circFrac (float): The fraction of the full circle to use for the layout. Default is 1.0.
        inwardSpace (float): Amount of space to leave in the middle of the tree (can be negative for inward-facing trees). Default is 0.0.
        normaliseHeight (function or None): A function to normalize the x-coordinates. Default is None, creates a normalisation that returns 0.0 at root and 1.0 at the most diverged tip.
        zorder (int or None): The z-order for the points. Default is None, which sets the z-order to 3.
        outline (bool or None): If True, adds an outline to the points. Default is None, which sets the outline to True.
        outline_size (int or function or None): The size of the outline. Default is None, which sets the outline size to twice the size of the points.
        outline_colour (str or function or None): The color of the outline. Default is None, which sets the outline color to 'k' (black).
        **kwargs: Additional keyword arguments to pass to the `ax.scatter` method.

        Returns:
        matplotlib.axes.Axes: The axes with the points added.

        Example:
        >>> tree.plotCircularPoints(ax)

        Docstring generated with ChatGPT 4o.
        """
        if outline_size is None:

            def func1(k):
                return size(k) * 2 if callable(size) else size * 2

            outline_size = func1

        if inwardSpace < 0:
            inwardSpace -= self.treeHeight

        circ_s = circStart * math.pi * 2
        circ = circFrac * math.pi * 2

        allXs = list(map(x_attr, self.Objects))
        if normaliseHeight is None:

            def func(value):
                return (value - min(allXs)) / (max(allXs) - min(allXs))

            normaliseHeight = func

        def linspace(start, stop, n):
            return list(start + ((stop - start) / (n - 1)) * i for i in range(n)) if n > 1 else stop

        xs = []
        ys = []
        colours = []
        sizes = []

        outline_xs = []
        outline_ys = []
        outline_colours = []
        outline_sizes = []
        for k in filter(target, self.Objects):
            x = normaliseHeight(x_attr(k) + inwardSpace)  ## find normalised x position along circle's radius
            y = circ_s + circ * y_attr(k) / self.ySpan  ## get y position along circle's perimeter
            X = math.sin(y) * x  ## transform
            Y = math.cos(y) * x  ## transform

            xs.append(X)
            ys.append(Y)
            colours.append(colour(k)) if callable(colour) else colours.append(colour)
            sizes.append(size(k)) if callable(size) else sizes.append(size)

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outline_colour(k)) if callable(outline_colour) else outline_colours.append(
                    outline_colour
                )
                outline_sizes.append(outline_size(k)) if callable(outline_size) else outline_sizes.append(outline_size)

        ax.scatter(
            xs,
            ys,
            s=sizes,
            facecolor=colours,
            edgecolor="none",
            zorder=zorder,
            **kwargs,
        )  ## put a circle at each tip
        if outline:
            ax.scatter(
                outline_xs,
                outline_ys,
                s=outline_sizes,
                facecolor=outline_colours,
                edgecolor="none",
                zorder=zorder - 1,
                **kwargs,
            )  ## put a circle at each tip

        return ax


def untangle(
    trees: list, cost_function: Callable[..., float] | None = None, iterations: int | None = None, verbose: bool = False
):
    """
    Minimise y-axis discrepancies between tips of trees in a list.
    Only the tangling of adjacent trees in the list is minimised, so the order of trees matters.
    Trees do not need to have the same number of tips but tip names should match.

    Parameters:
    trees (list): A list of tree objects to untangle.
    cost_function (function or None): A function to calculate the cost of y-axis discrepancies between tips.
                                      Default is None, which uses the squared difference between y axis positions.
    iterations (int or None): The number of iterations to perform. Default is None, which sets the iterations to 3.
    verbose (bool): If True, prints verbose output during the process. Default is False.

    Returns:
    list: The list of untangled tree objects.

    Example:
    >>> untangled_trees = untangle(list_of_trees, iterations=5, verbose=True)

    Docstring generated with ChatGPT 4o.
    """
    from itertools import permutations

    if iterations is None:
        iterations = 3
    if cost_function is None:

        def func(pair):
            return math.pow(abs(pair[0] - pair[1]), 2)

        cost_function = func

    y_positions = {
        T: {k.name: k.y for k in T.getExternal()} for T in trees
    }  ## get y positions of all the tips in every tree

    for iteration in range(iterations):
        if verbose:
            print("Untangling iteration %d" % (iteration + 1))
        first_trees = list(range(len(trees) - 1)) + [-1]  ## trees up to next-to-last + last
        next_trees = list(range(1, len(trees))) + [0]  ## trees from second + first
        for cur, nex in zip(first_trees, next_trees):  ## adjacent pairs
            tree1 = trees[cur]  ## fetch current tree
            tree2 = trees[nex]  ## fetch next tree
            if verbose:
                print("%d vs %d" % (cur, nex))
            for k in sorted(
                tree2.getInternal(), key=lambda branch: branch.height
            ):  ## iterate through nodes of next tree by height (start from root)
                clade_y_positions = sorted(
                    [y_positions[tree2][tip] for tip in k.leaves]
                )  ## sorted list of available y coordinates for node
                costs = {}  ## will store cost of all children permutations
                if len(k.children) >= 10:
                    raise RuntimeWarning("Node is too polytomic and untangling will take an astronomically long time")
                if verbose:
                    print(len(k.children))
                for permutation in permutations(k.children):  ## iterate over permutations of node's children
                    clade_order = sum(
                        [[child.name] if child.is_leaf() else list(child.leaves) for child in permutation],
                        [],
                    )  ## flat list of tip names as they would appear in permutation order
                    new_y_positions = {
                        clade_order[i]: clade_y_positions[i] for i in range(len(clade_y_positions))
                    }  ## assign available y positions in order

                    tip_costs = list(
                        map(
                            cost_function,
                            [
                                (y_positions[tree1][tip], new_y_positions[tip])
                                for tip in clade_order
                                if tip in y_positions[tree1]
                            ],
                        )
                    )
                    costs[permutation] = sum(tip_costs) / len(
                        tip_costs
                    )  ## compute cost of this permutation in relation to next tree

                best = sorted(costs.keys(), key=lambda w: -costs[w])[0]  ## get tree with smallest cost
                k.children = list(best)  ## reorder children according to minimised cost

            tree2.drawTree()  ## compute new y coordinates for nodes
            for k in tree2.getExternal():  ## iterate over tips
                y_positions[tree2][k.name] = k.y  ## remember new coordinates

    return trees


def make_tree(data: str, ll: tree | None = None, verbose: bool = False):
    """
    Parse a tree string and create a tree object.

    Parameters:
    data (str): The tree string to be parsed.
    ll (tree or None): An instance of a tree object. If None, a new tree object is created. Default is None.
    verbose (bool): If True, prints verbose output during the process. Default is False.

    Returns:
    tree: The tree object created from the parsed tree string.

    Example:
    >>> tree_string = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    >>> tree = make_tree(tree_string)

    Docstring generated with ChatGPT 4o.
    """
    # Pattern to match tips in BEAST format (integers)
    BEAST_TIP_RE = re.compile(r"(\(|,)([0-9]+)(\[|\:)")
    # Pattern to match tips with unencoded names
    NON_BEAST_TIP_RE = re.compile(r"(\(|,)(\'|\")*([^\(\):\[\'\"#]+)(\'|\"|)*(\[)*")
    RETICULATION_BEGIN_RE = re.compile(r"[\(,](#[A-Za-z0-9]+)")
    RETICULATION_END_RE = re.compile(r"\)(#[A-Za-z0-9]+)")
    MCC_COMMENTS_RE = re.compile(r"(?:\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&]+)\]")

    # Add in some checks that the data are in correct format
    assert data.endswith(";"), "Improperly formatted string: must end in semicolon"
    assert data.count("(") == data.count(")"), "Improperly formatted string: must have matching parentheses"

    if ll is None:  ## calling without providing a tree object - create one
        ll = tree()
    i = 0  ## is an adjustable index along the tree string, it is incremented to advance through the string

    ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop
    stored_i = None

    while i < len(data):  ## while there's characters left in the tree string - loop away
        if stored_i == i and verbose:
            print(f"{i} >{data[i]}<")

        assert stored_i != i, (
            f"\nTree string unparseable\nStopped at >>{data[i]}<<\nstring region looks like this: {data[i : i + 5000]}"
        )  # Ensure that the index has advanced; if not, raise an error indicating an unparseable string
        stored_i = i  # Store the current index at the end of the loop to check for infinite loops

        if data[i] == "(":  ## look for new nodes
            if verbose:
                print(f"{i} adding node")
            ll.add_node(i)  ## add node to current node in tree ll
            i += 1  ## advance in tree string by one character

        match = BEAST_TIP_RE.match(data[i - 1 : i + 100])  ## look for tips in BEAST format (integers).
        if match:
            if verbose:
                print(f"{i} adding leaf (BEAST) {match[2]!r}")
            ll.add_leaf(i, match[2])  ## add tip
            i += len(match[2])  ## advance in tree string by however many characters the tip is encoded

        assert ll.cur_node is not None

        match = NON_BEAST_TIP_RE.match(
            data[i - 1 : i + 200]
        )  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if match:
            if verbose:
                print(f"{i} adding leaf (non-BEAST) {match[3]!r}")
            ll.add_leaf(i, match[3].strip('"').strip("'"))  ## add tip
            i += (
                len(match[3]) + match[0].count("'") + match[0].count('"')
            )  ## advance in tree string by however many characters the tip is encoded

        match = re.match(r"\)([0-9]+)\[", data[i - 1 : i + 100])  ## look for multitype tree singletons.
        if match:
            if verbose:
                print(f"{i} adding multitype node {match[1]!r}")
            i += len(match[1])

        match = RETICULATION_BEGIN_RE.match(data[i - 1 : i + 200])  ## look for beginning of reticulate branch
        if match:
            if verbose:
                print(f"{i} adding outgoing reticulation branch {match[1]}")
            ll.add_reticulation(match[1])  ## add reticulate branch

            destination = None
            for k in ll.Objects:  ## iterate over branches parsed so far
                if k.traits.get("label", None) == match[1]:  ## if there's a branch with a matching id
                    if destination is None:  ## not set destination before
                        destination = k  ## destination is matching node
                    else:  ## destination seen before - raise an error (indicates reticulate branch ids are not unique)
                        raise ValueError(f"Reticulate branch not unique: {match[1]} seen elsewhere in the tree")
            if destination:  ## identified destination of this branch
                if verbose:
                    print(f"identified {match[1]} destination")
                assert isinstance(ll.cur_node, reticulation)
                ll.cur_node.target = destination  ## set current node's target as the destination
                setattr(destination, "contribution", ll.cur_node)  ## add contributing edge to destination
            else:
                if verbose:
                    print(f"destination of {match[1]} not identified yet")
            i += len(match[0]) - 1

        match = RETICULATION_END_RE.match(data[i - 1 : i + 200])  ## look for landing point of reticulate branch
        if match:
            if verbose:
                print(f"{i} adding incoming reticulation branch {match[1]}")
            ll.cur_node.traits["label"] = match[1]  ## set node label

            origin = None  ## branch is landing, check if its origin was seen previously
            for k in ll.Objects:  ## iterate over currently existing branches
                if (
                    isinstance(k, reticulation) and k.name == match[1]
                ):  ## check if any reticulate branches match the origin
                    if origin is None:  ## origin not identified yet
                        origin = k  ## origin is reticulate branch with the correct name
                    else:  ## origin has been identified - shouldn't happen, implies that multiple reticulate branches exist with the same name
                        raise ValueError(f"Reticulate branch not unique: {match[1]!r} seen elsewhere in the tree")
            if origin:  ## identified origin
                if verbose:
                    print(f"identified {match[1]!r} origin")
                origin.target = ll.cur_node  ## set origin's landing at this node
                setattr(ll.cur_node, "contribution", origin)  ## add contributing edge to this node
            else:
                if verbose:
                    print(f"origin of {match[1]!r} not identified yet")
            i += len(match[0]) - 1

        match = MCC_COMMENTS_RE.match(data[i:])
        if match:
            if verbose:
                print(f"{i} comment: {match[1]!r}")
            comment = match[1]
            numerics = re.findall(
                r"[,&][A-Za-z\_\.0-9]+=[0-9\-Ee\.]+", comment
            )  ## find all entries that have values as floats
            strings = re.findall(
                r"[,&][A-Za-z\_\.0-9]+=[\"|']*[A-Za-z\_0-9\.\+ :\/\(\)\&\-]+[\"|']*",
                comment,
            )  ## strings
            treelist = re.findall(
                r"[,&][A-Za-z\_\.0-9]+={[A-Za-z\_,{}0-9\. :\/\(\)\&]+}", comment
            )  ## complete history logged robust counting (MCMC trees)
            sets = re.findall(r'[,&][A-Za-z\_\.0-9\%]+={[A-Za-z\.\-0-9eE,"\_ :\/\(\)\&]+}', comment)  ## sets and ranges
            figtree = re.findall(r"\![A-Za-z]+=[A-Za-z0-9# :\/\(\)\&]+", comment)

            for vals in strings:
                tr, val = vals.split("=")
                tr = tr[1:]
                if "+" in val:
                    ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                    val = val.split("+")[0]
                ll.cur_node.traits[tr] = val.strip('"')

            for vals in numerics:  ## assign all parsed annotations to traits of current branch
                tr, val = vals.split("=")  ## split each value by =, left side is name, right side is value
                tr = tr[1:]
                if val.replace("E", "", 1).replace("e", "", 1).replace("-", "", 1).replace(".", "", 1).isdigit():
                    ll.cur_node.traits[tr] = float(val)

            for val in treelist:
                tr, val = val.split("=")
                tr = tr[1:]
                micromatch = []
                if val.count(",") == 2:
                    micromatch = re.findall(r"{([0-9\.\-e]+,[a-z_A-Z]+,[a-z_A-Z]+)}", val)
                elif val.count(",") == 3:
                    micromatch = re.findall(r"{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}", val)
                ll.cur_node.traits[tr] = []
                for val in micromatch:
                    val.split(",")
                    ll.cur_node.traits[tr].append(val.split(","))

            for vals in sets:
                tr, val = vals.split("=")
                tr = tr[1:]
                if "set" in tr:
                    ll.cur_node.traits[tr] = []
                    for v in val[1:-1].split(","):
                        if "set.prob" in tr:
                            ll.cur_node.traits[tr].append(float(v))
                        else:
                            ll.cur_node.traits[tr].append(v.strip('"'))
                else:
                    try:
                        ll.cur_node.traits[tr] = [float(v) for v in val[1:-1].split(",")]
                    except (ValueError, TypeError):
                        print(f"some other trait: {vals!r}")

            if len(figtree) > 0:
                print("FigTree comment found, ignoring")

            i += len(match[0])  ## advance in tree string by however many characters it took to encode labels

        match = re.match(r"([A-Za-z\_\-0-9\.]+)([:;\[])", data[i:])  ## look for old school node labels
        if match:
            if verbose:
                print(f"old school comment found: {match[1]}")
            ll.cur_node.traits["label"] = match[1]

            i += len(match[1])

        micromatch = re.match(r"(\:)*([0-9\.\-Ee]+)", data[i : i + 100])  ## look for branch lengths without comments
        if micromatch is not None:
            if verbose:
                print("adding branch length (%d) %.6f" % (i, float(micromatch.group(2))))
            ll.cur_node.length = float(micromatch.group(2))  ## set branch length of current node
            ## advance in tree string by however many characters it took to encode branch length
            i += len(micromatch[0])

        if data[i] == "," or data[i] == ")":  ## look for bifurcations or clade ends
            i += 1  ## advance in tree string
            ll.cur_node = ll.cur_node.parent

        if data[i] == ";":  ## look for string end
            return ll
            break  ## end loop


def make_treeJSON(JSONnode: dict, json_translation: dict, ll: tree | None = None, verbose: bool = False):
    """
    Parse an auspice JSON tree and create a baltic tree object.

    Parameters:
    JSONnode (dict): The JSON node to be parsed.
    json_translation (dict): A dictionary for translating JSON keys to tree attributes.
    ll (tree or None): An instance of a tree object. If None, a new tree object is created. Default is None.
    verbose (bool): If True, prints verbose output during the process. Default is False.

    Returns:
    tree: The tree object created from the parsed JSON.

    Docstring generated with ChatGPT 4o.
    """
    if "children" in JSONnode:  ## only nodes have children
        new_node = node()
    else:
        new_node = leaf()
        new_node.name = JSONnode[json_translation["name"]]  ## set leaf name to be the same

    if ll is None:
        ll = tree()
        ll.root = new_node
    if "attr" in JSONnode:
        attr = JSONnode.pop("attr")
        JSONnode.update(attr)

    new_node.parent = ll.cur_node  ## set parent-child relationships

    assert ll.cur_node is not None and is_node(ll.cur_node)
    ll.cur_node.children.append(new_node)
    new_node.index = JSONnode[json_translation["name"]]  ## indexing is based on name
    new_node.traits = {
        n: JSONnode[n] for n in list(JSONnode.keys()) if n != "children"
    }  ## set traits to non-children attributes
    ll.Objects.append(new_node)
    ll.cur_node = new_node

    if "children" in JSONnode:
        for child in JSONnode["children"]:
            make_treeJSON(child, json_translation, ll)
            assert ll.cur_node is not None
            ll.cur_node = ll.cur_node.parent
    return ll


def loadNewick(
    tree_path,
    tip_regex: str = r"\|([0-9]+\-[0-9]+\-[0-9]+)",
    date_fmt: str = "%Y-%m-%d",
    variableDate: bool = True,
    absoluteTime: bool = False,
    verbose: bool = False,
    sortBranches: bool = True,
):
    r"""
    Load a tree from a Newick file and process it.

    Parameters:
    tree_path (str or file-like object): The path to the Newick file or a file-like object containing the Newick formatted tree.
    tip_regex (str): A regular expression to extract dates from tip names. Default is '\|([0-9]+\-[0-9]+\-[0-9]+)'.
    date_fmt (str): The date format for the extracted dates. Default is '%Y-%m-%d'.
    variableDate (bool): If True, allows for variable date formats. Default is True.
    absoluteTime (bool): If True, converts the tree to absolute time using the tip dates encoded in tip names. Default is False.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    sortBranches (bool): If True, sorts the branches of the tree after loading. Default is True.

    Returns:
    tree: The tree object created from the Newick file.

    Raises:
    AssertionError: If the tree string cannot be found or if tip dates cannot be extracted when absoluteTime is True.

    Example:
    >>> tree = loadNewick("path/to/tree.newick", absoluteTime=False, verbose=True)

    Docstring generated with ChatGPT 4o.
    """
    ll = None

    handle = open(tree_path, "r") if isinstance(tree_path, str) else tree_path

    for line in handle:
        line = line.strip("\n")
        if "(" in line:
            treeString_start = line.index("(")
            ll = make_tree(line[treeString_start:], verbose=verbose)  ## send tree string to make_tree function
            if verbose:
                print("Identified tree string")

    assert ll, "Regular expression failed to find tree string"
    ll.traverse_tree(verbose=verbose)  ## traverse tree

    if sortBranches:
        ll.sortBranches()  ## traverses tree, sorts branches, draws tree

    if absoluteTime:
        tip_dates = []
        tip_names = []
        for k in ll.getExternal():
            tip_names.append(k.name)
            match = re.search(tip_regex, k.name)
            if match:
                tip_dates.append(decimalDate(match.group(1), fmt=date_fmt, variable=variableDate))
        assert len(tip_dates) > 0, (
            "Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s"
            % (tip_names[0], tip_regex, date_fmt)
        )
        ll.setAbsoluteTime(max(tip_dates))

    if isinstance(tree_path, str):
        handle.close()
    return ll


def loadNexus(
    tree_path,
    tip_regex=r"\|([0-9]+\-[0-9]+\-[0-9]+)",
    date_fmt="%Y-%m-%d",
    treestring_regex=r"tree [A-Za-z\_]+([0-9]+)",
    variableDate=True,
    absoluteTime=True,
    verbose=False,
    sortBranches=True,
):
    r"""
    Load a tree from a Nexus file and process it.

    Parameters:
    tree_path (str or file-like object): The path to the Nexus file or a file-like object containing the NEXUS formatted tree.
    tip_regex (str): A regular expression to extract dates from tip names. Default is '\|([0-9]+\-[0-9]+\-[0-9]+)'.
    date_fmt (str): The date format for the extracted dates. Default is '%Y-%m-%d'.
    treestring_regex (str): A regular expression to identify the tree string in the NEXUS file. Default is 'tree [A-Za-z\_]+([0-9]+)'.
    variableDate (bool): If True, allows for variable date formats. Default is True.
    absoluteTime (bool): If True, converts the tree to absolute time using the tip dates extracted from tip names. Default is True.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    sortBranches (bool): If True, sorts the branches of the tree after loading. Default is True.

    Returns:
    tree: The tree object created from the NEXUS file.

    Raises:
    AssertionError: If the tree string cannot be found or if tip dates cannot be extracted when absoluteTime is True.

    Example:
    >>> tree = loadNexus("path/to/tree.nexus", absoluteTime=True, verbose=True)

    Docstring generated with ChatGPT 4o.
    """
    tip_flag = False
    tips: dict[str, str] = {}
    tip_num = 0
    ll = None
    NTAXA_RE = re.compile(r"dimensions ntax=(\d+);", flags=re.I)
    TREE_RE = re.compile(treestring_regex, flags=re.I)
    TRANSLATE_RE = re.compile(r"([0-9]+) ([a-z0-9\-_/.'\" |?]+)", flags=re.I)

    if isinstance(tree_path, str):
        handle = Path(tree_path).read_text().splitlines()
    else:
        handle = tree_path.read().splitlines()

    for line in handle:
        match = NTAXA_RE.search(line)
        if match:
            tip_num = int(match.group(1))
            if verbose:
                print("File should contain %d taxa" % (tip_num))

        match = TREE_RE.search(line)
        if match:
            treeString_start = line.index("(")
            ll = make_tree(line[treeString_start:], verbose=verbose)  ## send tree string to make_tree function
            if verbose:
                print("Identified tree string")

        if tip_flag:
            match = TRANSLATE_RE.search(line)
            if match:
                tips[match.group(1)] = match.group(2).strip('"').strip("'")
                if verbose:
                    print("Identified tip translation %s: %s" % (match.group(1), tips[match.group(1)]))
            elif ";" not in line:
                print("tip not captured by regex:", line.replace("\t", ""))

        if "Translate" in line:
            tip_flag = True
        if ";" in line:
            tip_flag = False

    assert ll, "Failed to find tree string using regular expression"
    ll.traverse_tree()  ## traverse tree
    if sortBranches:
        ll.sortBranches()  ## traverses tree, sorts branches, draws tree
    if len(tips) > 0:
        ll.renameTips(tips)  ## renames tips from numbers to actual names
        ll.tipMap = tips
    if absoluteTime:
        tip_dates = []
        tip_names = []
        for k in ll.getExternal():
            tip_names.append(k.name)
            match = re.search(tip_regex, k.name)
            if match:
                tip_dates.append(decimalDate(match.group(1), fmt=date_fmt, variable=variableDate))

        assert len(tip_dates) > 0, (
            "Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s"
            % (tip_names[0], tip_regex, date_fmt)
        )
        ll.setAbsoluteTime(max(tip_dates))

    return ll


def loadJSON(
    json_object,
    json_translation={"name": "name", "absoluteTime": "num_date"},
    verbose=False,
    sort=True,
    stats=True,
):
    """
    Load a Nextstrain JSON file and create a tree object.

    Parameters:
    json_object (str or dict): The path to the JSON file, a URL to a Nextstrain JSON, or a JSON object.
    json_translation (dict): A dictionary that translates JSON attributes to tree attributes (e.g., baltic branch attribute 'absoluteTime' is called 'num_date' in Nextstrain JSONs).
                             Default is {'name': 'name', 'absoluteTime': 'num_date'}.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    sort (bool): If True, sorts the branches of the tree after loading. Default is True.
    stats (bool): If True, calculates tree statistics after loading. Default is True.

    Returns:
    tuple: A tuple containing the tree object created from the JSON and the metadata from the JSON.

    Raises:
    AssertionError: If the json_translation dictionary is missing the `name` attribute (crucial for tips) and least one branch length attribute (`absoluteTime`, `length` or `height`).
    KeyError: If a string attribute in json_translation is not found in the JSON data structure.
    AttributeError: If an attribute in json_translation is neither a string nor callable.

    Example:
    >>> tree, metadata = loadJSON("path/to/tree.json", verbose=True)

    Docstring generated with ChatGPT 4o.
    """
    required_keys = ["absoluteTime", "length", "height"]
    assert "name" in json_translation and any(key in json_translation for key in required_keys), (
        "JSON translation dictionary missing entries: %s"
        % (", ".join([entry for entry in ["name"] + required_keys if entry not in json_translation]))
    )
    if verbose:
        print("Reading JSON")

    if isinstance(json_object, str):  ## string provided - either nextstrain URL or local path
        if "nextstrain.org" in json_object:  ## nextsrain.org in URL - request it
            if verbose:
                print("Assume URL provided, loading JSON from nextstrain.org")
            from io import BytesIO as csio

            import requests

            auspice_json = json.load(csio(requests.get(json_object).content))
        else:  ## not nextstrain.org URL - assume local path to auspice v2 json
            if verbose:
                print("Loading JSON from local path")
            with open(json_object) as json_data:
                auspice_json = json.load(json_data)
    else:  ## not string, assume auspice v2 json object given
        if verbose:
            print("Loading JSON from object given")
        auspice_json = json_object

    json_meta = auspice_json["meta"]
    json_tree = auspice_json["tree"]
    ll = make_treeJSON(json_tree, json_translation, verbose=verbose)

    assert (
        "absoluteTime" in json_translation and ("length" not in json_translation or "height" not in json_translation)
    ) or ("absoluteTime" not in json_translation and ("length" in json_translation or "height" in json_translation)), (
        "Cannot use both absolute time and branch length, include only one in json_translation dictionary."
    )

    if verbose:
        print("Setting baltic traits from JSON")
    for k in ll.Objects:  ## make node attributes easier to access
        for key in k.traits["node_attrs"]:
            if isinstance(k.traits["node_attrs"][key], dict):
                if "value" in k.traits["node_attrs"][key]:
                    k.traits[key] = k.traits["node_attrs"][key]["value"]
                if "confidence" in k.traits["node_attrs"][key]:
                    k.traits["%s_confidence" % (key)] = k.traits["node_attrs"][key]["confidence"]
            elif key == "div":
                k.traits["divergence"] = k.traits["node_attrs"][key]

    for attr in json_translation:  ## iterate through attributes in json_translation
        for k in ll.Objects:  ## for every branch
            if isinstance(json_translation[attr], str):
                if json_translation[attr] in k.traits:
                    setattr(k, attr, k.traits[json_translation[attr]])  ## set attribute value for branch
                elif "node_attrs" in k.traits and json_translation[attr] in k.traits["node_attrs"]:
                    setattr(k, attr, k.traits["node_attrs"][json_translation[attr]])
                elif "branch_attrs" in k.traits and json_translation[attr] in k.traits["branch_attrs"]:
                    setattr(k, attr, k.traits["branch_attrs"][json_translation[attr]])
                else:
                    raise KeyError("String attribute %s not found in JSON" % (json_translation[attr]))
            elif callable(json_translation[attr]):
                setattr(k, attr, json_translation[attr](k))  ## set attribute value with a function for branch
            else:
                raise AttributeError("Attribute %s neither string nor callable" % (json_translation[attr]))

    for branch_unit in [
        "height",
        "absoluteTime",
    ]:  ## iterate between divergence and absolute time
        if branch_unit in json_translation:  ## it's available in tree
            for k in ll.Objects:  ## iterate over all branches
                cur_branch = getattr(k, branch_unit)  ## get parameter for this branch
                par_branch = getattr(k.parent, branch_unit)  ## get parameter for parental branch
                k.length = (
                    cur_branch - par_branch if cur_branch and par_branch else 0.0
                )  ## difference between current and parent is branch length (or, if parent unavailabel it's 0)

    if verbose:
        print("Traversing and drawing tree")

    ll.traverse_tree(verbose=verbose)
    ll.drawTree()
    if stats:
        ll.treeStats()  ## initial traversal, checks for stats
    if sort:
        ll.sortBranches()  ## traverses tree, sorts branches, draws tree

    cmap = {}
    for colouring in json_meta["colorings"]:
        if colouring["type"] == "categorical" and "scale" in colouring:
            cmap[colouring["key"]] = {}
            for entry in colouring["scale"]:
                key, value = entry
                cmap[colouring["key"]][key] = value
    setattr(ll, "cmap", cmap)

    return ll, json_meta


if __name__ == "__main__":
    import sys

    ll = make_tree(sys.argv[1])
    assert ll is not None
    ll.traverse_tree()
    sys.stdout.write("%s\n" % (ll.treeHeight))
