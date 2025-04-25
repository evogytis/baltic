import copy
import math
import logging
import random
import string
from functools import reduce
import numpy as np
from matplotlib.collections import LineCollection
from .node import Node
from .leaf import Leaf
from .clade import Clade
from .reticulation import Reticulation
from .bt_utils import root_to_tip

logger = logging.getLogger("baltic.Tree")


class Tree:  ## tree class

    def __init__(self):
        self.curNode = Node()  ## current node is a new instance of a node class
        self.curNode.index = "Root"  ## first object in the tree is the root
        self.curNode.length = 0.0  ## startind node branch length is 0
        self.curNode.height = 0.0  ## starting node height is 0
        self.root = None  # self.cur_node ## root of the tree is current node
        self.Objects = []  ## tree objects have a flat list of all branches in them
        self.tipMap = None
        self.treeHeight = (
            0  ## tree height is the distance between the root and the most recent tip
        )
        self.mostRecent = None
        self.ySpan = 0.0
        self.treeType = None  # indicates if the tree is a time or divergence tree

    def add_reticulation(self, name):
        logger.info(f"Creating new reticulation: {name}.")
        ret = Reticulation(name)
        ret.index = name
        ret.parent = self.curNode
        self.curNode.children.append(ret)
        self.Objects.append(ret)
        self.curNode = ret

    def add_node(self, i):
        logger.info(f"Creating new node: {i}.")
        newNode = Node()  ## new node instance
        newNode.index = i  ## new node's index is the position along the tree string
        if self.root is None:
            self.root = newNode
            self.root.length = 0.0

        newNode.parent = self.curNode  ## new node's parent is current node
        if not self.curNode.is_node():
            logger.error("Attempted to add a child to a non-node object.")
            logger.error(
                "Check if tip names have illegal characters like parentheses or commas."
            )
            raise TypeError()
        self.curNode.children.append(newNode)  ## new node is a child of current node
        self.curNode = newNode  ## current node is now new node
        self.Objects.append(
            self.curNode
        )  ## add new node to list of objects in the tree

    def add_leaf(self, i, name):
        logger.info(f"Creating new leaf: {name}.")
        newLeaf = Leaf(name)  ## new instance of leaf object
        newLeaf.index = i  ## index is position along tree string
        if self.root is None:
            self.root = newLeaf

        newLeaf.parent = self.curNode  ## leaf's parent is current node
        if not self.curNode.is_node():
            logger.error("Attempted to add a child to a non-node object.")
            logger.error(
                "Check if tip names have illegal characters like parentheses or commas."
            )
            raise TypeError()
        self.curNode.children.append(newLeaf)  ## assign leaf to parent's children
        self.curNode = newLeaf  ## current node is now new leaf
        self.Objects.append(self.curNode)  ## add leaf to all objects in the tree

    def subtree(self, startingNode=None, traverseCondition=None, stem=True):
        logger.info("Generating subtree.")
        if startingNode is None:
            logger.debug("No startingNode given, using root as default.")
            startingNode = self.root
        if traverseCondition is None:
            logger.debug("No traverseCondition given, using all branches.")
            traverseCondition = lambda k: True

        node = (
            startingNode.parent if stem else startingNode
        )  ## move up a node if we want the stem

        subtreeBranches = self.traverse_tree(
            node, includeCondition=lambda k: True, traverseCondition=traverseCondition
        )
        subtreeBranches = copy.deepcopy(subtreeBranches)

        logger.debug("Using stem.")
        if stem:  ## using stem - need to prune subtrees from root now
            unwantedBranches = []
            for child in node.children:  ## iterate over parent's children
                if (
                    child.index != startingNode.index
                ):  ## not at focal branch (unwanted sibling)
                    unwantedBranches += self.traverse_tree(
                        child, includeCondition=lambda w: True
                    )  ## add all branches resulting from traversals of unwanted siblings

            ## iterate over subtree branches, remember those that belong to unwanted subtrees
            remove = []
            for k in subtreeBranches:
                if k.index in [ub.index for ub in unwantedBranches]:
                    remove.append(k)
            for r in remove:  ## iterate over branches belong to unwanted subtrees
                subtreeBranches.remove(r)  ## remove from list

        ## nothing found or no leaf objects in traversal
        if subtreeBranches is None or not any(
            node.is_leaf() for node in subtreeBranches
        ):
            logger.error("No branches found in subtree traversal. Exiting.")
            return

        logger.debug("Creting new Tree object for subtree and assigning new branches.")
        localTree = Tree()  ## create a new tree object where the subtree will be
        localTree.Objects = subtreeBranches  ## assign branches to new tree object
        localTree.root = subtreeBranches[0]  ## root is the beginning of the traversal
        localTree.root.parent = None  ## tree begins strictly at node

        subtreeSet = set(
            subtreeBranches
        )  ## turn branches into set for quicker look up later

        logger.debug("Didn't use default traversal. ")
        ## didn't use default traverse condition,
        ## might need to deal with hanging nodes and prune children
        if traverseCondition is not None:
            for nd in localTree.get_internal():  ## iterate over nodes
                ## only keep children seen in traversal
                nd.children = [child for child in nd.children if child in subtreeSet]
            localTree.fix_hanging_nodes()

        logger.debug("Copying tipMap from original tree to subtree.")
        if self.tipMap:  ## if original tree has a tipMap dictionary
            localTipMap = {}
            ## copy over the relevant tip translations
            for tipNum in self.tipMap:
                if self.tipMap[tipNum] in [w.name for w in localTree.get_external()]:
                    localTipMap[tipNum] = self.tipMap[tipNum]
            localTree.tipMap = localTipMap
            # TODO: Make sure I broke this down correctly
            # localTree.tipMap={tipNum: self.tipMap[tipNum] for tipNum in self.tipMap if self.tipMap[tipNum] in [w.name for w in localTree.get_external()]}

        return localTree

    def make_single_type(self):
        logger.info("Converting to single-type tree.")
        while True:
            multitypeNodes = self.get_internal(lambda k: len(k.children) == 1)

            logger.debug(f"Number of multitype nodes in tree: {len(multitypeNodes)}")
            if not multitypeNodes:
                break

            for k in sorted(multitypeNodes, key=lambda x: -x.height):

                child = k.children[0]  ## fetch child
                grandparent = (
                    k.parent if k.parent.index else self.root
                )  ## fetch grandparent
                logger.debug(
                    f"At multitype node {k.index} with child {child.index} and grandparent {grandparent.index}"
                )
                child.parent = grandparent  ## child's parent is now grandparent
                grandparent.children.append(
                    child
                )  ## add child to grandparent's children
                grandparent.children.remove(
                    k
                )  ## remove old parent from grandparent's children
                grandparent.children = list(set(grandparent.children))

                child.length += k.length  ## adjust child length

                self.Objects.remove(k)  ## remove old parent from all objects
        self.sort_branches()

    def set_absolute_time(self, mostRecentSamplingDate):
        logger.debug("Setting absoluteTime for all branches.")
        logger.debug(f"MRSD: {mostRecentSamplingDate}")
        for k in self.Objects:  ## iterate over all objects
            k.absoluteTime = (
                mostRecentSamplingDate - self.treeHeight + k.height
            )  ## heights are in units of time from the root
        self.mostRecent = max(k.absoluteTime for k in self.Objects)

    def rescale(self, factor):
        logger.debug(f"Rescaling tree by factor of {factor}.")
        for k in self.Objects:
            k.length = k.length * factor
        self.traverse_tree()

    def treeStats(self):
        stats = self._calculate_tree_stats()
        print(
            "\nTree height: %.6f\nTree length: %.6f"
            % (stats["treeHeight"], stats["treeLength"])
        )  ## report the height and length of tree

        if stats["strictlyBifurcating"]:
            print("strictly bifurcating tree")
        if stats["multitypeTree"]:
            print("multitype tree")
        if stats["singletonTree"]:
            print("singleton tree")
        if stats["hasTraits"]:
            print("annotations present")

        print(
            "\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n"
            % (stats["numObjects"], stats["numNodes"], stats["numLeaves"])
        )  ## report numbers of different objects in the tree

    def treeStatsDict(self):
        return self._calculate_tree_stats()

    def _calculate_tree_stats(self):
        logger.debug("Calculating tree statistics.")
        stats = {}

        self.traverse_tree()  ## traverse the tree
        obs = self.Objects  ## convenient list of all objects in the tree
        stats["treeHeight"] = self.treeHeight
        stats["treeLength"] = sum([x.length for x in obs])

        nodes = self.get_internal()  ## get all nodes
        strictlyBifurcating = all(
            len(x.children) == 2 for x in nodes
        )  ## assume tree is not strictly bifurcating
        multiType = any(len(x.children) == 1 for x in nodes)
        singleton = len(nodes) == 0

        hasTraits = False  ## assume tree has no annotations
        max_annotations = max(
            len(k.traits) for k in obs
        )  ## check the largest number of annotations any branch has
        hasTraits = max_annotations > 0

        stats["strictlyBifurating"] = bool(strictlyBifurcating)
        stats["multitypeTree"] = bool(multiType)
        stats["singletonTree"] = bool(singleton)
        stats["hasTraits"] = bool(hasTraits)

        stats["numObjects"] = len(obs)
        stats["numNodes"] = len(nodes)
        stats["numLeaves"] = len(self.get_external())

        return stats

    def traverse_tree(
        self, curNode=None, includeCondition=None, traverseCondition=None, collect=None
    ):
        logger.info("Beginning tree traversal.")
        if curNode is None:  ## if no starting point defined - start from root
            logger.debug("Initiated traversal from root.")
            curNode = self.root

            if (
                traverseCondition is None and includeCondition is None
            ):  ## reset heights if traversing from scratch
                for k in self.Objects:  ## reset various parameters
                    if k.is_node():
                        k.leaves = set()
                        k.childHeight = None
                    k.height = None

        if traverseCondition is None:
            logger.debug("No traverseCondition given, using all branches.")
            traverseCondition = lambda k: True
        if includeCondition is None:
            logger.debug("No includeCondition given, using all leaflike.")
            includeCondition = lambda k: k.is_leaflike()

        if collect is None:  ## initiate collect list if not initiated
            collect = []

        if (
            curNode.parent and curNode.height is None
        ):  ## cur_node has a parent - set height if it doesn't have it already
            curNode.height = curNode.length + curNode.parent.height
        elif (
            curNode.height is None
        ):  ## cur_node does not have a parent (root), if height not set before it's zero
            curNode.height = 0.0

        # TODO: deprecate branchType
        logger.debug(f"At {curNode.index} ({curNode.branchType})")

        if includeCondition(curNode):  ## test if interested in cur_node
            collect.append(curNode)  ## add to collect list for reporting later

        if (
            curNode.is_leaf() and self.root != curNode
        ):  ## cur_node is a tip (and tree is not single tip)
            curNode.parent.leaves.add(curNode.name)  ## add to parent's list of tips

        elif curNode.is_node():  ## cur_node is node
            for child in filter(
                traverseCondition, curNode.children
            ):  ## only traverse through children we're interested
                logger.debug(f"Visiting child {child.index}.")
                self.traverse_tree(
                    curNode=child,
                    includeCondition=includeCondition,
                    traverseCondition=traverseCondition,
                    collect=collect,
                )  ## recurse through children
                logger.debug(f"child {child.index} done.")
            if len(curNode.children) < 1:
                logger.error(
                    f"Tried traversing through hanging node without children. Index: {curNode.index}."
                )
                raise AttributeError
            curNode.childHeight = max(
                [
                    child.childHeight if child.is_node() else child.height
                    for child in curNode.children
                ]
            )

            if curNode.parent:
                curNode.parent.leaves = curNode.parent.leaves.union(
                    curNode.leaves
                )  ## pass tips seen during traversal to parent
            self.treeHeight = (
                curNode.childHeight
            )  ## it's the highest child of the starting node
        return collect

    def rename_tips(self, tipNameMap=None):
        if tipNameMap is None:
            if self.tipMap is not None:
                logger.debug("No tipNameMap given, using tree.tipMap.")
                tipNameMap = self.tipMap
            else:
                raise ValueError("No dictionary provided for renaming tips.")

        for k in self.get_external():  ## iterate through leaf objects in tree
            k.name = tipNameMap[k.name]  ## change its name

    def reroot(self, branch=None, branchFrac=0.5, fixSingletons=True):
        if self.treeType != "time":
            logger.error("Cannon reroot on a time-calibrated tree.")
            raise AttributeError

        if branch == self.root:
            logger.warning("Rerooting attempted on existing root.")
            return self

        oldTreeLength = sum(self.get_parameter_list("length"))  ## get old tree length
        ###############
        if branch is None:  ## midpoint rooting
            logger.debug("No branch provided, rooting at midpoint")
            # Identify the largest pairwise distance
            maxDistance = 0.0

            for tip in self.get_external():  ## iterate over tips

                self.reroot(branch=tip, branchFrac=0.0)
                self.traverse_tree()  ## set heights

                highestTip = sorted(self.get_external(), key=lambda k: k.height)[-1]
                newMax = highestTip.height
                if newMax > maxDistance:  ## check if current height is bigger
                    tip1 = tip
                    tip2 = highestTip
                    maxDistance = newMax
                    # logger.info(f'setting {tip2.name} as current highest tip at height {max_distance}')

            # TODO: make sure that tip1 and tip2 always get assigned
            logger.debug(f"Midpoint rooting: Rooting on highest tip under current topology ({tip1.name}).")
            self.reroot(branch=tip1, branchFrac=0.0)

            # Depth to go from the ingroup tip toward the outgroup tip
            rootRemainder = 0.5 * maxDistance  # - (self.root.length or 0))
            assert rootRemainder >= 0
            # Identify the midpoint and reroot there.
            # Trace the path to the outgroup tip until all of the root depth has
            # been traveled/accounted for.
            path = tip2.get_path_to_root()[:-1]

            for node in path[::-1]:  ## iterate from old root to new
                rootRemainder -= node.length
                #             logger.info(f'iterating over path: {node.index} root remainder: {root_remainder}')
                if rootRemainder < 0:
                    outgroup_node = node
                    branchFrac = -rootRemainder / outgroup_node.length
                    break
            logger.debug(f"Midpoint rooting: rooting on node {outgroup_node.index} halfway from previous highest tip and current highest tip {tip2.name}.")
            self.reroot(
                branch=outgroup_node, branchFrac=branchFrac, fixSingletons=fixSingletons
            )
            logger.debug("Finished midpoint rooting.")
            return self

        ##############
        path = branch.get_path_to_root()[
            :-1
        ]  ## get path from new root to old root, ignore actual root node
        path = path[-2::-1]  ## invert

        ogLen = float(branch.length)  ## store branch length on the new root branch
        prevBranchLength = ogLen
        branch.length = ogLen * (
            branchFrac
        )  ## this modifies unchanged path branch length
        ##################

        characters = string.ascii_letters + string.digits  # A-Z, a-z, 0-9
        randomString = "".join(random.choices(characters, k=10))

        newRoot = Node()  ## create new root
        newRoot.index = f"new_root_{branch.index}_{branchFrac}_{randomString}"
        newRoot.length = 0.0
        newRoot.children.append(branch)

        logger.debug(f"Created new root {newRoot.index}")

        ######################
        if (
            len(path) == 1
        ):  ## new root is the old root (but maybe branch length is adjusted)
            newParent = newRoot
        else:  ##
            parent = path.pop(-2)  ## get previous branch
            parent.children.remove(branch)  ## remove child

            prevBranchLength, parent.length = parent.length, (
                prevBranchLength - branch.length
            )  ## store previous branch length, assign new branch length

            newRoot.children.append(parent)  ## add branch as child of new root

            parent.parent = newRoot  ## assign new parent
            newParent = parent  ## set new_parent variable

        #######################
        for parent in path[-2::-1]:  ## iterate from new root to old root
            parent.children.remove(newParent)  ## remove parent (used to be child)
            parent.parent = newParent  ## assign child as parent
            prevBranchLength, parent.length = (
                parent.length,
                prevBranchLength,
            )  ## pass branch length to the next node
            newParent.children.append(parent)  ## add former parent as child
            newParent = parent  ## move up to next node
        #################

        logger.debug("Cleaning up old root")

        oldRoot = self.root  ## get old root

        if (
            branch in oldRoot.children
        ):  ## remove child that was on the way to the new root from old root's children
            assert len(path) == 1
            oldRoot.children.remove(branch)
        else:
            oldRoot.children.remove(newParent)

        if (
            branch.parent == self.root
        ):  ## if we were rooting on a branch connected to the root
            oldRoot.length = prevBranchLength * (1 - branchFrac)
        else:
            oldRoot.length = prevBranchLength

        newParent.children.append(oldRoot)  ## add old root as child
        oldRoot.parent = newParent  ## set parent as old root
        ########

        logger.debug("Setting up root parent node (it's a baltic thing)")
        self.root = newRoot  ## set root to new root

        rootParent = (
            Node()
        )  ## baltic tree roots have an additional parent node to the root
        rootParent.index = "Root"
        rootParent.length = 0.0
        rootParent.height = 0.0

        rootParent.children.append(self.root)  ## set up basal baltic root parent node
        self.root.parent = rootParent
        self.Objects.append(self.root)

        branch.parent = newRoot
        #############
        if (
            fixSingletons
        ):  ## fixing singleton nodes (node's with 1 child, i.e. old root)
            logger.debug("Fixing singletons")
            self.make_single_type()  ## tree class can handle it

        ############
        newTreeLength = sum(
            [k.length for k in self.traverse_tree(includeCondition=lambda q: True)]
        )  ## get new tree length
        assert math.isclose(
            oldTreeLength, newTreeLength
        ), f"Tree length changed after rerooting - was {oldTreeLength}, now {newTreeLength}"  ## check if tree length hasn't changed after rerooting

        self.traverse_tree()

        return self

    def root_by_regression(self, stat="r^2", forcePositive=True):
        validRootingMethods = ["r^2", "correlation", "sum of squares"]
        assert (
            stat in validRootingMethods
        ),f"Invalid option for root-to-tip regression: {stat} (options are {validRootingMethods})"

        cll = copy.deepcopy(self)  ## deepcopy entire tree

        res = {}  ## will be used to track better regressions with new roots
        for k in self.Objects[1:]:  ## iterate over branches in tree (ignore root)
            cll = cll.reroot(
                [w for w in cll.Objects if w.index == k.index][0],
                fixSingletons=False,
            )  ## reroot on branch provided
            tips = cll.get_external()  ## get all tips

            xs, ys = zip(
                *[(w.absoluteTime, w.height) for w in tips]
            )  ## get collection dates
            res = root_to_tip(
                k, xs, ys, res, stat=stat, forcePositive=forcePositive
            )  ## check root-to-tip regression, update res if better
        ##############
        logger.debug(
                f"Rooting tree first time at node {res['root']} with regression: {res}"
            )
        self = self.reroot(
            res["root"]
        )  ## reroot on branch optimising stat
        ############
        if (
            len(self.root.children) == 2
        ):  ## if root is strictly bifurcating - also optimise branch fraction
            logger.debug("Bifurcating root, finding more precise rooting along new root")
            res["frac"] = 0.0  ## add frac argument
            left, right = self.root.children  ## get left and right children

            totalBranch = (
                left.length + right.length
            )  ## total amount of branch length adjustment available at root
            leftSubtree = self.traverse_tree(left)  ## get left subtree
            rightSubtree = self.traverse_tree(right)  ## get right subtree

            lxs = [
                k.absoluteTime for k in leftSubtree
            ]  ## get x coordinates for left subtree #TODO: these don't get used
            rxs = [
                k.absoluteTime for k in rightSubtree
            ]  ## get x coordinates for right subtree #TODO: see above

            for f in np.linspace(
                0, 1, 21
            ):  ## check root positions along best overall root
                adjustL = (
                    -left.length + totalBranch * f
                )  ## height adjustments for left subtree
                adjustR = -right.length + totalBranch * (
                    1 - f
                )  ## height adjustments for right subtree

                lys = [
                    k.height + adjustL for k in leftSubtree
                ]  ## adjust left subtree heights #TODO: see above
                rys = [
                    k.height + adjustR for k in rightSubtree
                ]  ## same for right #TODO: see above

                res = root_to_tip(
                    left, xs, ys, res, stat=stat, forcePositive=forcePositive, frac=f
                )  ## check regression

            self = self.reroot(
                branch=res["root"], branchFrac=res["frac"]
            )  ## reroot on branch optimising stat
        #############
        logger.info(f"Correlation coefficient: {res['correlation']}")
        logger.info(f"Sum of squares: {res['sum of squares']}")
        logger.info(f"Evolutionary rate: {res['slope']}")
        logger.info(f"Intercept (TMRCA): {max(xs) - res['intercept']}")
        logger.info(f"r^2: {res['r^2']}")

        return self

    def sort_branches(self, descending=True, sortFxn=None, sortByHeight=True):
        mod = -1 if descending else 1
        if sortFxn is None:
            sortFxn = lambda k: (
                (k.is_node(), -len(k.leaves) * mod, k.length * mod)
                if k.is_node()
                else (k.is_node(), k.length * mod)
            )
        if sortByHeight:  # Sort nodes by height and group nodes and leaves together
            """Sort descendants of each node."""

            for k in self.get_internal():  ## iterate over nodes
                k.children = sorted(k.children, key=sortFxn)
        else:  # Do not sort by height. Retain leaves at original positions. Only sort nodes
            for k in self.get_internal():
                leavesIdx = [
                    (i, ctr) for ctr, i in enumerate(k.children) if i.is_leaflike()
                ]  # Get original indices of leaves
                nodes = sorted(
                    [x for x in k.children if x.is_node()], key=sortFxn
                )  # Sort nodes only by number of descendants
                children = nodes
                for i in leavesIdx:  # Insert leaves back into same positions
                    children.insert(i[1], i[0])
                k.children = children
        self.assign_tree_coordinates()  ## update x and y positions of each branch, since y positions will have changed because of sorting

    def assign_tree_coordinates(
        self, order=None, widthFxn=None, padNodes=None
    ):
        """ formerly drawTree()
        """
        if order is None:
            order = self.traverse_tree(
                includeCondition=lambda k: k.is_leaflike()
            )  ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
            logger.debug("Drawing tree in pre-order")
        else:
            logger.debug("Drawing tree with provided order")

        nameOrder = {x.name: i for i, x in enumerate(order)}
        assert len(nameOrder) == len(order), "Non-unique names present in tree"
        if widthFxn is None:
            logger.debug("Drawing tree with default widths (1 unit for leaf objects, width+1 for clades)")
            skips = [1 if isinstance(x, Leaf) else x.width + 1 for x in order]
        else:
            skips = list(map(widthFxn, order))

        for k in self.Objects:  ## reset coordinates for all objects
            k.x = None
            k.y = None

        drawn = {}  ## drawn keeps track of what's been drawn
        for k in order:  ## iterate over tips
            x = k.height  ## x position is height
            y_idx = nameOrder[k.name]  ## assign y index
            y = (
                sum(skips[y_idx:]) - skips[y_idx] / 2.0
            )  ## sum across skips to find y position

            k.x = x  ## set x and y coordinates
            k.y = y
            drawn[k.index] = None  ## remember that this objects has been drawn

        if padNodes is not None:  ## will be padding nodes
            for n in padNodes:  ## iterate over nodes whose descendants will be padded
                idx = (
                    sorted([nameOrder[lf] for lf in n.leaves])
                    if n.is_node()
                    else [order.index(n)]
                )  ## indices of all tips to be padded
                for i, k in enumerate(order):  ## iterate over all tips

                    if i < idx[0]:  ## tip below clade
                        k.y += padNodes[n]  ## pad

                    if (i - 1) < idx[-1]:  ## tip above clade
                        k.y += padNodes[n]  ## pad again

            allYs = filter(
                None, self.get_parameter_list("y")
            )  ## get all y positions in tree that aren't None
            minY = min(allYs)  ## get min
            for k in self.get_external():  ## reset y positions so tree starts at y=0.5
                k.y -= minY - 0.5

        assert len([k for k in self.Objects if k.is_leaflike()]) == len(
            order
        ), "Number of tips in tree does not match number of unique tips, check if two or more collapsed clades were assigned the same name."
        storePlotted = 0

        while len(drawn) != len(
            self.Objects
        ):  # keep drawing the tree until everything is drawn
            logger.debug("Drawing iteration %d" % (len(drawn)))
            for k in filter(
                lambda w: w.index not in drawn, self.get_internal()
            ):  ## iterate through internal nodes that have not been drawn
                if len([q.y for q in k.children if q.y is not None]) == len(
                    k.children
                ):  ## all y coordinates of children known
                    logger.debug(f"Setting node {k.index} coordinates to"),
                    x = k.height  ## x position is height
                    childrenYCoords = [
                        q.y for q in k.children if q.y is not None
                    ]  ## get all existing y coordinates of the node
                    y = sum(childrenYCoords) / float(
                        len(childrenYCoords)
                    )  ## internal branch is in the middle of the vertical bar
                    k.x = x
                    k.y = y
                    drawn[k.index] = None  ## remember that this objects has been drawn
                    logger.debug(f"{k.y} ({len(drawn)} branches drawn)")
                    minYRange = min(
                        [
                            min(child.yRange) if child.is_node() else child.y
                            for child in k.children
                        ]
                    )  ## get lowest y coordinate across children
                    maxYRange = max(
                        [
                            max(child.yRange) if child.is_node() else child.y
                            for child in k.children
                        ]
                    )  ## get highest y coordinate across children
                    setattr(
                        k, "yRange", [minYRange, maxYRange]
                    )  ## assign the maximum extent of children's y coordinates

            if len(self.Objects) > len(drawn):
                assert len(drawn) > storePlotted, (
                    "Got stuck trying to find y positions of objects (%d branches drawn this iteration, %d branches during previous iteration out of %d total)"
                    % (len(drawn), storePlotted, len(self.Objects))
                )
            storePlotted = len(
                drawn
            )  ## remember how many branches were drawn this iteration

        yValues = [k.y for k in self.Objects]  ## all y values
        self.ySpan = (
            max(yValues) - min(yValues) + min(yValues) * 2
        )  ## determine appropriate y axis span of tree

        if self.root.is_node():
            self.root.x = min(
                [q.x - q.length for q in self.root.children if q.x is not None]
            )  ## set root x and y coordinates
            childrenYCoords = [q.y for q in self.root.children if q.y is not None]
            self.root.y = sum(childrenYCoords) / float(len(childrenYCoords))
        else:
            self.root.x = self.root.length

    def assign_unrooted_tree_coordinates(self, rotate=0.0, n=None, total=None):
        if n is None:
            total = sum(
                [1 if x.is_leaf() else x.width + 1 for x in self.get_external()]
            )
            n = self.root  # .children[0]
            for k in self.Objects:
                k._tau = 2 * math.pi * rotate
                k.x = 0.0
                k.y = 0.0

        w = (
            2 * math.pi * 1.0 / float(total)
            if n.is_leaf()
            else 2 * math.pi * len(n.leaves) / float(total)
        )

        if n.parent.x is None:
            n.parent.x = 0.0
            n.parent.y = 0.0

        n.x = n.parent.x + n.length * math.cos(n.traits["tau"] + w * 0.5)
        n.y = n.parent.y + n.length * math.sin(n.traits["tau"] + w * 0.5)
        eta = n._tau

        if n.is_node():
            for ch in n.children:
                w = (
                    2 * math.pi * 1.0 / float(total)
                    if ch.is_leaf()
                    else 2 * math.pi * len(ch.leaves) / float(total)
                )

                ch._tau = eta
                eta += w
                self.assign_unrooted_tree_coordinates(rotate, ch, total)

    def find_MRCA(self, descendants):
        assert (
            len(descendants) > 1
        ), "Not enough descendants to find common ancestor: %d" % (len(descendants))
        pathsToRoot = {
            k.index: set(k.get_path_to_root()) for k in descendants
        }  ## for every descendant create an empty set
        # for k in descendants: ## iterate through every descendant
        #     curNode=k ## start descent from descendant
        #     while curNode: ## while not at root
        #         pathsToRoot[k.index].add(curNode) ## remember every node visited along the way
        #         curNode=curNode.parent ## descend

        return sorted(
            reduce(set.intersection, pathsToRoot.values()),
            key=lambda k: (-len(k.leaves), k.height),
        )[
            -1
        ]  ## return the most recent branch that is shared across all paths to root

    def collapse_subtree_to_clade(
        self, cl, givenName, widthFunction=lambda k: len(k.leaves)
    ):
        assert cl.is_node(), "Cannot collapse non-node class"
        collapsedClade = Clade(givenName)
        collapsedClade.index = cl.index
        collapsedClade.leaves = cl.leaves
        collapsedClade.length = cl.length
        collapsedClade.height = cl.height
        collapsedClade.parent = cl.parent
        collapsedClade.absoluteTime = cl.absoluteTime
        collapsedClade.traits = cl.traits
        collapsedClade.width = widthFunction(cl)

        logger.debug(f"Replacing node {cl.index} (parent {cl.parent.index}) with a clade class")
        parent = cl.parent

        removeFromTree = self.traverse_tree(cl, includeCondition=lambda k: True)
        collapsedClade.subtree = removeFromTree
        assert len(removeFromTree) < len(
            self.Objects
        ), "Attempted collapse of entire tree"
        collapsedClade.lastHeight = max([x.height for x in removeFromTree])
        if [x.absoluteTime for x in removeFromTree].count(None) != len(removeFromTree):
            collapsedClade.lastAbsoluteTime = max(
                [x.absoluteTime for x in removeFromTree]
            )

        for k in removeFromTree:
            self.Objects.remove(k)

        parent.children.remove(cl)
        parent.children.append(collapsedClade)
        self.Objects.append(collapsedClade)
        collapsedClade.parent = parent
        if self.tipMap is not None:
            self.tipMap[givenName] = givenName

        self.traverse_tree()
        self.sort_branches()
        return collapsedClade

    def restore_all_collapsed_subtrees(self):
        while len([k for k in self.Objects if isinstance(k, Clade)]) > 0:
            clades = [k for k in self.Objects if isinstance(k, Clade)]
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

    def collapse_branches(
        self,
        collapseIfFxn=lambda x: x.traits["posterior"] <= 0.5,
        designatedNodes=[],
    ):
        """
        TODO: change desnignatedNode to None
        """
        newTree = copy.deepcopy(self)  ## work on a copy of the tree
        if (
            len(designatedNodes) == 0
        ):  ## no nodes were designated for deletion - relying on anonymous function to collapse nodes
            nodesToDelete = list(
                filter(
                    lambda n: n.is_node()
                    and collapseIfFxn(n) == True
                    and n != newTree.root,
                    newTree.Objects,
                )
            )  ## fetch a list of all nodes who are not the root and who satisfy the condition
        else:
            assert len([w for w in designatedNodes if w.is_node()]) == len(
                designatedNodes
            ), "Non-node class detected in list of nodes designated for deletion"
            assert (
                len([w for w in designatedNodes if w != newTree.root]) == 0
            ), "Root node was designated for deletion"

            nodesToDelete = list(
                filter(
                    lambda w: w.index in [q.index for q in designatedNodes],
                    newTree.Objects,
                )
            )  ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
        logger.debug(f"{len(nodesToDelete)} nodes set for collapsing: {[w.index for w in nodesToDelete]}")
        assert (
            len(nodesToDelete) < len(newTree.get_internal()) - 1
        ), "Chosen cutoff would remove all branches"
        while (
            len(nodesToDelete) > 0
        ):  ## as long as there are branches to be collapsed - keep reducing the tree

            logger.debug(f"Continuing collapse cycle, {len(nodesToDelete)} nodes left")
            for k in sorted(
                nodesToDelete, key=lambda x: -x.height
            ):  ## start with branches near the tips
                zeroNode = k.children  ## fetch the node's children
                k.parent.children += zeroNode  ## add them to the zero node's parent
                oldParent = k  ## node to be deleted is the old parent
                newParent = (
                    k.parent
                )  ## once node is deleted, the parent to all their children will be the parent of the deleted node
                if newParent is None:
                    newParent = self.root
                logger.debug(f"Removing node {oldParent.index}, attaching children {[w.index for w in k.children]} to node {newParent.index}")
                for (
                    w
                ) in (
                    newTree.Objects
                ):  ## assign the parent of deleted node as the parent to any children of deleted node
                    if w.parent == oldParent:
                        w.parent = newParent
                        w.length += oldParent.length
                        logger.debug(f"Fixing branch length for node {w.index}")
                k.parent.children.remove(
                    k
                )  ## remove traces of deleted node - it doesn't exist as a child, doesn't exist in the tree and doesn't exist in the nodes list
                newTree.Objects.remove(k)

                nodesToDelete.remove(k)  ## in fact, the node never existed

                if len(designatedNodes) == 0:
                    nodesToDelete == list(
                        filter(
                            lambda n: n.is_node()
                            and collapseIfFxn(n) == True
                            and n != newTree.root,
                            newTree.Objects,
                        )
                    )
                else:
                    assert len([w for w in designatedNodes if w.is_node()]) == len(
                        designatedNodes
                    ), "Non-node class detected in list of nodes designated for deletion"
                    assert (
                        len([w for w in designatedNodes if w != newTree.root]) == 0
                    ), "Root node was designated for deletion"
                    nodesToDelete = [
                        w
                        for w in newTree.Objects
                        if w.index in [q.index for q in designatedNodes]
                    ]

                logger.debug(f"Removing references to node {k.index}")
        newTree.sort_branches()  ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree  ## return collapsed tree

    def to_string(
        self,
        curNode=None,
        traits=None,
        nexus=False,
        stringFragment=None,
        traverseCondition=None,
        rename=None,
        quoteCharacter="'",
        json=False,
    ):
        if curNode is None:
            curNode = self.root  # .children[-1]
        if traits is None:
            traits = set(
                sum([list(k.traits.keys()) for k in self.Objects], [])
            )  ## fetch all trait keys
        if stringFragment is None:
            stringFragment = []
            if nexus:
                assert not json, "Nexus format not a valid option for JSON output"
                logger.debug("Exporting to NEXUS format")
                stringFragment.append("#NEXUS\nBegin trees;\ntree TREE1 = [&R] ")
        if traverseCondition is None:
            traverseCondition = lambda k: True

        comment = []  ## will hold comment
        if len(traits) > 0:  ## non-empty list of traits to output
            for tr in traits:  ## iterate through keys
                if tr in curNode.traits:  ## if key is available
                    logger.debug(f"trait {tr} available for {curNode.index} ({curNode.branchType}) type: {type(curNode.traits[tr])}")
                    if isinstance(curNode.traits[tr], str):  ## string value
                        comment.append(f'{tr}="{curNode.traits[tr]}"')
                        logger.debug(f"adding string comment {comment[-1]}")
                    elif isinstance(curNode.traits[tr], float) or isinstance(
                        curNode.traits[tr], int
                    ):  ## float or integer
                        comment.append(f"{tr}={curNode.traits[tr]}")
                        logger.debug(f"adding numeric comment {comment[-1]}")
                    elif isinstance(curNode.traits[tr], list):  ## lists
                        rangeComment = []
                        for val in curNode.traits[tr]:
                            if isinstance(val, str):  ## string
                                rangeComment.append(f'"{val}"')
                            elif isinstance(val, float) or isinstance(
                                val, int
                            ):  ## float or integer
                                rangeComment.append(f"{val}")
                            elif isinstance(
                                val, list
                            ):  ## list of lists, example complete history annotated on tree
                                rangeComment.append("{{{}}}".format(",".join(val)))
                        comment.append(f"{tr}={','.join(rangeComment)}")
                        logger.debug(f"adding range comment {comment[-1]}")
                else:
                    logger.debug(f"trait {tr} unavailable for {curNode.index} ({curNode.branchType})")

        if curNode.is_node():
            logger.debug(f"node: {curNode.index}")
            stringFragment.append("(")
            traverseChildren = list(filter(traverseCondition, curNode.children))
            assert (
                len(traverseChildren) > 0
            ), f"Node {curNode.index} does not have traversable children"
            for c, child in enumerate(
                traverseChildren
            ):  ## iterate through children of node if they satisfy traverse condition
                logger.debug(f"moving to child {child.index} of node {curNode.index}")
                self.to_string(
                    curNode=child,
                    traits=traits,
                    nexus=nexus,
                    stringFragment=stringFragment,
                    traverseCondition=traverseCondition,
                    rename=rename,
                    quoteCharacter=quoteCharacter,
                )
                if (c + 1) < len(
                    traverseChildren
                ):  ## not done with children, add comma for next iteration
                    stringFragment.append(",")
            stringFragment.append(")")  ## last child, node terminates

        elif curNode.is_leaf():
            if rename is None:
                treeName = curNode.name  ## designated numName
            else:
                assert isinstance(rename, dict), 'Variable "rename" is not a dictionary'
                assert (
                    curNode.name in rename
                ), f"Tip name {curNode.name} not in rename dictionary"
                treeName = rename[curNode.name]

            logger.debug(f"leaf: {curNode.index} ({treeName})")
            stringFragment.append(f"{quoteCharacter}{treeName}{quoteCharacter}")

        if len(comment) > 0:
            logger.debug(f"adding comment to {curNode.index}")
            comment = ",".join(comment)
            comment = "[&" + comment + "]"
            stringFragment.append(f"{comment}")  ## end of node, add annotations

        logger.debug(f"adding branch length to {curNode.index}")
        stringFragment.append(
            ":%8f" % (curNode.length)
        )  ## end of node, add branch length

        if curNode == self.root:  # .children[-1]:
            stringFragment.append(";")
            if nexus == True:
                stringFragment.append("\nEnd;")
            logger.debug("finished")
            return "".join(stringFragment)

    def get_all_tip_TMRCAs(self):
        assert (
            self.treeType == "time"
        ), "Error: can only calculate TMRCA matrix on time-calibrated tree."

        tipNames = [k.name for k in self.get_external()]
        tmrcaMatrix = {
            x: {y: None if x != y else 0.0 for y in tipNames} for x in tipNames
        }  ## pairwise matrix of tips

        for k in self.get_internal():  ## iterate over nodes
            all_children = list(k.leaves)  ## fetch all descendant tips of node

            for a, tipA in enumerate(all_children):
                for tipB in all_children[a + 1 :]:
                    if (
                        tmrcaMatrix[tipA][tipB] is None
                        or tmrcaMatrix[tipA][tipB] <= k.absoluteTime
                    ):  ## if node's time is more recent than previous entry - set new TMRCA value for pair of tips
                        tmrcaMatrix[tipA][tipB] = k.absoluteTime
                        tmrcaMatrix[tipB][tipA] = k.absoluteTime
        return tmrcaMatrix

    def reduce_tree(self, tipsToKeep):
        assert len(tipsToKeep) > 0, "No tips given to reduce the tree to."
        assert (
            len([k for k in tipsToKeep if not k.is_leaflike()]) == 0
        ), "Embedding contains %d branches that are not leaf-like." % (
            len([k for k in tipsToKeep if k.is_leaflike() == False])
        )
        logger.debug("Preparing branch hash for keeping %d branches" % (len(tipsToKeep)))
        branchHash = {k.index: k for k in tipsToKeep}
        embedding = []
        logger.debug("Deep copying tree")
        reducedTree = copy.deepcopy(self)  ## new tree object
        for k in reducedTree.Objects:  ## deep copy branches from current tree
            if k.index in branchHash:  ## if branch is designated as one to keep
                currentBranch = k
                logger.debug(f"Traversing to root from {currentBranch.index}")
                while currentBranch != reducedTree.root:  ## descend to root
                    logger.debug(f"at {currentBranch.index} root: {currentBranch == reducedTree.root}")
                    embedding.append(currentBranch)  ## keep track of the path to root
                    currentBranch = currentBranch.parent
        embedding.append(reducedTree.root)  ## add root to embedding
        logger.debug(f"Finished extracting embedding with {len(embedding)} branches ({len([w for w in embedding if w.is_leaf()])} tips, {len([w for w in embedding if w.is_node()])} nodes)")
        embedding = set(embedding)  ## prune down to only unique branches

        reducedTree.Objects = sorted(
            list(embedding), key=lambda x: x.height
        )  ## assign branches that are kept to new tree's Objects
        logger.debug("Pruning untraversed lineages")
        for k in reducedTree.get_internal():  ## iterate through reduced tree
            k.children = [
                c for c in k.children if c in embedding
            ]  ## only keep children that are present in lineage traceback
        reducedTree.root.children = [
            c for c in reducedTree.root.children if c in embedding
        ]  ## do the same for root

        reducedTree.fix_hanging_nodes()

        logger.debug("Last traversal and branch sorting")
        reducedTree.traverse_tree()  ## traverse
        reducedTree.sort_branches()  ## sort

        return reducedTree  ## return new tree

    def count_lineages_at_time(
        self, t, timeAttr="absoluteTime", inclusionConditionFxn=lambda x: True
    ):
        return len(
            [
                k
                for k in self.Objects
                if getattr(k.parent, timeAttr) is not None
                and getattr(k.parent, timeAttr) < t <= getattr(k, timeAttr)
                and inclusionConditionFxn(k)
            ]
        )

    def get_external(self, filterFxn=None, onlyLeaves=True):
        externals = list(
            filter(
                filterFxn,
                filter(
                    lambda k: k.is_leaf() if onlyLeaves else k.is_leaflike(),
                    self.Objects,
                ),
            )
        )
        return externals

    def get_internal(self, filterFxn=None):
        internals = list(filter(filterFxn, filter(lambda k: k.is_node(), self.Objects)))
        return internals

    def get_branches(self, filterFxn=lambda x: True, failIfNoResults=True):
        select = list(filter(filterFxn, self.Objects))

        if len(select) == 0:
            if failIfNoResults:
                raise Exception(
                    "No branches satisfying function were found amongst branches"
                )
            else:
                logger.warning(
                    "No results found matching the specified condition. Returning empty list."
                )
                return []
        # elif (
        #     len(select) == 1
        #     # TODO: dd a get_branch() method that returns a single branch (not a list of one)
        #     return select[-1]
        else:
            return select

    def get_parameter_list(self, statistic, useTraitsDict=False, filterFxn=None):
        if filterFxn is None:
            branches = self.Objects
        else:
            branches = filter(filterFxn, self.Objects)

        if useTraitsDict:
            params = [k.traits[statistic] for k in branches if statistic in k.traits]
        else:
            params = [getattr(k, statistic) for k in branches if hasattr(k, statistic)]

        return params

    def fix_hanging_nodes(self):
        while True:
            hangingNodes = [
                node for node in self.Objects if node.is_node() and not node.children
            ]  ## nodes without children (hanging nodes)
            if not hangingNodes:
                break

            for node in hangingNodes:
                node.parent.children.remove(node)
                self.Objects.remove(node)



    ############################################
    ##########   PLOTTING FUNCTIONS   ##########
    ############################################
    def plot_text(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        textContentFxn=None,
        zorder=4,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        if textContentFxn is None:
            textContentFxn = lambda k: k.name

        localKwargs = dict(kwargs)
        if "verticalalignment" not in localKwargs:
            localKwargs["verticalalignment"] = "center"

        for k in filter(targetFxn, self.Objects):
            x, y = xCoordinateFxn(k), yCoordinateFxn(k)
            z = zorder
            ax.text(x, y, textContentFxn(k), zorder=z, **localKwargs)
        return ax

    def plot_text_unrooted(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        textContentFxn=None,
        zorder=4,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        if textContentFxn is None:
            textContentFxn = lambda k: k.name

        for k in filter(targetFxn, self.Objects):
            localKwargs = dict(kwargs)

            x, y = xCoordinateFxn(k), yCoordinateFxn(k)
            z = zorder

            assert (
                "tau" in k.traits
            ), "Branch does not have angle tau computed by assign_unrooted_tree_coordinates()."

            rot = math.degrees(k.traits["tau"]) % 360

            if "horizontalalignment" not in localKwargs:
                localKwargs["horizontalalignment"] = (
                    "right" if 90 < rot < 270 else "left"
                )
            if "verticalalignment" not in localKwargs:
                localKwargs["verticalalignment"] = "center"

            rot = rot + 180 if 90 < rot < 270 else rot

            ax.text(
                x,
                y,
                textContentFxn(k),
                rotation=rot,
                rotation_mode="anchor",
                zorder=z,
                **localKwargs,
            )

        return ax

    def plot_text_annotations_circular(
        self,
        ax,
        targetFxn=None,
        textContentFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        circStart=0.0,
        circFrac=1.0,
        inwardSpace=0.0,
        normaliseHeight=None,
        zorder=4,
        **kwargs,
    ):
        if targetFxn is None:
            targetFxn = lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        if textContentFxn is None:
            textContentFxn = lambda k: k.name
        if zorder is None:
            zorder = 4

        circS = circStart * math.pi * 2
        circ = circFrac * math.pi * 2

        allXs = list(map(xCoordinateFxn, self.Objects))
        if normaliseHeight is None:
            normaliseHeight = lambda value: (value - min(allXs)) / (
                max(allXs) - min(allXs)
            )

        for k in filter(targetFxn, self.Objects):  ## iterate over branches
            localKwargs = dict(kwargs)  ## copy global kwargs into a local version

            x = normaliseHeight(
                xCoordinateFxn(k) + inwardSpace
            )  ## get branch x position
            y = yCoordinateFxn(k)  ## get y position

            y = circS + circ * y / self.ySpan
            X = math.sin(y)
            Y = math.cos(y)

            rot = math.degrees(y) % 360

            ## Set local keyword arguments that get passed on to ax.text()
            if "horizontalalignment" not in localKwargs:
                if 180 < rot < 360:
                    localKwargs["horizontalalignment"] = "right"
                else:
                    localKwargs["horizontalalignment"] = "left"
            if "verticalalignment" not in localKwargs:
                localKwargs["verticalalignment"] = "center"
            rot = 360 - rot - 90 if 180 < rot < 360 else 360 - rot + 90

            ax.text(
                X * x,
                Y * x,
                textContentFxn(k),
                rotation=rot,
                rotation_mode="anchor",
                zorder=zorder,
                **localKwargs,
            )

        return ax

    def plot_points(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        pointSize=None,
        pointSizeFxn=None,
        colour=None,
        colourFxn=None,
        outline=None,
        outlineSize=None,
        outlineSizeFxn=None,
        outlineColour=None,
        outlineColourFxn=None,
        zorder=3,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        # Point size logic
        if pointSize is not None and pointSizeFxn is not None:
            raise ValueError(
                "Cannot specify both pointSize and pointSizeFxn. Please use only one."
            )
        if pointSize is None and pointSizeFxn is None:
            pointSize = lambda k: 40
        else:
            pointSize = pointSize if pointSize is not None else pointSizeFxn
        # Outline size logic
        if outlineSize is not None and outlineSizeFxn is not None:
            raise ValueError(
                "Cannot specify both outlineSize and outlineSizeFxn. Please use only one."
            )
        if outlineSize is None and outlineSizeFxn is None:
            outlineSize = lambda k: pointSize(k) * 2
        else:
            outlineSize = outlineSize if outlineSize is not None else outlineSizeFxn
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            )
        if colour is None and colourFxn is None:
            colour = lambda k: "k"
        else:
            colour = colour if colour is not None else colourFxn
        # Outline colour logic
        if outlineColour is not None and outlineColourFxn is not None:
            raise ValueError(
                "Cannot specify both outlineColour and outlineColourFxn. Please use only one."
            )
        if outlineColour is None and outlineColourFxn is None:
            outlineColour = lambda k: "k"
        else:
            outlineColour = outlineColour if outlineColour is not None else outlineColourFxn

        xs = []
        ys = []
        colours = []
        sizes = []

        outline_xs = []
        outline_ys = []
        outline_colours = []
        outline_sizes = []
        for k in filter(targetFxn, self.Objects):
            xs.append(xCoordinateFxn(k))
            ys.append(yCoordinateFxn(k))
            colours.append(colour(k)) #  if callable(colour) else colours.append(colour)
            sizes.append(pointSize(k)) #  if callable(size) else sizes.append(size)

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outlineColour(k)) #  if callable(outlineColour) else outline_colours.append(outlineColour)
                outline_sizes.append(outlineSize(k)) #  if callable(outlineSize) else outline_sizes.append(outlineSize)

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

    def plot_tree(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        width=None,
        widthFxn=None,
        colour=None,
        colourFxn=None,
        connection_type=None,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: True
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        # Width logic
        if width is not None and widthFxn is not None:
            raise ValueError(
                "Cannot specify both width and widthFxn. Please use only one."
            )
        if width is None and widthFxn is None:
            width = lambda k: 2
        else:
            width = width if width is not None else widthFxn
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            )
        if colour is None and colourFxn is None:
            colour = lambda k: "k"
        else:
            colour = colour if colour is not None else colourFxn
        if connection_type is None:
            connection_type = "baltic"
        assert connection_type in [
            "baltic",
            "direct",
            "elbow",
        ], f'Unrecognised drawing type \"{connection_type}\"'

        branches = []
        colours = []
        linewidths = []
        for k in filter(targetFxn, self.Objects):  ## iterate over branches
            x = xCoordinateFxn(k)  ## get branch x position
            xp = xCoordinateFxn(k.parent) if k.parent else x  ## get parent x position
            y = yCoordinateFxn(k)  ## get y position

            try:
                colours.append(colour(k)) # if callable(colour) else colours.append(colour)
            except KeyError:
                colours.append(
                    (0.7, 0.7, 0.7)
                )  ## in case no colour available for branch set it to grey
            linewidths.append(width(k)) #  if callable(width) else linewidths.append(width)

            if (
                connection_type == "baltic"
            ):  ## each node has a single vertical line to which descendant branches are connected
                branches.append(((xp, y), (x, y)))
                if k.is_node():
                    yl, yr = yCoordinateFxn(k.children[0]), yCoordinateFxn(
                        k.children[-1]
                    )  ## y positions of first and last child
                    branches.append(((x, yl), (x, yr)))
                    linewidths.append(linewidths[-1])
                    colours.append(colours[-1])
            elif (
                connection_type == "elbow"
            ):  ## more standard connection where each branch connects to its parent via a right-angled line
                yp = yCoordinateFxn(k.parent) if k.parent else y  ## get parent x position
                branches.append(((xp, yp), (xp, y), (x, y)))
            elif (
                connection_type == "direct"
            ):  ## this gives triangular looking trees where descendants connect directly to their parents
                yp = yCoordinateFxn(k.parent)  ## get y position
                branches.append(((xp, yp), (x, y)))
            else:
                pass  ## for now

        if "capstyle" not in kwargs:
            kwargs["capstyle"] = "projecting"
        line_segments = LineCollection(branches, lw=linewidths, color=colours, **kwargs)
        ax.add_collection(line_segments)
        return ax

    def plot_circular_tree(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        width=None,
        widthFxn=None,
        colour=None,
        colourFxn=None,
        circStart=0.0,
        circFrac=1.0,
        inwardSpace=0.0,
        normaliseHeight=None,
        precision=15,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: True
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        # Width logic
        if width is not None and widthFxn is not None:
            raise ValueError(
                "Cannot specify both width and widthFxn. Please use only one."
            )
        if width is None and widthFxn is None:
            width = lambda k: 2
        else:
            width = width if width is not None else widthFxn
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            )
        if colour is None and colourFxn is None:
            colour = lambda k: "k"
        else:
            colour = colour if colour is not None else colourFxn

        if inwardSpace < 0:
            inwardSpace -= self.treeHeight

        branches = []
        colours = []
        linewidths = []

        circ_s = circStart * math.pi * 2
        circ = circFrac * math.pi * 2

        allXs = list(map(xCoordinateFxn, self.Objects))
        if normaliseHeight is None:
            normaliseHeight = lambda value: (value - min(allXs)) / (
                max(allXs) - min(allXs)
            )
        linspace = lambda start, stop, n: (
            list(start + ((stop - start) / (n - 1)) * i for i in range(n))
            if n > 1
            else stop
        )

        for k in filter(targetFxn, self.Objects):  ## iterate over branches
            x = normaliseHeight(xCoordinateFxn(k) + inwardSpace)  ## get branch x position
            xp = (
                normaliseHeight(xCoordinateFxn(k.parent) + inwardSpace)
                if k.parent.parent
                else x
            )  ## get parent x position
            y = yCoordinateFxn(k)  ## get y position

            try:
                (
                    colours.append(colour(k))
                    if callable(colour)
                    else colours.append(colour)
                )
            except KeyError:
                colours.append((0.7, 0.7, 0.7))
            linewidths.append(width(k)) #  if callable(width) else linewidths.append(width)

            y = circ_s + circ * y / self.ySpan
            X = math.sin(y)
            Y = math.cos(y)
            branches.append(((X * xp, Y * xp), (X * x, Y * x)))

            if k.is_node():
                yl, yr = yCoordinateFxn(k.children[0]), yCoordinateFxn(
                    k.children[-1]
                )  ## get leftmost and rightmost children's y coordinates
                yl = (
                    circ_s + circ * yl / self.ySpan
                )  ## transform y into a fraction of total y
                yr = circ_s + circ * yr / self.ySpan
                ybar = linspace(
                    yl, yr, precision
                )  ## what used to be vertical node bar is now a curved line

                xs = [
                    yx * x for yx in map(math.sin, ybar)
                ]  ## convert to polar coordinates
                ys = [yy * x for yy in map(math.cos, ybar)]

                branches += tuple(
                    zip(zip(xs, ys), zip(xs[1:], ys[1:]))
                )  ## add curved segment

                linewidths += [
                    linewidths[-1] for q in zip(ys, ys[1:])
                ]  ## repeat linewidths
                colours += [colours[-1] for q in zip(ys, ys[1:])]  ## repeat colours

        line_segments = LineCollection(
            branches,
            lw=linewidths,
            ls="-",
            color=colours,
            capstyle="projecting",
            zorder=1,
        )  ## create line segments
        ax.add_collection(line_segments)  ## add collection to axes
        return ax

    def plot_circular_points(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        yCoordinateFxn=None,
        pointSize=None,
        pointSizeFxn=None,
        colour=None,
        colourFxn=None,
        circStart=0.0,
        circFrac=1.0,
        inwardSpace=0.0,
        normaliseHeight=None,
        outlineColour=None,
        outlineColourFxn=None,
        outlineSize=None,
        outlineSizeFxn=None,
        zorder=3,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.x
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        # Point size logic
        if pointSize is not None and pointSizeFxn is not None:
            raise ValueError(
                "Cannot specify both pointSize and pointSizeFxn. Please use only one."
            )
        if pointSize is None and pointSizeFxn is None:
            pointSize = lambda k: 40
        else:
            pointSize = pointSize if pointSize is not None else pointSizeFxn
        # Outline size logic
        if outlineSize is not None and outlineSizeFxn is not None:
            raise ValueError(
                "Cannot specify both outlineSize and outlineSizeFxn. Please use only one."
            )
        if outlineSize is None and outlineSizeFxn is None:
            outlineSize = lambda k: pointSize(k) * 2
        else:
            outlineSize = outlineSize if outlineSize is not None else outlineSizeFxn
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            )
        if colour is None and colourFxn is None:
            colour = lambda k: "k"
        else:
            colour = colour if colour is not None else colourFxn
        # Outline colour logic
        if outlineColour is not None and outlineColourFxn is not None:
            raise ValueError(
                "Cannot specify both outlineColour and outlineColourFxn. Please use only one."
            )
        if outlineColour is None and outlineColourFxn is None:
            outlineColour = lambda k: "k"
        else:
            outlineColour = outlineColour if outlineColour is not None else outlineColourFxn

        if inwardSpace < 0:
            inwardSpace -= self.treeHeight

        circ_s = circStart * math.pi * 2
        circ = circFrac * math.pi * 2

        allXs = list(map(xCoordinateFxn, self.Objects))
        if normaliseHeight is None:
            normaliseHeight = lambda value: (value - min(allXs)) / (
                max(allXs) - min(allXs)
            )
        linspace = lambda start, stop, n: (  # ToDo: this never gets used
            list(start + ((stop - start) / (n - 1)) * i for i in range(n))
            if n > 1
            else stop
        )

        xs = []
        ys = []
        colours = []
        sizes = []

        outline_xs = []
        outline_ys = []
        outline_colours = []
        outline_sizes = []
        for k in filter(targetFxn, self.Objects):
            x = normaliseHeight(
                xCoordinateFxn(k) + inwardSpace
            )  ## find normalised x position along circle's radius
            y = (
                circ_s + circ * yCoordinateFxn(k) / self.ySpan
            )  ## get y position along circle's perimeter
            X = math.sin(y) * x  ## transform
            Y = math.cos(y) * x  ## transform

            xs.append(X)
            ys.append(Y)
            colours.append(colour(k)) #  if callable(colour) else colours.append(colour)
            sizes.append(pointSize(k)) #  if callable(pointSize) else sizes.append(pointSize)

            if outlineColour:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outlineSizeFxn(k)) #  if callable(outlineSizeFxn else outline_colours.append(outlineSizeFxn)
                outline_sizes.append(outlineSize(k)) #  if callable(outlineSize) else outline_sizes.append(outlineSize)

        ax.scatter(
            xs,
            ys,
            s=sizes,
            facecolor=colours,
            edgecolor="none",
            zorder=zorder,
            **kwargs,
        )  ## put a circle at each tip
        if outlineColour:
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
