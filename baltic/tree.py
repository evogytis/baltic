import copy
import math
import logging
import random
import string
import warnings
from functools import reduce
import numpy as np
from matplotlib.collections import LineCollection
from baltic.node import Node
from baltic.leaf import Leaf
from baltic.clade import Clade
from baltic.reticulation import Reticulation
from baltic.bt_utils import project_to_polar, project_polar_vector, unnest

logger = logging.getLogger("baltic.Tree")

class Tree: ## tree class
    def __init__(self, treeType):
        ## initial current node is a new instance of a node class which needs to be initialized
        self.curNode = Node() ## current node is a new instance of a node class
        self.curNode.index = 'Root' ## first object in the tree is the root to which the rest gets attached
        self.curNode.length = 0.0 ## startind node branch length is 0
        self.curNode.height = 0.0 ## starting node height is 0
        
        self.root = None #self.curNode ## root of the tree is current node
        self.Objects = [] ## tree objects have a flat list of all branches in them
        self.tipMap = None
        self.treeHeight = 0 ## tree height is the distance between the root and the most recent tip
        self.mostRecent = None
        self.ySpan = 0.0
        assert treeType in [
            "divergence", 
            "time"
        ], f"Failed to initialize tree with treeType {treeType}."
        self.treeType = treeType

    def add_reticulation(self, name):
        
        logger.info(f"Creating new reticulation: {name}.")
        ret = Reticulation(name)
        ret.index = name
        ret.parent = self.curNode
        self.curNode.children.append(ret)
        self.Objects.append(ret)
        self.curNode=ret

    def add_node(self, i):
        """Create a new internal node with appropriate parent-child links; add it to the tree.
        """
        logger.info(f"Creating new node: {i}.")
        new_node = Node() ## new node instance
        new_node.index = i ## new node's index is the position along the tree string
        if self.root is None:
            self.root = new_node
            self.root.length = 0.0

        new_node.parent = self.curNode ## new node's parent is current node
        if not self.curNode.is_node():
            logger.error("Attempted to add a child to a non-node object.")
            logger.error(
                "Check if tip names have illegal characters like parentheses or commas."
            )
            raise TypeError()

        self.curNode.children.append(new_node) ## new node is a child of current node
        self.curNode = new_node ## current node is now new node
        self.Objects.append(self.curNode) ## add new node to list of objects in the tree

    def add_leaf(self, i, name):
        """Create a new leaf with appropriate parent-child links; add it to the tree.
        """
        logger.info(f"Creating new leaf: {name}.")
        new_leaf = Leaf(name) ## new instance of leaf object
        new_leaf.index = i ## index is position along tree string
        if self.root is None:
            self.root = new_leaf

        new_leaf.parent = self.curNode ## leaf's parent is current node
        if not self.curNode.is_node():
            logger.error("Attempted to add a child to a non-node object.")
            logger.error(
                "Check if tip names have illegal characters like parentheses or commas."
            )
            raise TypeError()
        self.curNode.children.append(new_leaf) ## assign leaf to parent's children
        self.curNode = new_leaf ## current node is now new leaf
        self.Objects.append(self.curNode) ## add leaf to all objects in the tree

    def subtree(self,startingNode=None,traverseCondition=None,stem=True):
        """Generate a new subtree starting from a given root node according to a certain condition.
        """
        logger.info("Generating subtree.")
        if startingNode is None:
            logger.debug("No startingNode given, using root as default.")
            startingNode = self.root
        if traverseCondition is None:
            logger.debug("No traverseCondition given, using all branches.")
            traverseCondition = lambda k: True

                # move up a node if we include the stem
        # if stem and startingNode != self.root:
        #     node = startingNode.parent
        # else:
        node = startingNode

        subtreeBranches = self.traverse_tree(
            node, includeCondition=lambda k: True, traverseCondition=traverseCondition
        )
        subtreeBranches = copy.deepcopy(subtreeBranches)

        if stem and startingNode != self.root:  ## using stem - need to prune subtrees from root now
            logger.debug("Using stem.")
            stem_branch = copy.deepcopy([branch.parent for branch in subtreeBranches if branch.parent.index == startingNode.parent.index][-1]) ## identify parent of starting node
            stem_branch.children = [ch for ch in stem_branch.children if ch.index == startingNode.index] ## reduce children to just the node of interest
            stem_branch.children = [subtreeBranches[0]] ## stem has a single child
            subtreeBranches[0].parent = stem_branch ## assign stem as parent
            subtreeBranches.insert(0, stem_branch) ## add to list of subtree branches

            # ## add all branches resulting from traversals of unwanted siblings
            # unwantedBranches = []
            # for child in node.children:
            #     if child.index != startingNode.index:
            #         unwantedBranches.extend(self.traverse_tree(child, includeCondition=lambda w: True))

            # ## iterate over subtree branches, remember those that belong to unwanted subtrees
            # remove = []
            # for k in subtreeBranches:
            #     if k.index in [ub.index for ub in unwantedBranches]:
            #         remove.append(k)
            
            # print(len(subtreeBranches),len(unwantedBranches),len(remove))
            # for r in remove:
            #     subtreeBranches.remove(r)

        ## nothing found or no leaf objects in traversal
        if subtreeBranches is None or not any(node.is_leaf() for node in subtreeBranches):
            logger.error("No branches found in subtree traversal. Exiting.")
            return None

        ### Initialize the new tree
        logger.debug("Creating new Tree object for subtree and assigning new branches.")
        localTree = Tree(self.treeType)
        localTree.Objects = subtreeBranches
        localTree.root = subtreeBranches[0]

        superRoot = Node()
        superRoot.index = 'Root'
        superRoot.length = 0.0
        superRoot.height = 0.0
        superRoot.children.append(localTree.root)
        localTree.root.parent = superRoot

        ## turn branches into set for quicker look up later
        subtreeSet = set(subtreeBranches)

        if traverseCondition is not None:
            logger.debug("Didn't use default traversal.")
            ## didn't use default traverse condition,
            ## might need to deal with hanging nodes and prune children
            for nd in localTree.get_internal():
                ## only keep children seen in traversal
                nd.children = [child for child in nd.children if child in subtreeSet]
            logger.debug("Removing (hanging) nodes without children.")
            localTree.fix_hanging_nodes()

        if self.tipMap:
            logger.debug("Copying tipMap from original tree to subtree.")
            localTipMap = {}
            ## copy over the relevant tip translations
            for tipNum in self.tipMap:
                if self.tipMap[tipNum] in [w.name for w in localTree.get_external()]:
                    localTipMap[tipNum] = self.tipMap[tipNum]
            localTree.tipMap = localTipMap
            # TODO: Make sure I broke this down correctly
            # localTree.tipMap={tipNum: self.tipMap[tipNum] for tipNum in self.tipMap if self.tipMap[tipNum] in [w.name for w in localTree.get_external()]}

        localTree.traverse_tree()

        return localTree

    def make_single_type(self):
        logger.info("Converting to single-type tree.")
        while True:
            multitypeNodes = self.get_internal(lambda k: len(k.children) == 1)

            logger.debug(f"Number of multitype nodes in tree: {len(multitypeNodes)}")
            if not multitypeNodes:
                break

            for k in sorted(multitypeNodes, key=lambda x: -x.height):

                child = k.children[0]
                if k.parent is not None:
                    grandparent = k.parent
                else:
                    grandparent = self.root
                logger.debug(f"At parent multitype node {k.index} with child {child.index} and grandparent {grandparent.index}. Grandparent children: {[ch.index for ch in grandparent.children]}.")
                child.parent = grandparent
                grandparent.children.append(child)
                grandparent.children.remove(k)
                grandparent.children = list(set(grandparent.children))

                child.length += k.length  ## adjust child length

                self.Objects.remove(k)  ## remove old parent from all objects
        self.sort_branches()


    def set_absolute_time(self, mostRecentSamplingDate, justLeaves = False):
        logger.debug("Setting absoluteTime for all branches.")
        logger.debug(f"MRSD: {mostRecentSamplingDate}")
        for k in self.Objects:  ## iterate over all objects
            if justLeaves and k.is_node():
                continue
            k.absoluteTime = (
                mostRecentSamplingDate - self.treeHeight + k.height
            )  ## heights are in units of time from the root
        self.mostRecent = max(k.absoluteTime for k in self.Objects if k.absoluteTime)

    def _assign_date_uncertainty(self, dateUncertainties):
        for k in self.get_external():
            date, uncertainty = dateUncertainties[k.name]
            k.absoluteTimeRange = uncertainty

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

        stats["strictlyBifurcating"] = bool(strictlyBifurcating)
        stats["multitypeTree"] = bool(multiType)
        stats["singletonTree"] = bool(singleton)
        stats["hasTraits"] = bool(hasTraits)

        stats["numObjects"] = len(obs)
        stats["numNodes"] = len(nodes)
        stats["numLeaves"] = len(self.get_external())

        return stats

    def _partition_tree(self, partitionFxn, node=None, partitionLabel=None):
        """
        Assign a label of 10 random characters to branches depending on partitionFxn.
        If partitionFxn evaluates to True a new label is generated for passing on to descendants, otherwise the label of current node is passed on.
        Use case example: label parts of a tree that are introduced into a new location and stay in that location, generating a new label for exports to other countries.
        """
        if node is None: node = self.root
        
        if partitionLabel is None:
            import random,string
            characters = string.ascii_letters + string.digits  # A-Z, a-z, 0-9
            partitionLabel = ''.join(random.choices(characters, k = 10))
        
        generateLabelFxn = lambda k: k == self.root or partitionFxn(k) ## root of the tree also counts as entering a new state
        
        if generateLabelFxn(node): ## if True - generate new label
            import random,string
            characters = string.ascii_letters + string.digits  # A-Z, a-z, 0-9
            partitionLabel = ''.join(random.choices(characters, k = 10))
        
        node.traits['partition'] = partitionLabel
        
        if node.is_node():
            for ch in node.children:
                self._partition_tree(partitionFxn = generateLabelFxn, node=ch, partitionLabel=partitionLabel)
        
        return self

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

    def midpoint_root(self, fixSingletons=True):
        logger.debug("No branch provided, rooting at midpoint")
        # Identify the largest pairwise distance
        maxDistance = 0.0

        for tip in self.get_external():  ## iterate over tips

            self.reroot(branch=tip, branchFrac=0.0, fixSingletons=fixSingletons)
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
        self = self.reroot(branch=tip1, branchFrac=0.0, fixSingletons=fixSingletons)

        # Depth to go from the ingroup tip toward the outgroup tip
        rootRemainder = 0.5 * maxDistance  # - (self.root.length or 0))
        assert rootRemainder >= 0
        # Identify the midpoint and reroot there.
        # Trace the path to the outgroup tip until all of the root depth has
        # been traveled/accounted for.
        path = tip2.get_path_to_root()

        for node in path[::-1]:  ## iterate from old root to new
            rootRemainder -= node.length
            #             logger.info(f'iterating over path: {node.index} root remainder: {root_remainder}')
            if rootRemainder < 0:
                outgroup_node = node
                branchFrac = -rootRemainder / outgroup_node.length
                break
        logger.debug(f"Midpoint rooting: rooting on node {outgroup_node.index} halfway from previous highest tip and current highest tip {tip2.name}.")
        self = self.reroot(
            branch=outgroup_node, branchFrac=branchFrac, fixSingletons=fixSingletons
        )
        logger.debug("Finished midpoint rooting.")
        return self

    def reroot(self, branch=None, branchFrac=0.5, fixSingletons=True):
        if self.treeType == "time":
            logger.error("Cannot reroot a time-calibrated tree.")
            raise AttributeError

        if branch == self.root:
            logger.warning("Rerooting attempted on existing root.")
            return self

        if self.root.length > 0:
            logger.warning("Current root has non-zero branch length likely set up for visual purposes. This branch length will be suppressed for rerooting.")
            self.root.length = 0.0

        oldTreeLength = sum(self.get_parameter_list("length"))  ## get old tree length
        ###############
        if branch is None:  ## midpoint rooting
            self = self.midpoint_root(fixSingletons=fixSingletons)

            return self

        ##############
        path = branch.get_path_to_root()  ## get path from new root to old root, ignore actual root node
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

    def root_by_regression(
        self,
        stat: str = "r^2",
        forcePositive: bool = True,
        nJobs: int | None = None,
        baseRefineIters: int = 20,
        refineItersPerTip: int = 10,
        maxRefineIters: int = 400,
        returnRefinedDates: bool = True
    ):
        """
        Root the tree by evaluating root-to-tip regression for every possible root.
        Uncertain tip dates are refined deterministically (not Monte Carlo).

        Parameters:
            stat              - "r^2", "correlation", or "sum of squares"
            forcePositive     - enforce positive evolutionary rate
            nJobs             - parallel worker count
            baseRefineIters   - baseline refinement iterations
            refineItersPerTip - extra iterations in proportion to uncertain tips
            maxRefineIters    - hard cap on refinement iterations
            returnRefinedDates - return the final refined tip dates

        Returns:
            rerooted tree and (optionally) refined dates for uncertain tips.
        """
        from baltic.bt_utils import _rtt_worker, _root_to_tip, decimal_to_calendar_date
        from concurrent.futures import ProcessPoolExecutor

        validRootingMethods = ["r^2", "correlation", "sum of squares"]
        assert stat in validRootingMethods, (
            f"Invalid root-to-tip metric: {stat} (choose from {validRootingMethods})"
        )

        # Logging setup
        func_logger = logging.getLogger("baltic.tree.root_by_regression")
        func_logger.setLevel(logging.INFO)
        func_logger.propagate = False
        if not func_logger.handlers:
            import sys
            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(logging.Formatter("[root_by_regression] %(message)s"))
            handler.setLevel(logging.INFO)
            func_logger.addHandler(handler)

        func_logger.info(f"Starting root-to-tip regression, metric = {stat}")

        # Candidate roots
        candidates = self.Objects[1:]       # all non-root nodes
        if not candidates:
            func_logger.warning("Tree has no internal nodes; returning unchanged.")
            return self

        # Tip dates
        tips = self.get_external()
        tipDates = {k.name: k.absoluteTime for k in tips}

        assert [k.absoluteTimeRange for k in tips].count(None) != len(tips), f"Tips don't have absoluteTimeRange attribute. Check whether tree was imported with absoluteTime set to True"

        # Identify uncertain date tips (date ranges > tiny threshold)
        uncertainTips = [k for k in tips
                         if np.diff(k.absoluteTimeRange)[0] > 1e-12]
        uncertainRanges = {k.name: k.absoluteTimeRange for k in uncertainTips}

        nUncertain = len(uncertainTips)

        # Decide refinement iterations
        if nUncertain == 0:
            refineIters = 1
            func_logger.info(
                "All tip dates precise; performing a single regression per root candidate."
            )
        else:
            refineIters = (
                baseRefineIters +
                refineItersPerTip * nUncertain
            )
            refineIters = min(refineIters, maxRefineIters)

            func_logger.info(
                f"{nUncertain} uncertain tip dates → "
                f"{refineIters} refinement iterations per candidate root."
            )

        # Build worker tasks
        tasks = []
        for k in candidates:
            tasks.append(
                (
                    self,               # picklable tree
                    k.index,            # candidate root
                    copy.deepcopy(tipDates),
                    copy.deepcopy(uncertainRanges),
                    refineIters,        # max refinement steps
                    stat,
                    forcePositive
                )
            )

        # Determine number of jobs
        if nJobs is None or nJobs <= 0:
            from os import cpu_count
            nJobs = cpu_count() or 1

        func_logger.info(
            f"Evaluating {len(candidates)} candidate roots using {nJobs} workers..."
        )

        # Dispatch workers
        if len(tasks) == 1:
            results = [_rtt_worker(tasks[0])]
        else:
            with ProcessPoolExecutor(max_workers=nJobs) as pool:
                results = list(pool.map(_rtt_worker, tasks))

        # Determine the best root by score
        best_res = max(results, key=lambda r: r["score"] if r["score"] is not None else -np.inf)
        best_root_index = best_res["root_index"]
        refined_dates = best_res.get("assigned_uncertain_dates", {})  # still using old field name

        for k in uncertainTips:
            k.absoluteTime = refined_dates[k.name] ## assign new dates to uncertain tips

        func_logger.info(f"Best root index: {best_root_index}, score = {best_res['score']}")

        # Re-root the tree on the best node
        best_node = next(obj for obj in self.Objects if obj.index == best_root_index)
        self = self.reroot(best_node)

        if len(self.root.children) == 2:
            from baltic.bt_utils import _adjust_tip_dates_by_regression

            func_logger.info("Refining root position along the bifurcating root branch...")
            left, right = self.root.children

            totalBranch = left.length + right.length
            leftNodes = self.traverse_tree(left)
            rightNodes = self.traverse_tree(right)

            tips_current = self.get_external()

            tip_is_left = {k.name: False for k in tips_current}
            tip_is_right = {k.name: False for k in tips_current}

            for k in leftNodes:
                tip_is_left[k.name] = True
            for k in rightNodes:
                tip_is_right[k.name] = True

            xs = []
            for k in tips_current:
                if k.name in refined_dates:
                    xs.append(refined_dates[k.name])
                else:
                    xs.append(k.absoluteTime)

            res_frac = {"frac": 0.0}
            for f in np.linspace(0, 1, 21):
                adjustL = -left.length + totalBranch * f
                adjustR = -right.length + totalBranch * (1 - f)

                ys_f = []
                for k in tips_current:
                    if tip_is_left[k.name]:
                        ys_f.append(k.height + adjustL)
                    elif tip_is_right[k.name]:
                        ys_f.append(k.height + adjustR)
                    else:
                        ys_f.append(k.height)

                res_frac = _root_to_tip(
                    left,
                    xs,
                    ys_f,
                    res_frac,
                    stat=stat,
                    forcePositive=forcePositive,
                    frac=f
                )

            self = self.reroot(branch=res_frac["root"], branchFrac=res_frac["frac"])

            best_res = res_frac

            uncertainTips = [k for k in tips
                         if np.diff(k.absoluteTimeRange)[0] > 1e-12]
            uncertainRanges = {k.name: k.absoluteTimeRange for k in uncertainTips}
            _adjust_tip_dates_by_regression(uncertainTips, best_res['slope'], best_res['intercept'])
        
        slope = best_res.get("slope", None)
        intercept = best_res.get("intercept", None)
        r2 = best_res.get("r^2", None)
        corr = best_res.get("correlation", None)
        sse = best_res.get("sum of squares", None)

        if corr is not None:
            func_logger.info(f"Correlation: {corr:.4g}")
        if sse is not None:
            func_logger.info(f"Sum of squares: {sse:.4g}")
        if r2 is not None:
            func_logger.info(f"r^2: {r2:.4g}")
        if slope is not None:
            func_logger.info(f"Slope (rate): {slope:.4e}")
        if slope and intercept:
            tmrca = -intercept / slope
            func_logger.info(
                f"Intercept: {intercept:.4g} TMRCA = {tmrca:.3f} "
                f"({decimal_to_calendar_date(tmrca, fmt='%Y-%b-%d')})"
            )
        
        if returnRefinedDates:
            return self, refined_dates
        else:
            return self


    def sort_branches(self, descending=True, sortFxn=None, operationFxn=None):
        mod = 1 if descending else -1
        if sortFxn is None and operationFxn is None:
            logger.debug("Using default sort function.")
            # TODO This is really hard to parse. It is a variable mapped to a lambda function
            # that sometimes returns a 2-tuple and sometimes returns a 3-tuple, depending on
            # whether the argument to the function is a node or not
            # This should be rewritten as a normal function, and we should make sure
            # the different lengths of the tuple get handled sensibly downstream

            def sortFxn(node):

                if node.is_node():
                    return (mod, len(node.leaves) * mod, node.length)

                elif node.is_leaflike():
                    if node.is_leaf() or isinstance(node, Reticulation):
                        return (mod * -1, 1 * mod, node.length)
                    elif isinstance(node, Clade):
                        return (mod * -1, len(node.leaves) * mod, node.length)
        elif sortFxn is not None and operationFxn is not None:
            raise Exception(
                    "Cannot specify both sortFxn and operationFxn, pick one."
                )

        if sortFxn:
            for k in self.get_internal():
                k.children = sorted(k.children, key = sortFxn)
        else:
            for k in self.get_internal():
                k.children = operationFxn(k.children)

        self._assign_tree_coordinates()  ## update x and y positions of each branch, since y positions will have changed because of sorting

    def _assign_tree_coordinates(self, order=None, padNodes=None):
        """ formerly drawTree()
        """
        if order is None:
            ## order is a list of tips recovered from a tree traversal
            # to make sure they're plotted in the correct order along the primary dimension
            order = self.traverse_tree(includeCondition=lambda k: k.is_leaflike())
            logger.debug("Drawing tree in pre-order")
        else:
            logger.debug("Drawing tree with provided order")

        nameOrder = {x.name: i for i, x in enumerate(order)}
        assert len(nameOrder) == len(order), "Non-unique names present in tree"

        logger.debug("Drawing tree with default widths (1 unit for leaf objects, width+1 for clades)")
        skips = [1 if isinstance(x, Leaf) else x.width + 1 for x in order]


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
            if isinstance(k,Clade): ## need to reset clade subtree attributes too since they can be accessed for "skewed" clade style
                for collapsedBranch in k.subtree:
                    collapsedBranch.x = collapsedBranch.height
                    collapsedBranch.y = y

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
            for k in self.get_external(onlyLeaves=False):  ## reset y positions so tree starts at y=0.5
                k.y -= minY - 0.5

        assert len([k for k in self.Objects if k.is_leaflike()]) == len(
            order
        ), "Number of tips in tree does not match number of unique tips, check if two or more collapsed clades were assigned the same name."
        storePlotted = 0

        while len(drawn) != len(
            self.Objects
        ):  # keep drawing the tree until everything is drawn
            logger.debug(f"Drawing iteration {len(drawn)}")
            for k in filter(
                lambda w: w.index not in drawn, self.get_internal()
            ):  ## iterate through internal nodes that have not been drawn
                if len([q.y for q in k.children if q.y is not None]) == len(
                    k.children
                ):  ## all y coordinates of children known
                    logger.debug(f"Setting node {k.index} coordinates to XXX")  # TODO: XXX=k.x? {k.x}?
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
                    f"Got stuck trying to find y positions of objects ({len(drawn)} branches drawn this iteration, {storePlotted} branches during previous iteration out of {len(self.Objects)} total)"
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

    def _assign_unrooted_tree_coordinates(self, circStart=0.0, node=None, total=None, padNodes=None):
        
        if padNodes is None: padNodes = {}

        if node is None:
            # total = sum(
            #     [1 if x.is_leaf() else x.width for x in self.get_external()]
            # )
            self._assign_tree_coordinates(padNodes=padNodes)
            total = self.ySpan
            node = self.root  # .children[0]
            for k in self.Objects:
                k._tau = 2 * math.pi * circStart
                k.x = 0.0
                k.y = 0.0

        if node.is_leaf(): ## Leaf - default width is 1
            w = 2*math.pi * 1/total

        elif node.is_node(): ## Node - default width is number of tips + width of descendant clades
            descendants=self.traverse_tree(node,includeCondition=lambda w: True)
            w = 2*math.pi * sum([1 if d.is_leaf() else d.width for d in descendants if d.is_leaflike()])/total

            for desc in descendants:
                if desc in padNodes:
                    w += 2*math.pi * (padNodes[desc]*2)/total

        elif isinstance(node,Clade): ## Clade - width is whatever was assigned during collapse
            w = 2*math.pi * (node.width + 1)/total
        

        if node in padNodes: ## at a branch that's designated for additional padding
            w += 2*math.pi * (padNodes[node]*2)/total
        

        if node.parent.x is None:
            node.parent.x = 0.0
            node.parent.y = 0.0

        node.x, node.y = project_polar_vector(node.parent.x, node.parent.y, node._tau+w*0.5, node.length)
        eta = node._tau

        if node.is_node():
            for ch in node.children:
                if ch.is_leaf():
                    w = 2*math.pi * 1/total

                elif ch.is_node():
                    descendants=self.traverse_tree(ch,includeCondition=lambda w: True)
                    w = 2*math.pi * sum([1 if d.is_leaf() else d.width + 1 for d in descendants if d.is_leaflike()])/total

                    for desc in descendants:
                        if desc in padNodes:
                            w += 2*math.pi * (padNodes[desc]*2)/total
                    
                elif isinstance(ch,Clade):
                    w = 2*math.pi * (ch.width + 1)/total
                

                if ch in padNodes:
                    w += 2*math.pi * (padNodes[ch]*2)/total

                ch._tau = eta
                eta += w
                self._assign_unrooted_tree_coordinates(circStart, ch, total, padNodes)

    # def find_MRCA(self, descendants):
    #     assert (
    #         len(descendants) > 1
    #     ), f"Not enough descendants to find common ancestor: {len(descendants)}"
    #     pathsToRoot = {
    #         k.index: set(k.get_path_to_root()) for k in descendants
    #     }  ## for every descendant create an empty set
    #     # for k in descendants: ## iterate through every descendant
    #     #     curNode=k ## start descent from descendant
    #     #     while curNode: ## while not at root
    #     #         pathsToRoot[k.index].add(curNode) ## remember every node visited along the way
    #     #         curNode=curNode.parent ## descend

    #     # return the most recent branch that is shared across all paths to root
    #     return sorted(reduce(set.intersection, pathsToRoot.values()), key=lambda k: (-len(k.leaves), k.height))[-1]
    def find_MRCA(self, *descendants):
        """
        This is a ChatGPT suggested function that's meant to be much faster.
        """
        if len(descendants) == 1 and isinstance(descendants[0], list):
            descendants = descendants[0]

        strInput = all([isinstance(k, str) for k in descendants]) ## check if a number of strings were provided corresponding to tip names

        if strInput:
            tipBranches = self.get_external(lambda k: k.name in descendants)
            assert len(tipBranches) == len(descendants), f"Could not find tips amongst branches: {', '.join([k for k in descendants if k not in [k.name for k in self.Objects if k.is_leaf()]])}"
            descendants = tipBranches
        paths = [k.get_path_to_root()[::-1] for k in descendants]
        mrca = None

        # zip stops when any list ends
        for nodes in zip(*paths):
            # if all same object, continue
            if all(n is nodes[0] for n in nodes):
                mrca = nodes[0]
            else:
                break

        return mrca


    def collapse_subtree_to_clade(self,
                                    cl,
                                    givenName,
                                    widthFunction=lambda k: len(k.leaves)):
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
                    and bool(collapseIfFxn(n))
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
                    nodesToDelete = list(
                        filter(
                            lambda n: n.is_node()
                            and bool(collapseIfFxn(n))
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
                        comment.append(f"{tr}={{{','.join(rangeComment)}}}")
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
            comment = f"[&{comment}]"
            stringFragment.append(comment)  ## end of node, add annotations

        logger.debug(f"adding branch length to {curNode.index}")
        stringFragment.append(f":{curNode.length:.15f}")  ## end of node, add branch length

        if curNode == self.root:  # .children[-1]:
            stringFragment.append(";")
            if bool(nexus):
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
            len([k for k in tipsToKeep if not k.is_leaflike()])
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

    def get_leaf(self, tipName):
        assert isinstance(tipName, str), f"tipName should be str, not {type(tipName)}"

        matchingTip = [k for k in self.get_external() if k.name == tipName]
        assert len(matchingTip) == 1, f"Found {len(matchingTip)} matches"

        return matchingTip[0]

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

    def explode_tree(self, trait=None, customFxn=None, stem=True):
        
        if trait is not None and customFxn is not None:
            raise ValueError(
                "Cannot specify both a trait and a custom function to split the tree by. Please use only one."
            )
        elif trait is None and customFxn is None:
            raise ValueError(
                "Must specify either a trait or a custom function to split the tree by."
            )
        if trait is not None:
            traceFunction = lambda k: k == self.root or k.traits[trait] != k.parent.traits[trait] ## current branch trait not the same as parent or at root
        elif customFxn:
            traceFunction = lambda k: k == self.root or customFxn(k)
        
        subtrees=[]
        
        for k in self.Objects: ## iterate over every branch
            if traceFunction(k): ## check if branch satisfies splitting function
                subtree = self.subtree(startingNode = k, traverseCondition = lambda ch: not traceFunction(ch), stem=stem) ## extract subtree
                if subtree:
                    subtrees.append(subtree)
                else:
                    logger.error(f"Failed to extract valid subtree at branch {k.index}.")
                
        return subtrees

    def condense_tree(self, cutoffs = None, protectedTips = [], widthFxn = None):
        ## cutoffs is a tuple of values that define the minimum and maximum sizes of clades designated for collapsing
        ## protectedTips is a list of tips names as str that will remain uncollapsed
        ## widthFxn will determine how much y-axis space collapsed clades will occupy (default is same amount the uncollapsed subtree would)
        if cutoffs is None:
            N = len(self.get_external())
            minClade, maxClade = 3, int(N * 0.2) ## default is collapse any node with between 3 and up to 20% of the tree's leaves as descendants
        else:
            minClade, maxClade = cutoffs

        if widthFxn is None:
            widthFxn = lambda k: len(k.leaves)

        collapseCandidates = self.get_internal(lambda k: len(k.leaves.intersection(set(protectedTips))) == 0 and (minClade <= len(k.leaves) <= maxClade))
        nodesToCollapse = unnest(collapseCandidates, towardsRoot = True)

        for i,node in enumerate(nodesToCollapse):
            self.collapse_subtree_to_clade(node,f'collapsed clade {i}', widthFunction = widthFxn)

        return self
    ############################################
    ##########   PLOTTING FUNCTIONS   ##########
    ############################################

    def _plot_unrooted_clades(self,ax,xCoordinateFxn,yCoordinateFxn,endAttrFxn,targetFxn,colour,precision,style,cladeBaseWidth,**kwargs):
        
        valid_styles=['equal','skewed']
        assert style in valid_styles, f"Style {style} not recognised. Options are {valid_styles}"
        
        from matplotlib.patches import Polygon
        
        # if target==None: target=lambda k:True
        mainTargetFxn = lambda k: isinstance(k,Clade) and targetFxn(k)
        
        # if xCoordinateFxn==None: xCoordinateFxn=lambda k:k.x
        # if yCoordinateFxn==None: yCoordinateFxn=lambda k:k.y
        # if end_attr==None: end_attr=lambda k:k.lastHeight
        # if pad_nodes is None: pad_nodes = {}
        # if colour==None: colour=lambda k:(0.7,0.7,0.7) ## grey if nothing given
        
        localKwargs=dict(kwargs)
        #########################
        # self._drawTree(pad_node=pad_nodes)
        total_y=self.ySpan
        
        w=self.treeHeight*cladeBaseWidth ## we turn the clade triangle into a trapezoid to avoid a very acute angle where the clade joins the tree, so we turn it to a trapezoid whose edge adjacent to the branch is a tiny line
        
        for k in filter(mainTargetFxn,self.Objects):
            x0,y0=xCoordinateFxn(k),yCoordinateFxn(k) ## get child coordinate
            ###############
            tau=k._tau
            p_tau=k.parent._tau
            clade_length=k.lastHeight-k.height
            clade_width=k.width/total_y * math.pi*2
            narrow=0.5
            
            ################ trapezoid base (prevents clades from forming sharp triangles where they join the tree)
            x1a,y1a = project_polar_vector(x0, y0, p_tau, w)
            x1b,y1b = project_polar_vector(x0, y0, p_tau-math.radians(180), w)

            ###################
            ybar=np.linspace(tau,tau+clade_width,precision) ## arc starts at tau and finishes with tau+width
            
            if style=='equal':
                arc = [project_polar_vector(x0, y0, yb, clade_length) for yb in ybar]
            
            if style=='skewed':
                early_end = min([xCoordinateFxn(q) for q in k.subtree if q.is_leaf()]) ## earliest leaf height
                short_clade_length = early_end-k.height ## distance between earliest tip and beginning of clade
                
                arc = [project_polar_vector(x0, y0, yb, xdist) for xdist,yb in zip(np.linspace(short_clade_length,clade_length,precision),ybar)]
                
            coords=[(x1b,y1b),(x1a,y1a)] + arc ## trapezoid base + arc
            ############
            fc=colour(k) if callable(colour) else colour
            ec='k'

            if 'zorder' not in localKwargs: localKwargs['zorder']=100
            if 'capstyle' not in localKwargs: localKwargs['capstyle']='round'
            if 'edgecolor' not in localKwargs: localKwargs['edgecolor']=ec
            if 'facecolor' not in localKwargs: localKwargs['facecolor']=fc
            poly=Polygon(coords,**localKwargs)
            ax.add_patch(poly)
        
        return ax

    def plot_text(
        self,
        ax,
        targetFxn=None,
        xCoordinateFxn=None,
        xSpace=0.005,
        yCoordinateFxn=None,
        ySpace=0.0,
        textContentFxn=None,
        colour=None,
        colourFxn=None,
        treeType='rectangular',
        orientation='horizontal',
        padNodes=None,
        circStart=0.0,
        circFrac=1.0,
        inwardSpace=0.0,
        normaliseHeight=None,
        cladeEndAttrFxn=None,
        recomputeCoordinates=True,
        **kwargs,
    ):
        valid_treeTypes = ['rectangular','circular','unrooted']
        assert treeType in valid_treeTypes, f"Tree type {treeType} not recognised. Options are {valid_treeTypes}"

        valid_orientations = ['vertical','horizontal']
        assert orientation in valid_orientations, f"Tree orientation {orientation} not recognised. Options are {valid_orientations}"
        
        ## set defaults
        if colour is None and colourFxn is None:
            colourFxn = lambda k: "k"
        elif colourFxn is None:
            colourFxn = lambda k: colour
        if targetFxn is None: 
            targetFxn = lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.height + xSpace*self.treeHeight if self.treeType=='divergence' else k.absoluteTime + xSpace*self.treeHeight
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y + ySpace*self.treeHeight
        if textContentFxn is None:
            textContentFxn = lambda k: k.name
        if padNodes is None:
            padNodes = {}
        if cladeEndAttrFxn is None:
            cladeEndAttrFxn = lambda k: max([xCoordinateFxn(desc) for desc in k.subtree])

        if recomputeCoordinates:
            self._assign_tree_coordinates(padNodes=padNodes)

        if treeType == 'rectangular':
            ax = self._plot_rectangular_text(ax=ax,targetFxn=targetFxn,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,textContentFxn=textContentFxn,orientation=orientation,cladeEndAttrFxn=cladeEndAttrFxn,colourFxn=colourFxn,**kwargs)

        elif treeType == 'circular':
            if inwardSpace<0: inwardSpace-=self.treeHeight

            allXs = list(map(xCoordinateFxn, self.Objects))
            if normaliseHeight is None:
                normaliseHeight = lambda value: (value - min(allXs)) / (
                    max(allXs) - min(allXs)
                )

            ax = self._plot_circular_text(ax=ax,targetFxn=targetFxn,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,textContentFxn=textContentFxn,circStart=circStart,circFrac=circFrac,inwardSpace=inwardSpace,normaliseHeight=normaliseHeight,cladeEndAttrFxn=cladeEndAttrFxn,colourFxn=colourFxn,**kwargs)
        
        elif treeType == 'unrooted':
            self._assign_unrooted_tree_coordinates(circStart=circStart,padNodes=padNodes) ## compute unrooted coordinates
            xCoordinateFxn = lambda k: k.x ## use projected x coordinate now

            ax = self._plot_unrooted_text(ax=ax,targetFxn=targetFxn,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,textContentFxn=textContentFxn,cladeEndAttrFxn=cladeEndAttrFxn,colourFxn=colourFxn,**kwargs)

        return ax

    def _plot_rectangular_text(
        self,
        ax,
        targetFxn,
        xCoordinateFxn,
        yCoordinateFxn,
        textContentFxn,
        orientation,
        cladeEndAttrFxn,
        colourFxn,
        **kwargs,
    ):
        localKwargs = dict(kwargs)

        if orientation == 'vertical':
            if 'rotation' not in localKwargs:
                localKwargs['rotation']=90

        if 'verticalalignment' not in localKwargs:
                localKwargs['verticalalignment']='center'
        if 'horizontalalignment' not in localKwargs:
            localKwargs['horizontalalignment']='left'
        if 'rotation_mode' not in localKwargs:
            localKwargs['rotation_mode']='anchor'
        if 'zorder' not in localKwargs:
            localKwargs['zorder'] = 4


        for k in filter(targetFxn, self.Objects):
            x, y = xCoordinateFxn(k), yCoordinateFxn(k)
            if orientation == 'vertical':
                x, y = y, x
            colour = colourFxn(k)

            ax.text(x, y, textContentFxn(k), color=colour, **localKwargs)
        return ax

    def _plot_unrooted_text(
        self,
        ax,
        targetFxn,
        xCoordinateFxn,
        yCoordinateFxn,
        textContentFxn,
        cladeEndAttrFxn,
        colourFxn,
        **kwargs,
    ):
        localKwargs = dict(kwargs)

        if 'verticalalignment' not in localKwargs:
            localKwargs["verticalalignment"] = "center"
        if 'rotation_mode' not in localKwargs:
            localKwargs['rotation_mode']='anchor'
        if "zorder" not in localKwargs:
            localKwargs["zorder"] = 4

        for k in filter(targetFxn, self.Objects):    

            x, y = xCoordinateFxn(k), yCoordinateFxn(k)
            xp, yp = xCoordinateFxn(k.parent), yCoordinateFxn(k.parent)

            # assert (
            #     hasattr(k,"_tau")
            # ), "Branch does not have angle _tau computed by _assign_unrooted_tree_coordinates()."


            textRotation = math.degrees(math.atan2(yp - y, xp - x)) % 360
            ha = 'left' if 90 <= textRotation % 360 <= 270 else 'right' ## rotate labels to aid readability
            textRotation = (textRotation - 180) % 360 if 90 <= textRotation % 360 <= 270 else textRotation
            colour = colourFxn(k)

            ax.text(
                x,
                y,
                textContentFxn(k),
                rotation=textRotation,
                ha=ha,
                color=colour,
                **localKwargs,
            )

        return ax

    def _plot_circular_text(
        self,
        ax,
        targetFxn,
        textContentFxn,
        xCoordinateFxn,
        yCoordinateFxn,
        circStart,
        circFrac,
        inwardSpace,
        normaliseHeight,
        cladeEndAttrFxn,
        colourFxn,
        **kwargs,
    ):
        assert circFrac > 0.0, f"Circular tree layout not given any space (circFrac == {circFrac})"

        localKwargs = dict(kwargs)  ## copy global kwargs into a local version

        # if cladeEndAttrFxn is None: cladeEndAttrFxn = lambda k: max([xCoordinateFxn(desc) for desc in k.subtree])

        # self._assign_tree_coordinates(padNodes=padNodes) ## computes y coordinates with padding
        total_y = self.ySpan ## plot_text already made a call to compute coordinates

        if 'verticalalignment' not in localKwargs and 'va' not in localKwargs:
            localKwargs['verticalalignment'] = "center"
        if 'rotation_mode' not in localKwargs:
            localKwargs['rotation_mode'] = 'anchor'
        if 'zorder' not in localKwargs:
            localKwargs['zorder'] = 4

        for k in filter(targetFxn, self.Objects):  ## iterate over branches
            x = normaliseHeight(
                xCoordinateFxn(k) + inwardSpace
            )  ## get branch x position
            y = yCoordinateFxn(k)  ## get y position

            rotationRadians = (circStart + (circFrac * y/total_y)) * 2*math.pi
            textRotation = (90 - math.degrees(rotationRadians))%360
            if 'horizontalalignment' not in localKwargs and 'ha' not in localKwargs:
                ha = 'right' if 90 <= textRotation%360 <= 270 else 'left' ## rotate labels to aid readability
            else:
                ha = localKwargs['horizontalalignment'] if 'horizontalalignment' in localKwargs else localKwargs['ha']

            textRotation = (textRotation - 180)%360 if 90 <= textRotation%360 <= 270 else textRotation

            colour = colourFxn(k)

            x,y = project_to_polar(x=x,y=y,yRange=total_y,circleStart=circStart,circleFraction=circFrac)

            ax.text(
                x,
                y,
                textContentFxn(k),
                rotation=textRotation,
                ha=ha,
                color=colour,
                **localKwargs,
            )

        return ax

    def plot_aligned_tip_labels(self, ax, xSpace=0.005, connectingLines=True, **kwargs):
        ## need to implement orientation
        localKwargs=dict(kwargs)
        if 'treeType' in localKwargs and localKwargs['treeType'] == 'unrooted': warnings.warn("Unrooted tree layout cannot accommodate aligned text labels, ignoring and plotting tip labels at the tips.")
        if 'xCoordinateFxn' in localKwargs:
            warnings.warn("xCoordinateFxn provided but needs to be reassigned for aligned text, ignoring.")
            del localKwargs['xCoordinateFxn']
        if 'normaliseHeight' in localKwargs:
            warnings.warn("normaliseHeight provided but needs to be reassigned for aligned text, ignoring.")
            del localKwargs['normaliseHeight']

        ## xSpace is fraction of tree height after which tip labels are put (0 == all tip labels are placed at exactly the x coordinate of the last tip, 1 == all tip labels are placed a full tree height away from the last tip)
        xCoordinateFxn = lambda k: self.treeHeight * (1 + xSpace) if self.treeType == 'divergence' else self.root.absoluteTime + self.treeHeight * (1 + xSpace)
        
        ## this is for handling circular layouts (prevents automatic normaliseHeight from doing division by 0 because xCoordinateFxn returns a constant value)
        allXs = self.get_parameter_list('absoluteTime') if self.treeType == 'time' else self.get_parameter_list('height')
        normaliseHeight = lambda value: (value - min(allXs)) / (max(allXs) - min(allXs))
        
        self.plot_text(ax, xCoordinateFxn=xCoordinateFxn, normaliseHeight=normaliseHeight, **localKwargs)

        ## add connecting lines
        if connectingLines:
            lines = []
            colours = []
            
            if 'colour' in localKwargs:
                colourFxn = lambda k: localKwargs['colour']
            elif 'colourFxn' in localKwargs:
                colourFxn = localKwargs['colourFxn']
            else:
                colourFxn = lambda k: 'k'
                
            if 'treeType' in localKwargs and localKwargs['treeType'] == 'circular': ## set defaults for circular layout if they haven't been set before
                circStart = localKwargs['circStart'] if 'circStart' in localKwargs else 0.0
                circFrac = localKwargs['circFrac'] if 'circFrac' in localKwargs else 1.0
                inwardSpace = localKwargs['inwardSpace'] if 'inwardSpace' in localKwargs else 0.0
                if inwardSpace < 0: inwardSpace-=self.treeHeight
            elif 'treeType' in localKwargs and localKwargs['treeType'] == 'unrooted': ## ignore when dealing with unrooted trees, but warn
                warnings.warn("Unrooted tree layout cannot accommodate aligned text labels, ignoring and not plotting connecting lines.")

            ## iterate over tips, add lines
            for k in self.get_external():
                x = k.height if self.treeType == 'divergence' else k.absoluteTime
                
                if 'treeType' in localKwargs and localKwargs['treeType'] == 'circular': ## circular layout
                    x1, y1 = project_to_polar(normaliseHeight(x+inwardSpace), k.y, self.ySpan, circleStart=circStart, circleFraction=circFrac)
                    x2, y2 = project_to_polar(normaliseHeight(xCoordinateFxn(k)+inwardSpace), k.y, self.ySpan, circleStart=circStart, circleFraction=circFrac)
                elif 'treeType' in localKwargs and localKwargs['treeType'] == 'unrooted': ## ignore when dealing with unrooted trees
                    pass
                else: ## rectangular layout
                    x1, y1 = x, k.y
                    x2, y2 = xCoordinateFxn(k), k.y
                
                lines.append([(x1, y1), 
                              (x2, y2)])
                colours.append(colourFxn(k))

            lines = LineCollection(lines, lw=1, ls=':', colors=colours, clip_on=False)
            ax.add_collection(lines)
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
        outline=True,
        outlineSize=None,
        outlineSizeFxn=None,
        outlineColour=None,
        outlineColourFxn=None,
        padNodes=None,
        treeType='rectangular',
        orientation='horizontal',
        circStart=None,
        circFrac=None,
        inwardSpace=None,
        normaliseHeight=None,
        recomputeCoordinates=True,
        **kwargs,
    ):
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: k.is_leaf()
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.height if self.treeType == 'divergence' else k.absoluteTime
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y
        # Point size logic
        if pointSize is not None and pointSizeFxn is not None:
            raise ValueError(
                "Cannot specify both pointSize and pointSizeFxn. Please use only one."
            )
        if pointSize is None and pointSizeFxn is None:
            pointSizeFxn = lambda k: 40
        elif pointSizeFxn is None:
            pointSizeFxn = lambda k: pointSize
        # Outline size logic
        if outlineSize is not None and outlineSizeFxn is not None:
            raise ValueError(
                "Cannot specify both outlineSize and outlineSizeFxn. Please use only one."
            )
        if outlineSize is None and outlineSizeFxn is None:
            outlineSizeFxn = lambda k: ((2+math.sqrt(pointSizeFxn(k)/math.pi))**2)*math.pi
        elif outlineSizeFxn is None:
            outlineSizeFxn = lambda k: outlineSize
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            ) ## should be a warning, since this eventuality is handled in the next line
        if colour is None and colourFxn is None:
            colourFxn = lambda k: "k"
        elif colourFxn is None:
            colourFxn = lambda k: colour
        # Outline colour logic
        if outlineColour is not None and outlineColourFxn is not None:
            raise ValueError(
                "Cannot specify both outlineColour and outlineColourFxn. Please use only one."
            )
        if outlineColour is None and outlineColourFxn is None:
            outlineColourFxn = lambda k: "k"
        elif outlineColourFxn is None:
            outlineColourFxn = lambda k: outlineColour

        if padNodes is None:
            padNodes = {}

        assert treeType in [
            'rectangular', 
            'circular', 
            'unrooted'
        ], f"Unrecognised tree layout {treeType}"
        assert orientation in [
            'horizontal', 
            'vertical'
        ], f"Unrecognised orientation {orientation}"

        localKwargs=dict(kwargs)

        xs=[]
        ys=[]
        colours=[]
        sizes=[]

        outline_xs=[]
        outline_ys=[]
        outline_colours=[]
        outline_sizes=[]

        if recomputeCoordinates:  ## compute branch coordinates
            if treeType == 'unrooted':
                self._assign_unrooted_tree_coordinates(circStart=circStart, padNodes=padNodes)
            else:
                self._assign_tree_coordinates(padNodes=padNodes)
        
        ########## warnings for parameters that were set that do not belong to chosen tree layout
        if treeType=='rectangular':
            if circStart is not None: warnings.warn('Rectangular trees do not have a circStart parameter, ignoring')
            if circFrac is not None: warnings.warn('Rectangular trees do not have a circFrac parameter, ignoring')
            if inwardSpace is not None: warnings.warn('Rectangular trees do not have an inwardSpace parameter, ignoring')
            if normaliseHeight is not None: warnings.warn('Rectangular trees do not have a normaliseHeight parameter, ignoring')
            
        elif treeType=='circular':
            if inwardSpace is None: inwardSpace = 0.0
            if inwardSpace<0: inwardSpace-=self.treeHeight

            if circStart is None: circStart = 0.0
            if circFrac is None: circFrac = 1.0

            allXs=list(map(xCoordinateFxn,self.Objects)) + sum([list(map(xCoordinateFxn,w.subtree)) for w in self.Objects if isinstance(w,Clade)],[]) ## look at both branch heights and last heights of clades
            if normaliseHeight is None: normaliseHeight=lambda value: (value-min(allXs))/(max(allXs)-min(allXs))
        
        elif treeType=='unrooted':
            if circStart is None: circStart = 0.0
            if circFrac is not None: warnings.warn('Unrooted trees do not have a circFrac parameter, ignoring')
            if inwardSpace is not None: warnings.warn('Unrooted trees do not have an inwardSpace parameter, ignoring')
            if normaliseHeight is not None: warnings.warn('Unrooted trees do not have a normaliseHeight parameter, ignoring')
                
        ########
        for k in filter(targetFxn, self.Objects):
            if treeType=='rectangular':
                
                xs.append(xCoordinateFxn(k))
                ys.append(yCoordinateFxn(k))
            
            elif treeType=='circular':

                x=normaliseHeight(xCoordinateFxn(k)+inwardSpace) ## find normalised x position along circle's radius
                y=yCoordinateFxn(k) ## get y position along circle's perimeter

                px,py=project_to_polar(x=x,y=y,yRange=self.ySpan,circleStart=circStart,circleFraction=circFrac) ## polar x and polar y coordinates for current branch

                xs.append(px)
                ys.append(py)
            
            elif treeType=='unrooted':
                xs.append(k.x) ## unrooted trees should already have x and y coordinates computed 
                ys.append(k.y)
            colours.append(colourFxn(k))
            sizes.append(pointSizeFxn(k))

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outlineColourFxn(k))
                outline_sizes.append(outlineSizeFxn(k))

        if orientation=='vertical': ## flip x and y coordinates if plotting horizontally
            if treeType == 'rectangular':
                xs, ys = ys, xs
                outline_xs, outline_ys = outline_ys, outline_xs
            else:
                warnings.warn('Trees with %s layout do not have %s orientation, ignoring'%(treeType,orientation))

        if 'zorder' not in localKwargs:
            localKwargs['zorder'] = 3

        ax.scatter(
            xs,
            ys,
            s=sizes,
            facecolor=colours,
            edgecolor="none",
            **localKwargs,
        )  ## put a circle at each tip
        if outline:
            localKwargs['zorder'] -= 1 ## decrease zorder for outlines
            ax.scatter(
                outline_xs,
                outline_ys,
                s=outline_sizes,
                facecolor=outline_colours,
                edgecolor="none",
                **localKwargs,
            )  ## put a circle at each tip

        return ax

    def _plot_rectangular_clades(self,ax,xCoordinateFxn,yCoordinateFxn,endAttrFxn,targetFxn,colour,orientation,style,cladeBaseWidth,**kwargs):
        """
        Adds triangles for plotting collapsed clades.
        """

        valid_styles=['equal', 'skewed']
        assert style in valid_styles, f"Style {style} not recognised. Options are {valid_styles}"
        
        from matplotlib.patches import Polygon
        
        mainTargetFxn = lambda k: isinstance(k,Clade) and targetFxn(k) ## called on clade objects only, additional options allowed

        localKwargs = dict(kwargs)

        ###########
        w=self.ySpan*cladeBaseWidth
        
        if 'zorder' not in localKwargs: localKwargs['zorder'] = 2

        for k in filter(mainTargetFxn,self.Objects):
            lower_left=(xCoordinateFxn(k),yCoordinateFxn(k)-w)
            upper_left=(xCoordinateFxn(k),yCoordinateFxn(k)+w)
            
            upper_right=(endAttrFxn(k),yCoordinateFxn(k)+k.width/2)
            if style=='equal':
                lower_right=(endAttrFxn(k),yCoordinateFxn(k)-k.width/2)
            elif style=='skewed':
                early_end=min([xCoordinateFxn(q) for q in k.subtree if q.is_leaf()])
                lower_right=(early_end,yCoordinateFxn(k)-k.width/2)
            
            coords=[upper_left,upper_right,lower_right,lower_left]
            if orientation=='vertical':
                coords=[c[::-1] for c in coords]
                
            fc=colour(k) if callable(colour) else colour
            ec='k'

            if 'facecolor' not in localKwargs and 'fc' not in localKwargs:
                localKwargs['facecolor'] = fc
            if 'edgecolor' not in localKwargs and 'ec' not in localKwargs:
                localKwargs['edgecolor'] = ec

            poly=Polygon(coords,**localKwargs)
            ax.add_patch(poly)
            
        return ax

    def _plot_circular_clades(self,ax,xCoordinateFxn,yCoordinateFxn,endAttrFxn,targetFxn,colour,widthFxn,circStart,circFrac,inwardSpace,normaliseHeight,precision, 
        shape,style,cladeBaseWidth,**kwargs):
        """
        Adds triangles for plotting collapsed clades.
        """
        valid_styles=['equal', 'skewed']
        assert style in valid_styles, f"Style {style} not recognised. Options are {valid_styles}"
        
        valid_shapes=['triangle', 'rectangle']
        assert shape in valid_shapes, f"Clade shape {shape} not recognised. Options are {valid_shapes}"
        
        from matplotlib.patches import Polygon
        
        mainTargetFxn = lambda k: isinstance(k,Clade) and targetFxn(k)
        
        total_y = self.ySpan

        localKwargs=dict(kwargs)
        
        linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

        ##############
        w=self.ySpan*cladeBaseWidth
        
        for k in filter(mainTargetFxn,self.Objects):
            
            ###############
            x = normaliseHeight(xCoordinateFxn(k)+inwardSpace) ## x coordinate where clade begins
            y = yCoordinateFxn(k)

            ylo_right = y-k.width/2 ## lower y coordinate for clade end
            yhi_right = y+k.width/2 ## upper y coordinate for clade end

            ylo_left = y-w/2 ## trapezoid_base width / 2
            yhi_left = y+w/2

            clade_base = [project_to_polar(x=x,y=ylo_left,yRange=total_y,circleStart=circStart,circleFraction=circFrac),project_to_polar(x=x,y=yhi_left,yRange=total_y,circleStart=circStart,circleFraction=circFrac)] ## line defining the beginning of the clade

            end_x=normaliseHeight(endAttrFxn(k)+inwardSpace) ## x coordinate end where clade finishes

            hi_arc = [project_to_polar(x=end_x,y=arc_y,yRange=total_y,circleStart=circStart,circleFraction=circFrac) for arc_y in linspace(ylo_right,yhi_right,precision)] ## arc at the end of the clade
            lo_arc = [project_to_polar(x=x,y=arc_y,yRange=total_y,circleStart=circStart,circleFraction=circFrac) for arc_y in linspace(ylo_right,yhi_right,precision)] ## arc at the beginning of the clade

            #####################################
            if shape=='triangle': ## triangles start at a sharp-ish point
                coords = list(clade_base)

            elif shape=='rectangle': ## rectangles start with an arc
                coords = list(lo_arc)
            ######################################
            if style=='equal': ## if clade style is equal add an arc at the highest point
                coords += hi_arc[::-1]

            elif style=='skewed': ## if clade style is skewed need to interpolate smooth curve
                early_end=min([xCoordinateFxn(q) for q in k.subtree if q.is_leaf()])
                early_x=normaliseHeight(early_end+inwardSpace)

                ############
                blended_arc = [project_to_polar(x=arc_x,y=arc_y,yRange=total_y,circleStart=circStart,circleFraction=circFrac) for arc_x,arc_y in zip(linspace(end_x,early_x,precision),linspace(ylo_right,yhi_right,precision))] ## arc at the end of the clade
                
                coords += blended_arc[::-1]
            ######################################

            fc=colour(k) if callable(colour) else colour
            ec='k'
            
            if 'zorder' not in localKwargs: localKwargs['zorder']=10
            if 'capstyle' not in localKwargs: localKwargs['capstyle']='round'
            if 'edgecolor' not in localKwargs: localKwargs['edgecolor']=ec
            if 'facecolor' not in localKwargs: localKwargs['facecolor']=fc
            if widthFxn is not None: localKwargs['linewidth']=widthFxn(k) if callable(widthFxn) else widthFxn

            poly=Polygon(coords,**localKwargs)
            ax.add_patch(poly)
            
        return ax

    def _plot_rectangular_tree(self,ax,connectionType,xCoordinateFxn,yCoordinateFxn,targetFxn,widthFxn,colourFxn,orientation, 
                                plotClades,cladeColour,cladeEndAttrFxn,cladeStyle,cladeBaseWidth,**kwargs):


        localKwargs = dict(kwargs)

        branches=[]
        colours=[]
        linewidths=[]

        for k in filter(targetFxn,self.Objects): ## iterate over branches
            x = xCoordinateFxn(k) ## get branch x position
            xp = xCoordinateFxn(k.parent) if k.parent.parent else x ## get parent x position
            y = yCoordinateFxn(k) ## get y position

            try:
                colours.append(colourFxn(k)) # if callable(colour) else colours.append(colour)
            except KeyError: ## not guaranteed to be a KeyError
                colours.append((0.7,0.7,0.7)) ## in case no colour available for branch set it to grey
            
            linewidths.append(widthFxn(k)) if callable(widthFxn) else linewidths.append(widthFxn)

            coords=[]
            if connectionType=='baltic': ## each node has a single vertical line to which descendant branches are connected
                coords.append([(xp,y),(x,y)])
                if k.is_node():
                    yl,yr=yCoordinateFxn(k.children[0]),yCoordinateFxn(k.children[-1]) ## y positions of first and last child

                    coords.append([(x,yl),(x,yr)])
                    linewidths.append(linewidths[-1])
                    colours.append(colours[-1])
            
            elif connectionType=='elbow': ## more standard connection where each branch connects to its parent via a right-angled line
                yp = yCoordinateFxn(k.parent) if k.parent.parent else y ## get parent x position
            
                coords.append([(xp,yp),(xp,y),(x,y)])
            
            elif connectionType=='direct': ## this gives triangular looking trees where descendants connect directly to their parents
                yp = yCoordinateFxn(k.parent) if k.parent.parent else y ## get y position
            
                coords.append([(xp,yp),(x,y)])
            
            else:
                pass  ## for now

            if orientation=='vertical': ## flip x and y coordinates if plotting horizontally
                coords=[[coord[::-1] for coord in branch_coords] for branch_coords in coords] ## invert coordinates if plotting horizontally
            
            branches+=coords ## store branch for plotting

            if isinstance(k,Clade) and plotClades == True:
                ax = self._plot_rectangular_clades(ax,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,endAttrFxn=cladeEndAttrFxn,targetFxn=targetFxn,colour=cladeColour,orientation=orientation,style=cladeStyle,cladeBaseWidth=cladeBaseWidth,**localKwargs)

        if 'capstyle' not in localKwargs:
            localKwargs['capstyle']='projecting'
        if 'zorder' not in localKwargs:
            localKwargs['zorder'] = 1
        line_segments = LineCollection(branches,lw=linewidths,color=colours,**localKwargs)
        ax.add_collection(line_segments)

        return ax

    

    def _plot_circular_tree(self,ax,targetFxn,xCoordinateFxn,yCoordinateFxn,widthFxn,colourFxn,
                            circStart,circFrac,inwardSpace,normaliseHeight,connectionType,padNodes,precision,plotClades,cladeColour, 
                            cladeEndAttrFxn,cladeStyle,cladeShape,cladeBaseWidth,**kwargs):
        if cladeColour is None: cladeColour = (0.7,0.7,0.7)
        if cladeShape is None: cladeShape = 'triangle'
        if cladeEndAttrFxn is None: cladeEndAttrFxn = lambda k: max([xCoordinateFxn(desc) for desc in k.subtree])

        assert connectionType in [
            'baltic', 
            'direct', 
            'elbow'
        ], f'Unrecognised connection type "{connectionType}"'

        assert cladeShape in [
            'triangle',
            'rectangle'
        ], f'Unrecognised clade shape "{cladeShape}"'

        assert circFrac>0.0,'Circular tree layout not given any space (circFrac == %s)'%(circFrac)
        # if inwardSpace<0: inwardSpace-=self.treeHeight

        localKwargs=dict(kwargs)
        
        # self._drawTree(pad_nodes=pad_nodes) ## computes y coordinates with padding
        total_y = self.ySpan

        allXs=list(map(xCoordinateFxn,self.Objects)) + sum([list(map(xCoordinateFxn,w.subtree)) for w in self.Objects if isinstance(w,Clade)],[]) ## look at both branch heights and last heights of clades
        if normaliseHeight is None: normaliseHeight=lambda value: (value-min(allXs))/(max(allXs)-min(allXs))
        linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

        branches=[]
        colours=[]
        linewidths=[]

        for k in filter(targetFxn,self.Objects): ## iterate over branches

            x = normaliseHeight(xCoordinateFxn(k)+inwardSpace) ## get branch x position
            xp = normaliseHeight(xCoordinateFxn(k.parent)+inwardSpace) if k.parent.parent else x ## get parent x position
            y = yCoordinateFxn(k) ## get y position
            yp = yCoordinateFxn(k.parent) if k.parent.parent else y ## get parent y position


            try:
                colours.append(colourFxn(k))
            except KeyError:
                colours.append((0.7,0.7,0.7))
            linewidths.append(widthFxn(k)) if callable (widthFxn) else linewidths.append(widthFxn)

            px,py = project_to_polar(x=x,y=y,yRange=total_y,circleStart=circStart,circleFraction=circFrac) ## polar x and polar y coordinates for current branch
            pxp,pyp = project_to_polar(x=xp,y=yp,yRange=total_y,circleStart=circStart,circleFraction=circFrac) ## polar x and polar y coordinates for branch parent

            if connectionType=='baltic':

                branches.append((project_to_polar(x=xp,y=y,yRange=total_y,circleStart=circStart,circleFraction=circFrac),(px,py))) ## horizontal line from parent to current branch

                if k.is_node():
                    yl,yr=yCoordinateFxn(k.children[0]),yCoordinateFxn(k.children[-1]) ## get leftmost and rightmost children's y coordinates

                    coordinates = tuple([project_to_polar(x=x,y=arc_y,yRange=total_y,circleStart=circStart,circleFraction=circFrac) for arc_y in linspace(yl,yr,precision)]) ## compute polar coordinates for arc of node

                    branches.append(coordinates) ## add curved segment at node

                    linewidths+=[linewidths[-1] for q in coordinates] ## repeat linewidths
                    colours+=[colours[-1] for q in coordinates] ## repeat colours

            elif connectionType=='direct':

                branches.append(((pxp,pyp),(px,py))) ## straight line from parent to current branch
                linewidths+=[linewidths[-1]] ## repeat linewidths
                colours+=[colours[-1]] ## repeat colours

            elif connectionType=='elbow':

                coordinates = tuple([(px,py)]+[project_to_polar(x=xp,y=arc_y,yRange=total_y,circleStart=circStart,circleFraction=circFrac) for arc_y in linspace(y,yp,precision)]) ## compute polar coordinates for arc at parent + horizontal line to current branch

                branches.append(coordinates) ## straight line from parent to current branch
                linewidths+=[linewidths[-1] for q in coordinates] ## repeat linewidths
                colours+=[colours[-1] for q in coordinates] ## repeat colours

            if isinstance(k,Clade) and plotClades==True:
                # self,ax,xCoordinateFxn=None,yCoordinateFxn=None,end_attr=None,width=None,target=None,colour=None,circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,precision=15, 
        # clade_shape='triangle',style='equal',trapezoid_base=0.003,pad_nodes=None,**kwargs
                ax = self._plot_circular_clades(ax,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,endAttrFxn=cladeEndAttrFxn,targetFxn=targetFxn,colour=cladeColour,widthFxn=widthFxn,circStart=circStart,circFrac=circFrac, 
                                                inwardSpace=inwardSpace,normaliseHeight=normaliseHeight,precision=precision,shape=cladeShape,style=cladeStyle,cladeBaseWidth=cladeBaseWidth,**localKwargs)

        if 'capstyle' not in localKwargs: localKwargs['capstyle'] = 'projecting'
        if 'zorder' not in localKwargs: localKwargs['zorder'] = 1
        line_segments = LineCollection(branches,lw=linewidths,color=colours,**localKwargs) ## create line segments
        ax.add_collection(line_segments) ## add collection to axes
        return ax

    def plot_tree(self,ax,targetFxn=None,
                xCoordinateFxn=None,yCoordinateFxn=None,width=None,widthFxn=None,connectionType=None,
                colour=None,colourFxn=None,orientation=None,padNodes=None,treeType='rectangular', 
                circStart=None,circFrac=None,inwardSpace=None,normaliseHeight=None,precision=None, 
                plotClades=True,cladeColour=None,cladeEndAttrFxn=None,cladeStyle='equal',cladeShape=None,cladeBaseWidth=0.001,recomputeCoordinates=True,autoSort=True,**kwargs):
        ### Set default values ###
        if targetFxn is None:
            targetFxn=lambda k: True
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.height if self.treeType == 'divergence' else k.absoluteTime
        if yCoordinateFxn is None:
            yCoordinateFxn = lambda k: k.y

        if padNodes is None:
            padNodes = {}

        # Width logic
        if width is not None and widthFxn is not None:
            raise ValueError(
                "Cannot specify both width and widthFxn. Please use only one."
            )
        if width is None and widthFxn is None:
            widthFxn = lambda k: 2
        elif widthFxn is None:
            widthFxn = lambda k: width
        
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            )
        if colour is None and colourFxn is None:
            colourFxn = lambda k: "k"
        elif colourFxn is None:
            colourFxn = lambda k: colour

        if cladeColour is None: cladeColour=(0.7,0.7,0.7) ## calling with defaults, clade colour will be light gray
        if cladeEndAttrFxn is None: cladeEndAttrFxn = lambda k: max([xCoordinateFxn(desc) for desc in k.subtree])

        assert treeType in [
            'rectangular', 
            'circular', 
            'unrooted'
        ], f'Unrecognised tree type: "%s"'%(treeType)
        
        localKwargs = dict(kwargs)

        if autoSort == True:
            self.sort_branches()
        
        if recomputeCoordinates:
            self._assign_tree_coordinates(padNodes=padNodes)

        if treeType != 'unrooted':
            if connectionType is None:
                connectionType = 'baltic'

            assert connectionType in [
                    'baltic', 
                    'direct', 
                    'elbow'
            ], f'Unrecognised drawing type "%s"'%(connectionType)
        
        ################################
        if treeType=='rectangular':
            if orientation is None:
                orientation = 'horizontal'
            else:
                assert orientation in [
                    'horizontal', 
                    'vertical'
                ], f'Unrecognised tree orientation "%s"'%(orientation)

            if precision is not None: warnings.warn('Rectangular trees do not have a precision parameter, ignoring')
            if circStart is not None: warnings.warn('Rectangular trees do not have a circStart parameter, ignoring')
            if circFrac is not None: warnings.warn('Rectangular trees do not have a circFrac parameter, ignoring')
            if inwardSpace is not None: warnings.warn('Rectangular trees do not have a inwardSpace parameter, ignoring')
            if normaliseHeight is not None: warnings.warn('Rectangular trees do not have a normaliseHeight parameter, ignoring')
            if cladeShape is not None: warnings.warn('Rectangular trees do not support alternative clade shapes, ignoring') ## does this make sense?

            ax = self._plot_rectangular_tree(ax=ax,connectionType=connectionType,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,targetFxn=targetFxn,
                                                widthFxn=widthFxn, colourFxn=colourFxn,orientation=orientation,
                                                cladeColour=cladeColour,cladeEndAttrFxn=cladeEndAttrFxn,plotClades=plotClades, 
                                                cladeStyle=cladeStyle,cladeBaseWidth=cladeBaseWidth,**localKwargs)
            
        ###############################
        elif treeType=='circular':
            if circStart is None: circStart = 0.0
            if circFrac is None: circFrac = 1.0

            if inwardSpace is None: inwardSpace = 0.0
            if inwardSpace<0: inwardSpace-=self.treeHeight

            if precision is None: precision = 15
            if orientation is not None: warnings.warn('Circular trees do not have an orientation parameter, ignoring')

            ax = self._plot_circular_tree(ax=ax,connectionType=connectionType,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,targetFxn=targetFxn,widthFxn=widthFxn, colourFxn=colourFxn, 
                                            circStart=circStart, circFrac=circFrac, inwardSpace=inwardSpace, normaliseHeight=normaliseHeight, precision=precision, 
                                            padNodes=padNodes, plotClades=plotClades,
                                            cladeColour=cladeColour, cladeEndAttrFxn=cladeEndAttrFxn,cladeStyle=cladeStyle,cladeShape=cladeShape,cladeBaseWidth=cladeBaseWidth,**localKwargs)

            ax.set_aspect(1)
        ###########################
        elif treeType=='unrooted':
            if circStart is None: circStart = 0.0
            if precision is None: precision = 15
            if orientation is not None: warnings.warn('Unrooted trees do not have an orientation parameter, ignoring')
            if circFrac is not None: warnings.warn('Unrooted trees do not have a circFrac parameter, ignoring')
            if inwardSpace is not None: warnings.warn('Unrooted trees do not have a inwardSpace parameter, ignoring')
            if normaliseHeight is not None: warnings.warn('Unrooted trees do not have a normaliseHeight parameter, ignoring')
            if cladeShape is not None: warnings.warn('Unrooted trees do not support alternative clade shapes, ignoring')

            if connectionType is None:
                connectionType = 'direct'
            else:
                warnings.warn('Unrooted trees use a direct connectionType, ignoring')
            
            if recomputeCoordinates:
                self._assign_unrooted_tree_coordinates(circStart=circStart,padNodes=padNodes)
            
            xCoordinateFxn = lambda k: k.x ## using projected x coordinate now

            ax = self._plot_rectangular_tree(ax=ax,connectionType='direct',xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,targetFxn=targetFxn,orientation=None, ## only allow .x and .y here?
                                                widthFxn=widthFxn, colourFxn=colourFxn,
                                                cladeColour=cladeColour,cladeEndAttrFxn=cladeEndAttrFxn,cladeBaseWidth=cladeBaseWidth,cladeStyle=cladeStyle,plotClades=False,**localKwargs)

            if plotClades==True:
                ax = self._plot_unrooted_clades(ax=ax,xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,endAttrFxn=cladeEndAttrFxn,targetFxn=targetFxn,colour=cladeColour,precision=precision,style=cladeStyle,cladeBaseWidth=cladeBaseWidth,**localKwargs)

            ax.set_aspect(1)

        ax.autoscale() ## makes sure plot limits are auto-adjusted to include collections
        return ax

    def plot_exploded_tree(self,ax,trait=None,customFxn=None,stem=True,subtreeSortFxn=None,
                       targetFxn=None,
                       xCoordinateFxn=None,colour=None,colourFxn=None,tipPoints=True,originPoint=True,
                       verticalSpace=2,
                       width=None,widthFxn=None,
                       pointSize=None,pointSizeFxn=None,
                       outline=True,outlineSize=None,outlineSizeFxn=None,
                       outlineColour=None,outlineColourFxn=None,
                       padNodes=None,orientation='horizontal',connectionType='baltic',**kwargs):
        #### TODO: support for collapsed clades
        
        ### Set default values ###
        if targetFxn is None:
            targetFxn = lambda k: True
        if xCoordinateFxn is None:
            xCoordinateFxn = lambda k: k.height if self.treeType == 'divergence' else k.absoluteTime

        # Width logic
        if width is not None and widthFxn is not None:
            raise ValueError(
                "Cannot specify both width and widthFxn. Please use only one."
            )
        if width is None and widthFxn is None:
            widthFxn = lambda k: 2
        elif widthFxn is None:
            widthFxn = lambda k: width
        
        # Point size logic
        if pointSize is not None and pointSizeFxn is not None:
            raise ValueError(
                "Cannot specify both pointSize and pointSizeFxn. Please use only one."
            )
        if pointSize is None and pointSizeFxn is None:
            pointSizeFxn = lambda k: 40
        elif pointSizeFxn is None:
            pointSizeFxn = lambda k: pointSize
        # Outline size logic
        if outlineSize is not None and outlineSizeFxn is not None:
            raise ValueError(
                "Cannot specify both outlineSize and outlineSizeFxn. Please use only one."
            )
        if outlineSize is None and outlineSizeFxn is None:
            outlineSizeFxn = lambda k: ((2+math.sqrt(pointSizeFxn(k)/math.pi))**2)*math.pi
        elif outlineSizeFxn is None:
            outlineSizeFxn = lambda k: outlineSize
        # Colour logic
        if colour is not None and colourFxn is not None:
            raise ValueError(
                "Cannot specify both colour and colourFxn. Please use only one."
            ) ## should be a warning, since this eventuality is handled in the next line
        if colour is None and colourFxn is None:
            colourFxn = lambda k: "k"
            warnings.warn('Without an explicitly specified colourFxn each exploded subtree will be shown in the same colour.')
        elif colourFxn is None:
            colourFxn = lambda k: colour
            warnings.warn('Without an explicitly specified colourFxn each exploded subtree will be shown in the specified colour.')
        # Outline colour logic
        if outlineColour is not None and outlineColourFxn is not None:
            raise ValueError(
                "Cannot specify both outlineColour and outlineColourFxn. Please use only one."
            )
        if outlineColour is None and outlineColourFxn is None:
            outlineColourFxn = lambda k: "k"
        elif outlineColourFxn is None:
            outlineColourFxn = lambda k: outlineColour

        if padNodes is None:
            padNodes = {}
        
        assert orientation in [
            'horizontal', 
            'vertical'
        ], f"Unrecognised orientation {orientation}"
        
        assert connectionType in [
                'baltic', 
                'direct', 
                'elbow'
        ], f'Unrecognised drawing type "%s"'%(connectionType)
        
        if trait is not None and customFxn is not None:
            raise ValueError(
                "Cannot specify both a trait and a custom function to split the tree by. Please use only one."
            )
        elif trait is None and customFxn is None:
            raise ValueError(
                "Must specify either a trait or a custom function to split the tree by."
            )
        if trait is not None:
            traceFunction = lambda k: k == self.root or k.traits[trait] != k.parent.traits[trait] ## current branch trait not the same as parent or at root
        elif customFxn:
            traceFunction = customFxn

        localKwargs=dict(kwargs)
            
        treeList = self.explode_tree(customFxn = traceFunction, stem=stem)
        
        ######### plotting
        cumulative_y = 0

        if subtreeSortFxn is None:
            subtreeSortFxn = lambda subtree: -subtree.root.absoluteTime if subtree.treeType == 'time' else -subtree.root.height ## sort trees by their root time

        yCoordinateFxn = lambda k: k.y + cumulative_y

        if orientation == 'vertical':
            xCoordinateFxn, yCoordinateFxn = yCoordinateFxn, xCoordinateFxn
        
        for subtree in sorted(treeList, key = subtreeSortFxn):
            subtree.plot_tree(ax,targetFxn=targetFxn,
                              xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,
                              widthFxn=widthFxn,connectionType=connectionType,colourFxn=colourFxn,orientation=orientation,padNodes=padNodes)
            if tipPoints:
                subtree.plot_points(ax,targetFxn = lambda k: k.is_leaf(),
                                    xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,
                                    colourFxn=colourFxn,
                                    pointSizeFxn=pointSizeFxn,
                                    outline=outline,
                                    outlineSizeFxn=outlineSizeFxn,
                                    outlineColourFxn=outlineColourFxn,
                                    padNodes=padNodes,orientation=orientation)
            if originPoint:
                subtree.plot_points(ax,targetFxn = lambda k: k == subtree.root,
                                    xCoordinateFxn=xCoordinateFxn,yCoordinateFxn=yCoordinateFxn,
                                    colourFxn=colourFxn,
                                    pointSizeFxn=pointSizeFxn,
                                    outline=outline,
                                    outlineSizeFxn=outlineSizeFxn,
                                    outlineColourFxn=outlineColourFxn,
                                    padNodes=padNodes,orientation=orientation)
            
            cumulative_y += subtree.ySpan + verticalSpace

        ax.autoscale()
        return ax

    

    def to_auspice_json(self, traits=None, mostRecentDate=None):
        from baltic.bt_utils import branch_to_json

        if traits is None:
            traits = set([item for k in self.Objects for item in list(k.traits.keys())]) ## extract all traits dict keys across tree
            logger.warning(f"Traits identified in baltic object that will be included in the auspice JSON: {', '.join(traits)}")
        else:
            logger.warning(f"User-provided trait list: {', '.join(traits)}\nNote that uncertainties (e.g. parameter ranges or state probabilities) might not be parsed if you don't include '_95%_HPD' and '.set.prob' traits too.")

        if mostRecentDate is None and self.treeType == 'time':
            mostRecentDate = self.mostRecent
            assert mostRecentDate, f"mostRecentDate could not be recovered from time tree, provide it manually"

        assert isinstance(mostRecentDate, float), f"mostRecentDate is of type {type(mostRecentDate)}, not a float"
        if isinstance(mostRecentDate, float) and self.treeType == 'divergence':
            logger.warning(f"mostRecentDate set to {mostRecentDate} but treeType is set to {treeType}")
        
        #####
        auspiceJSON = {}
        auspiceJSON['version'] = 'v2'
        auspiceJSON['meta'] = {}
        auspiceJSON['meta']['colorings'] = []
        
        #### process traits - create meta field for colouring, reduce trait set to those that can be extracted from the traits dict of branches
        traitsMeta = {}
        for trait in traits:
            if trait.endswith('_95%_HPD'):
                traitName = trait.split('_95%_HPD')[0]
                traitsMeta[traitName] = {'key': traitName, 'title': traitName, 'type': 'continuous'}

            elif trait.endswith('.set.prob'):
                traitName = trait.split('.set.prob')[0]
                traitsMeta[traitName] = {'key': traitName, 'title': traitName, 'type': 'categorical'}

            elif trait.endswith('_median'):
                traitsMeta[trait] = {'key': trait, 'title': trait, 'type': 'continuous'}

            elif trait in ['posterior']:
                traitsMeta[trait] = {'key': trait, 'title': trait, 'type': 'continuous'}

        treeTraits = [traitsMeta[t]['key'] for t in traitsMeta if (traitsMeta[t]['key'].endswith('_95%_HPD') == False and traitsMeta[t]['key'].endswith('.set.prob') == False)]
        
        for trait in traitsMeta:
            auspiceJSON['meta']['colorings'].append(traitsMeta[trait])
        ####
        from importlib.metadata import version
        balticVersion = version('baltic')
        auspiceJSON['meta']['data_provenance'] = [{'name': 'baltic', 'version': balticVersion, 'url': 'https://pypi.org/project/baltic/'}]
        ########
        nodeOrder = self.traverse_tree(includeCondition=lambda k: k.is_node())
        for k in self.get_internal():
            k.traits['node_idx'] = f"NODE_{nodeOrder.index(k) + 1:07d}"

        jsonTree = branch_to_json(curNode=self.root, treeType=self.treeType, traits=treeTraits, mostRecentDate=mostRecentDate)
        auspiceJSON['tree'] = jsonTree
        
        return auspiceJSON
        
