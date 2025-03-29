import copy,math
import numpy as np
from functools import reduce
from matplotlib.collections import LineCollection
from .node import Node
from .leaf import Leaf
from .clade import Clade
from .reticulation import Reticulation
from .bt_utils import root_to_tip

class Tree: ## tree class

    def __init__(self):
        self.curNode=Node() ## current node is a new instance of a node class
        self.curNode.index='Root' ## first object in the tree is the root to which the rest gets attached
        self.curNode.length=0.0 ## startind node branch length is 0
        self.curNode.height=0.0 ## starting node height is 0
        self.root=None #self.cur_node ## root of the tree is current node
        self.Objects=[] ## tree objects have a flat list of all branches in them
        self.tipMap=None
        self.treeHeight=0 ## tree height is the distance between the root and the most recent tip
        self.mostRecent=None
        self.ySpan=0.0
        self.treeType=None # indicates if the tree is a time or divergence tree

    def add_reticulation(self,name):
        ret=Reticulation(name)
        ret.index=name
        ret.parent=self.curNode
        self.curNode.children.append(ret)
        self.Objects.append(ret)
        self.curNode=ret

    def add_node(self,i):
        newNode=Node() ## new node instance
        newNode.index=i ## new node's index is the position along the tree string
        if self.root is None:
            self.root=newNode
            self.root.length=0.0

        newNode.parent=self.curNode ## new node's parent is current node
        assert self.curNode.is_node(), 'Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses or commas.'
        self.curNode.children.append(newNode) ## new node is a child of current node
        self.curNode=newNode ## current node is now new node
        self.Objects.append(self.curNode) ## add new node to list of objects in the tree

    def add_leaf(self,i,name):
        newLeaf=Leaf(name) ## new instance of leaf object
        newLeaf.index=i ## index is position along tree string
        if self.root is None: self.root=newLeaf

        newLeaf.parent=self.curNode ## leaf's parent is current node
        assert self.curNode.is_node(), 'Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses.'
        self.curNode.children.append(newLeaf) ## assign leaf to parent's children
        # new_leaf.name=name
        self.curNode=newLeaf ## current node is now new leaf
        self.Objects.append(self.curNode) ## add leaf to all objects in the tree

    def subtree(self,startingNode=None,traverseCondition=None,stem=True):
        if startingNode is None: startingNode = self.root
        if traverseCondition is None: traverseCondition = lambda k: True

        node = startingNode.parent if stem else startingNode ## move up a node if we want the stem

        subtreeBranches=self.traverse_tree(node,includeCondition=lambda k:True,traverseCondition=traverseCondition)
        subtreeBranches=copy.deepcopy(subtreeBranches)

        if stem: ## using stem - need to prune subtrees from root now
            unwantedBranches=[]
            for child in node.children: ## iterate over parent's children
                if child.index!=startingNode.index: ## not at focal branch (unwanted sibling)
                    unwantedBranches+=self.traverse_tree(child,includeCondition=lambda w: True) ## add all branches resulting from traversals of unwanted siblings

            remove=[k for k in subtreeBranches if k.index in [ub.index for ub in unwantedBranches]] ## iterate over subtree branches, remember those that belong to unwanted subtrees 
            for r in remove: ## iterate over branches belong to unwanted subtrees
                subtreeBranches.remove(r) ## remove from list

        if subtreeBranches is None or not any(node.is_leaf() for node in subtreeBranches): ## nothing found or no leaf objects in traversal
            return None

        localTree=Tree() ## create a new tree object where the subtree will be
        localTree.Objects=subtreeBranches ## assign branches to new tree object

        localTree.root=subtreeBranches[0] ## root is the beginning of the traversal
        localTree.root.parent=None ## tree begins strictly at node

        subtreeSet=set(subtreeBranches) ## turn branches into set for quicker look up later

        if traverseCondition is not None: ## didn't use default traverse condition, might need to deal with hanging nodes and prune children
            for nd in localTree.get_internal(): ## iterate over nodes
                nd.children=[child for child in nd.children if child in subtreeSet] ## only keep children seen in traversal
            localTree.fix_hanging_nodes()

        if self.tipMap: ## if original tree has a tipMap dictionary
            localTree.tipMap={tipNum: self.tipMap[tipNum] for tipNum in self.tipMap if self.tipMap[tipNum] in [w.name for w in localTree.get_external()]} ## copy over the relevant tip translations

        return localTree

    def make_single_type(self,verbose=False):
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
            multitypeNodes=self.get_internal(lambda k: len(k.children)==1)

            if verbose:
                print('Number of multitype nodes in tree: %d'%(len(multitypeNodes)))
            if not multitypeNodes:
                break

            for k in sorted(multitypeNodes,key=lambda x:-x.height):

                child=k.children[0] ## fetch child
                grandparent=k.parent if k.parent.index else self.root ## fetch grandparent
                if verbose:
                    print('At multitype node %s with child %s and grandparent %s'%(k.index,child.index,grandparent.index))
                child.parent=grandparent ## child's parent is now grandparent
                grandparent.children.append(child) ## add child to grandparent's children
                grandparent.children.remove(k) ## remove old parent from grandparent's children
                grandparent.children=list(set(grandparent.children))

                child.length+=k.length ## adjust child length

                self.Objects.remove(k) ## remove old parent from all objects
        self.sort_branches()

    def set_absolute_time(self,date):
        """
        Places all objects in absolute time by providing the date of the most recent tip.
        
        Parameters:
        date (float): The date of the most recent tip in the tree in decimal date format (see the `decimalDate()` function).
        
        This method calculates the absolute time for each object in the tree based on its height and
        the provided date of the most recent tip (most recent date - tree height + node/leaf height). The absolute time is stored in the `absoluteTime`
        attribute of each object.
        
        Example:
        >>> tree.set_absolute_time(2023.0)

        Docstring generated with ChatGPT 4o.
        """
        for k in self.Objects: ## iterate over all objects
            k.absoluteTime=date - self.treeHeight + k.height ## heights are in units of time from the root
        self.mostRecent=max(k.absoluteTime for k in self.Objects)

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
        stats = self._calculate_tree_stats()
        print('\nTree height: %.6f\nTree length: %.6f'%(stats["treeHeight"], stats["treeLength"])) ## report the height and length of tree

        if stats["strictlyBifurcating"]: print('strictly bifurcating tree')
        if stats["multitypeTree"]: print('multitype tree')
        if stats["singletonTree"]: print('singleton tree')
        if stats["hasTraits"]: print('annotations present')

        print('\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n'%(stats["numObjects"],stats["numNodes"],stats["numLeaves"])) ## report numbers of different objects in the tree

    def treeStatsDict(self):
        return self._calculate_tree_stats()

    def _calculate_tree_stats(self):
        """Return a dictionary of tree stats for being plotted in the treeStats() method
        """

        stats = {}

        self.traverse_tree() ## traverse the tree
        obs=self.Objects ## convenient list of all objects in the tree
        stats["treeHeight"] = self.treeHeight
        stats["treeLength"] = sum([x.length for x in obs])

        nodes=self.get_internal() ## get all nodes
        strictlyBifurcating=all(len(x.children) == 2 for x in nodes) ## assume tree is not strictly bifurcating
        multiType=any(len(x.children) == 1 for x in nodes)
        singleton=len(nodes) == 0

        hasTraits=False ## assume tree has no annotations
        max_annotations = max(len(k.traits) for k in obs) ## check the largest number of annotations any branch has
        hasTraits = max_annotations > 0

        stats["strictlyBifurating"] = True if strictlyBifurcating else False
        stats["multitypeTree"] = True if multiType else False
        stats["singletonTree"] = True if singleton else False
        stats["hasTraits"] = True if hasTraits else False

        stats["numObjects"] = len(obs)
        stats["numNodes"] = len(nodes)
        stats["numLeaves"] = len(self.get_external())

        return stats

    def traverse_tree(self,curNode=None,includeCondition=None,traverseCondition=None,collect=None,verbose=False):
        """
        Traverses the tree starting from a specified node and collects nodes based on conditions.
        
        Parameters:
        curNode (node or None): The starting node for traversal. If None, starts from the root. Default is None.
        includeCondition (function or None): A function that determines whether a node should be included in the `collect` list.
                                              Default is None, which includes leaf-like nodes.
        traverseCondition (function or None): A function that determines whether a child node should be traversed.
                                               Default is None, which traverses all nodes.
        collect (list): A list to collect nodes that meet the include condition. Default is None, which collects and returns only leaf and leaf-like (`Clade` and `Reticulation`) objects.
        verbose (bool): If True, prints verbose output during traversal. Default is False.
        
        Returns:
        list: A list of nodes that pass `includeCondition` that were encountered during the traversal that meets the `traverseCondition`.
        
        Example:
        >>> descendantTipList = tree.traverse_tree(includeCondition=lambda node: node.is_leaf()) ## will return a list of all leaf (but not leaf-like) objects encountered during the traversal
        
        Docstring generated with ChatGPT 4o.
        """
        if curNode==None: ## if no starting point defined - start from root
            if verbose==True: print('Initiated traversal from root')
            curNode=self.root

            if traverseCondition==None and includeCondition==None: ## reset heights if traversing from scratch
                for k in self.Objects: ## reset various parameters
                    if k.is_node():
                        k.leaves=set()
                        k.childHeight=None
                    k.height=None

        if traverseCondition==None: traverseCondition=lambda k: True
        if includeCondition==None: includeCondition=lambda k: k.is_leaflike()

        if collect==None: ## initiate collect list if not initiated
            collect=[]

        if curNode.parent and curNode.height==None: ## cur_node has a parent - set height if it doesn't have it already
            curNode.height=curNode.length+curNode.parent.height
        elif curNode.height==None: ## cur_node does not have a parent (root), if height not set before it's zero
            curNode.height=0.0

        if verbose==True: print('at %s (%s)'%(curNode.index,curNode.branchType))

        if includeCondition(curNode): ## test if interested in cur_node
            collect.append(curNode) ## add to collect list for reporting later

        if curNode.is_leaf() and self.root!=curNode: ## cur_node is a tip (and tree is not single tip)
            curNode.parent.leaves.add(curNode.name) ## add to parent's list of tips

        elif curNode.is_node(): ## cur_node is node
            for child in filter(traverseCondition,curNode.children): ## only traverse through children we're interested
                if verbose==True: print('visiting child %s'%(child.index))
                self.traverse_tree(curNode=child,includeCondition=includeCondition,traverseCondition=traverseCondition,verbose=verbose,collect=collect) ## recurse through children
                if verbose==True: print('child %s done'%(child.index))
            assert len(curNode.children)>0, 'Tried traversing through hanging node without children. Index: %s'%(curNode.index)
            curNode.childHeight=max([child.childHeight if child.is_node() else child.height for child in curNode.children])

            if curNode.parent:
                curNode.parent.leaves=curNode.parent.leaves.union(curNode.leaves) ## pass tips seen during traversal to parent
            self.treeHeight=curNode.childHeight ## it's the highest child of the starting node
        return collect


    def rename_tips(self,tipNameMap=None):
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
        if tipNameMap==None and self.tipMap!=None:
            tipNameMap=self.tipMap

        if tipNameMap is None: raise ValueError("No dictionary provided for renaming tips.")

        for k in self.get_external(): ## iterate through leaf objects in tree
            # k.name=d[k.numName] ## change its name
            k.name=tipNameMap[k.name] ## change its name

    def reroot(self, branch=None, branchFrac=0.5, fixSingletons=True, verbose=False):
        """
        Shouldn't be allowed for time trees, complex data structures or trees without branch lengths.
        Code translated from Biopython.
        If branch is None - root at midpoint.
        """
        assert self.treeType != "time", "Error: cannot reroot a time-calibrated tree."

        if branch==self.root:
            if verbose:
                print('Rerooting attempted on existing root')
            return self

        oldTreeLength=sum(self.get_parameter_list('length')) ## get old tree length
        ###############
        if branch==None: ## midpoint rooting
            if verbose:
                print('No branch provided, rooting at midpoint')
            # Identify the largest pairwise distance
            maxDistance = 0.0
            
            for tip in self.get_external(): ## iterate over tips
                
                self.reroot(branch = tip, branchFrac = 0.0)
                self.traverse_tree() ## set heights
                
                highestTip = sorted(self.get_external(),key=lambda k: k.height)[-1]
                newMax = highestTip.height
                if newMax > maxDistance: ## check if current height is bigger 
                    tip1 = tip
                    tip2 = highestTip
                    maxDistance = newMax
    #                 print('setting %s as current highest tip at height %s'%(tip2.name,max_distance))
            
            if verbose:
                print('First step of midpoint root - rooting on highest tip under current topology (%s)'%(tip1.name))
            self.reroot(branch = tip1, branchFrac = 0.0)
            
            # Depth to go from the ingroup tip toward the outgroup tip
            rootRemainder = 0.5 * maxDistance# - (self.root.length or 0))
            assert rootRemainder >= 0
            # Identify the midpoint and reroot there.
            # Trace the path to the outgroup tip until all of the root depth has
            # been traveled/accounted for.
            path=tip2.get_path_to_root()[:-1]
            
            for node in path[::-1]: ## iterate from old root to new
                rootRemainder -= node.length
    #             print('iterating over path: %s root remainder: %s'%(node.index,root_remainder))
                if rootRemainder < 0:
                    outgroup_node = node
                    branchFrac = -rootRemainder / outgroup_node.length
                    break
    #         print('rerooting on %s with %s frac'%(node.index,branch_frac))
            if verbose:
                print('Second step of midpoint root - rooting on node %s halfway from previous highest tip and current highest tip %s'%(outgroup_node.index,tip2.name))
            self.reroot(branch = outgroup_node, branchFrac = branchFrac, fixSingletons = fixSingletons)
            
            return self
        
        ##############
        path=branch.get_path_to_root()[:-1] ## get path from new root to old root, ignore actual root node
        path=path[-2::-1] ## invert
        
        ogLen=float(branch.length) ## store branch length on the new root branch
        prevBranchLength = ogLen
        branch.length = ogLen*(branchFrac) ## this modifies unchanged path branch length
        ##################

        import random,string
        characters = string.ascii_letters + string.digits  # A-Z, a-z, 0-9
        randomString=''.join(random.choices(characters, k=10))


        newRoot=Node() ## create new root
        newRoot.index='new_root_%s_%s_%s'%(branch.index,branchFrac,randomString)
        newRoot.length=0.0
        newRoot.children.append(branch)

        if verbose:
            print('Created new root %s'%(newRoot.index))
        
        ######################
        if len(path) == 1: ## new root is the old root (but maybe branch length is adjusted)
            newParent = newRoot    
        else: ## 
            parent = path.pop(-2) ## get previous branch
            parent.children.remove(branch) ## remove child
            
            prevBranchLength, parent.length = parent.length, (prevBranchLength - branch.length) ## store previous branch length, assign new branch length
            
            newRoot.children.append(parent) ## add branch as child of new root
            
            parent.parent=newRoot ## assign new parent
            newParent=parent ## set new_parent variable
        
        #######################
        for parent in path[-2::-1]: ## iterate from new root to old root
            parent.children.remove(newParent) ## remove parent (used to be child)
            parent.parent=newParent ## assign child as parent
            prevBranchLength, parent.length = parent.length, prevBranchLength ## pass branch length to the next node
            newParent.children.append(parent) ## add former parent as child
            newParent=parent ## move up to next node
        #################

        if verbose:
            print('Cleaning up old root')

        oldRoot = self.root ## get old root
        
        if branch in oldRoot.children: ## remove child that was on the way to the new root from old root's children
            assert len(path) == 1
            oldRoot.children.remove(branch)
        else:
            oldRoot.children.remove(newParent)
        
        if branch.parent==self.root: ## if we were rooting on a branch connected to the root
            oldRoot.length = prevBranchLength*(1-branchFrac)
        else:
            oldRoot.length = prevBranchLength
        
        newParent.children.append(oldRoot) ## add old root as child
        oldRoot.parent=newParent ## set parent as old root
        ########

        if verbose:
            print('Setting up root parent node (it\'s a baltic thing)')
        self.root=newRoot ## set root to new root
        
        rootParent=Node() ## baltic tree roots have an additional parent node to the root
        rootParent.index='Root'
        rootParent.length=0.0
        rootParent.height=0.0
        
        rootParent.children.append(self.root) ## set up basal baltic root parent node
        self.root.parent=rootParent
        self.Objects.append(self.root)
        
        branch.parent=newRoot
        #############
        if fixSingletons: ## fixing singleton nodes (node's with 1 child, i.e. old root)
            if verbose:
                print('Fixing singletons')
            self.make_single_type(verbose=verbose) ## tree class can handle it
        
        ############
        newTreeLength=sum([k.length for k in self.traverse_tree(includeCondition=lambda q: True)]) ## get new tree length
        assert math.isclose(oldTreeLength,newTreeLength),'Tree length changed after rerooting - was %s, now %s'%(oldTreeLength,newTreeLength) ## check if tree length hasn't changed after rerooting
        
        self.traverse_tree()
        
        return self

    def root_by_regression(self,stat='r^2',forcePositive=True,verbose=False):
        """
        Run root-to-tip regression to find optimal root based on minimising sum of squares or maximising r^2 or correlation coefficient.
        New root is placed in the middle of the branch.
        """
        validRootingMethods=['r^2','correlation','sum of squares']
        assert stat in validRootingMethods,'Invalid option for root-to-tip regression: %s (options are %s)'%(stat,validRootingMethods)
        
        cll=copy.deepcopy(self) ## deepcopy entire tree
        
        res={} ## will be used to track better regressions with new roots
        for k in self.Objects[1:]: ## iterate over branches in tree (ignore root)
            cll = cll.reroot([w for w in cll.Objects if w.index==k.index][0],fixSingletons=False,verbose=verbose) ## reroot on branch provided
            tips=cll.get_external() ## get all tips
            
            xs,ys=zip(*[(w.absoluteTime,w.height) for w in tips]) ## get collection dates
            res = root_to_tip(k,xs,ys,res,stat=stat,forcePositive=forcePositive) ## check root-to-tip regression, update res if better
        ##############
        if verbose:
            print('Rooting tree first time at node %s with regression: %s'%(res['root'],res))
        self = self.reroot(res['root'],verbose=verbose) ## reroot on branch optimising stat
        ############
        if len(self.root.children) == 2: ## if root is strictly bifurcating - also optimise branch fraction
            if verbose:
                print('Bifurcating root, finding more precise rooting along new root')
            res['frac'] = 0.0 ## add frac argument
            left,right = self.root.children ## get left and right children
            
            totalBranch = left.length + right.length ## total amount of branch length adjustment available at root
            leftSubtree = self.traverse_tree(left) ## get left subtree
            rightSubtree = self.traverse_tree(right) ## get right subtree
            
            lxs = [k.absoluteTime for k in leftSubtree] ## get x coordinates for left subtree #TODO: these don't get used
            rxs = [k.absoluteTime for k in rightSubtree] ## get x coordinates for right subtree #TODO: see above
            
            for f in np.linspace(0,1,21): ## check root positions along best overall root
                adjustL = -left.length + totalBranch * f ## height adjustments for left subtree
                adjustR = -right.length + totalBranch * (1-f) ## height adjustments for right subtree
                
                lys = [k.height+adjustL for k in leftSubtree] ## adjust left subtree heights #TODO: see above
                rys = [k.height+adjustR for k in rightSubtree] ## same for right #TODO: see above
                
                res = root_to_tip(left,xs,ys,res,stat=stat,forcePositive=forcePositive,frac=f) ## check regression
            
            self = self.reroot(branch = res['root'], branchFrac = res['frac'],verbose=verbose) ## reroot on branch optimising stat
        #############
        print('Correlation coefficient: %s\nSum of squares: %s\nEvolutionary rate: %s\nIntercept (TMRCA): %s\nr^2: %s'%(res['correlation'],res['sum of squares'],res['slope'],max(xs)-res['intercept'],res['r^2']))
        
        return self


    def sort_branches(self,descending=True,sortFxn=None,sortByHeight=True):
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

        mod=-1 if descending else 1
        if sortFxn==None: sortFxn=lambda k: (k.is_node(),-len(k.leaves)*mod,k.length*mod) if k.is_node() else (k.is_node(),k.length*mod)
        if sortByHeight: # Sort nodes by height and group nodes and leaves together
            """ Sort descendants of each node. """

            for k in self.get_internal(): ## iterate over nodes
                k.children=sorted(k.children,key=sortFxn)
        else: # Do not sort by height. Retain leaves at original positions. Only sort nodes
            for k in self.get_internal():
                leavesIdx = [(i,ctr) for ctr, i in enumerate(k.children) if i.is_leaflike()] # Get original indices of leaves
                nodes=sorted([x for x in k.children if x.is_node()],key=sortFxn) # Sort nodes only by number of descendants
                children = nodes
                for i in leavesIdx: # Insert leaves back into same positions
                    children.insert(i[1], i[0])
                k.children = children
        self.assign_tree_coordinates() ## update x and y positions of each branch, since y positions will have changed because of sorting

    def assign_tree_coordinates(self,order=None,widthFxn=None,padNodes=None,verbose=False):
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
        if order==None:
            order=self.traverse_tree(includeCondition=lambda k: k.is_leaflike()) ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
            if verbose==True: print('Drawing tree in pre-order')
        else:
            if verbose==True: print('Drawing tree with provided order')

        nameOrder={x.name: i for i,x in enumerate(order)}
        assert len(nameOrder)==len(order), 'Non-unique names present in tree'
        if widthFxn==None:
            if verbose==True:
                print('Drawing tree with default widths (1 unit for leaf objects, width+1 for clades)')
            skips=[1 if isinstance(x,Leaf) else x.width+1 for x in order]
        else:
            skips=list(map(widthFxn,order))

        for k in self.Objects: ## reset coordinates for all objects
            k.x=None
            k.y=None

        drawn={} ## drawn keeps track of what's been drawn
        for k in order: ## iterate over tips
            x=k.height ## x position is height
            y_idx=nameOrder[k.name] ## assign y index
            y=sum(skips[y_idx:])-skips[y_idx]/2.0 ## sum across skips to find y position

            k.x=x ## set x and y coordinates
            k.y=y
            drawn[k.index]=None ## remember that this objects has been drawn

        if padNodes!=None: ## will be padding nodes
            for n in padNodes: ## iterate over nodes whose descendants will be padded
                idx=sorted([nameOrder[lf] for lf in n.leaves]) if n.is_node() else [order.index(n)] ## indices of all tips to be padded
                for i,k in enumerate(order): ## iterate over all tips

                    if i<idx[0]: ## tip below clade
                        k.y+=padNodes[n] ## pad

                    if (i-1)<idx[-1]: ## tip above clade
                        k.y+=padNodes[n] ## pad again

            allYs=filter(None,self.get_parameter_list('y')) ## get all y positions in tree that aren't None
            minY=min(allYs) ## get min
            for k in self.get_external(): ## reset y positions so tree starts at y=0.5
                k.y-=minY-0.5

        assert len([k for k in self.Objects if k.is_leaflike()])==len(order),'Number of tips in tree does not match number of unique tips, check if two or more collapsed clades were assigned the same name.'
        storePlotted=0

        while len(drawn)!=len(self.Objects): # keep drawing the tree until everything is drawn
            if verbose==True: print('Drawing iteration %d'%(len(drawn)))
            for k in filter(lambda w:w.index not in drawn,self.get_internal()): ## iterate through internal nodes that have not been drawn
                if len([q.y for q in k.children if q.y!=None])==len(k.children): ## all y coordinates of children known
                    if verbose==True: print('Setting node %s coordinates to'%(k.index)),
                    x=k.height ## x position is height
                    childrenYCoords=[q.y for q in k.children if q.y!=None] ## get all existing y coordinates of the node
                    y=sum(childrenYCoords)/float(len(childrenYCoords)) ## internal branch is in the middle of the vertical bar
                    k.x=x
                    k.y=y
                    drawn[k.index]=None ## remember that this objects has been drawn
                    if verbose==True: print('%s (%s branches drawn)'%(k.y,len(drawn)))
                    minYRange=min([min(child.yRange) if child.is_node() else child.y for child in k.children]) ## get lowest y coordinate across children
                    maxYRange=max([max(child.yRange) if child.is_node() else child.y for child in k.children]) ## get highest y coordinate across children
                    setattr(k,'yRange',[minYRange,maxYRange]) ## assign the maximum extent of children's y coordinates

            if len(self.Objects)>len(drawn):
                assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects (%d branches drawn this iteration, %d branches during previous iteration out of %d total)'%(len(drawn),storePlotted,len(self.Objects))
            storePlotted=len(drawn) ## remember how many branches were drawn this iteration

        yValues=[k.y for k in self.Objects] ## all y values
        self.ySpan=max(yValues)-min(yValues)+min(yValues)*2 ## determine appropriate y axis span of tree

        if self.root.is_node():
            self.root.x=min([q.x-q.length for q in self.root.children if q.x!=None]) ## set root x and y coordinates
            childrenYCoords=[q.y for q in self.root.children if q.y!=None]
            self.root.y=sum(childrenYCoords)/float(len(childrenYCoords))
        else:
            self.root.x=self.root.length

    def assign_unrooted_tree_coordinates(self,rotate=0.0,n=None,total=None):
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
        >>> tree.assign_unrooted_tree_coordinates(rotate=0.1)
        
        Docstring generated with ChatGPT 4o.
        """

        if n==None:
            total=sum([1 if x.is_leaf() else x.width+1 for x in self.get_external()])
            n=self.root#.children[0]
            for k in self.Objects:
                k._tau=2*math.pi*rotate
                k.x=0.0
                k.y=0.0

        w=2*math.pi*1.0/float(total) if n.is_leaf() else 2*math.pi*len(n.leaves)/float(total)

        if n.parent.x==None:
            n.parent.x=0.0
            n.parent.y=0.0

        n.x = n.parent.x + n.length * math.cos(n.traits['tau'] + w*0.5)
        n.y = n.parent.y + n.length * math.sin(n.traits['tau'] + w*0.5)
        eta=n._tau

        if n.is_node():
            for ch in n.children:
                w=2*math.pi*1.0/float(total) if ch.is_leaf() else 2*math.pi*len(ch.leaves)/float(total)

                ch._tau = eta
                eta += w
                self.assign_unrooted_tree_coordinates(rotate,ch,total)

    def find_MRCA(self,descendants):
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
        assert len(descendants)>1,'Not enough descendants to find common ancestor: %d'%(len(descendants))
        pathsToRoot={k.index: set(k.get_path_to_root()) for k in descendants} ## for every descendant create an empty set
        # for k in descendants: ## iterate through every descendant
        #     curNode=k ## start descent from descendant
        #     while curNode: ## while not at root
        #         pathsToRoot[k.index].add(curNode) ## remember every node visited along the way
        #         curNode=curNode.parent ## descend

        return sorted(reduce(set.intersection,pathsToRoot.values()),key=lambda k: (-len(k.leaves), k.height))[-1] ## return the most recent branch that is shared across all paths to root

    def collapse_subtree_to_clade(self,cl,givenName,verbose=False,widthFunction=lambda k:len(k.leaves)):
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
        assert cl.is_node(),'Cannot collapse non-node class'
        collapsedClade=Clade(givenName)
        collapsedClade.index=cl.index
        collapsedClade.leaves=cl.leaves
        collapsedClade.length=cl.length
        collapsedClade.height=cl.height
        collapsedClade.parent=cl.parent
        collapsedClade.absoluteTime=cl.absoluteTime
        collapsedClade.traits=cl.traits
        collapsedClade.width=widthFunction(cl)

        if verbose==True: print('Replacing node %s (parent %s) with a clade class'%(cl.index,cl.parent.index))
        parent=cl.parent

        removeFromTree=self.traverse_tree(cl,includeCondition=lambda k: True)
        collapsedClade.subtree=removeFromTree
        assert len(removeFromTree)<len(self.Objects),'Attempted collapse of entire tree'
        collapsedClade.lastHeight=max([x.height for x in removeFromTree])
        if [x.absoluteTime for x in removeFromTree].count(None)!=len(removeFromTree):
            collapsedClade.lastAbsoluteTime=max([x.absoluteTime for x in removeFromTree])

        for k in removeFromTree:
            self.Objects.remove(k)

        parent.children.remove(cl)
        parent.children.append(collapsedClade)
        self.Objects.append(collapsedClade)
        collapsedClade.parent=parent
        if self.tipMap!=None: self.tipMap[givenName]=givenName

        self.traverse_tree()
        self.sort_branches()
        return collapsedClade

    def restore_all_collapsed_subtrees(self):
        """
        Uncollapse all collapsed subtrees in the tree.
        
        This method restores all previously collapsed clades back to their original subtree structures.
        It iterates through all objects in the tree, identifies clades, and replaces each clade with its
        corresponding subtree that was stored in the `clade` class.
        
        Example:
        >>> tree.restore_all_collapsed_subtrees()
        
        Docstring generated with ChatGPT 4o.
        """
        while len([k for k in self.Objects if isinstance(k,Clade)])>0:
            clades=[k for k in self.Objects if isinstance(k,Clade)]
            for cl in clades:
                parent=cl.parent
                subtree=cl.subtree
                parent.children.remove(cl)
                parent.children.append(subtree[0])
                self.Objects+=subtree
                self.Objects.remove(cl)
                if self.tipMap!=None:
                    self.tipMap.pop(cl.name,None)
        self.traverse_tree()

    def collapse_branches(self,collapseIfFxn=lambda x:x.traits['posterior']<=0.5,designatedNodes=[],verbose=False):
        """
        TODO: change desnignatedNode to None
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
        >>> newTree = tree.collapseBranches()
        
        Docstring generated with ChatGPT 4o.
        """
        newTree=copy.deepcopy(self) ## work on a copy of the tree
        if len(designatedNodes)==0: ## no nodes were designated for deletion - relying on anonymous function to collapse nodes
            nodesToDelete=list(filter(lambda n: n.is_node() and collapseIfFxn(n)==True and n!=newTree.root, newTree.Objects)) ## fetch a list of all nodes who are not the root and who satisfy the condition
        else:
            assert len([w for w in designatedNodes if w.is_node()])==len(designatedNodes),'Non-node class detected in list of nodes designated for deletion'
            assert len([w for w in designatedNodes if w!=newTree.root])==0,'Root node was designated for deletion'

            nodesToDelete=list(filter(lambda w: w.index in [q.index for q in designatedNodes], newTree.Objects)) ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
        if verbose==True: print('%s nodes set for collapsing: %s'%(len(nodesToDelete),[w.index for w in nodesToDelete]))
        assert len(nodesToDelete)<len(newTree.get_internal())-1,'Chosen cutoff would remove all branches'
        while len(nodesToDelete)>0: ## as long as there are branches to be collapsed - keep reducing the tree

            if verbose==True: print('Continuing collapse cycle, %s nodes left'%(len(nodesToDelete)))
            for k in sorted(nodesToDelete,key=lambda x:-x.height): ## start with branches near the tips
                zeroNode=k.children ## fetch the node's children
                k.parent.children+=zeroNode ## add them to the zero node's parent
                oldParent=k ## node to be deleted is the old parent
                newParent=k.parent ## once node is deleted, the parent to all their children will be the parent of the deleted node
                if newParent==None:
                    newParent=self.root
                if verbose==True: print('Removing node %s, attaching children %s to node %s'%(oldParent.index,[w.index for w in k.children],newParent.index))
                for w in newTree.Objects: ## assign the parent of deleted node as the parent to any children of deleted node
                    if w.parent==oldParent:
                        w.parent=newParent
                        w.length+=oldParent.length
                        if verbose==True: print('Fixing branch length for node %s'%(w.index))
                k.parent.children.remove(k) ## remove traces of deleted node - it doesn't exist as a child, doesn't exist in the tree and doesn't exist in the nodes list
                newTree.Objects.remove(k)

                nodesToDelete.remove(k) ## in fact, the node never existed

                if len(designatedNodes)==0:
                    nodesToDelete==list(filter(lambda n: n.is_node() and collapseIfFxn(n)==True and n!=newTree.root, newTree.Objects))
                else:
                    assert len([w for w in designatedNodes if w.is_node()])==len(designatedNodes),'Non-node class detected in list of nodes designated for deletion'
                    assert len([w for w in designatedNodes if w!=newTree.root])==0,'Root node was designated for deletion'
                    nodesToDelete=[w for w in newTree.Objects if w.index in [q.index for q in designatedNodes]]

                if verbose==True: print('Removing references to node %s'%(k.index))
        newTree.sort_branches() ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree ## return collapsed tree

    def to_string(self,curNode=None,traits=None,verbose=False,nexus=False,stringFragment=None,traverseCondition=None,rename=None,quoteCharacter="'",json=False):
        """
        Output the topology of the tree with branch lengths and comments to a string.
        
        Parameters:
        curNode (node or None): The starting point for traversal. Default is None, which starts at the root.
        traits (list or None): A list of keys to output entries in the traits dictionary of each branch. Default is all available traits.
        verbose (bool): If True, prints verbose output during the process. Default is False.
        nexus (bool): If True, outputs in NEXUS format. Default is False, which outputs in Newick format.
        stringFragment (list or None): A list of characters that comprise the tree string. Default is None.
        traverseCondition (function or None): A function that determines whether a child branch should be traversed. Default is None which traverses all children.
        rename (dict or None): A dictionary to rename tip names. Default is None.
        quoteCharacter (str): The character to use for quoting tip names. Default is "'".
        json (bool): If True, outputs in auspice JSON format (somewhat experimental). Default is False.
        
        Returns:
        str: The tree string in the specified format.
        
        Example:
        >>> tree_string = tree.toString()
        
        Docstring generated with ChatGPT 4o.
        """
        if curNode==None: curNode=self.root#.children[-1]
        if traits==None: traits=set(sum([list(k.traits.keys()) for k in self.Objects],[])) ## fetch all trait keys
        if stringFragment==None:
            stringFragment=[]
            if nexus:
                assert not json,'Nexus format not a valid option for JSON output'
                if verbose==True: print('Exporting to NEXUS format')
                stringFragment.append('#NEXUS\nBegin trees;\ntree TREE1 = [&R] ')
        if traverseCondition==None: traverseCondition=lambda k: True

        comment=[] ## will hold comment
        if len(traits)>0: ## non-empty list of traits to output
            for tr in traits: ## iterate through keys
                if tr in curNode.traits: ## if key is available
                    if verbose==True: print('trait %s available for %s (%s) type: %s'%(tr,curNode.index,curNode.branchType,type(curNode.traits[tr])))
                    if isinstance(curNode.traits[tr],str): ## string value
                        comment.append('%s="%s"'%(tr,curNode.traits[tr]))
                        if verbose==True: print('adding string comment %s'%(comment[-1]))
                    elif isinstance(curNode.traits[tr],float) or isinstance(curNode.traits[tr],int): ## float or integer
                        comment.append('%s=%s'%(tr,curNode.traits[tr]))
                        if verbose==True: print('adding numeric comment %s'%(comment[-1]))
                    elif isinstance(curNode.traits[tr],list): ## lists
                        rangeComment=[]
                        for val in curNode.traits[tr]:
                            if isinstance(val,str): ## string
                                rangeComment.append('"%s"'%(val))
                            elif isinstance(val,float) or isinstance(val,int): ## float or integer
                                rangeComment.append('%s'%(val))
                            elif isinstance(val, list): ## list of lists, example complete history annotated on tree
                                rangeComment.append("{{{}}}".format(",".join(val)))
                        comment.append('%s={%s}'%(tr,','.join(rangeComment)))
                        if verbose==True: print('adding range comment %s'%(comment[-1]))
                elif verbose==True: print('trait %s unavailable for %s (%s)'%(tr,curNode.index,curNode.branchType))

        if curNode.is_node():
            if verbose==True: print('node: %s'%(curNode.index))
            stringFragment.append('(')
            traverseChildren=list(filter(traverseCondition,curNode.children))
            assert len(traverseChildren)>0,'Node %s does not have traversable children'%(curNode.index)
            for c,child in enumerate(traverseChildren): ## iterate through children of node if they satisfy traverse condition
                if verbose==True: print('moving to child %s of node %s'%(child.index,curNode.index))
                self.toString(curNode=child,traits=traits,verbose=verbose,nexus=nexus,stringFragment=stringFragment,traverseCondition=traverseCondition,rename=rename,quoteCharacter=quoteCharacter)
                if (c+1)<len(traverseChildren): ## not done with children, add comma for next iteration
                    stringFragment.append(',')
            stringFragment.append(')') ## last child, node terminates

        elif curNode.is_leaf():
            if rename==None:
                treeName=curNode.name ## designated numName
            else:
                assert isinstance(rename,dict), 'Variable "rename" is not a dictionary'
                assert curNode.name in rename, 'Tip name %s not in rename dictionary'%(curNode.name)
                treeName=rename[curNode.name]

            if verbose==True: print('leaf: %s (%s)'%(curNode.index,treeName))
            stringFragment.append("%s%s%s"%(quoteCharacter,treeName,quoteCharacter))

        if len(comment)>0:
            if verbose==True: print('adding comment to %s'%(curNode.index))
            comment=','.join(comment)
            comment='[&'+comment+']'
            stringFragment.append('%s'%(comment)) ## end of node, add annotations

        if verbose==True: print('adding branch length to %s'%(curNode.index))
        stringFragment.append(':%8f'%(curNode.length)) ## end of node, add branch length

        if curNode==self.root:#.children[-1]:
            stringFragment.append(';')
            if nexus==True:
                stringFragment.append('\nEnd;')
            if verbose==True: print('finished')
            return ''.join(stringFragment)

    def get_all_tip_TMRCAs(self):
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

        NOTE: This will fail if there are negative branch lengths for some reason. Though I suppose a lot of things will
        """
        assert self.treeType == "time", "Error: can only calculate TMRCA matrix on time-calibrated tree."

        tipNames=[k.name for k in self.get_external()]
        tmrcaMatrix={x:{y:None if x!=y else 0.0 for y in tipNames} for x in tipNames} ## pairwise matrix of tips

        for k in self.get_internal(): ## iterate over nodes
            all_children=list(k.leaves) ## fetch all descendant tips of node

            for a,tipA in enumerate(all_children):
                for tipB in all_children[a+1:]:
                    if tmrcaMatrix[tipA][tipB]==None or tmrcaMatrix[tipA][tipB]<=k.absoluteTime: ## if node's time is more recent than previous entry - set new TMRCA value for pair of tips
                        tmrcaMatrix[tipA][tipB]=k.absoluteTime
                        tmrcaMatrix[tipB][tipA]=k.absoluteTime
        return tmrcaMatrix

    def reduce_tree(self,tipsToKeep,verbose=False):
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
        assert len(tipsToKeep)>0,"No tips given to reduce the tree to."
        assert len([k for k in tipsToKeep if not k.is_leaflike()])==0, "Embedding contains %d branches that are not leaf-like."%(len([k for k in tipsToKeep if k.is_leaflike()==False]))
        if verbose==True: print("Preparing branch hash for keeping %d branches"%(len(tipsToKeep)))
        branchHash={k.index:k for k in tipsToKeep}
        embedding=[]
        if verbose==True: print("Deep copying tree")
        reducedTree=copy.deepcopy(self) ## new tree object
        for k in reducedTree.Objects: ## deep copy branches from current tree
            if k.index in branchHash: ## if branch is designated as one to keep
                currentBranch=k
                if verbose==True: print("Traversing to root from %s"%(currentBranch.index))
                while currentBranch!=reducedTree.root: ## descend to root
                    if verbose==True: print("at %s root: %s"%(currentBranch.index,currentBranch==reducedTree.root))
                    embedding.append(currentBranch) ## keep track of the path to root
                    currentBranch=currentBranch.parent
        embedding.append(reducedTree.root) ## add root to embedding
        if verbose==True: print("Finished extracting embedding with %s branches (%s tips, %s nodes)"%(len(embedding),len([w for w in embedding if w.is_leaf()]),len([w for w in embedding if w.is_node()])))
        embedding=set(embedding) ## prune down to only unique branches

        reducedTree.Objects=sorted(list(embedding),key=lambda x:x.height) ## assign branches that are kept to new tree's Objects
        if verbose==True: print("Pruning untraversed lineages")
        for k in reducedTree.get_internal(): ## iterate through reduced tree
            k.children = [c for c in k.children if c in embedding] ## only keep children that are present in lineage traceback
        reducedTree.root.children=[c for c in reducedTree.root.children if c in embedding] ## do the same for root

        reducedTree.fix_hanging_nodes()

        if verbose==True: print("Last traversal and branch sorting")
        reducedTree.traverse_tree() ## traverse
        reducedTree.sort_branches() ## sort

        return reducedTree ## return new tree

    def count_lineages_at_time(self,t,timeAttr='absoluteTime',inclusionConditionFxn=lambda x:True):
        """
        Count the number of lineages present at a specific time point.
        
        Parameters:
        t (float): The time point at which to count the lineages.
        timeAttr (str): The attribute used to determine the time of the nodes. Default is `absoluteTime`.
        condition (function): A function that determines whether a lineage should be included in the count. Default is a function that always returns True.
        
        Returns:
        int: The number of lineages present at the specified time point (branches whose time is above and parent is below the time point provided).
        
        Example:
        >>> num_lineages = tree.countLineages(2020.5)
        
        Docstring generated with ChatGPT 4o.
        """
        return len([k for k in self.Objects if getattr(k.parent,timeAttr)!=None and getattr(k.parent,timeAttr)<t<=getattr(k,timeAttr) and inclusionConditionFxn(k)])

    def get_external(self,filterFxn=None,onlyLeaves=True):
        """
        Get all leaf-like branches (`Leaf`, `Clade`, and `Reticulation` classes).
        
        Parameters:
        secondFilterFxn (function or None): An optional function to further filter the leaf branches based on an additional property. Default is None.
        
        Returns:
        list: A list of leaf branches that optionally satisfy the secondFilter condition.
        
        Example:
        >>> leaves = tree.get_external()
        >>> filteredLeaves = tree.get_external(lambda x: x.absoluteTime >= 2023.0)
        
        Docstring generated with ChatGPT 4o.
        """
        externals=list(filter(filterFxn,filter(lambda k: k.is_leaf() if onlyLeaves else k.is_leaflike(),self.Objects)))
        return externals

    def get_internal(self,filterFxn=None):
        """
        Get all branches belonging to the `node` class.
        
        Parameters:
        secondFilter (function or None): An optional function to further filter the internal nodes based on an additional property. Default is None.
        
        Returns:
        list: A list of node branches that optionally satisfy the secondFilter condition.
        
        Example:
        >>> nodes = tree.get_internal()
        >>> filteredNodes = tree.get_internal(lambda x: x.absoluteTime >= 2023.0)
        
        Docstring generated with ChatGPT 4o.
        """
        internals=list(filter(filterFxn,filter(lambda k: k.is_node(),self.Objects)))
        return internals

    def get_branches(self,filterFxn=lambda x:True,failIfNoResults=True):
        """
        Get branches that satisfy a specified condition.
        
        Parameters:
        filterFxn (function): A function that determines whether a branch should be included. Default is a function that always returns True.
        failIfNoResults (bool): If True, raises an exception if no branches satisfying the condition are found. Default is True.
        
        Returns:
        list or object: A list of branches that satisfy the condition (list is empty if warn is False). If only one branch satisfies the condition, returns that branch.
        
        Raises:
        Exception: If no branches satisfying the condition are found and warn is True.
        
        Example:
        >>> branches = tree.get_branches(lambda x: x.length > 0.5)
        >>> singleBranch = tree.get_branches(lambda x: x.index == 'node1', warn=False)
        
        Docstring generated with ChatGPT 4o.
        """
        select=list(filter(filterFxn,self.Objects))
        
        if len(select) == 0:
            if failIfNoResults:
                raise Exception('No branches satisfying function were found amongst branches')
            else:
                print("Warning, no results found matching the specified condition. Returning empty list.")
                return []
        elif len(select)==1: #TODO this if block should be removed because the function should always return something of type List, not sometimes a list and sometimes a branch
            #TODO: remove this bit and add a get_branch() method
            return select[-1]
        else:
            return select

    def get_parameter_list(self,statistic,useTraitsDict=False,filterFxn=None):
        """
        Return a list of either branch trait or attribute states across branches.
        
        Parameters:
        statistic (str): The name of the trait or attribute to retrieve.
        use_trait (bool): If True, retrieves the trait from the branch's traits dictionary. If False, retrieves the attribute directly from branch attributes. Default is False (retrieves attributes).
        filterFxn (function or None): A function that determines which branches to include. Default is None, which includes all branches in the tree.
        
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
        if filterFxn==None:
            branches=self.Objects
        else:
            branches=filter(filterFxn,self.Objects)

        if useTraitsDict:
            params=[k.traits[statistic] for k in branches if statistic in k.traits]
        else:
            params=[getattr(k,statistic) for k in branches if hasattr(k,statistic)]

        return params

    def fix_hanging_nodes(self):
        """
        Remove internal nodes without any children. Used in `reduceTree()` and `subtree()` functions internally.
        
        This method iterates over all objects in the tree and removes nodes that have no children. It continues to check
        for and remove hanging nodes until none are left.
        
        Example:
        >>> tree.fix_hanging_nodes()
        
        Docstring generated with ChatGPT 4o.
        """
        while True:
            hangingNodes = [node for node in self.Objects if node.is_node() and not node.children] ## nodes without children (hanging nodes)
            if not hangingNodes:
                break

            for node in hangingNodes:
                node.parent.children.remove(node)
                self.Objects.remove(node)

    def plot_text(self,
                    ax,
                    targetFxn=None,
                    xCoordinateFxn=None,
                    yCoordinateFxn=None,
                    textContentFxn=None,
                    zorder=4,
                    **kwargs):
        """
        Add text annotations to the tree plot.
        
        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the text to.
        target (function): A function to select which branches to annotate. Default selects all `leaf` nodes.
        x_attr (function or None): A function to determine the x-coordinate for the text. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the text. Default is None, which uses the branch's y attribute.
        text (function or None): A function to determine the text content. Default is None, which uses the `leaf` name attribute.
        zorder (int): The z-order for the text. Default is 4.
        **kwargs: Additional keyword arguments to pass to the `ax.text` method.
        
        Returns:
        matplotlib.axes.Axes: The axes with the text annotations added.
        
        Example:
        >>> tree.addText(ax, target=lambda node: node.is_node(), x_attr=lambda node: node.x - 7/365, y_attr=lambda node: node.y - 0.25, text=lambda node: node.traits['posterior'], ha='right', va='top') ## adds posterior values to the left and below internal nodes
        
        Docstring generated with ChatGPT 4o.
        """
        ### Set default values ###
        if targetFxn is None: lambda k: k.is_leaf()
        if xCoordinateFxn==None: xCoordinateFxn=lambda k: k.x
        if yCoordinateFxn==None: yCoordinateFxn=lambda k: k.y
        if textContentFxn==None: textContentFxn=lambda k: k.name

        localKwargs=dict(kwargs)
        if 'verticalalignment' not in localKwargs: localKwargs['verticalalignment']='center'

        for k in filter(targetFxn,self.Objects):
            x,y=xCoordinateFxn(k),yCoordinateFxn(k)
            z=zorder
            ax.text(x,y,textContentFxn(k),zorder=z,**localKwargs)
        return ax

    def plot_text_unrooted(self,
                            ax,
                            targetFxn=None,
                            xCoordinateFxn=None,
                            yCoordinateFxn=None,
                            textContentFxn=None,
                            zorder=4,
                            **kwargs):
        """
        Add text annotations to an unrooted tree plot.
        
        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the text to.
        target (function or None): A function to select which branches to annotate. Default is None, which selects all `leaf` nodes.
        x_attr (function or None): A function to determine the x-coordinate for the text. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the text. Default is None, which uses the branch's y attribute.
        text (function or None): A function to determine the text content. Default is None, which uses the branch's name attribute.
        zorder (int or None): The z-order for the text. Default is None, which sets the z-order to 4.
        **kwargs: Additional keyword arguments to pass to the `ax.text` method.
        
        Returns:
        matplotlib.axes.Axes: The axes with the text annotations added.
        
        Example:
        >>> tree.addTextUnrooted(ax) ## adds tip names to the tree
        
        Docstring generated with ChatGPT 4o.
        """
        if targetFxn==None: targetFxn=lambda k: k.is_leaf()
        if xCoordinateFxn==None: xCoordinateFxn=lambda k: k.x
        if yCoordinateFxn==None: yCoordinateFxn=lambda k: k.y
        if textContentFxn==None: textContentFxn=lambda k: k.name
        
        for k in filter(targetFxn,self.Objects):
            localKwargs=dict(kwargs)
            
            x,y=xCoordinateFxn(k),yCoordinateFxn(k)
            z=zorder
            
            assert 'tau' in k.traits, 'Branch does not have angle tau computed by assign_unrooted_tree_coordinates().'
            
            rot=math.degrees(k.traits['tau'])%360
            
            if 'horizontalalignment' not in localKwargs: localKwargs['horizontalalignment']='right' if 90<rot<270 else 'left'
            if 'verticalalignment' not in localKwargs: localKwargs['verticalalignment']='center'
            
            rot=rot+180 if 90<rot<270 else rot
            
            ax.text(x,y,textContentFxn(k),rotation=rot,rotation_mode='anchor',zorder=z,**localKwargs)
            
        return ax

    def plot_text_annotations_circular(self,
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
                                        **kwargs):
        """
        Add text annotations to a circular tree plot.

        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the text to.
        target (function or None): A function to select which branches to annotate. Default is None, which selects all `leaf` nodes.
        text (function or None): A function to determine the text content. Default is None, which uses the `leaf` name attribute.
        x_attr (function or None): A function to determine the x-coordinate for the text. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the text. Default is None, which uses the branch's y attribute.
        circStart (float): The starting angle (in fractions of 2*pi, i.e. radians) for the circular layout. Default is 0.0.
        circFrac (float): The fraction of the full circle to use for the layout. Default is 1.0.
        inwardSpace (float): Amount of space to leave in the middle of the tree (can be negative for inward-facing trees). Default is 0.0.
        normaliseHeight (function or None): A function to normalize the x-coordinates. Default is None, creates a normalisation that returns 0.0 at root and 1.0 at the most diverged tip.
        zorder (int or None): The z-order for the text. Default is None, which sets the z-order to 4.
        **kwargs: Additional keyword arguments to pass to the `ax.text` method.
        
        Returns:
        matplotlib.axes.Axes: The axes with the text annotations added.
        
        Example:
        >>> tree.addTextCircular(ax) ## adds `leaf` names to a circular tree plot
        
        Docstring generated with ChatGPT 4o.
        """

        if targetFxn==None: targetFxn=lambda k: k.is_leaf()
        if xCoordinateFxn==None: xCoordinateFxn=lambda k:k.x
        if yCoordinateFxn==None: yCoordinateFxn=lambda k:k.y
        if textContentFxn==None: textContentFxn=lambda k: k.name
        if zorder==None: zorder=4
        
        circS=circStart*math.pi*2
        circ=circFrac*math.pi*2
        
        allXs=list(map(xCoordinateFxn,self.Objects))
        if normaliseHeight==None: normaliseHeight=lambda value: (value-min(allXs))/(max(allXs)-min(allXs))
            
        for k in filter(targetFxn,self.Objects): ## iterate over branches
            localKwargs=dict(kwargs) ## copy global kwargs into a local version
            
            x=normaliseHeight(xCoordinateFxn(k)+inwardSpace) ## get branch x position
            y=yCoordinateFxn(k) ## get y position
            
            y=circS+circ*y/self.ySpan
            X=math.sin(y)
            Y=math.cos(y)
            
            rot=math.degrees(y)%360
            
            if 'horizontalalignment' not in localKwargs: localKwargs['horizontalalignment']='right' if 180<rot<360 else 'left' ## rotate labels to aid readability
            if 'verticalalignment' not in localKwargs: localKwargs['verticalalignment']='center'
            rot=360-rot-90 if 180<rot<360 else 360-rot+90
            
            ax.text(X*x,Y*x,textContentFxn(k),rotation=rot,rotation_mode='anchor',zorder=zorder,**localKwargs)
        
        return ax

    def plot_points(self,
                    ax,
                    x_attr=None,
                    y_attr=None,
                    target=None,
                    pointSize=None,
                    pointSizeFxn=None,
                    colour=None,
                    zorder=None,
                    outline=None,
                    outline_size=None,
                    outline_colour=None,
                    **kwargs):
        """
        Plot points on the tree plot.
        
        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to add the points to.
        x_attr (function or None): A function to determine the x-coordinate for the points. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the points. Default is None, which uses the branch's y attribute.
        target (function or None): A function to select which branches to annotate. Default is None, which selects all `leaf` objects.
        size (int or function or None): The size of the points. Default is None, which sets the size to 40.
        pointSizeFxn (function or None): A functioin to specify the size of the points. Default is none.
        colour (str or function or None): The color of the points. Default is None, which sets the color to 'k' (black).
        zorder (int): The z-order for the points. Default is 3.
        outline (bool or None): If True, adds an outline to the points. Default is None, which sets the outline to True.
        outline_size (int or function or None): The size of the outline. Default is None, which sets the outline size to twice the size of the points.
        outline_colour (str or function or None): The color of the outline. Default is None, which sets the outline color to 'k' (black).
        **kwargs: Additional keyword arguments to pass to the `ax.scatter` method.
        
        Returns:
        matplotlib.axes.Axes: The axes with the points added.
        
        Example:
        >>> tree.plotPoints(ax, target=lambda node: node.traits['posterior'] >= 0.95) ## will add circles at nodes with greater than 0.95 posterior support
        
        Docstring generated with ChatGPT 4o.
        """


        if target==None: target=lambda k: k.is_leaf()
        if x_attr==None: x_attr=lambda k:k.x
        if y_attr==None: y_attr=lambda k:k.y
        if size==None: size=40
        if colour==None: colour=lambda f:'k'
        if zorder==None: zorder=3

        if outline==None: outline=True
        if outline_size==None: outline_size=lambda k: size(k)*2 if callable(size) else size*2
        if outline_colour==None: outline_colour='k'

        xs=[]
        ys=[]
        colours=[]
        sizes=[]

        outline_xs=[]
        outline_ys=[]
        outline_colours=[]
        outline_sizes=[]
        for k in filter(target,self.Objects):
            xs.append(x_attr(k))
            ys.append(y_attr(k))
            colours.append(colour(k)) if callable(colour) else colours.append(colour)
            sizes.append(size(k)) if callable(size) else sizes.append(size)

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outline_colour(k)) if callable(outline_colour) else outline_colours.append(outline_colour)
                outline_sizes.append(outline_size(k)) if callable(outline_size) else outline_sizes.append(outline_size)

        ax.scatter(xs,ys,s=sizes,facecolor=colours,edgecolor='none',zorder=zorder,**kwargs) ## put a circle at each tip
        if outline:
            ax.scatter(outline_xs,outline_ys,s=outline_sizes,facecolor=outline_colours,edgecolor='none',zorder=zorder-1,**kwargs) ## put a circle at each tip

        return ax

    def plotTree(self,ax,connection_type=None,target=None,
             x_attr=None,y_attr=None,width=None,
             colour=None,**kwargs):
        """
        Plot the tree on a given matplotlib axes.
        
        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to plot the tree on.
        connection_type (str or None): The type of connection between nodes. Options are 'baltic' (parental branches are plotted as two straight lines - one horizontal, one vertical), 'direct' (diagonal line that directly connects parent and child branches), or 'elbow' (each child has its own angled branch connecting it to the parent). Default is 'baltic'.
        target (function or None): A function to select which branches to plot. Default is None, which selects all branches.
        x_attr (function or None): A function to determine the x-coordinate for the nodes. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the nodes. Default is None, which uses the branch's y attribute.
        width (int or function or None): The width of the lines. Default is None, which sets the width to 2.
        colour (str or function or None): The color of the lines. Default is None, which sets the color to 'k' (black).
        **kwargs: Additional keyword arguments to pass to the LineCollection.
        
        Returns:
        matplotlib.axes.Axes: The axes with the tree plot added.
        
        Example:
        >>> tree.plotTree(ax)
        
        Docstring generated with ChatGPT 4o.
        """
        if target==None: target=lambda k: True
        if x_attr==None: x_attr=lambda k: k.x
        if y_attr==None: y_attr=lambda k: k.y
        if width==None: width=2
        if colour==None: colour='k'
        if connection_type==None: connection_type='baltic'
        assert connection_type in ['baltic','direct','elbow'],'Unrecognised drawing type "%s"'%(tree_type)

        branches=[]
        colours=[]
        linewidths=[]
        for k in filter(target,self.Objects): ## iterate over branches
            x=x_attr(k) ## get branch x position
            xp=x_attr(k.parent) if k.parent else x ## get parent x position
            y=y_attr(k) ## get y position

            try:
                colours.append(colour(k)) if callable(colour) else colours.append(colour)
            except KeyError:
                colours.append((0.7,0.7,0.7)) ## in case no colour available for branch set it to grey
            linewidths.append(width(k)) if callable(width) else linewidths.append(width)

            if connection_type=='baltic': ## each node has a single vertical line to which descendant branches are connected
                branches.append(((xp,y),(x,y)))
                if k.is_node():
                    yl,yr=y_attr(k.children[0]),y_attr(k.children[-1]) ## y positions of first and last child
                    branches.append(((x,yl),(x,yr)))
                    linewidths.append(linewidths[-1])
                    colours.append(colours[-1])
            elif connection_type=='elbow': ## more standard connection where each branch connects to its parent via a right-angled line
                yp=y_attr(k.parent) if k.parent else y ## get parent x position
                branches.append(((xp,yp),(xp,y),(x,y)))
            elif connection_type=='direct': ## this gives triangular looking trees where descendants connect directly to their parents
                yp=y_attr(k.parent) ## get y position
                branches.append(((xp,yp),(x,y)))
            else:
                pass ## for now

        if 'capstyle' not in kwargs: kwargs['capstyle']='projecting'
        line_segments = LineCollection(branches,lw=linewidths,color=colours,**kwargs)
        ax.add_collection(line_segments)
        return ax

    def plotCircularTree(self,ax,target=None,x_attr=None,y_attr=None,width=None,colour=None,
                         circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,precision=15,**kwargs):
        """
        Plot the tree in a circular layout on a given matplotlib axes.
        
        Parameters:
        ax (matplotlib.axes.Axes): The matplotlib axes to plot the tree on.
        target (function or None): A function to select which branches to plot. Default is None, which selects all branches.
        x_attr (function or None): A function to determine the x-coordinate for the nodes. Default is None, which uses the branch's x attribute.
        y_attr (function or None): A function to determine the y-coordinate for the nodes. Default is None, which uses the branch's y attribute.
        width (int or function or None): The width of the lines. Default is None, which sets the width to 2.
        colour (str or function or None): The color of the lines. Default is None, which sets the color to 'k' (black).
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

        if target==None: target=lambda k: True
        if x_attr==None: x_attr=lambda k:k.x
        if y_attr==None: y_attr=lambda k:k.y
        if colour==None: colour='k'
        if width==None: width=2

        if inwardSpace<0: inwardSpace-=self.treeHeight

        branches=[]
        colours=[]
        linewidths=[]

        circ_s=circStart*math.pi*2
        circ=circFrac*math.pi*2

        allXs=list(map(x_attr,self.Objects))
        if normaliseHeight==None: normaliseHeight=lambda value: (value-min(allXs))/(max(allXs)-min(allXs))
        linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

        for k in filter(target,self.Objects): ## iterate over branches
            x=normaliseHeight(x_attr(k)+inwardSpace) ## get branch x position
            xp=normaliseHeight(x_attr(k.parent)+inwardSpace) if k.parent.parent else x ## get parent x position
            y=y_attr(k) ## get y position

            try:
                colours.append(colour(k)) if callable(colour) else colours.append(colour)
            except KeyError:
                colours.append((0.7,0.7,0.7))
            linewidths.append(width(k)) if callable (width) else linewidths.append(width)

            y=circ_s+circ*y/self.ySpan
            X=math.sin(y)
            Y=math.cos(y)
            branches.append(((X*xp,Y*xp),(X*x,Y*x)))

            if k.is_node():
                yl,yr=y_attr(k.children[0]),y_attr(k.children[-1]) ## get leftmost and rightmost children's y coordinates
                yl=circ_s+circ*yl/self.ySpan ## transform y into a fraction of total y
                yr=circ_s+circ*yr/self.ySpan
                ybar=linspace(yl,yr,precision) ## what used to be vertical node bar is now a curved line

                xs=[yx*x for yx in map(math.sin,ybar)] ## convert to polar coordinates
                ys=[yy*x for yy in map(math.cos,ybar)]

                branches+=tuple(zip(zip(xs,ys),zip(xs[1:],ys[1:]))) ## add curved segment

                linewidths+=[linewidths[-1] for q in zip(ys,ys[1:])] ## repeat linewidths
                colours+=[colours[-1] for q in zip(ys,ys[1:])] ## repeat colours

        line_segments = LineCollection(branches,lw=linewidths,ls='-',color=colours,capstyle='projecting',zorder=1) ## create line segments
        ax.add_collection(line_segments) ## add collection to axes
        return ax

    def plotCircularPoints(self,ax,x_attr=None,y_attr=None,target=None,size=None,colour=None,circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,
               zorder=None,outline=None,outline_size=None,outline_colour=None,**kwargs):
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

        if target==None: target=lambda k: k.is_leaf()
        if x_attr==None: x_attr=lambda k:k.x
        if y_attr==None: y_attr=lambda k:k.y
        if size==None: size=40
        if colour==None: colour='k'
        if zorder==None: zorder=3

        if outline==None: outline=True
        if outline_size==None: outline_size=lambda k: size(k)*2 if callable(size) else size*2
        if outline_colour==None: outline_colour=lambda k: 'k'

        if inwardSpace<0: inwardSpace-=self.treeHeight

        circ_s=circStart*math.pi*2
        circ=circFrac*math.pi*2

        allXs=list(map(x_attr,self.Objects))
        if normaliseHeight==None: normaliseHeight=lambda value: (value-min(allXs))/(max(allXs)-min(allXs))
        linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

        xs=[]
        ys=[]
        colours=[]
        sizes=[]

        outline_xs=[]
        outline_ys=[]
        outline_colours=[]
        outline_sizes=[]
        for k in filter(target,self.Objects):
            x=normaliseHeight(x_attr(k)+inwardSpace) ## find normalised x position along circle's radius
            y=circ_s+circ*y_attr(k)/self.ySpan ## get y position along circle's perimeter
            X=math.sin(y)*x ## transform
            Y=math.cos(y)*x ## transform

            xs.append(X)
            ys.append(Y)
            colours.append(colour(k)) if callable(colour) else colours.append(colour)
            sizes.append(size(k)) if callable(size) else sizes.append(size)

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outline_colour(k)) if callable(outline_colour) else outline_colours.append(outline_colour)
                outline_sizes.append(outline_size(k)) if callable(outline_size) else outline_sizes.append(outline_size)

        ax.scatter(xs,ys,s=sizes,facecolor=colours,edgecolor='none',zorder=zorder,**kwargs) ## put a circle at each tip
        if outline:
            ax.scatter(outline_xs,outline_ys,s=outline_sizes,facecolor=outline_colours,edgecolor='none',zorder=zorder-1,**kwargs) ## put a circle at each tip

        return ax
