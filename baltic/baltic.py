from matplotlib.collections import LineCollection
import re,copy,math,json,sys
import datetime as dt
from functools import reduce

__all__ = ['decimalDate', 'convertDate', 'reticulation', # make from baltic import * safe
           'clade', 'node', 'tree',
           'make_tree', 'make_treeJSON', 'loadJSON', 'loadNexus', 'loadNewick']

sys.setrecursionlimit(9001)

def decimalDate(date,fmt="%Y-%m-%d",variable=False):
    """ Converts calendar dates in specified format to decimal date. """
    delimiter=re.search('[^0-9A-Za-z%]',fmt) ## search for non-alphanumeric symbols in fmt (should be field delimiter)
    delimit=None
    if delimiter is not None:
        delimit=delimiter.group()

    if variable==True: ## if date is variable - extract what is available
        if delimit is not None:
            dateL=len(date.split(delimit)) ## split date based on symbol
        else:
            dateL=1 ## no non-alphanumeric characters in date, assume dealing with an imprecise date (something like just year)

        if dateL==2:
            fmt=delimit.join(fmt.split(delimit)[:-1]) ## reduce fmt down to what's available
        elif dateL==1:
            fmt=delimit.join(fmt.split(delimit)[:-2])

    adatetime=dt.datetime.strptime(date,fmt) ## convert to datetime object
    year = adatetime.year ## get year
    boy = dt.datetime(year, 1, 1) ## get beginning of the year
    eoy = dt.datetime(year + 1, 1, 1) ## get beginning of next year
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds())) ## return fractional year

def convertDate(x,start,end):
    """ Converts calendar dates between given formats """
    return dt.datetime.strftime(dt.datetime.strptime(x,start),end)

class reticulation: ## reticulation class (recombination, conversion, reassortment)
    def __init__(self,name):
        self.branchType='leaf'
        self.length=0.0
        self.height=0.0
        self.absoluteTime=None
        self.parent=None
        self.traits={}
        self.index=None
        self.name=name
        # self.numName=name
        self.x=None
        self.y=None
        self.width=0.5
        self.target=None

class clade: ## clade class
    def __init__(self,givenName):
        self.branchType='leaf' ## clade class poses as a leaf
        self.subtree=None ## subtree will contain all the branches that were collapsed
        self.leaves=None
        self.length=0.0
        self.height=None
        self.absoluteTime=None
        self.parent=None
        self.traits={}
        self.index=None
        self.name=givenName ## the pretend tip name for the clade
        # self.numName=givenName
        self.x=None
        self.y=None
        self.lastHeight=None ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime=None ## refers to the absolute time of the highest tip in the collapsed clade
        self.width=1

class node: ## node class
    def __init__(self):
        self.branchType='node'
        self.length=0.0 ## branch length, recovered from string
        self.height=None ## height, set by traversing the tree, which adds up branch lengths along the way
        self.absoluteTime=None ## branch end point in absolute time, once calibrations are done
        self.parent=None ## reference to parent node of the node
        self.children=[] ## a list of descendent branches of this node
        self.traits={} ## dictionary that will contain annotations from the tree string, e.g. {'posterior':1.0}
        self.index=None ## index of the character designating this object in the tree string, it's a unique identifier for every object in the tree
        self.childHeight=None ## the youngest descendant tip of this node
        self.x=None ## X and Y coordinates of this node, once drawTree() is called
        self.y=None
        ## contains references to all tips of this node
        self.leaves=set() ## is a set of tips that are descended from it

class leaf: ## leaf class
    def __init__(self):
        self.branchType='leaf'
        self.name=None ## name of tip after translation, since BEAST trees will generally have numbers for taxa but will provide a map at the beginning of the file
        # self.numName=None ## the original name of the taxon, would be an integer if coming from BEAST, otherwise can be actual name
        self.index=None ## index of the character that defines this object, will be a unique ID for each object in the tree
        self.length=None ## branch length
        self.absoluteTime=None ## position of tip in absolute time
        self.height=None ## height of tip
        self.parent=None ## parent
        self.traits={} ## trait dictionary
        self.x=None ## position of tip on x axis if the tip were to be plotted
        self.y=None ## position of tip on y axis if the tip were to be plotted

class tree: ## tree class
    def __init__(self):
        self.cur_node=node() ## current node is a new instance of a node class
        self.cur_node.index='Root' ## first object in the tree is the root to which the rest gets attached
        self.cur_node.length=0.0 ## startind node branch length is 0
        self.cur_node.height=0.0 ## starting node height is 0
        self.root=None #self.cur_node ## root of the tree is current node
        self.Objects=[] ## tree objects have a flat list of all branches in them
        self.tipMap=None
        self.treeHeight=0 ## tree height is the distance between the root and the most recent tip
        self.mostRecent=None
        self.ySpan=0.0

    def add_reticulation(self,name):
        """ Adds a reticulate branch. """
        ret=reticulation(name)
        ret.index=name
        ret.parent=self.cur_node
        self.cur_node.children.append(ret)
        self.Objects.append(ret)
        self.cur_node=ret

    def add_node(self,i):
        """ Attaches a new node to current node. """
        new_node=node() ## new node instance
        new_node.index=i ## new node's index is the position along the tree string
        if self.root is None:
            self.root=new_node

        new_node.parent=self.cur_node ## new node's parent is current node
        assert self.cur_node.branchType=='node', 'Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses.'
        self.cur_node.children.append(new_node) ## new node is a child of current node
        self.cur_node=new_node ## current node is now new node
        self.Objects.append(self.cur_node) ## add new node to list of objects in the tree

    def add_leaf(self,i,name):
        """ Attach a new leaf (tip) to current node. """
        new_leaf=leaf() ## new instance of leaf object
        new_leaf.index=i ## index is position along tree string
        if self.root is None: self.root=new_leaf

        new_leaf.parent=self.cur_node ## leaf's parent is current node
        self.cur_node.children.append(new_leaf) ## assign leaf to parent's children
        # new_leaf.numName=name ## numName is the name tip has inside tree string, BEAST trees usually have numbers for tip names
        new_leaf.name=name
        self.cur_node=new_leaf ## current node is now new leaf
        self.Objects.append(self.cur_node) ## add leaf to all objects in the tree

    def subtree(self,k=None,traverse_condition=None,stem=True):
        """ Generate a subtree (as a baltic tree object) from a traversal.
            k is the starting branch for traversal (default: root).
            traverse_condition is a function that determines whether a child branch should be visited (default: always true).
            Returns a new baltic tree instance.
            Note - custom traversal functions can result in multitype trees.
            If this is undesired call singleType() on the resulting subtree afterwards. """
        subtree=copy.deepcopy(self.traverse_tree(k,include_condition=lambda k:True,traverse_condition=traverse_condition))

        if subtree is None or len([k for k in subtree if k.branchType=='leaf'])==0:
            return None
        else:
            local_tree=tree() ## create a new tree object where the subtree will be
            local_tree.Objects=subtree ## assign branches to new tree object

            local_tree.root=subtree[0] ## root is the beginning of the traversal

            if stem==True: ## we want the stem
                local_tree.root.parent=copy.deepcopy(k.parent) ## means assigning an invisible parent to root
                local_tree.root.parent.height=0.0 ## set height to 0.0 so heights can be set
                if local_tree.root.parent.parent: local_tree.root.parent.parent=None ## remove reference to the invisible parent's parent
            else: ## don't want stem
                local_tree.root.parent=None ## tree begins strictly at node

            subtree_set=set(subtree) ## turn branches into set for quicker look up later

            if traverse_condition is not None: ## didn't use default traverse condition, might need to deal with hanging nodes and prune children
                for nd in local_tree.getInternal(): ## iterate over nodes
                    nd.children=list(filter(lambda k:k in subtree_set,nd.children)) ## only keep children seen in traversal
                local_tree.fixHangingNodes()

            if self.tipMap: ## if original tree has a tipMap dictionary
                # local_tree.tipMap={tipNum: self.tipMap[tipNum] for tipNum in filter(lambda w: w in {r.numName:None for r in local_tree.getExternal()},self.tipMap)} ## copy over the relevant tip translations
                local_tree.tipMap={tipNum: self.tipMap[tipNum] for tipNum in self.tipMap if self.tipMap[tipNum] in [w.name for w in local_tree.getExternal()]} ## copy over the relevant tip translations

            return local_tree

    def singleType(self):
        """ Removes any branches with a single child (multitype nodes). """
        multiTypeNodes=[k for k in self.Objects if k.branchType=='node' and len(k.children)==1]
        while len(multiTypeNodes)>0:
            multiTypeNodes=[k for k in self.Objects if k.branchType=='node' and len(k.children)==1]

            if self.root in multiTypeNodes:
                multiTypeNodes.remove(self.root)
                new_root=self.root.children[0]
                new_root.parent=None
                self.root=new_root

            for k in sorted(multiTypeNodes,key=lambda x:-x.height):
                child=k.children[0] ## fetch child
                grandparent=k.parent ## fetch grandparent

                child.parent=grandparent ## child's parent is now grandparent

                grandparent.children.append(child) ## add child to grandparent's children
                grandparent.children.remove(k) ## remove old parent from grandparent's children
                grandparent.children=list(set(grandparent.children))
                child.length+=k.length ## adjust child length

                multiTypeNodes.remove(k) ## remove old parent from multitype nodes
                self.Objects.remove(k) ## remove old parent from all objects
        self.sortBranches()

    def setAbsoluteTime(self,date):
        """ place all objects in absolute time by providing the date of the most recent tip """
        for i in self.Objects: ## iterate over all objects
            i.absoluteTime=date-self.treeHeight+i.height ## heights are in units of time from the root
        self.mostRecent=max([k.absoluteTime for k in self.Objects])

    def treeStats(self):
        """ provide information about the tree """
        self.traverse_tree() ## traverse the tree
        obs=self.Objects ## convenient list of all objects in the tree
        print('\nTree height: %.6f\nTree length: %.6f'%(self.treeHeight,sum([x.length for x in obs]))) ## report the height and length of tree

        nodes=self.getInternal() ## get all nodes
        strictlyBifurcating=False ## assume tree is not strictly bifurcating
        multiType=False
        singleton=False

        N_children=[len(x.children) for x in nodes]
        if len(N_children)==0:
            singleton=True
        else:
            minChildren,maxChildren=min(N_children),max(N_children) ## get the largest number of descendant branches of any node
            if maxChildren==2 and minChildren==2: ## if every node has at most two children branches
                strictlyBifurcating=True ## it's strictly bifurcating
            if minChildren==1:
                multiType=True

        hasTraits=False ## assume tree has no annotations
        maxAnnotations=max([len(x.traits) for x in obs]) ## check the largest number of annotations any branch has
        if maxAnnotations>0: ## if it's more than 0
            hasTraits=True ## there are annotations

        if strictlyBifurcating: print('strictly bifurcating tree') ## report
        if multiType: print('multitype tree') ## report
        if singleton: print('singleton tree')
        if hasTraits: print('annotations present') ## report

        print('\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n'%(len(obs),len(nodes),len(self.getExternal()))) ## report numbers of different objects in the tree

    def traverse_tree(self,cur_node=None,include_condition=lambda k:k.branchType=='leaf',traverse_condition=lambda k:True,collect=None,verbose=False):
        if cur_node==None: ## if no starting point defined - start from root
            for k in self.Objects: ## reset various parameters
                if k.branchType=='node':
                    k.leaves=set()
                    k.childHeight=None
                k.height=None

            if verbose==True: print('Initiated traversal from root')
            cur_node=self.root#.children[-1]

        if collect==None: ## initiate collect list if not initiated
            collect=[]

        if cur_node.parent and cur_node.height==None: ## cur_node has a parent - set height if it doesn't have it already
            cur_node.height=cur_node.length+cur_node.parent.height
        elif cur_node.height==None: ## cur_node does not have a parent (root), if height not set before it's zero
            cur_node.height=0.0

        if verbose==True: print('at %s (%s)'%(cur_node.index,cur_node.branchType))

        if include_condition(cur_node): ## test if interested in cur_node
            collect.append(cur_node) ## add to collect list for reporting later

        if cur_node.branchType=='leaf' and self.root!=cur_node: ## cur_node is a tip (and tree is not single tip)
            # cur_node.parent.leaves.add(cur_node.numName) ## add to parent's list of tips
            cur_node.parent.leaves.add(cur_node.name) ## add to parent's list of tips

        elif cur_node.branchType=='node': ## cur_node is node
            for child in filter(traverse_condition,cur_node.children): ## only traverse through children we're interested
                if verbose==True: print('visiting child %s'%(child.index))
                self.traverse_tree(cur_node=child,include_condition=include_condition,traverse_condition=traverse_condition,verbose=verbose,collect=collect) ## recurse through children
                if verbose==True: print('child %s done'%(child.index))
            assert len(cur_node.children)>0, 'Tried traversing through hanging node without children. Index: %s'%(cur_node.index)
            if verbose==True: print(cur_node.index,[child.height for child in cur_node.children])
            cur_node.childHeight=max([child.childHeight if child.branchType=='node' else child.height for child in cur_node.children])

            if cur_node.parent:
                cur_node.parent.leaves=cur_node.parent.leaves.union(cur_node.leaves) ## pass tips seen during traversal to parent
            self.treeHeight=cur_node.childHeight ## it's the highest child of the starting node
        return collect

    def renameTips(self,d=None):
        """ Give each tip its correct label using a dictionary. """
        if d==None and self.tipMap!=None:
            d=self.tipMap
        for k in self.getExternal(): ## iterate through leaf objects in tree
            # k.name=d[k.numName] ## change its name
            k.name=d[k.name] ## change its name

    def sortBranches(self,descending=True):
        """ Sort descendants of each node. """
        if descending==True:
            modifier=-1 ## define the modifier for sorting function later
        elif descending==False:
            modifier=1

        for k in self.getInternal(): ## iterate over nodes
            ## split node's offspring into nodes and leaves, sort each list individually
            nodes=sorted([x for x in k.children if x.branchType=='node'],key=lambda q:(-len(q.leaves)*modifier,q.length*modifier))
            leaves=sorted([x for x in k.children if x.branchType=='leaf'],key=lambda q:q.length*modifier)

            if modifier==1: ## if sorting one way - nodes come first, leaves later
                k.children=nodes+leaves
            elif modifier==-1: ## otherwise sort the other way
                k.children=leaves+nodes
        self.drawTree() ## update x and y positions of each branch, since y positions will have changed because of sorting

    def drawTree(self,order=None,width_function=None,verbose=False):
        """ Find x and y coordinates of each branch. """
        if order==None:
            order=self.traverse_tree() ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
            if verbose==True: print('Drawing tree in pre-order')
        else:
            if verbose==True: print('Drawing tree with provided order')

        # name_order={x.numName: i for i,x in enumerate(order)}
        name_order={x.name: i for i,x in enumerate(order)}
        if width_function==None:
            if verbose==True:
                print('Drawing tree with default widths (1 unit for leaf objects, width+1 for clades)')
            skips=[1 if isinstance(x,leaf) else x.width+1 for x in order]
        else:
            skips=list(map(width_function,order))

        for k in self.Objects: ## reset coordinates for all objects
            k.x=None
            k.y=None

        assert len(self.getExternal())==len(order),'Number of tips in tree does not match number of unique tips, check if two or more collapsed clades were assigned the same name.'
        storePlotted=0
        drawn={} ## drawn keeps track of what's been drawn
        while len(drawn)!=len(self.Objects): # keep drawing the tree until everything is drawn
            if verbose==True: print('Drawing iteration %d'%(len(drawn)))
            for k in filter(lambda w:w.index not in drawn,self.getExternal()+self.getInternal()): ## iterate through objects that have not been drawn
                if k.branchType=='leaf': ## if leaf - get position of leaf, draw branch connecting tip to parent node
                    if verbose==True: print('Setting leaf %s y coordinate to'%(k.index)),
                    x=k.height ## x position is height
                    # y_idx=name_order[k.numName] ## y position of leaf is given by the order in which tips were visited during the traversal
                    y_idx=name_order[k.name]
                    y=sum(skips[y_idx:])-skips[y_idx]/2.0 ## sum across skips to find y position
                    k.x=x ## set x and y coordinates
                    k.y=y
                    drawn[k.index]=None ## remember that this objects has been drawn
                    if verbose==True: print('%s (%s branches drawn)'%(k.y,len(drawn)))
                    if k.parent and hasattr(k.parent,'yRange')==False: ## if parent doesn't have a maximum extent of its children's y coordinates
                        setattr(k.parent,'yRange',[k.y,k.y]) ## assign it

                if k.branchType=='node': ## if parent is non-root node and y positions of all its children are known
                    if len([q.y for q in k.children if q.y!=None])==len(k.children):
                        if verbose==True: print('Setting node %s coordinates to'%(k.index)),
                        x=k.height ## x position is height
                        children_y_coords=[q.y for q in k.children if q.y!=None] ## get all existing y coordinates of the node
                        y=sum(children_y_coords)/float(len(children_y_coords)) ## internal branch is in the middle of the vertical bar
                        k.x=x
                        k.y=y
                        drawn[k.index]=None ## remember that this objects has been drawn
                        if verbose==True: print('%s (%s branches drawn)'%(k.y,len(drawn)))
                        minYrange=min([min(child.yRange) if child.branchType=='node' else child.y for child in k.children]) ## get lowest y coordinate across children
                        maxYrange=max([max(child.yRange) if child.branchType=='node' else child.y for child in k.children]) ## get highest y coordinate across children
                        setattr(k,'yRange',[minYrange,maxYrange]) ## assign the maximum extent of children's y coordinates

            if len(self.Objects)>len(drawn):
                assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects (%d branches drawn this iteration, %d branches during previous iteration out of %d total)'%(len(drawn),storePlotted,len(self.Objects))
            storePlotted=len(drawn)
            self.ySpan=sum(skips)

        if self.root.branchType=='node':
            self.root.x=min([q.x-q.length for q in self.root.children if q.x!=None]) ## set root x and y coordinates
            children_y_coords=[q.y for q in self.root.children if q.y!=None]
            self.root.y=sum(children_y_coords)/float(len(children_y_coords))
        else:
            self.root.x=self.root.length

    def drawUnrooted(self,n=None,total=None):
        """
        Calculate x and y coordinates in an unrooted arrangement.
        Code translated from https://github.com/nextstrain/auspice/commit/fc50bbf5e1d09908be2209450c6c3264f298e98c, written by Richard Neher.
        """
        if n==None:
            total=sum([1 if isinstance(x,leaf) else x.width+1 for x in self.getExternal()])
            n=self.root#.children[0]
            for k in self.Objects:
                k.traits['tau']=0.0
                k.x=0.0
                k.y=0.0

        if n.branchType=='leaf':
            w=2*math.pi*1.0/float(total)
        else:
            w=2*math.pi*len(n.leaves)/float(total)

        if n.parent.x==None:
            n.parent.x=0.0
            n.parent.y=0.0

        n.x = n.parent.x + n.length * math.cos(n.traits['tau'] + w*0.5)
        n.y = n.parent.y + n.length * math.sin(n.traits['tau'] + w*0.5)
        eta=n.traits['tau']

        if n.branchType=='node':
            for ch in n.children:
                if ch.branchType=='leaf':
                    w=2*math.pi*1.0/float(total)
                else:
                    w=2*math.pi*len(ch.leaves)/float(total)
                ch.traits['tau'] = eta
                eta += w
                self.drawUnrooted(ch,total)

    def commonAncestor(self,descendants,strict=False):
        """
        Find the most recent node object that gave rise to a given list of descendant branches.
        """
        assert len(descendants)>1,'Not enough descendants to find common ancestor: %d'%(len(descendants))
        paths_to_root={k.index: set() for k in descendants} ## for every descendant create an empty set
        for k in descendants: ## iterate through every descendant
            cur_node=k ## start descent from descendant
            while cur_node: ## while not at root
                paths_to_root[k.index].add(cur_node) ## remember every node visited along the way
                cur_node=cur_node.parent ## descend

        return sorted(reduce(set.intersection,paths_to_root.values()),key=lambda k: k.height)[-1] ## return the most recent branch that is shared across all paths to root

    def collapseSubtree(self,cl,givenName,verbose=False,widthFunction=lambda k:len(k.leaves)):
        """ Collapse an entire subtree into a clade object. """
        assert cl.branchType=='node','Cannot collapse non-node class'
        collapsedClade=clade(givenName)
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

        remove_from_tree=self.traverse_tree(cl,include_condition=lambda k: True)
        collapsedClade.subtree=remove_from_tree
        assert len(remove_from_tree)<len(self.Objects),'Attempted collapse of entire tree'
        collapsedClade.lastHeight=max([x.height for x in remove_from_tree])
        if [x.absoluteTime for x in remove_from_tree].count(None)!=len(remove_from_tree):
            collapsedClade.lastAbsoluteTime=max([x.absoluteTime for x in remove_from_tree])

        for k in remove_from_tree:
            self.Objects.remove(k)

        parent.children.remove(cl)
        parent.children.append(collapsedClade)
        self.Objects.append(collapsedClade)
        collapsedClade.parent=parent
        if self.tipMap!=None: self.tipMap[givenName]=givenName

        self.traverse_tree()
        self.sortBranches()
        return collapsedClade

    def uncollapseSubtree(self):
        """ Uncollapse all collapsed subtrees. """
        while len([k for k in self.Objects if isinstance(k,clade)])>0:
            clades=[k for k in self.Objects if isinstance(k,clade)]
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

    def collapseBranches(self,collapseIf=lambda x:x.traits['posterior']<=0.5,designated_nodes=[],verbose=False):
        """ Collapse all branches that satisfy a function collapseIf (default is an anonymous function that returns true if posterior probability is <=0.5).
            Alternatively, a list of nodes can be supplied to the script.
            Returns a deep copied version of the tree.
        """
        newTree=copy.deepcopy(self) ## work on a copy of the tree
        if len(designated_nodes)==0: ## no nodes were designated for deletion - relying on anonymous function to collapse nodes
            nodes_to_delete=list(filter(lambda n: n.branchType=='node' and collapseIf(n)==True and n!=newTree.root, newTree.Objects)) ## fetch a list of all nodes who are not the root and who satisfy the condition
        else:
            assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
            assert len([w for w in designated_nodes if w!=newTree.root])==0,'Root node was designated for deletion'

            nodes_to_delete=list(filter(lambda w: w.index in [q.index for q in designated_nodes], newTree.Objects)) ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
        if verbose==True: print('%s nodes set for collapsing: %s'%(len(nodes_to_delete),[w.index for w in nodes_to_delete]))
        # assert len(nodes_to_delete)<len(newTree.getInternal())-1,'Chosen cutoff would remove all branches'
        while len(nodes_to_delete)>0: ## as long as there are branches to be collapsed - keep reducing the tree

            if verbose==True: print('Continuing collapse cycle, %s nodes left'%(len(nodes_to_delete)))
            for k in sorted(nodes_to_delete,key=lambda x:-x.height): ## start with branches near the tips
                zero_node=k.children ## fetch the node's children
                k.parent.children+=zero_node ## add them to the zero node's parent
                old_parent=k ## node to be deleted is the old parent
                new_parent=k.parent ## once node is deleted, the parent to all their children will be the parent of the deleted node
                if new_parent==None:
                    new_parent=self.root
                if verbose==True: print('Removing node %s, attaching children %s to node %s'%(old_parent.index,[w.index for w in k.children],new_parent.index))
                for w in newTree.Objects: ## assign the parent of deleted node as the parent to any children of deleted node
                    if w.parent==old_parent:
                        w.parent=new_parent
                        w.length+=old_parent.length
                        if verbose==True: print('Fixing branch length for node %s'%(w.index))
                k.parent.children.remove(k) ## remove traces of deleted node - it doesn't exist as a child, doesn't exist in the tree and doesn't exist in the nodes list
                newTree.Objects.remove(k)

                nodes_to_delete.remove(k) ## in fact, the node never existed

                if len(designated_nodes)==0:
                    nodes_to_delete==list(filter(lambda n: n.branchType=='node' and collapseIf(n)==True and n!=newTree.root, newTree.Objects))
                else:
                    assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
                    assert len([w for w in designated_nodes if w!=newTree.root])==0,'Root node was designated for deletion'
                    nodes_to_delete=[w for w in newTree.Objects if w.index in [q.index for q in designated_nodes]]

                if verbose==True: print('Removing references to node %s'%(k.index))
        newTree.sortBranches() ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree ## return collapsed tree

    def toString(self,cur_node=None,traits=None,verbose=False,nexus=False,string_fragment=None,traverse_condition=None,rename=None,quotechar="'",json=False):
        """ Output the topology of the tree with branch lengths and comments to stringself.
            cur_node: starting point (default: None, starts at root)
            traits: list of keys that will be used to output entries in traits dict of each branch (default: all traits)
            numName: boolean, whether encoded (True) or decoded (default: False) tip names will be output
            verbose: boolean, debug
            nexus: boolean, whether to output newick (default: False) or nexus (True) formatted tree
            string_fragment: list of characters that comprise the tree string
        """
        if cur_node==None: cur_node=self.root#.children[-1]
        if traits==None: traits=set(sum([list(k.traits.keys()) for k in self.Objects],[])) ## fetch all trait keys
        if string_fragment==None:
            string_fragment=[]
            if nexus==True:
                assert json==False,'Nexus format not a valid option for JSON output'
                if verbose==True: print('Exporting to Nexus format')
                string_fragment.append('#NEXUS\nBegin trees;\ntree TREE1 = [&R] ')
        if traverse_condition==None: traverse_condition=lambda k: True

        comment=[] ## will hold comment
        if len(traits)>0: ## non-empty list of traits to output
            for tr in traits: ## iterate through keys
                if tr in cur_node.traits: ## if key is available
                    if verbose==True: print('trait %s available for %s (%s) type: %s'%(tr,cur_node.index,cur_node.branchType,type(cur_node.traits[tr])))
                    if isinstance(cur_node.traits[tr],str): ## string value
                        comment.append('%s="%s"'%(tr,cur_node.traits[tr]))
                        if verbose==True: print('adding string comment %s'%(comment[-1]))
                    elif isinstance(cur_node.traits[tr],float) or isinstance(cur_node.traits[tr],int): ## float or integer
                        comment.append('%s=%s'%(tr,cur_node.traits[tr]))
                        if verbose==True: print('adding numeric comment %s'%(comment[-1]))
                    elif isinstance(cur_node.traits[tr],list): ## lists
                        rangeComment=[]
                        for val in cur_node.traits[tr]:
                            if isinstance(val,str): ## string
                                rangeComment.append('"%s"'%(val))
                            elif isinstance(val,float) or isinstance(val,int): ## float or integer
                                rangeComment.append('%s'%(val))
                        comment.append('%s={%s}'%(tr,','.join(rangeComment)))
                        if verbose==True: print('adding range comment %s'%(comment[-1]))
                elif verbose==True: print('trait %s unavailable for %s (%s)'%(tr,cur_node.index,cur_node.branchType))

        if cur_node.branchType=='node':
            if verbose==True: print('node: %s'%(cur_node.index))
            string_fragment.append('(')
            traverseChildren=list(filter(traverse_condition,cur_node.children))
            assert len(traverseChildren)>0,'Node %s does not have traversable children'%(cur_node.index)
            for c,child in enumerate(traverseChildren): ## iterate through children of node if they satisfy traverse condition
                if verbose==True: print('moving to child %s of node %s'%(child.index,cur_node.index))
                self.toString(cur_node=child,traits=traits,verbose=verbose,nexus=nexus,string_fragment=string_fragment,traverse_condition=traverse_condition,rename=rename,quotechar=quotechar)
                if (c+1)<len(traverseChildren): ## not done with children, add comma for next iteration
                    string_fragment.append(',')
            string_fragment.append(')') ## last child, node terminates

        elif cur_node.branchType=='leaf':
            # if numName==False: ## if real names wanted
            #     assert cur_node.name!=None,'Tip does not have converted name' ## assert they have been converted
            #     treeName=cur_node.name ## designate real name
            # elif numName==True: ## if number names wanted
                # treeName=cur_node.numName ## designated numName
            if rename==None:
                treeName=cur_node.name ## designated numName
            else:
                assert isinstance(rename,dict), 'Variable "rename" is not a dictionary'
                assert cur_node.name in rename, 'Tip name %s not in rename dictionary'%(cur_node.name)
                treeName=rename[cur_node.name]

            if verbose==True: print('leaf: %s (%s)'%(cur_node.index,treeName))
            string_fragment.append("%s%s%s"%(quotechar,treeName,quotechar))

        if len(comment)>0:
            if verbose==True: print('adding comment to %s'%(cur_node.index))
            comment=','.join(comment)
            comment='[&'+comment+']'
            string_fragment.append('%s'%(comment)) ## end of node, add annotations

        if verbose==True: print('adding branch length to %s'%(cur_node.index))
        string_fragment.append(':%8f'%(cur_node.length)) ## end of node, add branch length

        if cur_node==self.root:#.children[-1]:
            string_fragment.append(';')
            if nexus==True:
                string_fragment.append('\nEnd;')
            if verbose==True: print('finished')
            return ''.join(string_fragment)

    def allTMRCAs(self,numName=True):
        if numName==False:
            assert len(self.tipMap)>0,'Tree does not have a translation dict for tip names'
            # tip_names=[self.tipMap[k.numName] for k in self.Objects if isinstance(k,leaf)]
            tip_names=[k.name for k in self.Objects if isinstance(k,leaf)]
        else:
            # tip_names=[k.numName for k in self.Objects if isinstance(k,leaf)]
            tip_names=[k.name for k in self.Objects if isinstance(k,leaf)]
        tmrcaMatrix={x:{y:None if x!=y else 0.0 for y in tip_names} for x in tip_names} ## pairwise matrix of tips

        for k in self.getInternal():
            if numName==True:
                all_children=list(k.leaves) ## fetch all descendant tips of node
            else:
                all_children=[self.tipMap[lv] for lv in k.leaves]
            for x in range(0,len(all_children)-1): ## for all pairwise comparisons of tips
                for y in range(x+1,len(all_children)):
                    tipA=all_children[x]
                    tipB=all_children[y]
                    if tmrcaMatrix[tipA][tipB]==None or tmrcaMatrix[tipA][tipB]<=k.absoluteTime: ## if node's time is more recent than previous entry - set new TMRCA value for pair of tips
                        tmrcaMatrix[tipA][tipB]=k.absoluteTime
                        tmrcaMatrix[tipB][tipA]=k.absoluteTime
        return tmrcaMatrix

    def reduceTree(self,keep,verbose=False):
        """
        Reduce the tree to just those tracking a small number of tips.
        Returns a new baltic tree object.
        """
        assert len(keep)>0,"No tips given to reduce the tree to."
        assert len([k for k in keep if k.branchType!='leaf'])==0, "Embedding contains %d non-leaf branches."%(len([k for k in keep if k.branchType!='leaf']))
        if verbose==True: print("Preparing branch hash for keeping %d branches"%(len(keep)))
        branch_hash={k.index:k for k in keep}
        embedding=[]
        if verbose==True: print("Deep copying tree")
        reduced_tree=copy.deepcopy(self) ## new tree object
        for k in reduced_tree.Objects: ## deep copy branches from current tree
            if k.index in branch_hash: ## if branch is designated as one to keep
                cur_b=k
                if verbose==True: print("Traversing to root from %s"%(cur_b.index))
                while cur_b!=reduced_tree.root: ## descend to root
                    if verbose==True: print("at %s root: %s"%(cur_b.index,cur_b==reduced_tree.root))
                    embedding.append(cur_b) ## keep track of the path to root
                    cur_b=cur_b.parent
        embedding.append(reduced_tree.root) ## add root to embedding
        if verbose==True: print("Finished extracting embedding with %s branches (%s tips, %s nodes)"%(len(embedding),len([w for w in embedding if w.branchType=='leaf']),len([w for w in embedding if w.branchType=='node'])))
        embedding=set(embedding) ## prune down to only unique branches

        reduced_tree.Objects=sorted(list(embedding),key=lambda x:x.height) ## assign branches that are kept to new tree's Objects
        if verbose==True: print("Pruning untraversed lineages")
        for k in reduced_tree.getInternal(): ## iterate through reduced tree
            k.children = [c for c in k.children if c in embedding] ## only keep children that are present in lineage traceback
        reduced_tree.root.children=[c for c in reduced_tree.root.children if c in embedding] ## do the same for root

        reduced_tree.fixHangingNodes()

        if verbose==True: print("Last traversal and branch sorting")
        reduced_tree.traverse_tree() ## traverse
        reduced_tree.sortBranches() ## sort

        return reduced_tree ## return new tree

    def countLineages(self,t,condition=lambda x:True):
        return len([k for k in self.Objects if k.parent.absoluteTime<t<=k.absoluteTime and condition(k)])

    def getExternal(self,secondFilter=None):
        """
        Get all branches whose branchType is "leaf".
        A function can be provided to filter internal nodes according to an additional property.
        """
        externals=list(self.Objects)
        # if self.root.branchType=='leaf':
        #     externals.append(self.root)
        return list(filter(secondFilter,filter(lambda k:k.branchType=='leaf',externals)))

    def getInternal(self,secondFilter=None):
        """
        Get all branches whose branchType is "node".
        A function can be provided to filter internal nodes according to an additional property.
        """
        internals=list(self.Objects)
        # if self.root.branchType=='node':
        #     internals.append(self.root)
        return list(filter(secondFilter,filter(lambda k:k.branchType=='node',internals)))

    def getBranches(self,attrs=lambda x:True,warn=True):
        select=list(filter(attrs,self.Objects))

        if len(select)==0 and warn==True:
            raise Exception('No branches satisfying function were found amongst branches')
        elif len(select)==0 and warn==False:
            return []
        elif len(select)==1:
            return select[-1]
        else:
            return select

    def getParameter(self,statistic,use_trait=True,which_branches=None):
        """
        Return either branch trait or attribute (default: trait, to switch to attribute set use_trait parameter to False) statistic across branches determined by the which_branches function (default: all objects in the tree).
        Note - branches which do not have the trait or attribute are skipped.
        """
        if which_branches==None:
            branches=self.Objects#+[self.root]
        else:
            branches=filter(which_branches,self.Objects)#+[self.root])

        if use_trait==False:
            params=[getattr(k,statistic) for k in branches if hasattr(k,statistic)]
        elif use_trait==True:
            params=[k.traits[statistic] for k in branches if statistic in k.traits]
        return params

    def fixHangingNodes(self):
        """
        Remove internal nodes without any children.
        """
        hangingCondition=lambda k:k.branchType=='node' and len(k.children)==0
        hangingNodes=list(filter(hangingCondition,self.Objects)) ## check for nodes without any children (hanging nodes)
        while len(hangingNodes)>0:
            for h in sorted(hangingNodes,key=lambda x:-x.height):
                h.parent.children.remove(h) ## remove old parent from grandparent's children
                hangingNodes.remove(h) ## remove old parent from multitype nodes
                self.Objects.remove(h) ## remove old parent from all objects
            hangingNodes=list(filter(hangingCondition,self.Objects)) ## regenerate list

    def addText(self,ax,target=lambda k:k.branchType=='leaf',position=lambda k:(k.x*1.01,k.y),text=lambda k:k.name,zorder_function=lambda k: 101,**kwargs):
        for k in filter(target,self.Objects):
            x,y=position(k)
            z=zorder_function(k)
            ax.text(x,y,text(k),zorder=z,**kwargs)
        return ax

    def plotPoints(self,ax,x_attr=None,y_attr=None,target=None,size_function=None,colour_function=None,
               zorder=None,outline=None,outline_size=None,outline_colour=None,**kwargs):
        if target==None: target=lambda k: k.branchType=='leaf'
        if x_attr==None: x_attr=lambda k:k.x
        if y_attr==None: y_attr=lambda k:k.y
        if size_function==None: size_function=lambda k:40
        if colour_function==None: colour_function=lambda f:'k'
        if zorder==None: zorder=3

        if outline==None: outline=True
        if outline_size==None: outline_size=lambda k: size_function(k)*2
        if outline_colour==None: outline_colour=lambda k: 'k'

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
            colours.append(colour_function(k))
            sizes.append(size_function(k))

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outline_colour(k))
                outline_sizes.append(outline_size(k))

        ax.scatter(xs,ys,s=sizes,facecolor=colours,edgecolor='none',zorder=zorder,**kwargs) ## put a circle at each tip
        if outline:
            ax.scatter(outline_xs,outline_ys,s=outline_sizes,facecolor=outline_colours,edgecolor='none',zorder=zorder-1,**kwargs) ## put a circle at each tip

        return ax


    def plotTree(self,ax,tree_type='rectangular',target=lambda k: True,
             x_attr=lambda k:k.x,y_attr=lambda k:k.y,branchWidth=lambda k:2,
             colour_function=lambda f:'k',**kwargs):

        assert tree_type in ['rectangular','unrooted','non-baltic'],'Unrecognised drawing type "%s"'%(tree_type)

        branches=[]
        colours=[]
        linewidths=[]
        for k in filter(target,self.Objects): ## iterate over branches
            x=x_attr(k) ## get branch x position
            xp=x_attr(k.parent) if k.parent else x ## get parent x position
            y=y_attr(k) ## get y position

            try:
                colours.append(colour_function(k))
            except KeyError:
                colours.append((0.7,0.7,0.7))
            linewidths.append(branchWidth(k))

            if tree_type=='rectangular':
                branches.append(((xp,y),(x,y)))
                if k.branchType=='node' and tree_type=='rectangular':
                    yl,yr=y_attr(k.children[0]),y_attr(k.children[-1])
                    branches.append(((x,yl),(x,yr)))
                    linewidths.append(linewidths[-1])
                    colours.append(colours[-1])
            elif tree_type=='non-baltic':
                yp=y_attr(k.parent) if k.parent else y ## get parent x position
                branches.append(((xp,yp),(xp,y),(x,y)))
            elif tree_type=='unrooted':
                yp=y_attr(k.parent) ## get y position
                branches.append(((xp,yp),(x,y)))
            else:
                pass ## for now

        line_segments = LineCollection(branches,lw=linewidths,ls='-',color=colours,capstyle='projecting')
        ax.add_collection(line_segments)

        return ax

    def plotCircularTree(self,ax,target=None,x_attr=None,y_attr=None,branchWidth=None,colour_function=None,
                         circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,precision=15,**kwargs):

        if target==None: target=lambda k: True
        if x_attr==None: x_attr=lambda k:k.x
        if y_attr==None: y_attr=lambda k:k.y
        if colour_function==None: colour_function=lambda f:'k'
        if branchWidth==None: branchWidth=lambda k:2

        if inwardSpace<0: inwardSpace-=self.treeHeight

        branches=[]
        colours=[]
        linewidths=[]

        circ_s=circStart*math.pi*2
        circ=circFrac*math.pi*2

        if normaliseHeight==None: normaliseHeight=lambda value,values: (value-min(values))/(max(values)-min(values))
        linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

        value_range=list(map(x_attr,self.Objects))
        for k in filter(target,self.Objects): ## iterate over branches
            x=normaliseHeight(x_attr(k)+inwardSpace,value_range) ## get branch x position
            xp=normaliseHeight(x_attr(k.parent)+inwardSpace,value_range) if k.parent.parent else x ## get parent x position
            y=y_attr(k) ## get y position

            try:
                colours.append(colour_function(k))
            except KeyError:
                colours.append((0.7,0.7,0.7))
            linewidths.append(branchWidth(k))

            y=circ_s+circ*y/self.ySpan
            X=math.sin(y)
            Y=math.cos(y)
            branches.append(((X*xp,Y*xp),(X*x,Y*x)))

            if k.branchType=='node':
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

    def plotCircularPoints(self,ax,x_attr=None,y_attr=None,target=None,size_function=None,colour_function=None,circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,
               zorder=None,outline=None,outline_size=None,outline_colour=None,**kwargs):
        if target==None: target=lambda k: k.branchType=='leaf'
        if x_attr==None: x_attr=lambda k:k.x
        if y_attr==None: y_attr=lambda k:k.y
        if size_function==None: size_function=lambda k:40
        if colour_function==None: colour_function=lambda f:'k'
        if zorder==None: zorder=3

        if outline==None: outline=True
        if outline_size==None: outline_size=lambda k: size_function(k)*2
        if outline_colour==None: outline_colour=lambda k: 'k'

        if inwardSpace<0: inwardSpace-=self.treeHeight

        circ_s=circStart*math.pi*2
        circ=circFrac*math.pi*2

        if normaliseHeight==None: normaliseHeight=lambda value,values: (value-min(values))/(max(values)-min(values))
        linspace=lambda start,stop,n: list(start+((stop-start)/(n-1))*i for i in range(n)) if n>1 else stop

        value_range=list(map(x_attr,self.Objects))

        xs=[]
        ys=[]
        colours=[]
        sizes=[]

        outline_xs=[]
        outline_ys=[]
        outline_colours=[]
        outline_sizes=[]
        for k in filter(target,self.Objects):
            x=normaliseHeight(x_attr(k)+inwardSpace,value_range) ## find normalised x position along circle's radius
            y=circ_s+circ*y_attr(k)/self.ySpan ## get y position along circle's perimeter
            X=math.sin(y)*x ## transform
            Y=math.cos(y)*x ## transform

            xs.append(X)
            ys.append(Y)
            colours.append(colour_function(k))
            sizes.append(size_function(k))

            if outline:
                outline_xs.append(xs[-1])
                outline_ys.append(ys[-1])
                outline_colours.append(outline_colour(k))
                outline_sizes.append(outline_size(k))

        ax.scatter(xs,ys,s=sizes,facecolor=colours,edgecolor='none',zorder=zorder,**kwargs) ## put a circle at each tip
        if outline:
            ax.scatter(outline_xs,outline_ys,s=outline_sizes,facecolor=outline_colours,edgecolor='none',zorder=zorder-1,**kwargs) ## put a circle at each tip

        return ax


def make_tree(data,ll=None,verbose=False):
    """
    data is a tree string, ll (LL) is an instance of a tree object
    """
    if isinstance(data,str)==False: ## tree string is not an instance of string (could be unicode) - convert
        data=str(data)

    if ll==None: ## calling without providing a tree object - create one
        ll=tree()
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    stored_i=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop

    while i < len(data): ## while there's characters left in the tree string - loop away
        if stored_i == i and verbose==True: print('%d >%s<'%(i,data[i]))

        assert (stored_i != i),'\nTree string unparseable\nStopped at >>%s<<\nstring region looks like this: %s'%(data[i],data[i:i+5000]) ## make sure that you've actually parsed something last time, if not - there's something unexpected in the tree string
        stored_i=i ## store i for later

        if data[i] == '(': ## look for new nodes
            if verbose==True: print('%d adding node'%(i))
            ll.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character

        cerberus=re.match('(\(|,)([0-9]+)(\[|\:)',data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if cerberus is not None:
            if verbose==True: print('%d adding leaf (BEAST) %s'%(i,cerberus.group(2)))
            ll.add_leaf(i,cerberus.group(2)) ## add tip
            i+=len(cerberus.group(2)) ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match('(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/ ]+)(\'|\"|)(\[)*',data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if cerberus is not None:
            if verbose==True: print('%d adding leaf (non-BEAST) %s'%(i,cerberus.group(3)))
            ll.add_leaf(i,cerberus.group(3).strip('"').strip("'"))  ## add tip
            i+=len(cerberus.group(3))+cerberus.group().count("'")+cerberus.group().count('"') ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match('\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if cerberus is not None:
            if verbose==True: print('%d adding multitype node %s'%(i,cerberus.group(1)))
            i+=len(cerberus.group(1))

        cerberus=re.match('[\(,](#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for beginning of reticulate branch
        if cerberus is not None:
            if verbose==True: print('%d adding outgoing reticulation branch %s'%(i,cerberus.group(1)))
            ll.add_reticulation(cerberus.group(1)) ## add reticulate branch

            destination=None
            for k in ll.Objects: ## iterate over branches parsed so far
                if 'label' in k.traits and k.traits['label']==cerberus.group(1): ## if there's a branch with a matching id
                    if destination==None: ## not set destination before
                        destination=k ## destination is matching node
                    else: ## destination seen before - raise an error (indicates reticulate branch ids are not unique)
                        raise Exception('Reticulate branch not unique: %s seen elsewhere in the tree'%(cerberus.group(1)))
            if destination: ## identified destination of this branch
                if verbose==True: print('identified %s destination'%(cerberus.group(1)))
                ll.cur_node.target=destination ## set current node's target as the destination
                setattr(destination,"contribution",ll.cur_node) ## add contributing edge to destination
            else:
                if verbose==True: print('destination of %s not identified yet'%(cerberus.group(1)))
            i+=len(cerberus.group())-1

        cerberus=re.match('\)(#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for landing point of reticulate branch
        if cerberus is not None:
            if verbose==True: print('%d adding incoming reticulation branch %s'%(i,cerberus.group(1)))
            ll.cur_node.traits['label']=cerberus.group(1) ## set node label

            origin=None ## branch is landing, check if its origin was seen previously
            for k in ll.Objects: ## iterate over currently existing branches
                if isinstance(k,reticulation) and k.name==cerberus.group(1): ## check if any reticulate branches match the origin
                    if origin == None: ## origin not identified yet
                        origin=k ## origin is reticulate branch with the correct name
                    else: ## origin has been identified - shouldn't happen, implies that multiple reticulate branches exist with the same name
                        raise Exception('Reticulate branch not unique: %s seen elsewhere in the tree'%(cerberus.group(1)))
            if origin: ## identified origin
                if verbose==True: print('identified %s origin'%(cerberus.group(1)))
                origin.target=ll.cur_node ## set origin's landing at this node
                setattr(ll.cur_node,"contribution",origin) ## add contributing edge to this node
            else:
                if verbose==True: print('origin of %s not identified yet'%(cerberus.group(1)))
            i+=len(cerberus.group())-1

        cerberus=re.match('(\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&]+)\]',data[i:])## look for MCC comments
        if cerberus is not None:
            if verbose==True: print('%d comment: %s'%(i,cerberus.group(2)))
            comment=cerberus.group(2)
            numerics=re.findall('[,&][A-Za-z\_\.0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall('[,&][A-Za-z\_\.0-9]+=["|\']*[A-Za-z\_0-9\.\+ :\/\(\)\&\-]+["|\']*',comment) ## strings
            treelist=re.findall('[,&][A-Za-z\_\.0-9]+={[A-Za-z\_,{}0-9\. :\/\(\)\&]+}',comment) ## complete history logged robust counting (MCMC trees)
            sets=re.findall('[,&][A-Za-z\_\.0-9\%]+={[A-Za-z\.\-0-9eE,\"\_ :\/\(\)\&]+}',comment) ## sets and ranges
            figtree=re.findall('\![A-Za-z]+=[A-Za-z0-9# :\/\(\)\&]+',comment)

            for vals in strings:
                tr,val=vals.split('=')
                tr=tr[1:]
                if '+' in val:
                    val=val.split('+')[0] ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                ll.cur_node.traits[tr]=val.strip('"')

            for vals in numerics: ## assign all parsed annotations to traits of current branch
                tr,val=vals.split('=') ## split each value by =, left side is name, right side is value
                tr=tr[1:]
                if val.replace('E','',1).replace('e','',1).replace('-','',1).replace('.','',1).isdigit():
                    ll.cur_node.traits[tr]=float(val)

            for val in treelist:
                tr,val=val.split('=')
                tr=tr[1:]
                microcerberus=re.findall('{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}',val)
                ll.cur_node.traits[tr]=[]
                for val in microcerberus:
                    codon,timing,start,end=val.split(',')
                    ll.cur_node.traits[tr].append((int(codon),float(timing),start,end))

            for vals in sets:
                tr,val=vals.split('=')
                tr=tr[1:]
                if 'set' in tr:
                    ll.cur_node.traits[tr]=[]
                    for v in val[1:-1].split(','):
                        if 'set.prob' in tr:
                            ll.cur_node.traits[tr].append(float(v))
                        else:
                            ll.cur_node.traits[tr].append(v.strip('"'))
                else:
                    try:
                        ll.cur_node.traits[tr]=list(map(float,val[1:-1].split(',')))
                    except:
                        print('some other trait: %s'%(vals))

            if len(figtree)>0:
                print('FigTree comment found, ignoring')

            i+=len(cerberus.group()) ## advance in tree string by however many characters it took to encode labels

        cerberus=re.match('([A-Za-z\_\-0-9\.]+)(\:|\;)',data[i:])## look for old school node labels
        if cerberus is not None:
            if verbose==True: print('old school comment found: %s'%(cerberus.group(1)))
            ll.cur_node.traits['label']=cerberus.group(1)

            i+=len(cerberus.group(1))

        microcerberus=re.match('(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
        if microcerberus is not None:
            if verbose==True: print('adding branch length (%d) %.6f'%(i,float(microcerberus.group(2))))
            ll.cur_node.length=float(microcerberus.group(2)) ## set branch length of current node
            i+=len(microcerberus.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            ll.cur_node=ll.cur_node.parent

        if data[i] == ';': ## look for string end
            return ll
            break ## end loop

def make_treeJSON(JSONnode,json_translation,ll=None,verbose=False):
    if 'children' in JSONnode: ## only nodes have children
        new_node=node()
    else:
        new_node=leaf()
        # new_node.numName=JSONnode[json_translation['name']] ## set leaf numName
        new_node.name=JSONnode[json_translation['name']] ## set leaf name to be the same

    if ll is None:
        ll=tree()
        ll.root=new_node
    if 'attr' in JSONnode:
        attr = JSONnode.pop('attr')
        JSONnode.update(attr)

    new_node.parent=ll.cur_node ## set parent-child relationships
    ll.cur_node.children.append(new_node)
    new_node.index=JSONnode[json_translation['name']] ## indexing is based on name
    new_node.traits={n:JSONnode[n] for n in list(JSONnode.keys()) if n!='children'} ## set traits to non-children attributes
    ll.Objects.append(new_node)
    ll.cur_node=new_node

    if 'children' in JSONnode:
        for child in JSONnode['children']:
            make_treeJSON(child,json_translation,ll)
            ll.cur_node=ll.cur_node.parent
    return ll

def loadNewick(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',variableDate=True,absoluteTime=False,verbose=False):
    ll=None
    if isinstance(tree_path,str):
        handle=open(tree_path,'r')
    else:
        handle=tree_path

    for line in handle:
        l=line.strip('\n')
        if '(' in l:
            treeString_start=l.index('(')
            ll=make_tree(l[treeString_start:],verbose=verbose) ## send tree string to make_tree function
            if verbose==True: print('Identified tree string')

    assert ll,'Regular expression failed to find tree string'
    ll.traverse_tree(verbose=verbose) ## traverse tree
    ll.sortBranches() ## traverses tree, sorts branches, draws tree

    if absoluteTime==True:
        tipDates=[]
        tipNames=[]
        for k in ll.getExternal():
            # n=k.numName
            n=k.name
            # k.name=k.numName
            tipNames.append(n)
            cerberus=re.search(tip_regex,n)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))
        assert len(tipDates)>0,'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s'%(tipNames[0],tip_regex,date_fmt)
        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)
    return ll

def loadNexus(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',treestring_regex='tree [A-Za-z\_]+([0-9]+)',variableDate=True,absoluteTime=True,verbose=False):
    tipFlag=False
    tips={}
    tipNum=0
    ll=None
    if isinstance(tree_path,str):
        handle=open(tree_path,'r')
    else:
        handle=tree_path

    for line in handle:
        l=line.strip('\n')

        cerberus=re.search('dimensions ntax=([0-9]+);',l.lower())
        if cerberus is not None:
            tipNum=int(cerberus.group(1))
            if verbose==True: print('File should contain %d taxa'%(tipNum))

        cerberus=re.search(treestring_regex,l)
        if cerberus is not None:
            treeString_start=l.index('(')
            ll=make_tree(l[treeString_start:]) ## send tree string to make_tree function
            if verbose==True: print('Identified tree string')

        if tipFlag==True:
            cerberus=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if cerberus is not None:
                tips[cerberus.group(1)]=cerberus.group(2).strip('"').strip("'")
                if verbose==True: print('Identified tip translation %s: %s'%(cerberus.group(1),tips[cerberus.group(1)]))
            elif ';' not in l:
                print('tip not captured by regex:',l.replace('\t',''))

        if 'translate' in l.lower():
            tipFlag=True
        if ';' in l:
            tipFlag=False

    assert ll,'Regular expression failed to find tree string'
    ll.traverse_tree() ## traverse tree
    ll.sortBranches() ## traverses tree, sorts branches, draws tree
    if len(tips)>0:
        ll.renameTips(tips) ## renames tips from numbers to actual names
        ll.tipMap=tips
    if absoluteTime==True:
        tipDates=[]
        tipNames=[]
        for k in ll.getExternal():
            # if len(tips)>0:
            #     n=k.name
            # else:
            #     n=k.numName
            tipNames.append(k.name)
            cerberus=re.search(tip_regex,k.name)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))

        assert len(tipDates)>0,'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s'%(tipNames[0],tip_regex,date_fmt)
        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)

    return ll

def loadJSON(json_object,json_translation={'name':'name','absoluteTime':'num_date'},verbose=False,sort=True,stats=True):
    """
    Load a nextstrain JSON by providing either the path to JSON or a file handle.
    json_translation is a dictionary that translates JSON attributes to baltic branch attributes (e.g. 'absoluteTime' is called 'num_date' in nextstrain JSONs).
    Note that to avoid conflicts in setting node heights you can either define the absolute time of each node or branch lengths (e.g. if you want a substitution tree).
    """
    assert 'name' in json_translation and ('absoluteTime' in json_translation or 'length' in json_translation or 'height' in json_translation),'JSON translation dictionary missing entries: %s'%(', '.join([entry for entry in ['name','height','absoluteTime','length'] if (entry in json_translation)==False]))
    if verbose==True: print('Reading JSON')

    if isinstance(json_object,str): ## string provided - either nextstrain URL or local path
        if 'nextstrain.org' in json_object: ## nextsrain.org in URL - request it
            if verbose==True: print('Assume URL provided, loading JSON from nextstrain.org')
            import requests
            from io import BytesIO as csio
            auspice_json=json.load(csio(requests.get(json_object).content))
        else: ## not nextstrain.org URL - assume local path to auspice v2 json
            if verbose==True: print('Loading JSON from local path')
            with open(json_object) as json_data:
                auspice_json = json.load(json_data)
    else: ## not string, assume auspice v2 json object given
        if verbose==True: print('Loading JSON from object given')
        auspice_json=json_object

    json_meta=auspice_json['meta']
    json_tree=auspice_json['tree']
    ll=make_treeJSON(json_tree,json_translation,verbose=verbose)

    assert ('absoluteTime' in json_translation and ('length' not in json_translation or 'height' not in json_translation)) or ('absoluteTime' not in json_translation and ('length' in json_translation or 'height' in json_translation)),'Cannot use both absolute time and branch length, include only one in json_translation dictionary.'

    if verbose==True: print('Setting baltic traits from JSON')
    for k in ll.Objects: ## make node attributes easier to access
        for key in k.traits['node_attrs']:
            if isinstance(k.traits['node_attrs'][key],dict):
                if 'value' in k.traits['node_attrs'][key]:
                    k.traits[key]=k.traits['node_attrs'][key]['value']
                if 'confidence' in k.traits['node_attrs'][key]:
                    k.traits['%s_confidence'%(key)]=k.traits['node_attrs'][key]['confidence']
            elif key=='div':
                k.traits['divergence']=k.traits['node_attrs'][key]

    for attr in json_translation: ## iterate through attributes in json_translation
        for k in ll.Objects: ## for every branch
            if isinstance(json_translation[attr],str):
                setattr(k,attr,k.traits[json_translation[attr]]) ## set attribute value for branch
            else:
                setattr(k,attr,json_translation[attr](k)) ## set attribute value with a function for branch

    if 'absoluteTime' in json_translation: ## if using absoluteTime need to set branch lengths for traversals
        for k in ll.Objects:
            if k.absoluteTime and k.parent.absoluteTime:
                k.length=k.absoluteTime-k.parent.absoluteTime
            else:
                k.length=0.0

    if verbose==True: print('Traversing and drawing tree')

    ll.traverse_tree(verbose=verbose)
    ll.drawTree()
    if stats==True:
        ll.treeStats() ## initial traversal, checks for stats
    if sort==True:
        ll.sortBranches() ## traverses tree, sorts branches, draws tree

    cmap={}
    for colouring in json_meta['colorings']:
        if colouring['type']=='categorical' and 'scale' in colouring:
            cmap[colouring['key']]={}
            for entry in colouring['scale']:
                key,value=entry
                cmap[colouring['key']][key]=value
    setattr(ll,'cmap',cmap)

    return ll,json_meta

if __name__ == '__main__':
    import sys
    ll=make_tree(sys.argv[1],ll)
    ll.traverse_tree()
    sys.stdout.write('%s\n'%(ll.treeHeight))
