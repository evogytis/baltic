from matplotlib.collections import LineCollection
import re,copy,math,json,sys
import datetime as dt
from functools import reduce
from matplotlib.collections import LineCollection

__all__ = ['decimalDate', 'convertDate', 'calendarDate', 'reticulation', # make from baltic import * safe
           'clade', 'leaf', 'node', 'tree',
           'make_tree', 'make_treeJSON', 'loadJSON', 'loadNexus', 'loadNewick', 'untangle']

sys.setrecursionlimit(9001)

def decimalDate(date,fmt="%Y-%m-%d",variable=False):
    """ Converts calendar dates in specified format to decimal date. """
    if fmt == "":
        return date
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

def calendarDate(timepoint,fmt='%Y-%m-%d'):
    """ Converts decimal dates to a specified calendar date format. """
    year = int(timepoint)
    rem = timepoint - year

    base = dt.datetime(year, 1, 1)
    result = base + dt.timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    return dt.datetime.strftime(result,fmt)

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
        self.x=None
        self.y=None
        self.width=0.5
        self.target=None

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return False

    def is_node(self):
        return False

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
        self.x=None
        self.y=None
        self.lastHeight=None ## refers to the height of the highest tip in the collapsed clade
        self.lastAbsoluteTime=None ## refers to the absolute time of the highest tip in the collapsed clade
        self.width=1

    def is_leaflike(self):
        return True

    def is_leaf(self):
        return False

    def is_node(self):
        return False

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

    def is_leaflike(self):
        return False

    def is_leaf(self):
        return False

    def is_node(self):
        return True

class leaf: ## leaf class
    def __init__(self):
        self.branchType='leaf'
        self.name=None ## name of tip after translation, since BEAST trees will generally have numbers for taxa but will provide a map at the beginning of the file
        self.index=None ## index of the character that defines this object, will be a unique ID for each object in the tree
        self.length=None ## branch length
        self.absoluteTime=None ## position of tip in absolute time
        self.height=None ## height of tip
        self.parent=None ## parent
        self.traits={} ## trait dictionary
        self.x=None ## position of tip on x axis if the tip were to be plotted
        self.y=None ## position of tip on y axis if the tip were to be plotted

    def is_leaflike(self):
        return False

    def is_leaf(self):
        return True

    def is_node(self):
        return False

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
        assert self.cur_node.is_node(), 'Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses.'
        self.cur_node.children.append(new_node) ## new node is a child of current node
        self.cur_node=new_node ## current node is now new node
        self.Objects.append(self.cur_node) ## add new node to list of objects in the tree

    def add_leaf(self,i,name):
        """ Attach a new leaf (tip) to current node. """
        new_leaf=leaf() ## new instance of leaf object
        new_leaf.index=i ## index is position along tree string
        if self.root is None: self.root=new_leaf

        new_leaf.parent=self.cur_node ## leaf's parent is current node
        assert self.cur_node.is_node(), 'Attempted to add a child to a non-node object. Check if tip names have illegal characters like parentheses.'
        self.cur_node.children.append(new_leaf) ## assign leaf to parent's children
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

        if subtree is None or len([k for k in subtree if k.is_leaf()])==0:
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
                local_tree.tipMap={tipNum: self.tipMap[tipNum] for tipNum in self.tipMap if self.tipMap[tipNum] in [w.name for w in local_tree.getExternal()]} ## copy over the relevant tip translations

            return local_tree

    def singleType(self):
        """ Removes any branches with a single child (multitype nodes). """
        multiTypeNodes=[k for k in self.Objects if k.is_node() and len(k.children)==1]
        while len(multiTypeNodes)>0:
            multiTypeNodes=[k for k in self.Objects if k.is_node() and len(k.children)==1]

            for k in sorted(multiTypeNodes,key=lambda x:-x.height):
                child=k.children[0] ## fetch child
                grandparent=k.parent if k.parent.index else self.root ## fetch grandparent

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

    def traverse_tree(self,cur_node=None,include_condition=None,traverse_condition=None,collect=None,verbose=False):
        if cur_node==None: ## if no starting point defined - start from root
            if verbose==True: print('Initiated traversal from root')
            cur_node=self.root

            if traverse_condition==None and include_condition==None: ## reset heights if traversing from scratch
                for k in self.Objects: ## reset various parameters
                    if k.is_node():
                        k.leaves=set()
                        k.childHeight=None
                    k.height=None

        if traverse_condition==None: traverse_condition=lambda k: True
        if include_condition==None: include_condition=lambda k: k.is_leaf()

        if collect==None: ## initiate collect list if not initiated
            collect=[]

        if cur_node.parent and cur_node.height==None: ## cur_node has a parent - set height if it doesn't have it already
            cur_node.height=cur_node.length+cur_node.parent.height
        elif cur_node.height==None: ## cur_node does not have a parent (root), if height not set before it's zero
            cur_node.height=0.0

        if verbose==True: print('at %s (%s)'%(cur_node.index,cur_node.branchType))

        if include_condition(cur_node): ## test if interested in cur_node
            collect.append(cur_node) ## add to collect list for reporting later

        if cur_node.is_leaf() and self.root!=cur_node: ## cur_node is a tip (and tree is not single tip)
            cur_node.parent.leaves.add(cur_node.name) ## add to parent's list of tips

        elif cur_node.is_node(): ## cur_node is node
            for child in filter(traverse_condition,cur_node.children): ## only traverse through children we're interested
                if verbose==True: print('visiting child %s'%(child.index))
                self.traverse_tree(cur_node=child,include_condition=include_condition,traverse_condition=traverse_condition,verbose=verbose,collect=collect) ## recurse through children
                if verbose==True: print('child %s done'%(child.index))
            assert len(cur_node.children)>0, 'Tried traversing through hanging node without children. Index: %s'%(cur_node.index)
            cur_node.childHeight=max([child.childHeight if child.is_node() else child.height for child in cur_node.children])

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

    def sortBranches(self,descending=True,sort_function=None,sortByHeight=True):
        mod=-1 if descending else 0
        if sortByHeight: # Sort nodes by height and group nodes and leaves together
            """ Sort descendants of each node. """                
            if sort_function==None: sort_function=lambda k: (k.is_node(),-len(k.leaves)*mod,k.length*mod) if k.is_node() else (k.is_node(),k.length*mod)

            for k in self.getInternal(): ## iterate over nodes
                k.children=sorted(k.children,key=sort_function)
        else: # Do not sort by height. Retain leaves at original positions. Only sort nodes
            for k in self.getInternal():
                leavesIdx = [(i,ctr) for ctr, i in enumerate(k.children) if i.branchType=="leaf"] # Get original indices of leaves
                nodes=sorted([x for x in k.children if x.branchType=='node'],key=lambda q:-len(q.leaves)*mod) # Sort nodes only by number of descendants
                children = nodes
                for i in leavesIdx: # Insert leaves back into same positions
                    children.insert(i[1], i[0])
                k.children = children
        self.drawTree() ## update x and y positions of each branch, since y positions will have changed because of sorting

    def drawTree(self,order=None,width_function=None,pad_nodes=None,verbose=False):
        """ Find x and y coordinates of each branch. """
        if order==None:
            order=self.traverse_tree() ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
            if verbose==True: print('Drawing tree in pre-order')
        else:
            if verbose==True: print('Drawing tree with provided order')

        name_order={x.name: i for i,x in enumerate(order)}
        assert len(name_order)==len(order), 'Non-unique names present in tree'
        if width_function==None:
            if verbose==True:
                print('Drawing tree with default widths (1 unit for leaf objects, width+1 for clades)')
            skips=[1 if isinstance(x,leaf) else x.width+1 for x in order]
        else:
            skips=list(map(width_function,order))

        for k in self.Objects: ## reset coordinates for all objects
            k.x=None
            k.y=None

        drawn={} ## drawn keeps track of what's been drawn
        for k in order: ## iterate over tips
            x=k.height ## x position is height
            y_idx=name_order[k.name] ## assign y index
            y=sum(skips[y_idx:])-skips[y_idx]/2.0 ## sum across skips to find y position

            k.x=x ## set x and y coordinates
            k.y=y
            drawn[k.index]=None ## remember that this objects has been drawn

        if pad_nodes!=None: ## will be padding nodes
            for n in pad_nodes: ## iterate over nodes whose descendants will be padded
                idx=sorted([name_order[lf] for lf in n.leaves]) if n.is_node() else [order.index(n)] ## indices of all tips to be padded
                for i,k in enumerate(order): ## iterate over all tips

                    if i<idx[0]: ## tip below clade
                        k.y+=pad_nodes[n] ## pad

                    if (i-1)<idx[-1]: ## tip above clade
                        k.y+=pad_nodes[n] ## pad again

            all_ys=filter(None,self.getParameter('y')) ## get all y positions in tree that aren't None
            minY=min(all_ys) ## get min
            for k in self.getExternal(): ## reset y positions so tree starts at y=0.5
                k.y-=minY-0.5

        assert len(self.getExternal())==len(order),'Number of tips in tree does not match number of unique tips, check if two or more collapsed clades were assigned the same name.'
        storePlotted=0

        while len(drawn)!=len(self.Objects): # keep drawing the tree until everything is drawn
            if verbose==True: print('Drawing iteration %d'%(len(drawn)))
            for k in filter(lambda w:w.index not in drawn,self.getInternal()): ## iterate through internal nodes that have not been drawn
                if len([q.y for q in k.children if q.y!=None])==len(k.children): ## all y coordinates of children known
                    if verbose==True: print('Setting node %s coordinates to'%(k.index)),
                    x=k.height ## x position is height
                    children_y_coords=[q.y for q in k.children if q.y!=None] ## get all existing y coordinates of the node
                    y=sum(children_y_coords)/float(len(children_y_coords)) ## internal branch is in the middle of the vertical bar
                    k.x=x
                    k.y=y
                    drawn[k.index]=None ## remember that this objects has been drawn
                    if verbose==True: print('%s (%s branches drawn)'%(k.y,len(drawn)))
                    minYrange=min([min(child.yRange) if child.is_node() else child.y for child in k.children]) ## get lowest y coordinate across children
                    maxYrange=max([max(child.yRange) if child.is_node() else child.y for child in k.children]) ## get highest y coordinate across children
                    setattr(k,'yRange',[minYrange,maxYrange]) ## assign the maximum extent of children's y coordinates

            if len(self.Objects)>len(drawn):
                assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects (%d branches drawn this iteration, %d branches during previous iteration out of %d total)'%(len(drawn),storePlotted,len(tree.Objects))
            storePlotted=len(drawn) ## remember how many branches were drawn this iteration

        yvalues=[k.y for k in self.Objects] ## all y values
        self.ySpan=max(yvalues)-min(yvalues)+min(yvalues)*2 ## determine appropriate y axis span of tree

        if self.root.is_node():
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

        if n.is_leaf():
            w=2*math.pi*1.0/float(total)
        else:
            w=2*math.pi*len(n.leaves)/float(total)

        if n.parent.x==None:
            n.parent.x=0.0
            n.parent.y=0.0

        n.x = n.parent.x + n.length * math.cos(n.traits['tau'] + w*0.5)
        n.y = n.parent.y + n.length * math.sin(n.traits['tau'] + w*0.5)
        eta=n.traits['tau']

        if n.is_node():
            for ch in n.children:
                if ch.is_leaf():
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
        assert cl.is_node(),'Cannot collapse non-node class'
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
            nodes_to_delete=list(filter(lambda n: n.is_node() and collapseIf(n)==True and n!=newTree.root, newTree.Objects)) ## fetch a list of all nodes who are not the root and who satisfy the condition
        else:
            assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
            assert len([w for w in designated_nodes if w!=newTree.root])==0,'Root node was designated for deletion'

            nodes_to_delete=list(filter(lambda w: w.index in [q.index for q in designated_nodes], newTree.Objects)) ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
        if verbose==True: print('%s nodes set for collapsing: %s'%(len(nodes_to_delete),[w.index for w in nodes_to_delete]))
        assert len(nodes_to_delete)<len(newTree.getInternal())-1,'Chosen cutoff would remove all branches'
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
                    nodes_to_delete==list(filter(lambda n: n.is_node() and collapseIf(n)==True and n!=newTree.root, newTree.Objects))
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

        if cur_node.is_node():
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

        elif cur_node.is_leaf():
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

    def allTMRCAs(self):
        tip_names=[k.name for k in self.getExternal()]
        tmrcaMatrix={x:{y:None if x!=y else 0.0 for y in tip_names} for x in tip_names} ## pairwise matrix of tips

        for k in self.getInternal(): ## iterate over nodes
            all_children=list(k.leaves) ## fetch all descendant tips of node

            for a,tipA in enumerate(all_children):
                for tipB in all_children[a+1:]:
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
        assert len([k for k in keep if k.is_leaf()])==0, "Embedding contains %d non-leaf branches."%(len([k for k in keep if k.is_leaf()==False]))
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
        if verbose==True: print("Finished extracting embedding with %s branches (%s tips, %s nodes)"%(len(embedding),len([w for w in embedding if w.is_leaf()]),len([w for w in embedding if w.is_node()])))
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

    def countLineages(self,t,attr='absoluteTime',condition=lambda x:True):
        return len([k for k in self.Objects if getattr(k.parent,attr)!=None and getattr(k.parent,attr)<t<=getattr(k,attr) and condition(k)])

    def getExternal(self,secondFilter=None):
        """
        Get all branches whose branchType is "leaf".
        A function can be provided to filter internal nodes according to an additional property.
        """
        externals=list(filter(secondFilter,filter(lambda k: k.is_leaf(),self.Objects)))
        return externals

    def getInternal(self,secondFilter=None):
        """
        Get all branches whose branchType is "node".
        A function can be provided to filter internal nodes according to an additional property.
        """
        internals=list(filter(secondFilter,filter(lambda k: k.is_node(),self.Objects)))
        return internals

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

    def getParameter(self,statistic,use_trait=False,which=None):
        """
        Return either branch trait or attribute (default: trait, to switch to attribute set use_trait parameter to False) statistic across branches determined by the which_branches function (default: all objects in the tree).
        Note - branches which do not have the trait or attribute are skipped.
        """
        if which==None:
            branches=self.Objects
        else:
            branches=filter(which,self.Objects)

        if use_trait==False:
            params=[getattr(k,statistic) for k in branches if hasattr(k,statistic)]
        elif use_trait==True:
            params=[k.traits[statistic] for k in branches if statistic in k.traits]
        return params

    def fixHangingNodes(self):
        """
        Remove internal nodes without any children.
        """
        hangingCondition=lambda k: k.is_node() and len(k.children)==0
        hangingNodes=list(filter(hangingCondition,self.Objects)) ## check for nodes without any children (hanging nodes)
        while len(hangingNodes)>0:
            for h in sorted(hangingNodes,key=lambda x:-x.height):
                h.parent.children.remove(h) ## remove old parent from grandparent's children
                hangingNodes.remove(h) ## remove old parent from multitype nodes
                self.Objects.remove(h) ## remove old parent from all objects
            hangingNodes=list(filter(hangingCondition,self.Objects)) ## regenerate list

    def addText(self,ax,target=None,x_attr=None,y_attr=None,text=None,zorder=None,**kwargs):
        if target==None: target=lambda k: k.is_leaf()
        if x_attr==None: x_attr=lambda k: k.x
        if y_attr==None: y_attr=lambda k: k.y
        if text==None: text=lambda k: k.name
        if zorder==None: zorder=4
        for k in filter(target,self.Objects):
            x,y=x_attr(k),y_attr(k)
            z=zorder
            ax.text(x,y,text(k),zorder=z,**kwargs)
        return ax

    def plotPoints(self,ax,x_attr=None,y_attr=None,target=None,size=None,colour=None,
               zorder=None,outline=None,outline_size=None,outline_colour=None,**kwargs):
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
                colours.append((0.7,0.7,0.7))
            linewidths.append(width(k)) if callable(width) else linewidths.append(width)

            if connection_type=='baltic':
                branches.append(((xp,y),(x,y)))
                if k.is_node():
                    yl,yr=y_attr(k.children[0]),y_attr(k.children[-1])
                    branches.append(((x,yl),(x,yr)))
                    linewidths.append(linewidths[-1])
                    colours.append(colours[-1])
            elif connection_type=='elbow':
                yp=y_attr(k.parent) if k.parent else y ## get parent x position
                branches.append(((xp,yp),(xp,y),(x,y)))
            elif connection_type=='direct':
                yp=y_attr(k.parent) ## get y position
                branches.append(((xp,yp),(x,y)))
            else:
                pass ## for now

        line_segments = LineCollection(branches,lw=linewidths,color=colours,capstyle='projecting',**kwargs)
        ax.add_collection(line_segments)
        return ax

    def plotCircularTree(self,ax,target=None,x_attr=None,y_attr=None,width=None,colour=None,
                         circStart=0.0,circFrac=1.0,inwardSpace=0.0,normaliseHeight=None,precision=15,**kwargs):

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

def untangle(trees,cost_function=None,iterations=None,verbose=False):
    """
    Minimise y-axis discrepancies between tips of trees in a list.
    Only the tangling of adjacent trees in the list is minimised, so the order of trees matters.
    Trees do not need to have the same number of tips but tip names should match.
    """
    from itertools import permutations

    if iterations==None: iterations=3
    if cost_function==None: cost_function=lambda pair: math.pow(abs(pair[0]-pair[1]),2)

    y_positions={T: {k.name: k.y for k in T.getExternal()} for T in trees} ## get y positions of all the tips in every tree

    for iteration in range(iterations):
        if verbose: print('Untangling iteration %d'%(iteration+1))
        first_trees=list(range(len(trees)-1))+[-1] ## trees up to next-to-last + last
        next_trees=list(range(1,len(trees)))+[0] ## trees from second + first
        for cur,nex in zip(first_trees,next_trees): ## adjacent pairs
            tree1=trees[cur] ## fetch current tree
            tree2=trees[nex] ## fetch next tree
            if verbose: print('%d vs %d'%(cur,nex))
            for k in sorted(tree2.getInternal(),key=lambda branch: branch.height): ## iterate through nodes of next tree by height (start from root)
                clade_y_positions=sorted([y_positions[tree2][tip] for tip in k.leaves]) ## sorted list of available y coordinates for node
                costs={} ## will store cost of all children permutations
                if len(k.children)>=10: raise RuntimeWarning('Node is too polytomic and untangling will take an astronomically long time')
                if verbose==True: print(len(k.children))
                for permutation in permutations(k.children): ## iterate over permutations of node's children
                    clade_order=sum([[child.name] if child.is_leaf() else list(child.leaves) for child in permutation],[]) ## flat list of tip names as they would appear in permutation order
                    new_y_positions={clade_order[i]: clade_y_positions[i] for i in range(len(clade_y_positions))} ## assign available y positions in order

                    tip_costs=list(map(cost_function,[(y_positions[tree1][tip],new_y_positions[tip]) for tip in clade_order if tip in y_positions[tree1]]))
                    costs[permutation]=sum(tip_costs)/len(tip_costs) ## compute cost of this permutation in relation to next tree

                best=sorted(costs.keys(),key=lambda w: -costs[w])[0] ## get tree with smallest cost
                k.children=list(best) ## reorder children according to minimised cost

            tree2.drawTree() ## compute new y coordinates for nodes
            for k in tree2.getExternal(): ## iterate over tips
                y_positions[tree2][k.name]=k.y ## remember new coordinates

    return trees

def make_tree(data,ll=None,verbose=False):
    """
    data is a tree string, ll (LL) is an instance of a tree object
    """
    patterns = {
        'beast_tip': r'(\(|,)([0-9]+)(\[|\:)',
        'non_beast_tip': r'(\(|,)(\'|\")*([^\(\):\[\'\"#]+)(\'|\"|)*(\[)*'
    }
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

        cerberus=re.match(patterns['beast_tip'],data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if cerberus is not None:
            if verbose==True: print('%d adding leaf (BEAST) %s'%(i,cerberus.group(2)))
            ll.add_leaf(i,cerberus.group(2)) ## add tip
            i+=len(cerberus.group(2)) ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match(patterns['non_beast_tip'],data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if cerberus is not None:
            if verbose==True: print('%d adding leaf (non-BEAST) %s'%(i,cerberus.group(3)))
            ll.add_leaf(i,cerberus.group(3).strip('"').strip("'"))  ## add tip
            i+=len(cerberus.group(3))+cerberus.group().count("'")+cerberus.group().count('"') ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match(r'\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if cerberus is not None:
            if verbose==True: print('%d adding multitype node %s'%(i,cerberus.group(1)))
            i+=len(cerberus.group(1))

        cerberus=re.match(r'[\(,](#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for beginning of reticulate branch
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

        cerberus=re.match(r'\)(#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for landing point of reticulate branch
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

        cerberus=re.match(r'(\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&]+)\]',data[i:])## look for MCC comments
        if cerberus is not None:
            if verbose==True: print('%d comment: %s'%(i,cerberus.group(2)))
            comment=cerberus.group(2)
            numerics=re.findall('[,&][A-Za-z\_\.0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall('[,&][A-Za-z\_\.0-9]+=["|\']*[A-Za-z\_0-9\.\+ :\/\(\)\&\-]+[\"|\']*',comment) ## strings
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
                microcerberus = []
                if val.count(",") == 2:
                    microcerberus=re.findall(r'{([0-9\.\-e]+,[a-z_A-Z]+,[a-z_A-Z]+)}',val)
                elif val.count(",") == 3:
                    microcerberus=re.findall(r'{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}',val)
                ll.cur_node.traits[tr]=[]
                for val in microcerberus:
                    val_split = val.split(',')
                    ll.cur_node.traits[tr].append(val.split(","))

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

        cerberus=re.match(r'([A-Za-z\_\-0-9\.]+)(\:|\;)',data[i:])## look for old school node labels
        if cerberus is not None:
            if verbose==True: print('old school comment found: %s'%(cerberus.group(1)))
            ll.cur_node.traits['label']=cerberus.group(1)

            i+=len(cerberus.group(1))

        microcerberus=re.match(r'(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
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

def loadNewick(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',variableDate=True,absoluteTime=False,verbose=False, sortBranches = True):
    """
    Load newick file
    """
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
    if sortBranches: 
        ll.sortBranches() ## traverses tree, sorts branches, draws tree

    if absoluteTime==True:
        tipDates=[]
        tipNames=[]
        for k in ll.getExternal():
            n=k.name
            tipNames.append(n)
            cerberus=re.search(tip_regex,n)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))
        assert len(tipDates)>0,'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s'%(tipNames[0],tip_regex,date_fmt)
        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)
    if isinstance(tree_path,str):
        handle.close()
    return ll

def loadNexus(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',treestring_regex='tree [A-Za-z\_]+([0-9]+)',variableDate=True,absoluteTime=True,verbose=False, sortBranches=True):
    """
    Load nexus file
    """
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

        cerberus=re.search(treestring_regex,l.lower())
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
    if sortBranches:
        ll.sortBranches() ## traverses tree, sorts branches, draws tree
    if len(tips)>0:
        ll.renameTips(tips) ## renames tips from numbers to actual names
        ll.tipMap=tips
    if absoluteTime==True:
        tipDates=[]
        tipNames=[]
        for k in ll.getExternal():
            tipNames.append(k.name)
            cerberus=re.search(tip_regex,k.name)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))

        assert len(tipDates)>0,'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s'%(tipNames[0],tip_regex,date_fmt)
        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)

    if isinstance(tree_path,str):
        handle.close()
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
                if json_translation[attr] in k.traits:
                    setattr(k,attr,k.traits[json_translation[attr]]) ## set attribute value for branch
                elif 'node_attrs' in k.traits and json_translation[attr] in k.traits['node_attrs']:
                    setattr(k,attr,k.traits['node_attrs'][json_translation[attr]])
                elif 'branch_attrs' in k.traits and json_translation[attr] in k.traits['branch_attrs']:
                    setattr(k,attr,k.traits['branch_attrs'][json_translation[attr]])
                else:
                    raise KeyError('String attribute %s not found in JSON'%(json_translation[attr]))
            elif callable(json_translation[attr]):
                setattr(k,attr,json_translation[attr](k)) ## set attribute value with a function for branch
            else:
                raise AttributeError('Attribute %s neither string nor callable'%(json_translation[attr]))

    for branch_unit in ['height','absoluteTime']: ## iterate between divergence and absolute time
        if branch_unit in json_translation: ## it's available in tree
            for k in ll.Objects: ## iterate over all branches
                cur_branch=getattr(k,branch_unit) ## get parameter for this branch
                par_branch=getattr(k.parent,branch_unit) ## get parameter for parental branch
                k.length=cur_branch-par_branch if cur_branch and par_branch else 0.0 ## difference between current and parent is branch length (or, if parent unavailabel it's 0)

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
