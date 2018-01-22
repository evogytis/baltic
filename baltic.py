import re
import copy
import math
import datetime as dt
import json

def unique(o, idfun=repr):
    """Reduce a list down to its unique elements."""
    seen = {}
    return [seen.setdefault(idfun(e),e) for e in o if idfun(e) not in seen]

def decimalDate(date,fmt="%Y-%m-%d",variable=False,dateSplitter='-'):
    """ Converts calendar dates in specified format to decimal date. """
    if variable==True: ## if date is variable - extract what is available
        dateL=len(date.split(dateSplitter))
        if dateL==2:
            fmt=dateSplitter.join(fmt.split(dateSplitter)[:-1])
        elif dateL==1:
            fmt=dateSplitter.join(fmt.split(dateSplitter)[:-2])

    adatetime=dt.datetime.strptime(date,fmt) ## convert to datetime object
    year = adatetime.year ## get year
    boy = dt.datetime(year, 1, 1) ## get beginning of the year
    eoy = dt.datetime(year + 1, 1, 1) ## get beginning of next year
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds())) ## return fractional year

def convertDate(x,start,end):
    """ Converts calendar dates between given formats """
    return dt.datetime.strftime(dt.datetime.strptime(x,start),end)

class clade: ## clade class
    def __init__(self,givenName):
        self.branchType='leaf' ## clade class poses as a leaf
        self.subtree=None ## subtree will contain all the branches that were collapsed
        self.length=0.0
        self.height=None
        self.absoluteTime=None
        self.parent=None
        self.traits={}
        self.index=None
        self.name=givenName ## the pretend tip name for the clade
        self.numName=givenName
        self.leaves=[]
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
        self.numChildren=0 ## number of tips that are descended from this node
        self.x=None ## X and Y coordinates of this node, once drawTree() is called
        self.y=None
        ## contains references to all tips of this node
        self.leaves=[] ## is a sorted list of tips that are descended from it

class leaf: ## leaf class
    def __init__(self):
        self.branchType='leaf'
        self.name=None ## name of tip after translation, since BEAST trees will generally have numbers for taxa but will provide a map at the beginning of the file
        self.numName=None ## the original name of the taxon, would be an integer if coming from BEAST, otherwise can be actual name
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
        self.root=self.cur_node ## root of the tree is current node
        self.Objects=[] ## tree objects have a flat list of all branches in them
        self.nodes=[] ## nodes is a list of node objects in tree
        self.leaves=[] ## leaves is a list of leaf objects in tree
        self.tipMap=None
        self.treeHeight=0 ## tree height is the distance between the root and the most recent tip
        self.ySpan=0.0

    def add_node(self,i):
        """ Attaches a new node to current node. """
        new_node=node() ## new node instance
        new_node.index=i ## new node's index is the position along the tree string
        new_node.parent=self.cur_node ## new node's parent is current node
        self.cur_node.children.append(new_node) ## new node is a child of current node
        self.cur_node=new_node ## current node is now new node
        self.Objects.append(self.cur_node) ## add new node to list of objects in the tree
        self.nodes.append(self.cur_node)

    def add_leaf(self,i,name):
        """ Attach a new leaf (tip) to current node. """
        new_leaf=leaf() ## new instance of leaf object
        new_leaf.index=i ## index is position along tree string
        new_leaf.numName=name ## numName is the name tip has inside tree string, BEAST trees usually have numbers for tip names
        new_leaf.parent=self.cur_node ## leaf's parent is current node
        self.cur_node.children.append(new_leaf) ## assign leaf to parent's children
        self.cur_node=new_leaf ## current node is now new leaf
        self.Objects.append(self.cur_node) ## add leaf to all objects in the tree
        self.leaves.append(self.cur_node)

    def subtree(self,k=None,subtree=[],traitName=None,converterDict=None):
        """ Generate a subtree (as a baltic tree object) from a traversal.
            If a trait name is provided the traversal occurs within the trait value of the starting node.
            Note - trait-specific traversal can result in multitype trees.
            If this is undesired call singleType() on the resulting subtree afterwards. """
        if len(subtree)==0:
            if traitName:
                subtree=copy.deepcopy(self.traverseWithinTrait(k,traitName,converterDict))
            else:
                subtree=copy.deepcopy(self.traverse_tree(k,include_all=True))
        else:
            subtree=copy.deepcopy(subtree)

        if subtree is None or [w.branchType=='leaf' for w in subtree].count(True)==0:
            return None
        else:
            local_tree=tree() ## create a new tree object where the subtree will be
            local_tree.Objects=subtree ## assign branches to new tree object
            local_tree.root.children.append(subtree[0]) ## connect tree object's root with subtree
            subtree[0].parent=local_tree.root ## subtree's root's parent is tree object's root
            local_tree.root.absoluteTime=subtree[0].absoluteTime-subtree[0].length ## root's absolute time is subtree's root time

            if traitName: ## relying on within-trait traversal
                for nd in local_tree.Objects:
                    if nd.branchType=='node':
                        nd.children=[ch for ch in nd.children if ch in subtree] ## prune each node's children down to branches that were present in the traversal
                hangingNodes=[h for h in local_tree.Objects if h.branchType=='node' and len(h.children)==0] ## check for nodes without any children (hanging nodes)
                while len(hangingNodes)>0:
                    for h in sorted(hangingNodes,key=lambda x:-x.height):
                        h.parent.children.remove(h) ## remove old parent from grandparent's children
                        hangingNodes.remove(h) ## remove old parent from multitype nodes
                        local_tree.Objects.remove(h) ## remove old parent from all objects
                    hangingNodes=[h for h in local_tree.Objects if h.branchType=='node' and len(h.children)==0]

            local_tree.sortBranches() ## sort branches, draw small tree
            return local_tree

    def singleType(self):
        """ Removes any branches with a single child (multitype nodes). """
        multiTypeNodes=[k for k in self.Objects if k.branchType=='node' and len(k.children)==1]
        while len(multiTypeNodes)>0:
            multiTypeNodes=[k for k in self.Objects if k.branchType=='node' and len(k.children)==1]
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

    def treeStats(self):
        """ provide information about the tree """
        self.traverse_tree() ## traverse the tree
        obs=self.Objects ## convenient list of all objects in the tree
        print('\nTree height: %.6f\nTree length: %.6f'%(self.treeHeight,sum([x.length for x in obs]))) ## report the height and length of tree

        nodes=[k for k in obs if k.branchType=='node'] ## get all nodes
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

        if strictlyBifurcating:
            print('strictly bifurcating tree') ## report
        if multiType:
            print('multitype tree') ## report
        if singleton:
            print('singleton tree')
        if hasTraits:
            print('annotations present') ## report


        print('\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n'%(len(obs),len(nodes),len(obs)-len(nodes))) ## report numbers of different objects in the tree

    def traverse_tree(self,startNode=None,include_all=False,verbose=False):
        """ Traverses tree from root. If a starting node is not defined begin traversal from root. By default returns a list of leaf objects that have been visited, optionally returns a list of all objects in the tree. """
        if startNode==None: ## if no starting point defined - start from root
            cur_node=self.root
            start=None ## remember that traversal was started from scratch
            startNode=cur_node
        elif startNode.branchType=='leaf':
            start=startNode
            if include_all==True:
                return [startNode]
            else:
                return [startNode.numName]
        else:
            start=startNode
            cur_node=startNode ## otherwise start from starting node

        self.leaves=[k for k in self.Objects if isinstance(k,leaf)]
        self.nodes=[k for k in self.Objects if isinstance(k,node)]

        if verbose==True:
            print('Verbose traversal initiated')

        for k in self.Objects: ## reset leaves and number of children for every node
            if isinstance(k,node):
                k.leaves=[]
                k.numChildren=0

        seen=[] ## remember what's been visited
        collected=[] ## collect leaf objects along the way
        maxHeight=0 ## check what the maximum distance between the root and the most recent tip is
        if startNode!=None and startNode.height!=None:
            height=startNode.height
        else:
            height=0.0 ## begin at height 0.0

        root=False ## root becomes true once you're unable to go back any further

        while root==False: ## cycle indefinitely as long as there's more of the tree to explore
            # if seen.count(cur_node.index)>3:
            #     if verbose==True:
            #         print 'node %s seen too many times, breaking'%(cur_node.index)
            #     root=True
            #     break

            if cur_node.branchType=='node': ## if currently dealing with a node
                if verbose==True:
                    print('encountered node %s'%(cur_node.index))
                if include_all==True: ## if node hasn't been collected and we want to collect all objects - add it for reporting later
                    collected.append(cur_node)
                ## check whether all children of the node have been seen or whether we're not currently at the root
                ## as long as all children have been seen will descend downwards
                while sum([1 if x.index in seen else 0 for x in cur_node.children])==len(cur_node.children) and cur_node!=startNode.parent:
                    if verbose==True:
                        print('seen all children of node %s'%(cur_node.index))
                    ## check wheteher current node's most recent child is higher than the known highest point in the tree
                    if cur_node.childHeight <= highestTip:
                        cur_node.childHeight = highestTip
                    elif cur_node.childHeight > highestTip:
                        highestTip = cur_node.childHeight

                    if cur_node.index==startNode.index: ## if currently at root - set the root flag to True and break the while loop
                        if verbose==True:
                            print('reached starting point %s'%(cur_node.index))
                        cur_node.height=height
                        root=True
                        break
                    else: ## otherwise...
                        if verbose==True:
                            print('heading to parent of %s'%(cur_node.index))
                        cur_node.parent.numChildren+=cur_node.numChildren ## add the number of current node's children to its parent
                        cur_node.parent.leaves+=cur_node.leaves ## add the list of tip names descended from current node to its parent
                        cur_node.parent.leaves=unique(cur_node.parent.leaves) ## reduce to only unique names
                        cur_node.parent.leaves=sorted(cur_node.parent.leaves) ## sort children
                        cur_node.height=height ## set height
                        height-=float(cur_node.length) ## prepare height value for the eventual descent downwards in the tree
                        cur_node=cur_node.parent ## current node is now current node's parent

                last_seen=[1 if x.index in seen else 0 for x in cur_node.children] ## put 1 for every child of the node that you have seen, 0 for every one that hasn't been seen

                try: ## try finding the first unvisited child of a node
                    idx_to_visit=last_seen.index(0)
                except ValueError:
                    idx_to_visit=0 ## no children seen yet

                if verbose==True:
                    print('visiting %s next (child of %s)'%(cur_node.children[idx_to_visit].index,cur_node.index))

                cur_node.height=height ## set height of current node
                height+=float(cur_node.children[idx_to_visit].length) ## prepare for heading towards the previously unvisited child branch
                cur_node=cur_node.children[idx_to_visit] ## set current node to unvisited child
                seen.append(cur_node.index) ## remember that the child has now been visited

            elif cur_node.branchType=='leaf': ## node dealing with node, are we dealing with a leaf (tip)?
                if verbose==True:
                    print('encountered leaf %s (%s or %s)'%(cur_node.index,cur_node.numName,cur_node.name))
                cur_node.parent.numChildren+=1 ## parent has one more new child (congratulations!)
                cur_node.parent.leaves.append(cur_node.numName) ## add the name of leaf to its parent's list of children names
                cur_node.parent.leaves=sorted(cur_node.parent.leaves) ## sort parent's children list
                seen.append(cur_node.index) ## leaf has now officially been seen
                highestTip=float(height) ## set the height of the leaf as being potentially the highest tip
                cur_node.height=height ## set leaf's height

                #if cur_node not in collected: ## if leaf hasn't been collected - add it for reporting later
                collected.append(cur_node)
                if maxHeight<=float(cur_node.height): ## is this the highest point we've seen in the tree so far?
                    maxHeight=float(cur_node.height)

                height-=float(cur_node.length) ## prepare for heading back
                cur_node=cur_node.parent ## current node is now leaf's parent

        self.treeHeight=float(maxHeight) ## tree height of this tree is the height of the highest tip
        return unique(collected) ## return a list of collected leaf objects

    def renameTips(self,d):
        """ Give each tip its correct label using a dictionary. """
        if self.tipMap!=None:
            d=self.tipMap
        for k in self.leaves: ## iterate through leaf objects in tree
            k.name=d[k.numName] ## change its name

    def sortBranches(self,descending=True):
        """ Sort descendants of each node. """
        if descending==True:
            modifier=-1 ## define the modifier for sorting function later
        elif descending==False:
            modifier=1

        for k in self.Objects: ## iterate over nodes
            if k.branchType=='node':
                ## split node's offspring into nodes and leaves, sort each list individually
                nodes=sorted([x for x in k.children if x.branchType=='node'],key=lambda q:(-len(q.leaves)*modifier,q.length*modifier))
                leaves=sorted([x for x in k.children if x.branchType=='leaf'],key=lambda q:q.length*modifier)

                if modifier==1: ## if sorting one way - nodes come first, leaves later
                    k.children=nodes+leaves
                elif modifier==-1: ## otherwise sort the other way
                    k.children=leaves+nodes
        self.nodes=[k for k in self.Objects if k.branchType=='node']
        self.drawTree() ## update x and y positions of each branch, since y positions will have changed because of sorting

    def drawTree(self,order=None,verbose=False):
        """ Find x and y coordinates of each branch. """
        if order==None:
            order=[x for x in self.traverse_tree() if x.branchType=='leaf'] ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
            if verbose==True:
                print('Drawing tree in pre-order')
        else:
            if verbose==True:
                print('Drawing tree with provided order')

        name_order=[x.numName for x in order]
        skips=[1 if isinstance(x,leaf) else x.width+1 for x in order]

        for k in self.Objects: ## reset coordinates for all objects
            k.x=None
            k.y=None

        storePlotted=0
        drawn=[] ## drawn keeps track of what's been drawn
        while len(drawn)!=len(self.Objects): # keep drawing the tree until everything is drawn
            if verbose==True:
                print('Drawing iteration %d'%(len(drawn)))
            for k in [x for x in self.Objects if x.index not in drawn]: ## iterate through objects that have not been drawn
                if k.branchType=='leaf': ## if leaf - get position of leaf, draw branch connecting tip to parent node
                    if verbose==True:
                        print('Setting leaf %s y coordinate to'%(k.index), end=' ')
                    x=k.height ## x position is height
                    y_idx=name_order.index(k.numName) ## y position of leaf is given by the order in which tips were visited during the traversal
                    y=sum(skips[y_idx:]) ## sum across skips to find y position
                    if verbose==True:
                        print('%s'%(y))
                    if isinstance(k,clade) and skips[y_idx]>1: ## if dealing with collapsed clade - adjust y position to be in the middle of the skip
                        y-=skips[y_idx]/2.0
                        if verbose==True:
                            print('adjusting clade y position to %s'%(y))
                    k.x=x ## set x and y coordinates
                    k.y=y
                    drawn.append(k.index) ## remember that this objects has been drawn

                if k.branchType=='node': ## if parent is non-root node and y positions of all its children are known
                    if len([q.y for q in k.children if q.y!=None])==len(k.children):
                        if verbose==True:
                            print('Setting node %s coordinates'%(k.index))
                        x=k.height ## x position is height
                        children_y_coords=[q.y for q in k.children if q.y!=None] ## get all existing y coordinates of the node
                        y=sum(children_y_coords)/float(len(children_y_coords)) ## internal branch is in the middle of the vertical bar
                        k.x=x
                        k.y=y
                        drawn.append(k.index) ## remember that this objects has been drawn

            assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects'
            storePlotted=len(drawn)
            self.ySpan=sum(skips)

    def drawUnrooted(self,n=None,total=None):
        """
        Calculate x and y coordinates in an unrooted arrangement.
        Code translated from https://github.com/nextstrain/auspice/commit/fc50bbf5e1d09908be2209450c6c3264f298e98c, written by Richard Neher.
        """
        if n==None:
            total=sum([1 if isinstance(x,leaf) else x.width+1 for x in [w for w in self.Objects if w.branchType=='leaf']])
            n=self.root.children[0]
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


    def traverseWithinTrait(self,startNode,traitName,converterDict=None):
        """ Traverse the tree staying within the trait value of the node provided.
            Can also use a converterDict if the trait value needs to be converted.
            Returns list of branches encountered in the traversal. """
        cur_node=startNode
        seen=[]
        collected=[]

        if converterDict==None: ## no dictionary to convert trait values, trying to stay inside trait value of starting node
            stayWithin=startNode.traits[traitName]
        else: ## there's a dictionary, trying to stay within same dictionary value
            stayWithin=converterDict[startNode.traits[traitName]]

        root=False

        if isinstance(startNode,leaf): ## quit immediately if starting from a leaf - nowhere to go
            collected.append(startNode)
            return collected

        else:
            valid=0
            for child in startNode.children: ## iterate over children of the starting node
                if converterDict==None: ## no dictionary to convert trait values, trying to stay inside trait value of starting node
                    chTrait=child.traits[traitName]
                else: ## there's a dictionary, trying to stay within same dictionary value
                    chTrait=converterDict[child.traits[traitName]]

                if chTrait==stayWithin: ## if child is under the same trait value as parent
                    valid+=1 ## add one valid child
            if valid==0: ## need at least one valid child, otherwise return None
                return None

        while root==False:
            if isinstance(cur_node,node):
                ## if all children have been seen and not at root
                while sum([1 if child.index in seen else 0 for child in cur_node.children])==len(cur_node.children) and cur_node.parent!=self.root:

                    if converterDict==None: ## if about to traverse into different trait state - break while loop
                        curTrait=cur_node.traits[traitName]
                    else:
                        curTrait=converterDict[cur_node.traits[traitName]]

                    if curTrait!=stayWithin or cur_node.index==startNode.parent.index or cur_node.index==startNode.index or cur_node.parent.index=='Root':
                        root=True
                        return collected
                    else: ## otherwise keep heading backwards
                        cur_node=cur_node.parent

                seen.append(cur_node.index) ## current node seen

                if cur_node not in collected: ## add current node to collection
                    collected.append(cur_node)

                for child in cur_node.children: ## iterate over children
                    if converterDict==None:
                        childTrait=child.traits[traitName]
                    else:
                        childTrait=converterDict[child.traits[traitName]]

                    if childTrait!=stayWithin: ## if child trait not what is wanted - pretend it's been visited
                        seen.append(child.index)

                last_seen=[1 if child.index in seen else 0 for child in cur_node.children] ## 1 if child was seen, otherwise 0
                if 0 in last_seen: ## if some of the children unseen - go to them
                    idx_to_visit=last_seen.index(0)
                    cur_node=cur_node.children[idx_to_visit]
                else: ## otherwise head back, nothing to see any more
                    cur_node=cur_node.parent

            elif isinstance(cur_node,leaf):
                seen.append(cur_node.index)
                if cur_node not in collected:
                    collected.append(cur_node)
                cur_node=cur_node.parent

    def commonAncestor(self,descendants,numName=False,strict=False):
        types=[desc.__class__ for desc in descendants]
        assert len(set(types))==1,'More than one type of data detected in descendants list'
        if numName==False:
            assert sum([1 if k in [w.name for w in self.Objects if w.branchType=='leaf'] else 0 for k in descendants])==len(descendants),'Not all specified descendants are in tree: %s'%(descendants)
        else:
            assert sum([1 if k in [w.numName for w in self.Objects if w.branchType=='leaf'] else 0 for k in descendants])==len(descendants),'Not all specified descendants are in tree: %s'%(descendants)
        dtype=list(set(types))[0]
        allAncestors=sorted([k for k in self.Objects if k.branchType=='node' and len(k.leaves)>=len(descendants)],key=lambda x:x.height)
        if numName==False:
            ancestor=[k for k in allAncestors if sum([[self.tipMap[w] for w in k.leaves].count(l) for l in descendants])==len(descendants)][-1]
        else:
            ancestor=[k for k in allAncestors if sum([[w for w in k.leaves].count(l) for l in descendants])==len(descendants)][-1]

        if strict==False:
            return ancestor
        elif strict==True and len(ancestor.leaves)==len(descendants):
            return ancestor
        elif strict==True and len(ancestor.leaves)>len(descendants):
            return None

    def collapseSubtree(self,cl,givenName,verbose=False,widthFunction=lambda x:x):
        """ Collapse an entire subtree into a clade object. """
        assert cl.branchType=='node','Cannot collapse non-node class'
        collapsedClade=clade(givenName)
        collapsedClade.index=cl.index
        collapsedClade.length=cl.length
        collapsedClade.height=cl.height
        collapsedClade.parent=cl.parent
        collapsedClade.absoluteTime=cl.absoluteTime
        collapsedClade.traits=cl.traits
        collapsedClade.width=widthFunction(widthFunction(cl.leaves))

        if verbose==True:
            print('Replacing node %s (parent %s) with a clade class'%(cl.index,cl.parent.index))
        parent=cl.parent
        #collapsedClade.subtree=cl

        remove_from_tree=self.traverse_tree(cl,include_all=True)
        collapsedClade.subtree=remove_from_tree
        assert len(remove_from_tree)<len(self.Objects),'Attempted collapse of entire tree'
        collapsedClade.lastHeight=max([x.height for x in remove_from_tree])
        collapsedClade.lastAbsoluteTime=max([x.absoluteTime for x in remove_from_tree])

        for k in remove_from_tree:
            self.Objects.remove(k)

        parent.children.remove(cl)
        parent.children.append(collapsedClade)
        self.Objects.append(collapsedClade)
        collapsedClade.parent=parent
        if self.tipMap!=None:
            self.tipMap[givenName]=givenName

        self.traverse_tree()
        self.sortBranches()

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
            nodes_to_delete=[n for n in newTree.Objects if n.branchType=='node' and collapseIf(n)==True and n.parent!=newTree.root] ## fetch a list of all nodes who are not the root and who satisfy the condition
        else:
            assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
            assert len([w for w in designated_nodes if w.parent.index=='Root'])==0,'Root node was designated for deletion'
            nodes_to_delete=[w for w in newTree.Objects if w.index in [q.index for q in designated_nodes]] ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
        if verbose==True:
            print('%s nodes set for collapsing: %s'%(len(nodes_to_delete),[w.index for w in nodes_to_delete]))
        #assert len(nodes_to_delete)<len(newTree.nodes)-1,'Chosen cutoff would remove all branches'
        while len(nodes_to_delete)>0: ## as long as there are branches to be collapsed - keep reducing the tree

            if verbose==True:
                print('Continuing collapse cycle, %s nodes left'%(len(nodes_to_delete)))
            for k in sorted(nodes_to_delete,key=lambda x:-x.height): ## start with branches near the tips
                zero_node=k.children ## fetch the node's children
                k.parent.children+=zero_node ## add them to the zero node's parent
                old_parent=k ## node to be deleted is the old parent
                new_parent=k.parent ## once node is deleted, the parent to all their children will be the parent of the deleted node
                if new_parent==None:
                    new_parent=self.root
                if verbose==True:
                    print('Removing node %s, attaching children %s to node %s'%(old_parent.index,[w.index for w in k.children],new_parent.index))
                for w in newTree.Objects: ## assign the parent of deleted node as the parent to any children of deleted node
                    if w.parent==old_parent:
                        w.parent=new_parent
                        w.length+=old_parent.length
                        if verbose==True:
                            print('Fixing branch length for node %s'%(w.index))
                k.parent.children.remove(k) ## remove traces of deleted node - it doesn't exist as a child, doesn't exist in the tree and doesn't exist in the nodes list
                newTree.Objects.remove(k)
                newTree.nodes.remove(k)

                nodes_to_delete.remove(k) ## in fact, the node never existed

                if len(designated_nodes)==0:
                    nodes_to_delete=[n for n in newTree.Objects if n.branchType=='node' and collapseIf(n)==True and n.parent!=newTree.root]
                else:
                    assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
                    assert len([w for w in designated_nodes if w.parent.index=='Root'])==0,'Root node was designated for deletion'
                    nodes_to_delete=[w for w in newTree.Objects if w.index in [q.index for q in designated_nodes]]


                if verbose==True:
                    print('Removing references to node %s'%(k.index))
        newTree.sortBranches() ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree ## return collapsed tree

    def toString(self,traits=[],numName=False,verbose=False,nexus=False):
        """ Output the topology of the tree with branch lengths to string """
        cur_node=self.root.children[-1]
        seen=[]
        tree_string=[]
        root=False

        if cur_node.branchType=='node':
            tree_string.append('(')

        while root==False:
            if verbose==True:
                print('Branch this iteration: %s'%(cur_node.index))
            if cur_node.branchType=='node':
                if verbose==True:
                    print('Encountered node %s'%(cur_node.index))
                ## if all children have been seen and not at root
                while sum([1 if x.index in seen else 0 for x in cur_node.children])==len(cur_node.children):
                    if verbose==True:
                        print('Recursing downwards, currently at %s'%(cur_node.index))
                        print('%s %s'%(len(set(seen)),len(cur_node.leaves)))
                    if cur_node.index=='Root':
                        root=True
                        break
                    else:
                        comment=[]
                        if len(traits)>0:
                            for tr in traits:
                                if tr in cur_node.traits:
                                    if isinstance(cur_node.traits[tr],str):
                                        comment.append('%s="%s"'%(tr,cur_node.traits[tr]))
                                    elif isinstance(cur_node.traits[tr],float):
                                        comment.append('%s=%s'%(tr,cur_node.traits[tr]))
                                    elif isinstance(cur_node.traits[tr],list):
                                        rangeComment=[]
                                        for val in cur_node.traits[tr]:
                                            if isinstance(val,str):
                                                rangeComment.append('"%s"'%(val))
                                            elif isinstance(val,float):
                                                rangeComment.append('%s'%(val))
                                        comment.append('%s={%s}'%(tr,','.join(rangeComment)))
                        if len(comment)>0:
                            comment=','.join(comment)
                            comment='[&'+comment+']'
                            tree_string.append(')%s:%8f'%(comment,cur_node.length)) ## end of node, add branch length with annotations
                        else:
                            tree_string.append('):%8f'%(cur_node.length)) ## end of node, add branch length
                        cur_node=cur_node.parent ## go back

                seen.append(cur_node.index)
                last_seen=[1 if x.index in seen else 0 for x in cur_node.children]

                if 0 in last_seen: ## there's unvisited children
                    idx_to_visit=last_seen.index(0)
                    if idx_to_visit>0: ## first child definitely seen, so visiting the second child at the very least
                        tree_string.append(',') ## add comma to indicate bifurcation
                    cur_node=cur_node.children[idx_to_visit]

                    if cur_node.branchType=='node':
                        tree_string.append('(') ## dealing with new node, add (
                    else:
                        comment=[]
                        if len(traits)>0: ## if a designated trait list exists
                            for tr in traits: ## iterate over traits
                                if tr in cur_node.traits:
                                    if isinstance(cur_node.traits[tr],str):
                                        comment.append('%s="%s"'%(tr,cur_node.traits[tr]))
                                    elif isinstance(cur_node.traits[tr],float):
                                        comment.append('%s=%s'%(tr,cur_node.traits[tr]))
                                    elif isinstance(cur_node.traits[tr],list):
                                        rangeComment=[]
                                        for val in cur_node.traits[tr]:
                                            if isinstance(val,str):
                                                rangeComment.append('"%s"'%(val))
                                            elif isinstance(val,float):
                                                rangeComment.append('%s'%(val))
                                        comment.append('%s={%s}'%(tr,','.join(rangeComment)))

                        if numName==False: ## if real names wanted
                            assert cur_node.name!=None,'Tip does not have converted names' ## assert they have been converted
                            treeName=cur_node.name ## designate real name
                        elif numName==True: ## if number names wanted
                            treeName=cur_node.numName ## designated numName

                        if len(comment)>0:
                            comment=','.join(comment) ## join up traits with commas
                            comment='[&'+comment+']' ## add starting and end bracket
                            tree_string.append('\'%s\'%s:%8f'%(treeName,comment,cur_node.length)) ## dealing with tip, write out name, add branch length with annotation
                        else:
                            tree_string.append('\'%s\':%8f'%(treeName,cur_node.length)) ## dealing with tip, write out name, add branch length

                else: ## all children seen, clade's end
                    comment=[]
                    if len(traits)>0:
                        for tr in traits:
                            if tr in cur_node.traits:
                                    if isinstance(cur_node.traits[tr],str):
                                        comment.append('%s="%s"'%(tr,cur_node.traits[tr]))
                                    elif isinstance(cur_node.traits[tr],float):
                                        comment.append('%s=%s'%(tr,cur_node.traits[tr]))
                                    elif isinstance(cur_node.traits[tr],list):
                                        rangeComment=[]
                                        for val in cur_node.traits[tr]:
                                            if isinstance(val,str):
                                                rangeComment.append('"%s"'%(val))
                                            elif isinstance(val,float):
                                                rangeComment.append('%s'%(val))
                                        comment.append('%s={%s}'%(tr,','.join(rangeComment)))
                    if len(comment)>0:
                        comment=','.join(comment)
                        comment='[&'+comment+']'
                        tree_string.append(')%s:%8f'%(comment,cur_node.length))
                    else:
                        tree_string.append('):%8f'%(cur_node.length))
                    cur_node=cur_node.parent

            elif cur_node.branchType=='leaf':
                if verbose==True:
                    print('Encountered leaf %s'%(cur_node.index))

                seen.append(cur_node.index)
                if cur_node.parent.index=='Root':
                    root=True
                    break
                else:
                    cur_node=cur_node.parent
        # if ''.join(tree_string).count(')')-1==''.join(tree_string).count('('):
        #     tree_string.insert(0,'(')

        if nexus==True:
            return '#NEXUS\nBegin trees;\ntree TREE1 = [&R] %s;\nEnd;'%(''.join(tree_string))
        else:
            return ''.join(tree_string)+';' ## the coup de grace

    def allTMRCAs(self):
        tip_names=[k.numName for k in self.Objects if isinstance(k,leaf)]
        tmrcaMatrix={x:{y:None for y in tip_names} for x in tip_names} ## pairwise matrix of tips

        for k in self.Objects:
            if isinstance(k,node):
                all_children=k.leaves ## fetch all descendant tips of node
                #print all_children
                for x in range(0,len(all_children)-1): ## for all pairwise comparisons of tips
                    for y in range(x+1,len(all_children)):
                        tipA=all_children[x]
                        tipB=all_children[y]
                        if tmrcaMatrix[tipA][tipB]==None or tmrcaMatrix[tipA][tipB]<=k.absoluteTime: ## if node's time is more recent than previous entry - set new TMRCA value for pair of tips
                            tmrcaMatrix[tipA][tipB]=k.absoluteTime
                            tmrcaMatrix[tipB][tipA]=k.absoluteTime
        return tmrcaMatrix

    def reduceTree(self,keep):
        """
        Reduce the tree to just those tracking a small number of tips.
        Returns a new baltic tree object.
        """
        assert len(keep)>0,"No tips given to reduce the tree to."
        assert len([k for k in keep if k.branchType!='leaf'])==0, "Embedding contains %d non-leaf branches."%(len([k for k in keep if k.branchType!='leaf']))
        embedding=[]
        reduced_tree=copy.deepcopy(self) ## new tree object
        for k in reduced_tree.Objects: ## deep copy branches from current tree
            if k.index in [q.index for q in keep]: ## if branch is designated as one to keep
                cur_b=k
                while cur_b: ## descend to root
                    embedding.append(cur_b) ## keep track of the path to root
                    cur_b=cur_b.parent
                embedding=list(set(embedding)) ## prune down to only unique branches

        reduced_tree.Objects=sorted([k for k in embedding if k.index!='Root'],key=lambda x:x.height) ## assign branches that are kept to new tree's Objects
        reduced_tree.root=[k for k in embedding if k.index=='Root'][0]
        # reduced_tree.tipMap=self.tipMap

        for k in reduced_tree.Objects: ## iterate through reduced tree
            if k.branchType=='node': ## node
                k.children = [c for c in k.children if c in embedding] ## only keep children that are present in lineage traceback

        reduced_tree.traverse_tree() ## traverse
        reduced_tree.sortBranches() ## sort

        return reduced_tree ## return new tree

    def countLineages(self,t,condition=lambda x:True):
        return len([k for k in self.Objects if k.parent.absoluteTime<t<=k.absoluteTime and condition(k)])
    
def make_tree(data,ll,verbose=False):
    """
    data is a tree string, ll (LL) is an instance of a tree object
    """
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    stored_i=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop

    while i < len(data): ## while there's characters left in the tree string - loop away
        if stored_i == i and verbose==True:
            print('%d >%s<'%(i,data[i]))

        assert (stored_i != i),'\nTree string unparseable\nStopped at >>%s<<\nstring region looks like this: %s'%(data[i],data[i:i+5000]) ## make sure that you've actually parsed something last time, if not - there's something unexpected in the tree string
        stored_i=i ## store i for later

        if data[i] == '(': ## look for new nodes
            if verbose==True:
                print('%d adding node'%(i))
            ll.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character

        cerberus=re.match('(\(|,)([0-9]+)(\[|\:)',data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if cerberus is not None:
            if verbose==True:
                print('%d adding leaf (BEAST) %s'%(i,cerberus.group(2)))
            ll.add_leaf(i,cerberus.group(2)) ## add tip
            i+=len(cerberus.group(2)) ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match('(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/]+)(\'|\"|)(\[)*',data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if cerberus is not None:
            if verbose==True:
                print('%d adding leaf (non-BEAST) %s'%(i,cerberus.group(3)))
            ll.add_leaf(i,cerberus.group(3).strip('"').strip("'"))  ## add tip
            i+=len(cerberus.group(3))+cerberus.group().count("'")+cerberus.group().count('"') ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match('\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if cerberus is not None:
            if verbose==True:
                print('%d adding multitype node %s'%(i,cerberus.group(1)))
            i+=len(cerberus.group(1))

        cerberus=re.match('(\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!#]+)\]',data[i:])## look for MCC comments
        #cerberus=re.match('\[&[A-Za-z\_\-{}\,0-9\.\%=\"\+]+\]',data[i:])## look for MCC comments
        if cerberus is not None:
            if verbose==True:
                print('%d comment: %s'%(i,cerberus.group(2)))
            comment=cerberus.group(2)
            numerics=re.findall('[,&][A-Za-z\_\.0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall('[,&][A-Za-z\_\.0-9]+=["|\']*[A-Za-z\_0-9\.\+]+["|\']*',comment) ## strings
            treelist=re.findall('[,&][A-Za-z\_\.0-9]+={[A-Za-z\_,{}0-9\.]+}',comment) ## complete history logged robust counting (MCMC trees)
            sets=re.findall('[,&][A-Za-z\_\.0-9\%]+={[A-Za-z\.\-0-9eE,\"\_]+}',comment) ## sets and ranges
            figtree=re.findall('\![A-Za-z]+=[A-Za-z0-9#]+',comment)

            for vals in strings:
                tr,val=vals.split('=')
                tr=tr[1:]
                if '+' in val:
                    val=val.split('+')[0] ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                ll.cur_node.traits[tr]=val.strip('"')

            for vals in numerics: ## assign all parsed annotations to traits of current branch
                tr,val=vals.split('=') ## split each value by =, left side is name, right side is value
                tr=tr[1:]
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
            if verbose==True:
                print('old school comment found: %s'%(cerberus.group(1)))
            ll.cur_node.traits['label']=cerberus.group(1)

            i+=len(cerberus.group(1))

        microcerberus=re.match('(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
        if microcerberus is not None:
            if verbose==True:
                print('adding branch length (%d) %.6f'%(i,float(microcerberus.group(2))))
            ll.cur_node.length=float(microcerberus.group(2)) ## set branch length of current node
            i+=len(microcerberus.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            ll.cur_node=ll.cur_node.parent

        if data[i] == ';': ## look for string end
            break ## end loop

def make_treeJSON(tree,parent,JSONnode,json_translation):
    if 'attr' in JSONnode:
        attr = JSONnode.pop('attr')
        JSONnode.update(attr)

    if 'children' in JSONnode: ## only nodes have children
        new_node=node()
    else:
        new_node=leaf()
        new_node.numName=JSONnode[json_translation['name']] ## set leaf numName
        new_node.name=JSONnode[json_translation['name']] ## set leaf name to be the same

    new_node.parent=tree.cur_node ## set parent-child relationships
    tree.cur_node.children.append(new_node)
    new_node.index=JSONnode[json_translation['name']] ## indexing is based on name
    new_node.traits={n:JSONnode[n] for n in list(JSONnode.keys()) if n!='children'} ## set traits to non-children attributes
    tree.Objects.append(new_node)
    tree.cur_node=new_node

    if isinstance(new_node, node):
        tree.nodes.append(new_node)
    else:
        tree.leaves.append(new_node)

    if 'children' in JSONnode:
        for child in JSONnode['children']:
            make_treeJSON(tree,JSONnode,child,json_translation)
            tree.cur_node=tree.cur_node.parent


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
            if verbose==True:
                print('File should contain %d taxa'%(tipNum))

        cerberus=re.search(treestring_regex,l)
        if cerberus is not None:
            treeString_start=l.index('(')
            ll=tree() ## new instance of tree
            make_tree(l[treeString_start:],ll) ## send tree string to make_tree function
            if verbose==True:
                print('Identified tree string')

        if tipFlag==True:
            cerberus=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if cerberus is not None:
                tips[cerberus.group(1)]=cerberus.group(2).strip('"').strip("'")
                if verbose==True:
                    print('Identified tip translation %s: %s'%(cerberus.group(1),tips[cerberus.group(1)]))
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
        for k in ll.leaves:
            if len(tips)>0:
                n=k.name
            else:
                n=k.numName

            cerberus=re.search(tip_regex,n)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))

        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)

    return ll


def loadJSON(tree_path,json_translation={'name':'strain','height':'tvalue'},json_meta=None,verbose=False,sort=True,stats=True):
    assert 'name' in json_translation and 'height' in json_translation,'JSON translation dictionary missing entries: %s'%(', '.join([entry for entry in ['name','height'] if (entry in json_translation)==False]))
    ll=tree()
    if verbose==True:
        print('Reading JSON')
    with open(tree_path) as json_data:
        d = json.load(json_data)
        make_treeJSON(ll,'root',d,json_translation)

    if 'absoluteTime' in json_translation:
        if verbose==True:
            print('Setting absolute time')
        for k in ll.Objects:
            setattr(k,'absoluteTime',k.traits[json_translation['absoluteTime']])

    if verbose==True:
        print('Setting heights')
    for k in ll.Objects:
        setattr(k,'height',k.traits[json_translation['height']])

    if verbose==True:
        print('Setting lengths')
    for k in ll.Objects:
        k.length=k.height-k.parent.height

    if verbose==True:
        print('Traversing and drawing tree')

    ll.traverse_tree()
    ll.drawTree()
    if stats==True:
        ll.treeStats() ## initial traversal, checks for stats
    if sort==True:
        ll.sortBranches() ## traverses tree, sorts branches, draws tree

    if json_meta:
        metadata=json.load(open(json_meta['file'],'r'))
        cmap=dict(metadata['color_options'][json_meta['traitName']]['color_map'])
        setattr(ll,'cmap',cmap)

    return ll

if __name__ == '__main__':
    import sys
    ll=tree()
    make_tree(sys.argv[1],ll)
    ll.traverse_tree()
    sys.stdout.write('%s\n'%(ll.treeHeight))
