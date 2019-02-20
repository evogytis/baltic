import re
import copy
import math
import datetime as dt
import json

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
        self.leaves=None
        self.length=0.0
        self.height=None
        self.absoluteTime=None
        self.parent=None
        self.traits={}
        self.index=None
        self.name=givenName ## the pretend tip name for the clade
        self.numName=givenName
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
        self.root=None #self.cur_node ## root of the tree is current node
        self.Objects=[] ## tree objects have a flat list of all branches in them
        self.tipMap=None
        self.treeHeight=0 ## tree height is the distance between the root and the most recent tip
        self.ySpan=0.0

    def add_node(self,i):
        """ Attaches a new node to current node. """
        new_node=node() ## new node instance
        new_node.index=i ## new node's index is the position along the tree string
        if self.root is None:
            self.root=new_node

        new_node.parent=self.cur_node ## new node's parent is current node
        self.cur_node.children.append(new_node) ## new node is a child of current node
        self.cur_node=new_node ## current node is now new node
        self.Objects.append(self.cur_node) ## add new node to list of objects in the tree

    def add_leaf(self,i,name):
        """ Attach a new leaf (tip) to current node. """
        new_leaf=leaf() ## new instance of leaf object
        new_leaf.index=i ## index is position along tree string
        if self.root is None:
            self.root=new_leaf

        new_leaf.parent=self.cur_node ## leaf's parent is current node
        self.cur_node.children.append(new_leaf) ## assign leaf to parent's children
        new_leaf.numName=name ## numName is the name tip has inside tree string, BEAST trees usually have numbers for tip names
        self.cur_node=new_leaf ## current node is now new leaf
        self.Objects.append(self.cur_node) ## add leaf to all objects in the tree

    def subtree(self,k=None,traverse_condition=None):
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
            local_tree.root=subtree[0] ## connect tree object's root with subtree

            subtree_set=set(subtree) ## turn branches into set for quicker look up later
            if traverse_condition is not None: ## didn't use default traverse condition, might need to deal with hanging nodes and prune children
                for nd in local_tree.getInternal(): ## iterate over nodes
                    nd.children=filter(lambda k:k in subtree_set,nd.children) ## only keep children seen in traversal

                local_tree.fixHangingNodes()

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

        if strictlyBifurcating:
            print('strictly bifurcating tree') ## report
        if multiType:
            print('multitype tree') ## report
        if singleton:
            print('singleton tree')
        if hasTraits:
            print('annotations present') ## report

        print('\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n'%(len(obs),len(nodes),len(obs)-len(nodes))) ## report numbers of different objects in the tree

    def traverse_tree(self,cur_node=None,include_condition=lambda k:k.branchType=='leaf',traverse_condition=lambda k:True,collect=None,verbose=False):
        if cur_node==None: ## if no starting point defined - start from root
            for k in self.Objects: ## reset various parameters
                if k.branchType=='node':
                    k.leaves=set()
                    k.childHeight=None
                k.height=None

            if verbose==True:
                print('Initiated traversal from root')
            cur_node=self.root#.children[-1]

        if collect==None: ## initiate collect list if not initiated
            collect=[]

        if cur_node.parent and cur_node.height==None: ## cur_node has a parent - set height if it doesn't already
            cur_node.height=cur_node.length+cur_node.parent.height
        elif cur_node.height==None: ## cur_node does not have a parent (root), if height not set before it's zero
            cur_node.height=0.0

        if verbose==True:
            print('at %s (%s)'%(cur_node.index,cur_node.branchType))

        if include_condition(cur_node): ## test if interested in cur_node
            collect.append(cur_node) ## add to collect list for reporting later

        if cur_node.branchType=='leaf': ## cur_node is a tip
            cur_node.parent.leaves.add(cur_node.numName) ## add to parent's list of tips

        elif cur_node.branchType=='node': ## cur_node is node
            for child in filter(traverse_condition,cur_node.children): ## only traverse through children we're interested
                if verbose==True:
                    print('visiting child %s'%(child.index))
                self.traverse_tree(cur_node=child,include_condition=include_condition,traverse_condition=traverse_condition,verbose=verbose,collect=collect) ## recurse through children
                if verbose==True:
                    print('child %s done'%(child.index))
            assert len(cur_node.children)>0, 'Tried traversing through hanging node without children. Index: %s'%(cur_node.index)
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
            k.name=d[k.numName] ## change its name

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

    def drawTree(self,order=None,verbose=False):
        """ Find x and y coordinates of each branch. """
        if order==None:
            order=self.traverse_tree() ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
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
        drawn={} ## drawn keeps track of what's been drawn
        while len(drawn)!=len(self.Objects): # keep drawing the tree until everything is drawn
            if verbose==True:
                print('Drawing iteration %d'%(len(drawn)))
            for k in filter(lambda w:w.index not in drawn,self.Objects): ## iterate through objects that have not been drawn
                if k.branchType=='leaf': ## if leaf - get position of leaf, draw branch connecting tip to parent node
                    if verbose==True:
                        print('Setting leaf %s y coordinate to'%(k.index))
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
                    drawn[k.index]=None ## remember that this objects has been drawn
                    if hasattr(k.parent,'yRange')==False: ## if parent doesn't have a maximum extent of its children's y coordinates
                        setattr(k.parent,'yRange',[k.y,k.y]) ## assign it

                if k.branchType=='node': ## if parent is non-root node and y positions of all its children are known
                    if len([q.y for q in k.children if q.y!=None])==len(k.children):
                        if verbose==True:
                            print('Setting node %s coordinates'%(k.index))
                        x=k.height ## x position is height
                        children_y_coords=[q.y for q in k.children if q.y!=None] ## get all existing y coordinates of the node
                        y=sum(children_y_coords)/float(len(children_y_coords)) ## internal branch is in the middle of the vertical bar
                        k.x=x
                        k.y=y
                        drawn[k.index]=None ## remember that this objects has been drawn
                        minYrange=min([min(child.yRange) if child.branchType=='node' else child.y for child in k.children]) ## get lowest y coordinate across children
                        maxYrange=max([max(child.yRange) if child.branchType=='node' else child.y for child in k.children]) ## get highest y coordinate across children
                        setattr(k,'yRange',[minYrange,maxYrange]) ## assign the maximum extent of children's y coordinates

            assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects'
            storePlotted=len(drawn)
            self.ySpan=sum(skips)

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

    def commonAncestor(self,descendants,numName=False,strict=False):
        types=[desc.__class__ for desc in descendants]
        assert len(set(types))==1,'More than one type of data detected in descendants list'
        if numName==False:
            assert sum([1 if k in [w.name for w in self.getExternal()] else 0 for k in descendants])==len(descendants),'Not all specified descendants are in tree: %s'%(descendants)
        else:
            assert sum([1 if k in [w.numName for w in self.getExternal()] else 0 for k in descendants])==len(descendants),'Not all specified descendants are in tree: %s'%(descendants)
        dtype=list(set(types))[0]
        allAncestors=sorted([k for k in self.Objects if (k.branchType=='node' or isinstance(k,clade)) and len(k.leaves)>=len(descendants)],key=lambda x:x.height)
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

        if verbose==True:
            print('Replacing node %s (parent %s) with a clade class'%(cl.index,cl.parent.index))
        parent=cl.parent

        remove_from_tree=self.traverse_tree(cl,include_condition=lambda k: True)
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
            nodes_to_delete=filter(lambda n: n.branchType=='node' and collapseIf(n)==True and n!=newTree.root, newTree.Objects) ## fetch a list of all nodes who are not the root and who satisfy the condition
        else:
            assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
            assert len([w for w in designated_nodes if w!=newTree.root])==0,'Root node was designated for deletion'

            nodes_to_delete=filter(lambda w: w.index in [q.index for q in designated_nodes], newTree.Objects) ## need to look up nodes designated for deletion by their indices, since the tree has been copied and nodes will have new memory addresses
        if verbose==True:
            print('%s nodes set for collapsing: %s'%(len(nodes_to_delete),[w.index for w in nodes_to_delete]))
        # assert len(nodes_to_delete)<len(newTree.getInternal())-1,'Chosen cutoff would remove all branches'
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

                nodes_to_delete.remove(k) ## in fact, the node never existed

                if len(designated_nodes)==0:
                    nodes_to_delete==filter(lambda n: n.branchType=='node' and collapseIf(n)==True and n!=newTree.root, newTree.Objects)
                else:
                    assert [w.branchType for w in designated_nodes].count('node')==len(designated_nodes),'Non-node class detected in list of nodes designated for deletion'
                    assert len([w for w in designated_nodes if w!=newTree.root])==0,'Root node was designated for deletion'
                    nodes_to_delete=[w for w in newTree.Objects if w.index in [q.index for q in designated_nodes]]

                if verbose==True:
                    print('Removing references to node %s'%(k.index))
        newTree.sortBranches() ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree ## return collapsed tree

    def toString(self,cur_node=None,traits=None,numName=False,verbose=False,nexus=False,string_fragment=None,traverse_condition=None,json=False):
        """ Output the topology of the tree with branch lengths and comments to stringself.
            cur_node: starting point (default: None, starts at root)
            traits: list of keys that will be used to output entries in traits dict of each branch (default: all traits)
            numName: boolean, whether encoded (True) or decoded (default: False) tip names will be output
            verbose: boolean, debug
            nexus: boolean, whether to output newick (default: False) or nexus (True) formatted tree
            string_fragment: list of characters that comprise the tree string
        """
        if cur_node==None:
            cur_node=self.root#.children[-1]
        if traits==None: ## if None
            traits=set(sum([k.traits.keys() for k in self.Objects],[])) ## fetch all trait keys
        if string_fragment==None:
            string_fragment=[]
            if nexus==True:
                assert json==False,'Nexus format not a valid option for JSON output'
                if verbose==True:
                    print('Exporting to Nexus format')
                string_fragment.append('#NEXUS\nBegin trees;\ntree TREE1 = [&R] ')
        if traverse_condition==None:
            traverse_condition=lambda k: True

        comment=[] ## will hold comment
        if len(traits)>0: ## non-empty list of traits to output
            for tr in traits: ## iterate through keys
                if tr in cur_node.traits: ## if key is available
                    if verbose==True:
                        print('trait %s available for %s (%s) type: %s'%(tr,cur_node.index,cur_node.branchType,type(cur_node.traits[tr])))
                    if isinstance(cur_node.traits[tr],str): ## string value
                        comment.append('%s="%s"'%(tr,cur_node.traits[tr]))
                        if verbose==True:
                            print('adding string comment %s'%(comment[-1]))
                    elif isinstance(cur_node.traits[tr],float) or isinstance(cur_node.traits[tr],int): ## float or integer
                        comment.append('%s=%s'%(tr,cur_node.traits[tr]))
                        if verbose==True:
                            print('adding numeric comment %s'%(comment[-1]))
                    elif isinstance(cur_node.traits[tr],list): ## lists
                        rangeComment=[]
                        for val in cur_node.traits[tr]:
                            if isinstance(val,str): ## string
                                rangeComment.append('"%s"'%(val))
                            elif isinstance(val,float) or isinstance(val,int): ## float or integer
                                rangeComment.append('%s'%(val))
                        comment.append('%s={%s}'%(tr,','.join(rangeComment)))
                        if verbose==True:
                            print('adding range comment %s'%(comment[-1]))
                elif verbose==True:
                    print('trait %s unavailable for %s (%s)'%(tr,cur_node.index,cur_node.branchType))

        if cur_node.branchType=='node':
            if verbose==True:
                print('node: %s'%(cur_node.index))
            string_fragment.append('(')
            traverseChildren=filter(traverse_condition,cur_node.children)
            assert len(traverseChildren)>0,'Node %s does not have traversable children'%(cur_node.index)
            for c,child in enumerate(traverseChildren): ## iterate through children of node if they satisfy traverse condition
                if verbose==True:
                    print('moving to child %s of node %s'%(child.index,cur_node.index))
                self.toString(cur_node=child,traits=traits,numName=numName,verbose=verbose,nexus=nexus,string_fragment=string_fragment,traverse_condition=traverse_condition)
                if (c+1)<len(traverseChildren): ## not done with children, add comma for next iteration
                    string_fragment.append(',')
            string_fragment.append(')') ## last child, node terminates

        elif cur_node.branchType=='leaf':
            if numName==False: ## if real names wanted
                assert cur_node.name!=None,'Tip does not have converted name' ## assert they have been converted
                treeName=cur_node.name ## designate real name
            elif numName==True: ## if number names wanted
                treeName=cur_node.numName ## designated numName
            if verbose==True:
                print('leaf: %s (%s)'%(cur_node.index,treeName))
            string_fragment.append("'%s'"%(treeName))

        if len(comment)>0:
            if verbose==True:
                print('adding comment to %s'%(cur_node.index))
            comment=','.join(comment)
            comment='[&'+comment+']'
            string_fragment.append('%s'%(comment)) ## end of node, add annotations

        if verbose==True:
            print('adding branch length to %s'%(cur_node.index))
        string_fragment.append(':%8f'%(cur_node.length)) ## end of node, add branch length

        if cur_node==self.root:#.children[-1]:
            string_fragment.append(';')
            if nexus==True:
                string_fragment.append('\nEnd;')
            if verbose==True:
                print('finished')
            return ''.join(string_fragment)

    def allTMRCAs(self,numName=True):
        if numName==False:
            assert len(self.tipMap)>0,'Tree does not have a translation dict for tip names'
            tip_names=[self.tipMap[k.numName] for k in self.Objects if isinstance(k,leaf)]
        else:
            tip_names=[k.numName for k in self.Objects if isinstance(k,leaf)]
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
        if verbose==True:
            print("Preparing branch hash for keeping %d branches"%(len(keep)))
        branch_hash={k.index:k for k in keep}
        embedding=[]
        if verbose==True:
            print("Deep copying tree")
        reduced_tree=copy.deepcopy(self) ## new tree object
        for k in reduced_tree.Objects: ## deep copy branches from current tree
            if k.index in branch_hash: ## if branch is designated as one to keep
                cur_b=k
                if verbose==True:
                    print("Traversing to root from %s"%(cur_b.index))
                while cur_b!=reduced_tree.root: ## descend to root
                    if verbose==True:
                        print("at %s root: %s"%(cur_b.index,cur_b==reduced_tree.root))
                    embedding.append(cur_b) ## keep track of the path to root
                    cur_b=cur_b.parent
        embedding.append(reduced_tree.root) ## add root to embedding
        if verbose==True:
            print("Finished extracting embedding")
        embedding=set(embedding) ## prune down to only unique branches

        reduced_tree.Objects=sorted(list(embedding),key=lambda x:x.height) ## assign branches that are kept to new tree's Objects
        if verbose==True:
            print("Pruning untraversed lineages")
        for k in reduced_tree.getInternal(): ## iterate through reduced tree
            k.children = [c for c in k.children if c in embedding] ## only keep children that are present in lineage traceback
        reduced_tree.root.children=[c for c in reduced_tree.root.children if c in embedding] ## do the same for root

        reduced_tree.fixHangingNodes()

        if verbose==True:
            print("Last traversal and branch sorting")
        reduced_tree.traverse_tree() ## traverse
        reduced_tree.sortBranches() ## sort

        return reduced_tree ## return new tree

    def countLineages(self,t,condition=lambda x:True):
        return len([k for k in self.Objects if k.parent.absoluteTime<t<=k.absoluteTime and condition(k)])

    def getExternal(self):
        return filter(lambda k:k.branchType=='leaf',self.Objects)

    def getInternal(self):
        return filter(lambda k:k.branchType=='node',self.Objects)

    def getBranches(self,attrs=lambda x:True):
        select=filter(attrs,self.Objects)

        if len(select)==0:
            raise Exception('No branches satisfying function were found amongst branches')
        elif len(select)==1:
            return select[-1]
        else:
            return select

    def fixHangingNodes(self):
        """
        Remove internal nodes without any children.
        """
        hangingCondition=lambda k:k.branchType=='node' and len(k.children)==0
        hangingNodes=filter(hangingCondition,self.Objects) ## check for nodes without any children (hanging nodes)
        while len(hangingNodes)>0:
            for h in sorted(hangingNodes,key=lambda x:-x.height):
                h.parent.children.remove(h) ## remove old parent from grandparent's children
                hangingNodes.remove(h) ## remove old parent from multitype nodes
                self.Objects.remove(h) ## remove old parent from all objects
            hangingNodes=filter(hangingCondition,self.Objects) ## regenerate list

    def addText(self,ax,target=lambda k:k.branchType=='leaf',position=lambda k:(k.x*1.01,k.y),text=lambda k:k.numName,zorder_function=lambda k: 101,**kwargs):
        for k in filter(target,self.Objects):
            x,y=position(k)
            z=zorder_function(k)
            ax.text(x,y,text(k),zorder=z,**kwargs)
        return ax

    def plotPoints(self,ax,x_attr=lambda k:k.height,y_attr=lambda k:k.y,target=lambda k:k.branchType=='leaf',size_function=lambda k:40,colour_function=lambda k:'k',zorder_function=lambda k: 100,**kwargs):
        for k in filter(target,self.Objects):
            y=y_attr(k) ## get y coordinates
            x=x_attr(k) ## x coordinate
            c=colour_function(k)
            size=size_function(k)
            z=zorder_function(k)
            ax.scatter(x,y,s=size,facecolor=c,edgecolor='none',zorder=z,**kwargs) ## put a circle at each tip
        return ax

    def plotTree(self,ax,type='rectangular',target=lambda k: True,x_attr=lambda k:k.height,y_attr=lambda k:k.y,branchWidth=lambda k:2,colour_function=lambda f:'k',zorder_function=lambda k: 98,**kwargs):
        assert type in ['rectangular','unrooted'],'Unrecognised drawing type "%s"'%(type)
        for k in filter(target,self.Objects): ## iterate over branches in the tree
            y=y_attr(k) ## get y coordinates
            x=x_attr(k) ## x coordinate
            xp=x_attr(k.parent) ## get parent's x
            if xp==None:
                xp=x
            c=colour_function(k)
            b=branchWidth(k)
            z=zorder_function(k)
            if type=='rectangular':
                if k.branchType=='node': ## if node...
                    yl=y_attr(k.children[0]) ## get y coordinates of first and last child
                    yr=y_attr(k.children[-1])
                    ax.plot([x,x],[yl,yr],color=c,lw=b,zorder=z,**kwargs) ## plot vertical bar connecting node to both its offspring

                ax.plot([x,xp],[y,y],color=c,lw=b,zorder=z,**kwargs) ## plot horizontal branch to parent
            elif type=='unrooted':
                yp=y_attr(k.parent)
                ax.plot([x,xp],[y,yp],color=c,lw=b,zorder=z,**kwargs)
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

        cerberus=re.match('(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/ ]+)(\'|\"|)(\[)*',data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
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

        cerberus=re.match('(\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&]+)\]',data[i:])## look for MCC comments
        if cerberus is not None:
            if verbose==True:
                print('%d comment: %s'%(i,cerberus.group(2)))
            comment=cerberus.group(2)
            numerics=re.findall('[,&][A-Za-z\_\.0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall('[,&][A-Za-z\_\.0-9]+=["|\']*[A-Za-z\_0-9\.\+ :\/\(\)\&]+["|\']*',comment) ## strings
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
            return ll
            break ## end loop

def make_treeJSON(JSONnode,json_translation,ll=None,verbose=False):
    if 'children' in JSONnode: ## only nodes have children
        new_node=node()
    else:
        new_node=leaf()
        new_node.numName=JSONnode[json_translation['name']] ## set leaf numName
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

def loadNewick(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',variableDate=True,absoluteTime=True,verbose=False):
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
            if verbose==True:
                print('Identified tree string')

    assert ll,'Regular expression failed to find tree string'
    ll.traverse_tree(verbose=verbose) ## traverse tree
    ll.sortBranches() ## traverses tree, sorts branches, draws tree

    if absoluteTime==True:
        tipDates=[]
        for k in ll.getExternal():
            n=k.numName
            k.name=k.numName
            cerberus=re.search(tip_regex,n)
            if cerberus is not None:
                tipDates.append(decimalDate(cerberus.group(1),fmt=date_fmt,variable=variableDate))

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
            if verbose==True:
                print('File should contain %d taxa'%(tipNum))

        cerberus=re.search(treestring_regex,l)
        if cerberus is not None:
            treeString_start=l.index('(')
            ll=make_tree(l[treeString_start:]) ## send tree string to make_tree function
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
        for k in ll.getExternal():
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

def loadJSON(tree_path,json_translation={'name':'strain','absoluteTime':'num_date'},json_meta=None,verbose=False,sort=True,stats=True):
    """
    Load a nextstrain JSON by providing either the path to JSON or a file handle.
    json_translation is a dictionary that translates JSON attributes to baltic branch attributes (e.g. 'absoluteTime' is called 'num_date' in nextstrain JSONs).
    Note that to avoid conflicts in setting node heights you can either define the absolute time of each node or branch lengths (e.g. if you want a substitution tree).
    """
    assert 'name' in json_translation and ('absoluteTime' in json_translation or 'length' in json_translation),'JSON translation dictionary missing entries: %s'%(', '.join([entry for entry in ['name','height','absoluteTime','length'] if (entry in json_translation)==False]))
    if verbose==True:
        print('Reading JSON')

    if isinstance(tree_path,str):
        with open(tree_path) as json_data:
            d = json.load(json_data)
            ll=make_treeJSON(d,json_translation,verbose=verbose)
    else:
        ll=make_treeJSON(json.load(tree_path),json_translation,verbose=verbose)

    assert ('absoluteTime' in json_translation and 'length' not in json_translation) or ('absoluteTime' not in json_translation and 'length' in json_translation),'Cannot use both absolute time and branch length, include only one in json_translation dictionary.'

    for attr in json_translation: ## iterate through attributes in json_translation
        for k in ll.Objects: ## for every branch
            setattr(k,attr,k.traits[json_translation[attr]]) ## set attribute value for branch

    if 'absoluteTime' in json_translation: ## if using absoluteTime need to set branch lengths for traversals
        for k in ll.Objects:
            if json_translation['absoluteTime'] in k.parent.traits:
                k.length=k.traits[json_translation['absoluteTime']]-k.parent.traits[json_translation['absoluteTime']]
            else:
                k.length=0.0

    if verbose==True:
        print('Traversing and drawing tree')

    ll.traverse_tree(verbose=verbose)
    ll.drawTree()
    if stats==True:
        ll.treeStats() ## initial traversal, checks for stats
    if sort==True:
        ll.sortBranches() ## traverses tree, sorts branches, draws tree

    if json_meta:
        if isinstance(json_meta,str):
            metadata=json.load(open(json_meta['file'],'r'))
        else:
            metadata=json.load(json_meta['file'])
        cmap=dict(metadata['color_options'][json_meta['traitName']]['color_map'])
        setattr(ll,'cmap',cmap)
    return ll

if __name__ == '__main__':
    import sys
    ll=make_tree(sys.argv[1],ll)
    ll.traverse_tree()
    sys.stdout.write('%s\n'%(ll.treeHeight))
