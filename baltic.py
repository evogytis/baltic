import re
import copy
import datetime as dt

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

    def add_node(self,i):
        """ attaches a new node to current node """
        new_node=node() ## new node instance
        new_node.index=i ## new node's index is the position along the tree string
        new_node.parent=self.cur_node ## new node's parent is current node
        self.cur_node.children.append(new_node) ## new node is a child of current node
        self.cur_node=new_node ## current node is now new node
        self.Objects.append(self.cur_node) ## add new node to list of objects in the tree
        self.nodes.append(self.cur_node)
        
    def add_leaf(self,i,name):
        """ attaches a new leaf (tip) to current node """
        new_leaf=leaf() ## new instance of leaf object
        new_leaf.index=i ## index is position along tree string
        new_leaf.numName=name ## numName is the name tip has inside tree string, BEAST trees usually have numbers for tip names
        new_leaf.parent=self.cur_node ## leaf's parent is current node
        self.cur_node.children.append(new_leaf) ## assign leaf to parent's children
        self.cur_node=new_leaf ## current node is now new leaf
        self.Objects.append(self.cur_node) ## add leaf to all objects in the tree
        self.leaves.append(self.cur_node)

    def setAbsoluteTime(self,date):
        """ place all objects in absolute time by providing the date of the most recent tip """
        
        for i in self.Objects: ## iterate over all objects
            i.absoluteTime=date-self.treeHeight+i.height ## heights are in units of time from the root

    def treeStats(self):
        """ provide information about the tree """
        self.traverse_tree() ## traverse the tree
        obs=self.Objects ## convenient list of all objects in the tree
        print '\nTree height: %.6f\nTree length: %.6f'%(self.treeHeight,sum([x.length for x in obs])) ## report the height and length of tree
        
        nodes=self.nodes ## get all nodes
        strictlyBifurcating=False ## assume tree is not strictly bifurcating
        multiType=False
        N_children=[len(x.children) for x in nodes]
        minChildren,maxChildren=min(N_children),max(N_children) ## get the largest number of descendant branches of any node
        if maxChildren==2: ## if every node has at most two children branches
            strictlyBifurcating=True ## it's strictly bifurcating
        if minChildren==1:
            multiType=True
        print '\nTree is strictly bifurcating = %s'%(strictlyBifurcating) ## report
        print '\nTree is multitype = %s'%(multiType) ## report
        
        hasTraits=False ## assume tree has no annotations
        maxAnnotations=max([len(x.traits) for x in obs]) ## check the largest number of annotations any branch has
        if maxAnnotations>0: ## if it's more than 0
            hasTraits=True ## there are annotations
        print '\nTree has annotations = %s'%(hasTraits) ## report
        
        print '\nNumbers of objects in tree: %d (%d nodes and %d leaves)\n'%(len(obs),len(nodes),len(obs)-len(nodes)) ## report numbers of different objects in the tree
            
    def traverse_tree(self,startNode=None,include_all=False,verbose=False):
        """ Traverses tree from root. If a starting node is not defined begin traversal from root. By default returns a list of leaf objects that have been visited, optionally returns a list of all objects in the tree. """
        if startNode==None: ## if no starting point defined - start from root
            cur_node=self.root
            startNode=cur_node      
        elif startNode.branchType=='leaf':
            return [startNode.numName]
        else:
            cur_node=startNode ## otherwise start from starting node
        
        if verbose==True:
            print 'Verbose traversal initiated'
        
        for k in self.Objects: ## reset leaves and number of children for every node
            if isinstance(k,node):
                k.leaves=[]
                k.numChildren=0
        
        seen=[] ## remember what's been visited
        collected=[] ## collect leaf objects along the way
        maxHeight=0 ## check what the maximum distance between the root and the most recent tip is
        height=0.0 ## begin at height 0.0
        root=False ## root becomes true once you're unable to go back any further

        while root==False: ## cycle indefinitely as long as there's more of the tree to explore
            if seen.count(cur_node.index)>3:
                if verbose==True:
                    print 'node %s seen too many times, breaking'%(cur_node.index)
                root=True
                break
                
            if isinstance(cur_node,node): ## if currently dealing with a node
                if verbose==True:
                    print 'encountered node %s'%(cur_node.index)
                if include_all==True: ## if node hasn't been collected and we want to collect all objects - add it for reporting later
                    collected.append(cur_node)
                ## check whether all children of the node have been seen or whether we're not currently at the root
                ## as long as all children have been seen will descend downwards
                while sum([1 if x.index in seen else 0 for x in cur_node.children])==len(cur_node.children) and cur_node!=startNode:
                    if verbose==True:
                        print 'seen all children of node %s'%(cur_node.index)
                    ## check wheteher current node's most recent child is higher than the known highest point in the tree
                    if cur_node.childHeight <= highestTip:
                        cur_node.childHeight = highestTip
                    elif cur_node.childHeight > highestTip:
                        highestTip = cur_node.childHeight
                    
                    if cur_node.parent.index==startNode.index: ## if currently at root - set the root flag to True and break the while loop
                        if verbose==True:
                            print 'reached starting point %s (parent %s)'%(cur_node.index,cur_node.parent.index)
                        cur_node.height=height
                        root=True
                        break
                    else: ## otherwise...
                        if verbose==True:
                            print 'heading to parent of %s'%(cur_node.index)
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
                    print 'visiting %s next (child of %s)'%(cur_node.children[idx_to_visit].index,cur_node.index)
                
                cur_node.height=height ## set height of current node
                height+=float(cur_node.children[idx_to_visit].length) ## prepare for heading towards the previously unvisited child branch
                cur_node=cur_node.children[idx_to_visit] ## set current node to unvisited child
                seen.append(cur_node.index) ## remember that the child has now been visited
                
            elif isinstance(cur_node,leaf): ## node dealing with node, are we dealing with a leaf (tip)?
                if verbose==True:
                    print 'encountered leaf %s (%s or %s)'%(cur_node.index,cur_node.numName,cur_node.name)
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
        """ give each tip its correct label using a dictionary """
        for k in self.leaves: ## iterate through leaf objects in tree
            k.name=d[k.numName] ## change its name
    
    def sortBranches(self,descending=True):
        """ sort descendants of each node """
        if descending==True:
            modifier=1 ## define the modifier for sorting function later
        elif descending==False:
            modifier=-1

        for k in self.nodes: ## iterate over nodes
            ## split node's offspring into nodes and leaves, sort each list individually
            nodes=sorted([x for x in k.children if isinstance(x,node)],key=lambda q:(-len(q.leaves)*modifier,q.length*modifier))
            leaves=sorted([x for x in k.children if isinstance(x,leaf)],key=lambda q:q.length*modifier)
            
            if modifier==1: ## if sorting one way - nodes come first, leaves later
                k.children=nodes+leaves
            elif modifier==-1: ## otherwise sort the other way
                k.children=leaves+nodes
                    
        self.drawTree() ## update x and y positions of each branch, since y positions will have changed because of sorting
        
    def drawTree(self):
        """ find x and y coordinates of each branch """
        order=[x.numName for x in self.traverse_tree() if x.branchType=='leaf'] ## order is a list of tips recovered from a tree traversal to make sure they're plotted in the correct order along the vertical tree dimension
        
        for k in self.Objects: ## reset coordinates for all objects
            k.x=None
            k.y=None
        
        storePlotted=0
        drawn=[] ## drawn keeps track of what's been drawn
        while len(drawn)!=len(self.Objects): # keep drawing the tree until everything is drawn
            for k in [x for x in self.Objects if x.index not in drawn]: ## iterate through objects that have not been drawn
                if isinstance(k,leaf): ## if leaf - get position of leaf, draw branch connecting tip to parent node
                    x=k.height ## x position is height
                    y=order.index(k.numName) ## y position of leaf is given by the order in which tips were visited during the traversal
                    k.x=x ## set x and y coordinates
                    k.y=y
                    drawn.append(k.index) ## remember that this objects has been drawn

                if isinstance(k,node): ## if parent is non-root node and y positions of all its children are known
                    if len([q.y for q in k.children if q.y!=None])==len(k.children):
                        x=k.height ## x position is height              
                        children_y_coords=[q.y for q in k.children if q.y!=None] ## get all existing y coordinates of the node
                        y=sum(children_y_coords)/float(len(children_y_coords)) ## internal branch is in the middle of the vertical bar
                        k.x=x
                        k.y=y
                        drawn.append(k.index) ## remember that this objects has been drawn
                        
            assert len(drawn)>storePlotted,'Got stuck trying to find y positions of objects'
            storePlotted=len(drawn)

    def traverseWithinTrait(self,startNode,traitName,converterDict=None):
        """ Traverse the tree staying within a particular trait value """
        cur_node=startNode
        seen=[]
        collected=[]

        if converterDict==None: ## no dictionary to convert trait values, trying to stay inside trait value of starting node
            stayWithin=startNode.traits[traitName]
        else: ## there's a dictionary, trying to stay within same dictionary value
            stayWithin=converterDict[startNode.traits[traitName]]
 
        root=False
        
        if isinstance(startNode,leaf): ## quite immediatelly if starting from a leaf - nowhere to go
            collected.append(startNode)
            return collected

        while root==False:
            if isinstance(cur_node,node):
                ## if all children have been seen and not at root
                while sum([1 if child.index in seen else 0 for child in cur_node.children])==len(cur_node.children) and cur_node.parent!=self.root:
                    if converterDict==None: ## if about to traverse into different trait state - break while loop
                        curTrait=cur_node.traits[traitName]
                    else:
                        curTrait=converterDict[cur_node.traits[traitName]]
                    
                    if curTrait!=stayWithin:
                        root=True
                        return collected

                    if cur_node.index==startNode.index or cur_node.index=='Root': ## if back at starting node or the root - break while loop
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
                if cur_node.index==startNode.index or cur_node.parent.index=='Root': ## if back at starting node or the root - break while loop
                    root=True
                    return collected
                else:
                    cur_node=cur_node.parent
            
    def collapseNodes(self,trait,cutoff):
        """ Collapse all nodes whose trait value is below the cutoff value """
        newTree=copy.deepcopy(self) ## work on a copy of the tree
        nodes=newTree.nodes ## fetch a list of all nodes
        zero_count=sum([1 if q.traits[trait]<cutoff else 0 for q in nodes]) ## count how many branches fail the test

        assert zero_count<len(nodes),'Chosen cutoff would remove all branches'
        
        while zero_count>0: ## as long as there are branches to be collapsed - keep reducing the tree
            for k in sorted(nodes,key=lambda x:x.height): ## start with branches near the tips
                if k.traits[trait]<=cutoff: ## if trait value is less than cutoff - remove it cleanly
                    zero_node=k.children ## fetch the node's children
                    k.parent.children+=zero_node ## add them to the zero node's parent
                    
                    old_parent=k ## node to be deleted is the old parent
                    new_parent=k.parent ## once node is deleted, the parent to all their children will be the parent of the deleted node
                    
                    for w in newTree.Objects: ## assign the parent of zero node as the parent to any children of zero node
                        if w.parent==old_parent:
                            w.parent=new_parent
                            w.length+=old_parent.length

                    k.parent.children.remove(k) ## remove traces of zero node - it doesn't exist as a child, doesn't exist in the tree and doesn't exist in the nodes list
                    newTree.Objects.remove(k)
                    nodes.remove(k)

                    zero_count-=1 ## one more zero node taken care of
                    
        newTree.sortBranches() ## sort the tree to traverse, draw and sort tree to adjust y coordinates
        return newTree ## return collapsed tree
    
    def toString(self):
        """ Output the topology of the tree with branch lengths to string """
        cur_node=self.root
        seen=[]
        tree_string=[]
        root=False

        while root==False:
            if isinstance(cur_node,node):
                ## if all children have been seen and not at root
                while sum([1 if x.index in seen else 0 for x in cur_node.children])==len(cur_node.children) and cur_node!=self.root:
                    if cur_node.parent.index=='Root':
                        root=True
                        break
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
                    
                    if isinstance(cur_node,node):
                        tree_string.append('(') ## dealing with new node, add (
                    else:
                        tree_string.append('\'%s\':%8f'%(cur_node.name,cur_node.length)) ## dealing with tip, write out name, add branch length
                    
                else: ## all children seen, clade's end
                    tree_string.append('):%8f'%(cur_node.length))
                    cur_node=cur_node.parent
                
            elif isinstance(cur_node,leaf):
                seen.append(cur_node.index)
                if cur_node.parent.index=='Root':
                    root=True
                    break
                else:
                    cur_node=cur_node.parent
                
        return ''.join(tree_string)+';' ## the coup de grace

    def allTMRCAs(self):
        tip_names=[k.numName for k in self.Objects if isinstance(k,leaf)]
        tmrcaMatrix={x:{y:0.0 for y in tip_names} for x in tip_names}
        #print tmrcaMatrix
        #print tip_names
        for k in self.Objects:
            if isinstance(k,node):
                all_children=k.leaves
                #print all_children
                for x in range(0,len(all_children)-1):
                    for y in range(x+1,len(all_children)):
                        tipA=all_children[x]
                        tipB=all_children[y]
                        if tmrcaMatrix[tipA][tipB]<=k.absoluteTime:
                            tmrcaMatrix[tipA][tipB]=k.absoluteTime
                            tmrcaMatrix[tipB][tipA]=k.absoluteTime
        return tmrcaMatrix

def make_tree(data,ll,verbose=False):
    """
    data is a tree string, ll (LL) is an instance of a tree object
    """
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    stored_i=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop
    
    while i < len(data): ## while there's characters left in the tree string - loop away
        if stored_i == i and verbose==True:
            print '%d >%s<'%(i,data[i])
        
        assert (stored_i != i),'\nTree string unparseable\nStopped at >>%s<<\nstring region looks like this: %s'%(data[i],data[i:i+5000]) ## make sure that you've actually parsed something last time, if not - there's something unexpected in the tree string
        stored_i=i ## store i for later
        
        if data[i] == '(': ## look for new nodes
            if verbose==True:
                print '%d adding node'%(i)
            ll.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character
            
        cerberus=re.match('(\(|,)([0-9]+)(\[|\:)',data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if cerberus is not None:
            if verbose==True:
                print '%d adding leaf (BEAST) %s'%(i,cerberus.group(2))
            ll.add_leaf(i,cerberus.group(2)) ## add tip
            i+=len(cerberus.group(2)) ## advance in tree string by however many characters the tip is encoded
            
        cerberus=re.match('(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/]+)(\'|\"|)(\[)*',data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if cerberus is not None:
            if verbose==True:
                print '%d adding leaf (non-BEAST) %s'%(i,cerberus.group(3))
            ll.add_leaf(i,cerberus.group(3).strip('"').strip("'"))  ## add tip
            i+=len(cerberus.group(3))+cerberus.group().count("'")+cerberus.group().count('"') ## advance in tree string by however many characters the tip is encoded

        cerberus=re.match('\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if cerberus is not None:
            if verbose==True:
                print '%d adding multitype node %s'%(i,cerberus.group(1))
            i+=len(cerberus.group(1))

        cerberus=re.match('(\:)*\[&([A-Za-z\_\-{}\,0-9\.\%=\"\+]+)\]',data[i:])## look for MCC comments
        #cerberus=re.match('\[&[A-Za-z\_\-{}\,0-9\.\%=\"\+]+\]',data[i:])## look for MCC comments
        if cerberus is not None:
            if verbose==True:
                print '%d comment: %s'%(i,cerberus.group(2))
            comment=cerberus.group(2)
            numerics=re.findall('[A-Za-z\_\.0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall('[A-Za-z\_\.0-9]+="[A-Za-z\_0-9\.\+]+"',comment) ## strings
            treelist=re.findall('[A-Za-z\_\.0-9]+={[A-Za-z\_,{}0-9\.]+}',comment) ## complete history logged robust counting (MCMC trees)
            sets=re.findall('[A-Za-z\_\.0-9\%]+={[A-Za-z\.\-0-9eE,\"\_]+}',comment) ## sets and ranges
            ## ranges not included - who needs those anyway?
            
            for vals in numerics: ## assign all parsed annotations to traits of current branch
                tr,val=vals.split('=') ## split each value by =, left side is name, right side is value
                ll.cur_node.traits[tr]=float(val)

            for vals in strings:
                tr,val=vals.split('=')
                if '+' in val:
                    val=val.split('+')[0] ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                ll.cur_node.traits[tr]=val.strip('"')

            for val in treelist:
                tr,val=val.split('=')
                microcerberus=re.findall('{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}',val)
                ll.cur_node.traits[tr]=[]
                for val in microcerberus:
                    codon,timing,start,end=val.split(',')
                    ll.cur_node.traits[tr].append((int(codon),float(timing),start,end))
                    
            for vals in sets:
                tr,val=vals.split('=')
                if 'set' in tr:
                    ll.cur_node.traits[tr]=[]
                    for v in val[1:-1].split(','):
                        if 'set.prob' in tr:
                            ll.cur_node.traits[tr].append(float(v))
                        else:
                            ll.cur_node.traits[tr].append(v.strip('"'))
                elif 'range' in tr or 'HPD' in tr:
                    ll.cur_node.traits[tr]=map(float,val[1:-1].split(','))
                else:
                    print 'some other trait: %s'%(vals)
            
            i+=len(cerberus.group()) ## advance in tree string by however many characters it took to encode labels
            
        microcerberus=re.match('(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
        if microcerberus is not None:
            if verbose==True:
                print 'adding branch length (%d) %.6f'%(i,float(microcerberus.group(2)))
            ll.cur_node.length=float(microcerberus.group(2)) ## set branch length of current node
            i+=len(microcerberus.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            ll.cur_node=ll.cur_node.parent

        if data[i] == ';': ## look for string end
            break ## end loop

def loadNexus(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',absoluteTime=True):
    tipFlag=False
    tips={}
    tipNum=0
    for line in open(tree_path,'r'):
        l=line.strip('\n')

        cerberus=re.search('dimensions ntax=([0-9]+);',l.lower())
        if cerberus is not None:
            tipNum=int(cerberus.group(1))

        cerberus=re.search('tree [A-Za-z\_]+([0-9]+) = \[&R\] ',l)
        if cerberus is not None:
            treeString_start=l.index('(')
            ll=tree() ## new instance of tree
            make_tree(l[treeString_start:],ll) ## send tree string to make_tree function

        if tipFlag==True:
            cerberus=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if cerberus is not None:
                tips[cerberus.group(1)]=cerberus.group(2).strip('"').strip("'")
            elif ';' not in l:
                print 'tip not captured by regex:',l.replace('\t','')

        if 'translate' in l.lower():
            tipFlag=True
        if ';' in l:
            tipFlag=False

    ll.treeStats() ## initial traversal, checks for stats
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
                tipDates.append(decimalDate(cerberus.group(1)))
    
        highestTip=max(tipDates)
        ll.setAbsoluteTime(highestTip)
    
    return ll

if __name__ == '__main__':
    import sys
    ll=tree()
    make_tree(sys.argv[1],ll)
    ll.traverse_tree()
    sys.stdout.write('%s\n'%(ll.treeHeight))