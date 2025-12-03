"""This module provides the core BALTIC functions ``make_tree`` and ``make_tree_JOSN``, the core parser functions which create BALTIC tree, node, and branchLike objects.

By convention, BALTIC (and this base module) are loaded under the name ``bt``.


Example
-------
After loading the module, a basic Newick formatted tree string can be read by:

>>> treeString='((A:1.0,B:2.0):1.0,C:3.0);'
>>> ll = bt.make_tree(treeString, treeType="divergence")


Notes
-----
This version of BALTIC (v0.1) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


Attributes
----------
logger : logging.logger
    Default logger which will be passed to other baltic functions.
"""
__all__ = ['make_tree', 'make_tree_JSON']
import re
import sys
import logging
from baltic.tree import Tree
from baltic.node import Node
from baltic.leaf import Leaf
from baltic.reticulation import Reticulation

logger = logging.getLogger("baltic")
sys.setrecursionlimit(9001)

def make_tree(data, treeType, tre=None):
    """Parses a tree string to create a baltic tree.

    Trees are parsed according to the following standard algorithm:

    * Every time an opening parenthesis (``(``) is encountered in the tree string a new instance of `Node` class is created. The new class' ``.index`` attribute is set to the index along the tree string where it was encountered, giving that particular class a unique identifier within the tree string. The ``.parent`` attribute is set to whatever the previous object encountered was, similarly, since the last encountered object could only be another node, the current node is added to its parents list of ``children``. Finally we set our new node as the 'current' node of the tree and append the node to the list of objects (``.Objects``, which are branches) contained in the tree.
    * Every time a string is encountered which may or may not be surrounded by quotation marks (``'`` or ``"``) or have the beginning of an annotation block (``[``) we create a new ``Leaf`` object. It also receives an ``.index`` identifier, like the ``Node`` object. Unlike the ``Node``, however, the ``.numName`` attribute is also set as the string that defined the tip. In BEAST X trees it will be the number that identifies the tip, but it could also be a regular string.
    * Next, baltic looks for annotations, which are the blocks in the format ``[&parameter1=1.0,parameter2=0.0]``. These are transformed into the ``.traits`` dictionary for the branch. In this example the branch being parsed would receive a dictionary with two keys: ``cur_branch.traits = {'parameter1' : 1.0, 'parameter2' : 0.0}``.
    * Annotations should be followed by branch lengths preceded by a colon (``:``). The branch length is assigned to the current branch's ``.length`` attribute.
    * Bifurcations (or multifurcations) in the tree string are notated by commas (``,``) and ends of clades are denoted by closing parentheses (``)``). Both mean that whatever comes next is in relation to the parent branch of whatever branch we were dealing with earlier.
    * Finally tree strings are finished with a semi colon (``;``).


    Parameters
    ----------
    data : str
        Input tree string to be parsed.
    treeType : {'time', 'divergence'}
        The type of the tree, which defines meaning of branch lengths (e.g. for a time-calibrated phylogeny whose branches are in units of subst./site/year, 'time' should be used.)
    tre : baltic.Tree, optional
        Reference to an existing tree object. If omitted, a new object will be created.

    Returns
    -------
    tre : baltic.Tree
        Output tree object.
    """
    patterns = {
        'beast_tip': r'(\(|,)([0-9]+)(\[|\:)', # Pattern to match tips in BEAST format (integers)
        'non_beast_tip': r'(\(|,)(\'|\")*([^\(\):\[\'\"#]+)(\'|\"|)*(\[)*' # Pattern to match tips with unencoded names
    }
    if not isinstance(data,str): ## tree string is not an instance of string (could be unicode) - convert
        data=str(data)

    # Add in some checks that the data are in correct format
    assert data.endswith(";"), "Improperly formatted string: must end in semicolon"
    assert data.count("(")==data.count(")"), "Improperly formatted string: must have matching parentheses"

    if tre is None: ## calling without providing a tree object - create one
        tre=Tree(treeType)
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    storedI=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop

    while i < len(data): ## while there's characters left in the tree string - loop away
        if storedI == i: logger.debug(f'{i} >{data[i]}<')

        assert (storedI != i),f'\nTree string unparseable\nStopped at >>{data[i]}<<\nstring region looks like this: {data[i:i+5000]}' # Ensure that the index has advanced; if not, raise an error indicating an unparseable string
        storedI=i # Store the current index at the end of the loop to check for infinite loops

        if data[i] == '(': ## look for new nodes
            logger.debug(f'{i} adding node')
            tre.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character

        match=re.match(patterns['beast_tip'],data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if match:
            logger.debug(f'{i} adding leaf (BEAST) {match.group(2)}')
            tre.add_leaf(i,match.group(2)) ## add tip
            i+=len(match.group(2)) ## advance in tree string by however many characters the tip is encoded

        match=re.match(patterns['non_beast_tip'],data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if match:
            logger.debug(f'{i} adding leaf (non-BEAST) {match.group(3)}')
            tre.add_leaf(i,match.group(3).strip('"').strip("'"))  ## add tip
            i+=len(match.group(3))+match.group().count("'")+match.group().count('"') ## advance in tree string by however many characters the tip is encoded

        match=re.match(r'\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if match:
            logger.debug(f'{i} adding multitype node {match.group(1)}')
            i+=len(match.group(1))

        match=re.match(r'[\(,](#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for beginning of reticulate branch
        if match:
            logger.debug(f'{i} adding outgoing reticulation branch {match.group(1)}')
            tre.add_reticulation(match.group(1)) ## add reticulate branch

            destination=None
            for k in tre.Objects: ## iterate over branches parsed so far
                if 'label' in k.traits and k.traits['label']==match.group(1): ## if there's a branch with a matching id
                    if destination is None: ## not set destination before
                        destination=k ## destination is matching node
                    else: ## destination seen before - raise an error (indicates reticulate branch ids are not unique)
                        logger.error(f'Reticulate branch not unique: {match.group(1)} seen elsewhere in the tree')
            if destination: ## identified destination of this branch
                logger.debug(f'identified {match.group(1)} destination')
                tre.curNode.target=destination ## set current node's target as the destination
                setattr(destination,"contribution",tre.curNode) ## add contributing edge to destination
            else:
                logger.debug(f'destination of {match.group(1)} not identified yet')
            i+=len(match.group())-1

        match=re.match(r'\)(#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for landing point of reticulate branch
        if match:
            logger.debug(f'{i} adding incoming reticulation branch {match.group(1)}')
            tre.curNode.traits['label']=match.group(1) ## set node label

            origin=None ## branch is landing, check if its origin was seen previously
            for k in tre.Objects: ## iterate over currently existing branches
                if isinstance(k,Reticulation) and k.name==match.group(1): ## check if any reticulate branches match the origin
                    if origin is None: ## origin not identified yet
                        origin=k ## origin is reticulate branch with the correct name
                    else: ## origin has been identified - shouldn't happen, implies that multiple reticulate branches exist with the same name
                        raise Exception(f'Reticulate branch not unique: {match.group(1)} seen elsewhere in the tree')
            if origin: ## identified origin
                logger.debug(f'identified {match.group(1)} origin')
                origin.target=tre.curNode ## set origin's landing at this node
                setattr(tre.curNode,"contribution",origin) ## add contributing edge to this node
            else:
                logger.debug(f'origin of {match.group(1)} not identified yet')
            i+=len(match.group())-1

        match=re.match(r'(\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&]+)\]',data[i:])## look for MCC comments
        if match:
            logger.debug(f'{i} comment: {match.group(2)}')
            comment=match.group(2)
            numerics=re.findall(r'[,&][A-Za-z\_\.\-0-9]+=[0-9\-Ee\.]+',comment) ## find all entries that have values as floats
            strings=re.findall(r'[,&][A-Za-z\_\.\-0-9]+=["|\']*[A-Za-z\_0-9\.\+ :\/\(\)\&\-]+[\"|\']*',comment) ## strings
            treelist=re.findall(r'[,&][A-Za-z\_\.\-0-9]+={[A-Za-z\_,{}0-9\. :\/\(\)\&]+}',comment) ## complete history logged robust counting (MCMC trees)
            sets=re.findall(r'[,&][A-Za-z\_\.\-0-9\%]+={[A-Za-z\.\-0-9eE,\"\_ :\/\(\)\&]+}',comment) ## sets and ranges
            figtree=re.findall(r'\![A-Za-z]+=[A-Za-z0-9# :\/\(\)\&]+',comment)

            for vals in strings:
                tr,val=vals.split('=') # tr := trait; val := value
                tr=tr[1:]
                if '+' in val:
                    val=val.split('+')[0] ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                tre.curNode.traits[tr]=val.strip('"')

            for vals in numerics: ## assign all parsed annotations to traits of current branch
                tr,val=vals.split('=') ## split each value by =, left side is name, right side is value
                tr=tr[1:]
                if val.replace('E','',1).replace('e','',1).replace('-','',1).replace('.','',1).isdigit():
                    tre.curNode.traits[tr]=float(val)

            for val in treelist:
                tr,val=val.split('=')
                tr=tr[1:]
                micromatch = []
                if val.count(",") == 2:
                    micromatch=re.findall(r'{([0-9\.\-e]+,[a-z_A-Z]+,[a-z_A-Z]+)}',val)
                elif val.count(",") == 3:
                    micromatch=re.findall(r'{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}',val)
                tre.curNode.traits[tr]=[]
                for val in micromatch:
                    val_split = val.split(',') #NOTE this variable is never used. remove?
                    tre.curNode.traits[tr].append(val.split(","))

            for vals in sets:
                tr,val=vals.split('=')
                tr=tr[1:]
                if 'set' in tr:
                    tre.curNode.traits[tr]=[]
                    for v in val[1:-1].split(','):
                        if 'set.prob' in tr:
                            tre.curNode.traits[tr].append(float(v))
                        else:
                            tre.curNode.traits[tr].append(v.strip('"'))
                else:
                    try:
                        tre.curNode.traits[tr]=list(map(float,val[1:-1].split(',')))
                    except:
                        logger.error(f'some other trait: {vals}')

            if len(figtree)>0:
                logger.info('FigTree comment found, ignoring')

            i+=len(match.group()) ## advance in tree string by however many characters it took to encode labels

        # match=re.match(r'([A-Za-z\_\-0-9\.]+)(\:|\;)',data[i:])## look for old school node labels
        match=re.match(r'([A-Za-z\_\-0-9\.\/]+)(\:|\;|\[)',data[i:])## look for old school node labels

        if match:
            logger.debug(f'old school comment found: {match.group(1)}')
            tre.curNode.traits['label']=match.group(1)

            i+=len(match.group(1))

        micromatch=re.match(r'(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
        if micromatch is not None:
            logger.debug(f'adding branch length ({i}) {float(micromatch.group(2))}')
            tre.curNode.length=float(micromatch.group(2)) ## set branch length of current node
            i+=len(micromatch.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            tre.curNode=tre.curNode.parent

        if data[i] == ';': ## look for string end
            return tre


def make_tree_JSON(jsonNode,jsonTranslationDict,treeType,tre=None,):
    """Creates a baltic tree object from a JSON.

    Parameters
    ----------
    jsonNode
        TODO: ADD THIS
    jsonTranslationDict
        TODO: ADD THIS
    treeType : {'time', 'divergence'}
        The type of the tree, which defines meaning of branch lengths (e.g. for a time-calibrated phylogeny whose branches are in units of subst./site/year, 'time' should be used.)
    tre : baltic.Tree, optional
        Reference to an existing tree object. If omitted, a new object will be created.

    Returns
    -------
    tre : baltic.Tree
        Output tree object.
    """
    if 'children' in jsonNode: ## only nodes have children
        newNode=Node()
    else:
        newNode=Leaf(name=jsonNode[jsonTranslationDict['name']])

    if tre is None:
        tre=Tree(treeType)
        tre.root=newNode
    if 'attr' in jsonNode:
        attr = jsonNode.pop('attr')
        jsonNode.update(attr)

    newNode.parent=tre.curNode ## set parent-child relationships
    tre.curNode.children.append(newNode)
    newNode.index=jsonNode[jsonTranslationDict['name']] ## indexing is based on name
    newNode.traits={n:jsonNode[n] for n in list(jsonNode.keys()) if n!='children'} ## set traits to non-children attributes
    tre.Objects.append(newNode)
    tre.curNode=newNode

    if 'children' in jsonNode:
        for child in jsonNode['children']:
            make_tree_JSON(child,jsonTranslationDict,tre)
            tre.curNode=tre.curNode.parent
    return tre

if __name__ == '__main__':
    pass
