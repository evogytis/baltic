__all__ = ['make_tree', 'make_treeJSON']

import re,sys
from .tree import Tree
from .node import Node
from .leaf import Leaf
from .reticulation import Reticulation

sys.setrecursionlimit(9001)

def make_tree(data,ll=None,verbose=False):
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
    patterns = {
        'beast_tip': r'(\(|,)([0-9]+)(\[|\:)', # Pattern to match tips in BEAST format (integers)
        'non_beast_tip': r'(\(|,)(\'|\")*([^\(\):\[\'\"#]+)(\'|\"|)*(\[)*' # Pattern to match tips with unencoded names
    }
    if isinstance(data,str)==False: ## tree string is not an instance of string (could be unicode) - convert
        data=str(data)

    # Add in some checks that the data are in correct format
    assert data.endswith(";"), "Improperly formatted string: must end in semicolon"
    assert data.count("(")==data.count(")"), "Improperly formatted string: must have matching parentheses"

    if ll==None: ## calling without providing a tree object - create one
        ll=Tree()
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    stored_i=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop

    while i < len(data): ## while there's characters left in the tree string - loop away
        if stored_i == i and verbose==True: print('%d >%s<'%(i,data[i]))

        assert (stored_i != i),'\nTree string unparseable\nStopped at >>%s<<\nstring region looks like this: %s'%(data[i],data[i:i+5000]) # Ensure that the index has advanced; if not, raise an error indicating an unparseable string
        stored_i=i # Store the current index at the end of the loop to check for infinite loops

        if data[i] == '(': ## look for new nodes
            if verbose==True: print('%d adding node'%(i))
            ll.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character

        match=re.match(patterns['beast_tip'],data[i-1:i+100]) ## look for tips in BEAST format (integers).
        if match:
            if verbose==True: print('%d adding leaf (BEAST) %s'%(i,match.group(2)))
            ll.add_leaf(i,match.group(2)) ## add tip
            i+=len(match.group(2)) ## advance in tree string by however many characters the tip is encoded

        match=re.match(patterns['non_beast_tip'],data[i-1:i+200])  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if match:
            if verbose==True: print('%d adding leaf (non-BEAST) %s'%(i,match.group(3)))
            ll.add_leaf(i,match.group(3).strip('"').strip("'"))  ## add tip
            i+=len(match.group(3))+match.group().count("'")+match.group().count('"') ## advance in tree string by however many characters the tip is encoded

        match=re.match(r'\)([0-9]+)\[',data[i-1:i+100]) ## look for multitype tree singletons.
        if match:
            if verbose==True: print('%d adding multitype node %s'%(i,match.group(1)))
            i+=len(match.group(1))

        match=re.match(r'[\(,](#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for beginning of reticulate branch
        if match:
            if verbose==True: print('%d adding outgoing reticulation branch %s'%(i,match.group(1)))
            ll.add_reticulation(match.group(1)) ## add reticulate branch

            destination=None
            for k in ll.Objects: ## iterate over branches parsed so far
                if 'label' in k.traits and k.traits['label']==match.group(1): ## if there's a branch with a matching id
                    if destination==None: ## not set destination before
                        destination=k ## destination is matching node
                    else: ## destination seen before - raise an error (indicates reticulate branch ids are not unique)
                        raise Exception('Reticulate branch not unique: %s seen elsewhere in the tree'%(match.group(1)))
            if destination: ## identified destination of this branch
                if verbose==True: print('identified %s destination'%(match.group(1)))
                ll.cur_node.target=destination ## set current node's target as the destination
                setattr(destination,"contribution",ll.cur_node) ## add contributing edge to destination
            else:
                if verbose==True: print('destination of %s not identified yet'%(match.group(1)))
            i+=len(match.group())-1

        match=re.match(r'\)(#[A-Za-z0-9]+)',data[i-1:i+200]) ## look for landing point of reticulate branch
        if match:
            if verbose==True: print('%d adding incoming reticulation branch %s'%(i,match.group(1)))
            ll.cur_node.traits['label']=match.group(1) ## set node label

            origin=None ## branch is landing, check if its origin was seen previously
            for k in ll.Objects: ## iterate over currently existing branches
                if isinstance(k,Reticulation) and k.name==match.group(1): ## check if any reticulate branches match the origin
                    if origin == None: ## origin not identified yet
                        origin=k ## origin is reticulate branch with the correct name
                    else: ## origin has been identified - shouldn't happen, implies that multiple reticulate branches exist with the same name
                        raise Exception('Reticulate branch not unique: %s seen elsewhere in the tree'%(match.group(1)))
            if origin: ## identified origin
                if verbose==True: print('identified %s origin'%(match.group(1)))
                origin.target=ll.cur_node ## set origin's landing at this node
                setattr(ll.cur_node,"contribution",origin) ## add contributing edge to this node
            else:
                if verbose==True: print('origin of %s not identified yet'%(match.group(1)))
            i+=len(match.group())-1

        match=re.match(r'(\:)*\[(&[A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&]+)\]',data[i:])## look for MCC comments
        if match:
            if verbose==True: print('%d comment: %s'%(i,match.group(2)))
            comment=match.group(2)
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
                micromatch = []
                if val.count(",") == 2:
                    micromatch=re.findall(r'{([0-9\.\-e]+,[a-z_A-Z]+,[a-z_A-Z]+)}',val)
                elif val.count(",") == 3:
                    micromatch=re.findall(r'{([0-9]+,[0-9\.\-e]+,[A-Z]+,[A-Z]+)}',val)
                ll.cur_node.traits[tr]=[]
                for val in micromatch:
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

            i+=len(match.group()) ## advance in tree string by however many characters it took to encode labels

        # match=re.match(r'([A-Za-z\_\-0-9\.]+)(\:|\;)',data[i:])## look for old school node labels
        match=re.match(r'([A-Za-z\_\-0-9\.]+)(\:|\;|\[)',data[i:])## look for old school node labels

        if match:
            if verbose==True: print('old school comment found: %s'%(match.group(1)))
            ll.cur_node.traits['label']=match.group(1)

            i+=len(match.group(1))

        micromatch=re.match(r'(\:)*([0-9\.\-Ee]+)',data[i:i+100]) ## look for branch lengths without comments
        if micromatch is not None:
            if verbose==True: print('adding branch length (%d) %.6f'%(i,float(micromatch.group(2))))
            ll.cur_node.length=float(micromatch.group(2)) ## set branch length of current node
            i+=len(micromatch.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            ll.cur_node=ll.cur_node.parent

        if data[i] == ';': ## look for string end
            return ll
            break ## end loop

def make_treeJSON(JSONnode,json_translation,ll=None,verbose=False):
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
    if 'children' in JSONnode: ## only nodes have children
        new_node=Node()
    else:
        new_node=Leaf()
        new_node.name=JSONnode[json_translation['name']] ## set leaf name to be the same

    if ll is None:
        ll=Tree()
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

if __name__ == '__main__':
    import sys
    ll=make_tree(sys.argv[1],ll)
    ll.traverse_tree()
    sys.stdout.write('%s\n'%(ll.treeHeight))
