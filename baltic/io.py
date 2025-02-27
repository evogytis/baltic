import re,json
from .baltic import make_tree, make_treeJSON
from .bt_utils import decimalDate
from .tree import Tree
from .node import Node
from .leaf import Leaf

def loadNewick(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',variableDate=True,absoluteTime=False,verbose=False, sortBranches = True):
    """
    Load a tree from a Newick file and process it.
    
    Parameters:
    tree_path (str or file-like object): The path to the Newick file or a file-like object containing the Newick formatted tree.
    tip_regex (str): A regular expression to extract dates from tip names. Default is '\|([0-9]+\-[0-9]+\-[0-9]+)'.
    date_fmt (str): The date format for the extracted dates. Default is '%Y-%m-%d'.
    variableDate (bool): If True, allows for variable date formats. Default is True.
    absoluteTime (bool): If True, converts the tree to absolute time using the tip dates encoded in tip names. Default is False.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    sortBranches (bool): If True, sorts the branches of the tree after loading. Default is True.
    
    Returns:
    tree: The tree object created from the Newick file.
    
    Raises:
    AssertionError: If the tree string cannot be found or if tip dates cannot be extracted when absoluteTime is True.
    
    Example:
    >>> tree = loadNewick("path/to/tree.newick", absoluteTime=False, verbose=True)
    
    Docstring generated with ChatGPT 4o.
    """
    ll=None

    handle = open(tree_path, 'r') if isinstance(tree_path, str) else tree_path

    for line in handle:
        l=line.strip('\n')
        if '(' in l:
            treeString_start=l.index('(')
            ll=make_tree(l[treeString_start:],verbose=verbose) ## send tree string to make_tree function
            if verbose==True: print('Identified tree string')

    assert ll,'Regular expression failed to find tree string'
    ll.traverse_tree(verbose=verbose) ## traverse tree
    
    if sortBranches: ll.sortBranches() ## traverses tree, sorts branches, draws tree

    if absoluteTime==True:
        tip_dates=[]
        tip_names=[]
        for k in ll.getExternal():
            tip_names.append(k.name)
            match=re.search(tip_regex,k.name)
            if match:
                tip_dates.append(decimalDate(match.group(1),fmt=date_fmt,variable=variableDate))
        assert len(tip_dates)>0,'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s'%(tip_names[0],tip_regex,date_fmt)
        ll.setAbsoluteTime(max(tip_dates))

    if isinstance(tree_path,str):
        handle.close()
    return ll

def loadNexus(tree_path,tip_regex='\|([0-9]+\-[0-9]+\-[0-9]+)',date_fmt='%Y-%m-%d',treestring_regex='tree [A-Za-z\_]+([0-9]+)',variableDate=True,absoluteTime=True,verbose=False, sortBranches=True):
    """
    Load a tree from a Nexus file and process it.
    
    Parameters:
    tree_path (str or file-like object): The path to the Nexus file or a file-like object containing the NEXUS formatted tree.
    tip_regex (str): A regular expression to extract dates from tip names. Default is '\|([0-9]+\-[0-9]+\-[0-9]+)'.
    date_fmt (str): The date format for the extracted dates. Default is '%Y-%m-%d'.
    treestring_regex (str): A regular expression to identify the tree string in the NEXUS file. Default is 'tree [A-Za-z\_]+([0-9]+)'.
    variableDate (bool): If True, allows for variable date formats. Default is True.
    absoluteTime (bool): If True, converts the tree to absolute time using the tip dates extracted from tip names. Default is True.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    sortBranches (bool): If True, sorts the branches of the tree after loading. Default is True.
    
    Returns:
    tree: The tree object created from the NEXUS file.
    
    Raises:
    AssertionError: If the tree string cannot be found or if tip dates cannot be extracted when absoluteTime is True.
    
    Example:
    >>> tree = loadNexus("path/to/tree.nexus", absoluteTime=True, verbose=True)
    
    Docstring generated with ChatGPT 4o.
    """
    tip_flag=False
    tips={}
    tip_num=0
    ll=None

    handle = open(tree_path, 'r') if isinstance(tree_path, str) else tree_path

    for line in handle:
        l=line.strip('\n')

        match=re.search('Dimensions ntax=([0-9]+);',l)
        if match:
            tip_num=int(match.group(1))
            if verbose==True: print('File should contain %d taxa'%(tip_num))

        match=re.search(treestring_regex,l)
        if match:
            treeString_start=l.index('(')
            ll=make_tree(l[treeString_start:],verbose=verbose) ## send tree string to make_tree function
            if verbose==True: print('Identified tree string')

        if tip_flag:
            match=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if match:
                tips[match.group(1)]=match.group(2).strip('"').strip("'")
                if verbose==True: print('Identified tip translation %s: %s'%(match.group(1),tips[match.group(1)]))
            elif ';' not in l:
                print('tip not captured by regex:',l.replace('\t',''))

        if 'Translate' in l:
            tip_flag=True
        if ';' in l:
            tip_flag=False

    assert ll,'Failed to find tree string using regular expression'
    ll.traverse_tree() ## traverse tree
    if sortBranches:
        ll.sortBranches() ## traverses tree, sorts branches, draws tree
    if len(tips)>0:
        ll.renameTips(tips) ## renames tips from numbers to actual names
        ll.tipMap=tips
    if absoluteTime==True:
        tip_dates=[]
        tip_names=[]
        for k in ll.getExternal():
            tip_names.append(k.name)
            match=re.search(tip_regex,k.name)
            if match:
                tip_dates.append(decimalDate(match.group(1),fmt=date_fmt,variable=variableDate))

        assert len(tip_dates)>0,'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: %s\nDate regex set to: %s\nExpected date format: %s'%(tip_names[0],tip_regex,date_fmt)
        ll.setAbsoluteTime(max(tip_dates))

    if isinstance(tree_path,str):
        handle.close()
    return ll

def loadJSON(json_object,json_translation={'name':'name','absoluteTime':'num_date'},verbose=False,sort=True,stats=True):
    """
    Load a Nextstrain JSON file and create a tree object.
    
    Parameters:
    json_object (str or dict): The path to the JSON file, a URL to a Nextstrain JSON, or a JSON object.
    json_translation (dict): A dictionary that translates JSON attributes to tree attributes (e.g., baltic branch attribute 'absoluteTime' is called 'num_date' in Nextstrain JSONs).
                             Default is {'name': 'name', 'absoluteTime': 'num_date'}.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    sort (bool): If True, sorts the branches of the tree after loading. Default is True.
    stats (bool): If True, calculates tree statistics after loading. Default is True.
    
    Returns:
    tuple: A tuple containing the tree object created from the JSON and the metadata from the JSON.
    
    Raises:
    AssertionError: If the json_translation dictionary is missing the `name` attribute (crucial for tips) and least one branch length attribute (`absoluteTime`, `length` or `height`).
    KeyError: If a string attribute in json_translation is not found in the JSON data structure.
    AttributeError: If an attribute in json_translation is neither a string nor callable.
    
    Example:
    >>> tree, metadata = loadJSON("path/to/tree.json", verbose=True)
    
    Docstring generated with ChatGPT 4o.
    """
    length_keys = ['absoluteTime', 'length', 'height']
    assert 'name' in json_translation and any(key in json_translation for key in required_keys),'JSON translation dictionary missing entries: %s'%(', '.join([entry for entry in ['name']+length_keys if (entry in json_translation)==False]))
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