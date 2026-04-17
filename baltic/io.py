"""This module provides the ``baltic`` input and output functions.

**Notes**

This version of ``baltic`` (v0.1) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""
import re
import json
from io import BytesIO as csio
import logging
import requests
from baltic.baltic import make_tree, make_tree_JSON
from baltic.bt_utils import calendar_to_decimal_date

logger = logging.getLogger("baltic.io")

def load_newick(treePath,
                treeType,
                tipRegex=r'\|([0-9\-]+)$',
                dateFmt='%Y-%m-%d',
                variableDate=True,
                absoluteTime=False,
                sortBranches=True,
                setNodes=False):
    """
    Load a tree from a Newick file or file-like object.

    **Parameters**

    treePath : str or file-like
        Path to a Newick file or an open handle containing tree text.

    treeType : {'divergence', 'time'}
        Interpretation of branch lengths in the parsed tree.

    tipRegex : str, optional
        Regular expression used to extract tip dates from leaf names.

    dateFmt : str, optional
        Date format used to parse the value captured by *tipRegex*.

    variableDate : bool, optional
        Whether partially specified tip dates should be interpreted as date
        ranges.

    absoluteTime : bool, optional
        If ``True``, assign absolute times from the tip labels.

    sortBranches : bool, optional
        If ``True``, sort branches after parsing.

    setNodes : bool, optional
        If ``True``, propagate absolute times to internal nodes when date
        information is available.

    **Returns**

    :class:`baltic.tree.Tree`
        Parsed tree object.

    **Examples**

    >>> import io
    >>> import baltic as bt
    >>> handle = io.StringIO("((A:1.0,B:1.5):0.5,C:2.0);")
    >>> ll = bt.io.load_newick(handle, treeType="divergence")
    >>> len(ll.get_external())
    3
    """
    ll = None

    handle = open(treePath, 'r') if isinstance(treePath, str) else treePath

    for line in handle:
        l = line.strip('\n')
        if '(' in l:
            treeStringStart = l.index('(')
            ll = make_tree(l[treeStringStart:], treeType) ## send tree string to make_tree function
            logger.debug('Identified tree string')

    assert ll, 'Regular expression failed to find tree string'

    if sortBranches: ll.sort_branches() ## traverses tree, sorts branches, draws tree

    if absoluteTime:
        process_tip_dates(tree=ll, tipRegex=tipRegex, dateFmt=dateFmt, variableDate=variableDate, setNodes=setNodes)

    if isinstance(treePath, str):
        handle.close()

    ll.traverse_tree() ## traverse tree
    return ll

def load_nexus(treePath,
                treeType,
                tipRegex=r'\|([0-9\-]+)$',
                dateFmt='%Y-%m-%d',
                treestringRegex=r'tree [A-Za-z\_]+([0-9]+)',
                variableDate=True,
                absoluteTime=True,
                sortBranches=True,
                setNodes=True):
    """
    Load a tree from a Nexus file or file-like object.

    **Parameters**

    treePath : str or file-like
        Path to a Nexus file or an open handle containing Nexus content.

    treeType : {'divergence', 'time'}
        Interpretation of branch lengths in the parsed tree.

    tipRegex : str, optional
        Regular expression used to extract tip dates from translated tip names.

    dateFmt : str, optional
        Date format used to parse the value captured by *tipRegex*.

    treestringRegex : str, optional
        Regular expression used to identify the tree line in the Nexus file.

    variableDate : bool, optional
        Whether partially specified tip dates should be interpreted as date
        ranges.

    absoluteTime : bool, optional
        If ``True``, assign absolute times from the tip labels.

    sortBranches : bool, optional
        If ``True``, sort branches after parsing.

    setNodes : bool, optional
        If ``True``, propagate absolute times to internal nodes when date
        information is available.

    **Returns**

    :class:`baltic.tree.Tree`
        Parsed tree object.

    **Examples**

    >>> import io
    >>> import baltic as bt
    >>> nexus = io.StringIO(
    ...     "#NEXUS\\n"
    ...     "Begin trees;\\n"
    ...     "Translate\\n"
    ...     "    1 A|2020-01-01,\\n"
    ...     "    2 B|2020-02-01;\\n"
    ...     "tree TREE1 = [&R] (1:0.1,2:0.2);\\n"
    ...     "End;\\n"
    ... )
    >>> ll = bt.io.load_nexus(nexus, treeType="time")
    >>> sorted(tip.name for tip in ll.get_external())
    ['A|2020-01-01', 'B|2020-02-01']
    """
    tip_flag = False
    tips = {}
    tip_num = 0
    ll = None

    handle = open(treePath, 'r') if isinstance(treePath, str) else treePath

    for line in handle:
        l = line.strip('\n')

        match = re.search('Dimensions ntax=([0-9]+);',l)
        if match:
            tip_num = int(match.group(1))
            logger.debug(f'File should contain {tip_num} taxa')

        match = re.search(treestringRegex,l)
        if match:
            treeString_start = l.index('(')
            ll = make_tree(l[treeString_start:], treeType) ## send tree string to make_tree function
            logger.debug('Identified tree string')

        if tip_flag:
            match = re.search(r'([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if match:
                tips[match.group(1)] = match.group(2).strip('"').strip("'")
                logger.debug(f'Identified tip translation {match.group(1)}: {tips[match.group(1)]}')
            elif ';' not in l:
                print('tip not captured by regex:', l.replace('\t',''))

        if 'Translate' in l:
            tip_flag = True
        if ';' in l:
            tip_flag = False

    assert ll,'Failed to find tree string using regular expression'

    if sortBranches: ll.sort_branches() ## traverses tree, sorts branches, draws tree

    if len(tips) > 0:
        ll.rename_tips(tips) ## renames tips from numbers to actual names
        ll.tipMap = tips
    if absoluteTime:
        process_tip_dates(tree = ll, tipRegex = tipRegex, dateFmt = dateFmt, variableDate = variableDate)

    if isinstance(treePath,str):
        handle.close()

    if any(k.name.endswith('_ancestor_taxon') for k in ll.get_external()):
        ancestorTaxa = [k for k in ll.get_external() if k.name.endswith('_ancestor_taxon')]
        regularTaxa = [k for k in ll.get_external() if k.name.endswith('_ancestor_taxon') == False]

        logger.warning(f"{len(ancestorTaxa)} ancestral taxa ('_ancestor_taxon' suffix) from travel-aware phylogeographic analysis detected. They will be removed whilst preserving phylogenetic relationships, resulting in a multitype tree.")

        ll = ll.reduce_tree(regularTaxa)

    ll.traverse_tree() ## traverse tree
    return ll

def load_JSON(jsonObject, treeType, jsonTranslation=None, sort=True, stats=True):
    """
    Load a ``baltic`` tree from an Auspice-style JSON source.

    **Parameters**

    jsonObject : str or dict
        Local path, Nextstrain URL, or already loaded JSON object.

    treeType : {'divergence', 'time'}
        Interpretation of branch lengths in the parsed tree.

    jsonTranslation : dict, optional
        Mapping from ``baltic`` attribute names to JSON keys or callables.

    sort : bool, optional
        If ``True``, sort branches after parsing.

    stats : bool, optional
        If ``True``, report tree statistics after parsing.

    **Returns**

    tuple[:class:`baltic.tree.Tree`, dict]
        Parsed tree and the associated metadata block from the JSON.

    **Examples**

    >>> import baltic as bt
    >>> auspice = {
    ...     "meta": {"colorings": []},
    ...     "tree": {
    ...         "name": "root",
    ...         "node_attrs": {"div": 0.0},
    ...         "children": [
    ...             {"name": "A", "node_attrs": {"div": 1.0}},
    ...             {"name": "B", "node_attrs": {"div": 1.2}},
    ...         ],
    ...     },
    ... }
    >>> ll, meta = bt.io.load_JSON(auspice, treeType="divergence", stats=False)
    >>> sorted(tip.name for tip in ll.get_external())
    ['A', 'B']
    >>> meta["colorings"]
    []
    """

    assert treeType in ['divergence', 'time'], f"Unrecognised treeType {treeType}. Must be 'divergence' or 'time'."

    if jsonTranslation is None:
        jsonTranslation = {'name': 'name', 'absoluteTime': 'num_date', 'height': 'div'}

    assert 'name' in jsonTranslation, f"jsonTranslation dict missing translation for baltic leaf attribute 'name' (default expectation: 'name': 'name')."

    if treeType == 'divergence':
        assert 'div' in jsonTranslation.values(), f"jsonTranslation dict missing translation for baltic branch attribute 'height' (default expectation: 'height': 'div')."
    elif treeType == 'time':
        assert 'num_date' in jsonTranslation.values(), f"jsonTranslation dict missing translation for baltic branch attribute 'absoluteTime' (default expectation: 'absoluteTime': 'num_date')."

    logger.debug('Reading JSON')

    if isinstance(jsonObject, str): ## string provided - either nextstrain URL or local path
        if 'nextstrain.org' in jsonObject: ## nextsrain.org in URL - request it
            logger.debug('Assume URL provided, loading JSON from nextstrain.org')

            auspice_json = json.load(csio(requests.get(jsonObject).content))
        else: ## not nextstrain.org URL - assume local path to auspice v2 json
            logger.debug('Loading JSON from local path')
            with open(jsonObject) as json_data:
                auspice_json = json.load(json_data)
    else: ## not string, assume auspice v2 json object given
        logger.debug('Loading JSON from object given')
        auspice_json = jsonObject

    json_meta = auspice_json['meta']
    if 'root_sequence' in auspice_json:
        json_meta['root_sequence'] = auspice_json['root_sequence'] ## transfer root sequence attribute to meta (by default it's a top level key in auspice json: {version, meta, tree, root_sequence})

    json_tree = auspice_json['tree']
    ll = make_tree_JSON(json_tree, jsonTranslation, treeType=treeType)

    logger.debug('Setting baltic traits from JSON')
    for k in ll.Objects: ## make node attributes easier to access
        for key in k.traits['node_attrs']:
            if isinstance(k.traits['node_attrs'][key], dict):
                if 'value' in k.traits['node_attrs'][key]:
                    k.traits[key] = k.traits['node_attrs'][key]['value']
                if 'confidence' in k.traits['node_attrs'][key]:
                    k.traits[f'{key}_confidence'] = k.traits['node_attrs'][key]['confidence']
            elif key == 'div':
                k.traits['divergence'] = k.traits['node_attrs'][key]

    for attr in jsonTranslation: ## iterate through attributes in json_translation
        for k in ll.Objects: ## for every branch
            if isinstance(jsonTranslation[attr], str):
                if jsonTranslation[attr] in k.traits:
                    setattr(k, attr,k.traits[jsonTranslation[attr]]) ## set attribute value for branch
                elif 'node_attrs' in k.traits and jsonTranslation[attr] in k.traits['node_attrs']:
                    setattr(k,attr, k.traits['node_attrs'][jsonTranslation[attr]])
                elif 'branch_attrs' in k.traits and jsonTranslation[attr] in k.traits['branch_attrs']:
                    setattr(k,attr, k.traits['branch_attrs'][jsonTranslation[attr]])
                else:
                    pass
                    # raise KeyError(f'String attribute {jsonTranslation[attr]} not found in JSON') ## since baltic objects are 'time' or 'divergence', baltic-exported JSONs will be missing one or the other and will raise an error here
            elif callable(jsonTranslation[attr]):
                setattr(k, attr, jsonTranslation[attr](k)) ## set attribute value with a function for branch
            else:
                raise AttributeError(f'Attribute {jsonTranslation[attr]} neither string nor callable')

    if treeType == 'divergence':
        branchUnit = 'height'
    elif treeType == 'time':
        branchUnit = 'absoluteTime'

    for k in ll.Objects: ## iterate over all branches
        curBranch = getattr(k, branchUnit) ## get parameter for this branch
        parBranch = getattr(k.parent, branchUnit) ## get parameter for parental branch
        k.length = curBranch - parBranch if curBranch and parBranch else 0.0 ## difference between current and parent is branch length (or, if parent unavailabel it's 0)

    logger.debug('Traversing and drawing tree')

    ll.traverse_tree()
    ll._assign_tree_coordinates()

    if stats:
        ll.treeStats() ## initial traversal, checks for stats
    if sort:
        ll.sort_branches() ## traverses tree, sorts branches, draws tree

    cmap={}
    for colouring in json_meta['colorings']:
        if colouring['type'] == 'categorical' and 'scale' in colouring:
            cmap[colouring['key']] = {}
            for entry in colouring['scale']:
                key, value = entry
                cmap[colouring['key']][key] = value
    setattr(ll, 'cmap', cmap)

    return ll, json_meta


def process_tip_dates(tree, tipRegex, dateFmt, variableDate, setNodes=True):
    """
    Extract sampling dates from tip labels and assign absolute times.

    **Parameters**

    tree : :class:`baltic.tree.Tree`
        Tree whose external branches should be inspected.

    tipRegex : str
        Regular expression used to capture the date token from each tip name.

    dateFmt : str
        Date format used to parse the captured token.

    variableDate : bool
        Whether partially specified dates should be interpreted with
        uncertainty ranges.

    setNodes : bool, optional
        If ``True``, assign absolute times to internal nodes as well as tips.

    **Examples**

    >>> import baltic as bt
    >>> ll = bt.make_tree("((A|2020-01-01:0.1,B|2020-02-01:0.2):0.3,C|2020-03-01:0.4);", treeType="divergence")
    >>> ll.sort_branches()
    >>> bt.io.process_tip_dates(ll, tipRegex=r"\\|([0-9\\-]+)$", dateFmt="%Y-%m-%d", variableDate=True, setNodes=False)
    >>> all(tip.absoluteTime is not None for tip in ll.get_external())
    True
    """
    tipDates = []
    tipNames = []
    tipDateUncertainties = {}
    for k in tree.get_external():
        tipNames.append(k.name)
        match = re.search(tipRegex, k.name)
        if match:
            date = calendar_to_decimal_date(match.group(1), fmt = dateFmt, variable = variableDate)
            if variableDate and (date[1][1] - date[1][0]) == 0: ## only keep dates that are certain for setting absolute time
                tipDates.append(date[0]) ## date will be tuple of (date, (dateMin,dateMax)) if variableDate == True
            elif variableDate == False:
                tipDates.append(date)

            tipDateUncertainties[k.name] = calendar_to_decimal_date(match.group(1), fmt = dateFmt, variable = True) ## logic around here might need more thinking
        else:
            logger.warning(f"Collection date of leaf {k.name} was not captured by regex {tipRegex}.")

    assert len(tipDates) > 0, f'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: {tipNames[0]}\nDate regex set to: {tipRegex}\nExpected date format: {dateFmt}'

    if setNodes:
        tree.set_absolute_time(max(tipDates), justLeaves=True if tree.treeType == 'divergence' else False)
    else:
        for k in tree.get_external():
            mean, interval = tipDateUncertainties[k.name]
            k.absoluteTime = mean

    if variableDate:
        tree._assign_date_uncertainty(tipDateUncertainties)
