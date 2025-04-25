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
                tipRegex='\|([0-9]+\-[0-9]+\-[0-9]+)',
                dateFmt='%Y-%m-%d',
                variableDate=True,
                absoluteTime=False,
                sortBranches=True):
    ll=None

    handle = open(treePath, 'r') if isinstance(treePath, str) else treePath

    for line in handle:
        l=line.strip('\n')
        if '(' in l:
            treeStringStart=l.index('(')
            ll=make_tree(l[treeStringStart:],treeType) ## send tree string to make_tree function
            logger.debug('Identified tree string')

    assert ll,'Regular expression failed to find tree string'
    ll.traverse_tree() ## traverse tree

    if sortBranches: ll.sort_branches() ## traverses tree, sorts branches, draws tree

    if absoluteTime:
        tip_dates=[]
        tip_names=[]
        for k in ll.get_external():
            tip_names.append(k.name)
            match=re.search(tipRegex,k.name)
            if match:
                tip_dates.append(calendar_to_decimal_date(match.group(1),fmt=dateFmt,variable=variableDate))
        assert len(tip_dates)>0, f'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: {tip_names[0]}\nDate regex set to: {tipRegex}\nExpected date format: {dateFmt}'
        ll.set_absolute_time(max(tip_dates))

    if isinstance(treePath,str):
        handle.close()
    return ll

def load_nexus(treePath,
                treeType,
                tipRegex='\|([0-9]+\-[0-9]+\-[0-9]+)',
                dateFmt='%Y-%m-%d',
                treestringRegex='tree [A-Za-z\_]+([0-9]+)',
                variableDate=True,
                absoluteTime=True,
                sortBranches=True):
    tip_flag=False
    tips={}
    tip_num=0
    ll=None

    handle = open(treePath, 'r') if isinstance(treePath, str) else treePath

    for line in handle:
        l=line.strip('\n')

        match=re.search('Dimensions ntax=([0-9]+);',l)
        if match:
            tip_num=int(match.group(1))
            logger.debug(f'File should contain {tip_num} taxa')

        match=re.search(treestringRegex,l)
        if match:
            treeString_start=l.index('(')
            ll=make_tree(l[treeString_start:],treeType) ## send tree string to make_tree function
            logger.debug('Identified tree string')

        if tip_flag:
            match=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l)
            if match:
                tips[match.group(1)]=match.group(2).strip('"').strip("'")
                logger.debug(f'Identified tip translation {match.group(1)}: {tips[match.group(1)]}')
            elif ';' not in l:
                print('tip not captured by regex:',l.replace('\t',''))

        if 'Translate' in l:
            tip_flag=True
        if ';' in l:
            tip_flag=False

    assert ll,'Failed to find tree string using regular expression'
    ll.traverse_tree() ## traverse tree
    if sortBranches:
        ll.sort_branches() ## traverses tree, sorts branches, draws tree
    if len(tips)>0:
        ll.renameTips(tips) ## renames tips from numbers to actual names
        ll.tipMap=tips
    if absoluteTime:
        tip_dates=[]
        tip_names=[]
        for k in ll.get_external():
            tip_names.append(k.name)
            match=re.search(tipRegex,k.name)
            if match:
                tip_dates.append(calendar_to_decimal_date(match.group(1),fmt=dateFmt,variable=variableDate))

        assert len(tip_dates)>0,f'Regular expression failed to find tip dates in tip names, review regex pattern or set absoluteTime option to False.\nFirst tip name encountered: {tip_names[0]}\nDate regex set to: {tipRegex}\nExpected date format: {dateFmt}'
        ll.set_absolute_time(max(tip_dates))

    if isinstance(treePath,str):
        handle.close()
    return ll

def load_JSON(jsonObject,
                treeType,
                jsonTranslation=None,
                sort=True,
                stats=True):
    if jsonTranslation is None:
        jsonTranslation={'name':'name','absoluteTime':'num_date'}


    lengthKeys = ['absoluteTime', 'length', 'height']
    # TODO the next line breaks because there is no `required_keys`
    assert 'name' in jsonTranslation and any(key in jsonTranslation for key in required_keys),f"JSON translation dictionary missing entries: {', '.join([entry for entry in ['name']+lengthKeys if not entry in jsonTranslation])}"
    logger.debug('Reading JSON')

    if isinstance(jsonObject,str): ## string provided - either nextstrain URL or local path
        if 'nextstrain.org' in jsonObject: ## nextsrain.org in URL - request it
            logger.debug('Assume URL provided, loading JSON from nextstrain.org')

            auspice_json=json.load(csio(requests.get(jsonObject).content))
        else: ## not nextstrain.org URL - assume local path to auspice v2 json
            logger.debug('Loading JSON from local path')
            with open(jsonObject) as json_data:
                auspice_json = json.load(json_data)
    else: ## not string, assume auspice v2 json object given
        logger.debug('Loading JSON from object given')
        auspice_json=jsonObject

    json_meta=auspice_json['meta']
    json_tree=auspice_json['tree']
    ll=make_tree_JSON(json_tree,jsonTranslation,treeType=treeType)

    assert ('absoluteTime' in jsonTranslation and ('length' not in jsonTranslation or 'height' not in jsonTranslation)) or ('absoluteTime' not in jsonTranslation and ('length' in jsonTranslation or 'height' in jsonTranslation)),'Cannot use both absolute time and branch length, include only one in json_translation dictionary.'

    logger.debug('Setting baltic traits from JSON')
    for k in ll.Objects: ## make node attributes easier to access
        for key in k.traits['node_attrs']:
            if isinstance(k.traits['node_attrs'][key],dict):
                if 'value' in k.traits['node_attrs'][key]:
                    k.traits[key]=k.traits['node_attrs'][key]['value']
                if 'confidence' in k.traits['node_attrs'][key]:
                    k.traits[f'{key}_confidence']=k.traits['node_attrs'][key]['confidence']
            elif key=='div':
                k.traits['divergence']=k.traits['node_attrs'][key]

    for attr in jsonTranslation: ## iterate through attributes in json_translation
        for k in ll.Objects: ## for every branch
            if isinstance(jsonTranslation[attr],str):
                if jsonTranslation[attr] in k.traits:
                    setattr(k,attr,k.traits[jsonTranslation[attr]]) ## set attribute value for branch
                elif 'node_attrs' in k.traits and jsonTranslation[attr] in k.traits['node_attrs']:
                    setattr(k,attr,k.traits['node_attrs'][jsonTranslation[attr]])
                elif 'branch_attrs' in k.traits and jsonTranslation[attr] in k.traits['branch_attrs']:
                    setattr(k,attr,k.traits['branch_attrs'][jsonTranslation[attr]])
                else:
                    raise KeyError(f'String attribute {jsonTranslation[attr]} not found in JSON')
            elif callable(jsonTranslation[attr]):
                setattr(k,attr,jsonTranslation[attr](k)) ## set attribute value with a function for branch
            else:
                raise AttributeError(f'Attribute {jsonTranslation[attr]} neither string nor callable')

    for branch_unit in ['height','absoluteTime']: ## iterate between divergence and absolute time
        if branch_unit in jsonTranslation: ## it's available in tree
            for k in ll.Objects: ## iterate over all branches
                cur_branch=getattr(k,branch_unit) ## get parameter for this branch
                par_branch=getattr(k.parent,branch_unit) ## get parameter for parental branch
                k.length=cur_branch-par_branch if cur_branch and par_branch else 0.0 ## difference between current and parent is branch length (or, if parent unavailabel it's 0)

    logger.debug('Traversing and drawing tree')

    ll.traverse_tree()
    ll.assign_tree_coordinates()
    if stats:
        ll.treeStats() ## initial traversal, checks for stats
    if sort:
        ll.sort_branches() ## traverses tree, sorts branches, draws tree

    cmap={}
    for colouring in json_meta['colorings']:
        if colouring['type']=='categorical' and 'scale' in colouring:
            cmap[colouring['key']]={}
            for entry in colouring['scale']:
                key,value=entry
                cmap[colouring['key']][key]=value
    setattr(ll,'cmap',cmap)

    return ll,json_meta
