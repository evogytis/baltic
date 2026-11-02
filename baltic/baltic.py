"""This module provides the core ``baltic`` functions :func:`make_tree` and :func:`make_tree_JSON`.

**Notes**

This version of ``baltic`` (v1.0 (Cedar)) contains many API changes from previous versions, and is not backwards-compatible. If you find pieces of documentation that refer to the old API, please let us know and we will try to update them with the next update.


**Attributes**

logger : ``logging.Logger``
    Default logger which will be passed to other ``baltic`` functions.
"""
__all__ = ['make_tree', 'make_tree_JSON', 'TreeParseError']
import re
import sys
import logging
from baltic.tree import Tree
from baltic.node import Node
from baltic.leaf import Leaf
from baltic.reticulation import Reticulation

logger = logging.getLogger("baltic")
sys.setrecursionlimit(9001)

## characters allowed inside [&...] annotations, as a regex character class body
_ANNOTATION_CHARS = r'A-Za-z\_\-{}\,0-9\.\%=\"\'\+!# :\/\(\)\&'
## the same without quote characters: quoted values are matched separately and may hold any character
_ANNOTATION_CHARS_UNQUOTED = r'A-Za-z\_\-{}\,0-9\.\%=\+!# :\/\(\)\&'

## a decimal number, optionally signed and in scientific notation (MrBayes writes 1.00000000e+00)
_NUMBER = r'[-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?'
_NUMBER_RE = re.compile(_NUMBER)
## annotation keys: BEAST (height_95%_HPD, location.set.prob), MrBayes (prob(percent), prob+-sd)
## and IQ-TREE concordance factors (gCF/gDF1/gDF2/gDFP)
_ANNOTATION_KEY_RE = re.compile(r'[A-Za-z0-9_.%()+/-]+')
## unquoted annotation values that are read as strings
_ANNOTATION_WORD_RE = re.compile(r'[A-Za-z0-9_.+ :/()&-]+')


class TreeParseError(ValueError, AssertionError):
    """
    Raised when :func:`make_tree` cannot interpret a tree string.

    The message gives the position in the string, a short snippet around it and,
    where possible, a hint about what to fix. It subclasses ``ValueError`` and,
    for backwards compatibility with code written against earlier versions (which
    raised assertion errors), ``AssertionError``.
    """


def _parse_error(data, pos, message, hint=None, context=30):
    """Build a :class:`TreeParseError` pointing at position ``pos`` of ``data``."""
    start, end = max(0, pos - context), min(len(data), pos + context)
    prefix = '...' if start > 0 else ''
    suffix = '...' if end < len(data) else ''
    snippet = data[start:end].replace('\n', ' ').replace('\t', ' ')
    lines = [f'{message} at position {pos}',
             f'  {prefix}{snippet}{suffix}',
             '  ' + ' ' * (len(prefix) + pos - start) + '^']
    if hint:
        lines.append(hint)
    return TreeParseError('\n'.join(lines))


def _stalled_parse_error(data, i):
    """Explain why parsing could not advance past position ``i``."""
    j = i + 1 if data[i] == ':' else i ## annotations may follow a colon
    if data.startswith('[&', j):
        if data.find(']', j) == -1:
            return _parse_error(data, j, 'Unterminated annotation (no closing "]")')
        allowed = re.compile(f'[{_ANNOTATION_CHARS_UNQUOTED}]')
        quote = None
        for pos in range(j + 2, len(data)): ## scan the annotation, skipping quoted values
            char = data[pos]
            if quote is not None:
                if char == quote:
                    quote = None
            elif char in ('"', "'"):
                quote, opened = char, pos
            elif char == ']':
                break
            elif not allowed.match(char):
                return _parse_error(data, pos,
                                    f'Unsupported character {char!r} in annotation',
                                    hint='Check the tree file; outside quoted values, annotations may only contain '
                                         'letters, digits, spaces and the characters _-{},.%=+!#:/()&')
        else:
            if quote is not None:
                return _parse_error(data, opened, 'Unterminated quote in annotation')
            return _parse_error(data, j, 'Unterminated annotation (no closing "]")')
    return _parse_error(data, i, f'Tree string unparseable: unexpected character {data[i]!r}')


def _split_annotation_fields(comment):
    """Split a BEAST annotation body outside quotes and nested braces."""
    start = 0
    depth = 0
    quote = None
    for pos, char in enumerate(comment):
        if quote is not None:
            if char == quote:
                quote = None
        elif char in ('"', "'"):
            quote = char
        elif char == '{':
            depth += 1
        elif char == '}':
            depth -= 1
            if depth < 0:
                raise ValueError("Unbalanced braces in BEAST annotation")
        elif char == ',' and depth == 0:
            yield comment[start:pos].strip()
            start = pos + 1
    if quote is not None or depth != 0:
        raise ValueError("Unterminated quote or braces in BEAST annotation")
    yield comment[start:].strip()


def _unquote_label(match):
    """Text of a matched quoted node label, with doubled quotes (the Newick escape) undone."""
    quote = match.group()[0]
    text = match.group(1) if quote == "'" else match.group(2)
    return text.replace(quote * 2, quote)


def _set_member_to_float(value):
    """Convert one member of a numeric annotation set, accepting only numbers (and NaN/Inf)."""
    value = value.strip()
    if _NUMBER_RE.fullmatch(value) or value.lower() in ('nan', 'inf', '-inf', '+inf'):
        return float(value)
    raise ValueError(value)


def _first_of_tie(value):
    """
    Reduce a tied discrete state such as ``"Asia+Europe"`` to its first state.

    Only values whose parts are all non-numeric words are split, so MrBayes values
    such as ``prob+-sd="100+-0"`` are kept whole. Returns ``(value, tied)``.
    """
    parts = value.split('+')
    if len(parts) > 1 and all(part and not _NUMBER_RE.fullmatch(part) for part in parts):
        return parts[0], True
    return value, False


def _read_annotation_field(field):
    """
    Read one ``key=value`` field of a ``[&...]`` annotation.

    Returns ``(key, value, tied)``, or ``None`` for FigTree display fields (``!key=value``),
    which are ignored. Raises ``ValueError(message, hint)`` if the field cannot be read,
    so that nothing is dropped silently.
    """
    key, separator, value = field.partition('=')
    key, value = key.strip(), value.strip()

    if key.startswith('!'): ## FigTree display settings, e.g. !color=#ff0000
        logger.info('FigTree comment found, ignoring')
        return None

    unreadable = (f'Could not read annotation field {field!r}',
                  'Check the tree file; the field is not in a format baltic recognises.')
    if not separator or not _ANNOTATION_KEY_RE.fullmatch(key):
        raise ValueError(*unreadable)

    tied = False
    if len(value) >= 2 and value[0] in ('"', "'") and value[-1] == value[0]: ## quoted string
        value, tied = _first_of_tie(value[1:-1])
        value = {'true': True, 'false': False}.get(value, value)

    elif value.startswith('{') and value.endswith('}'):
        body = value[1:-1].strip()
        if not body: ## empty set, e.g. segments={} on a reassortment network node
            value = []
        elif '{' in body: ## complete history logged by robust counting (MCMC trees)
            history = []
            for entry in re.findall(r'{([^{}]+)}', body):
                entryFields = [entryField.strip() for entryField in entry.split(',')]
                if len(entryFields) in [3,4]:
                    if key in ['all_S','all_N'] or key.endswith(('.all_S','.all_N')):
                        entryFields = [int(entryFields[0]), float(entryFields[1])] + entryFields[2:]
                    history.append(entryFields)
            value = history
        elif 'set' in key and 'set.prob' not in key: ## set of discrete states
            value = [member.strip().strip('"') for member in body.split(',')]
        else: ## set of state probabilities, or a range or HPD interval
            try:
                value = [_set_member_to_float(member) for member in body.split(',')]
            except ValueError:
                raise ValueError(f'Could not read annotation field {field!r}: expected a set of numbers',
                                 'Fields named without "set" must hold numbers, e.g. height_95%_HPD={0.1,0.2}.') from None

    elif _NUMBER_RE.fullmatch(value):
        value = float(value)

    elif _ANNOTATION_WORD_RE.fullmatch(value): ## unquoted string
        value, tied = _first_of_tie(value)
        value = {'true': True, 'false': False}.get(value, value)

    else:
        raise ValueError(*unreadable)

    return key, value, tied


def make_tree(data, treeType, tre=None):
    """
    Parse a Newick-like tree string into a :class:`baltic.tree.Tree`.

    Handles BEAST integer tips, ``[&...]`` annotation comments and ``#label``
    reticulations in addition to plain Newick. Node labels, quoted or not (e.g.
    IQ-TREE support ``)95.5/100:`` or ASTRAL ``)'[pp1=0.95;q1=0.6]':``), are stored
    as written in ``traits['label']``; when a name is followed directly by a quoted
    label (``N5'[pp1=1]'``), the name goes in ``label`` and the quoted text in
    ``quoted_label``. The returned tree has not been
    traversed; the loaders in :mod:`baltic.io` do that, so call
    :meth:`baltic.tree.Tree.traverse_tree` if using this directly.

    .. note::

       If every branch length parses as ``0.0`` -- as happens with a topology-only
       cladogram such as ``"((A,B),C);"`` -- branch lengths are **rewritten** to make
       the tree ultrametric with a height of 1.0, and a warning is logged. The
       lengths in the resulting tree are therefore invented, not parsed.

    **Parameters**

    data : str
        Tree string to parse.

    treeType : {'divergence', 'time'}
        Interpretation of branch lengths in the resulting tree.

    tre : :class:`baltic.tree.Tree`, optional
        Existing tree object to populate. If omitted, a new tree is created.

    **Returns**

    :class:`baltic.tree.Tree`
        Parsed tree object.

    **Raises**

    TreeParseError
        If the string does not end in a semicolon, if parentheses are unbalanced,
        if parsing stalls on a character it cannot interpret (for example an
        unsupported character inside an annotation), if the string ends before
        its final semicolon, if a branch length is malformed (e.g. ``1.2e``), if an
        annotation field cannot be read, or if a
        reticulation label is used more than once in the tree. The message gives
        the position and a snippet of the string. ``TreeParseError`` subclasses
        both ``ValueError`` and ``AssertionError``.

    **Examples**

    >>> import baltic as bt
    >>> ll = bt.make_tree("((A:1.0,B:1.5):0.5,C:2.0);", treeType="divergence")
    >>> [tip.name for tip in ll.get_external()]
    ['A', 'B', 'C']
    """
    patterns = {
        'beast_tip': re.compile(r'(\(|,)([0-9]+)(\[|\:)'), # Pattern to match tips in BEAST format (integers)
        'non_beast_tip': re.compile(r'(\(|,)(\'|\")*([^\(\):\[\'\"#,]+)(\'|\"|)*(\[)*'), # Pattern to match tips with unencoded names
        'multitype': re.compile(r'\)([0-9]+)\['), # Pattern to match multitype tree singletons
        'reticulation_start': re.compile(r'[\(,](#[A-Za-z0-9]+)'), # Pattern to match the beginning of a reticulate branch
        'reticulation_end': re.compile(r'\)(#[A-Za-z0-9]+)'), # Pattern to match the landing point of a reticulate branch
        'annotation': re.compile(rf'(\:)*\[(&(?:"[^"]*"|\'[^\']*\'|[{_ANNOTATION_CHARS_UNQUOTED}])+)\]'), # Pattern to match MCC comments; quoted values may hold any character, e.g. label="[pp1=0.95;q1=0.6]"
        'node_label': re.compile(r'([A-Za-z\_\-0-9\.\/]+)(\:|\;|\[|,|\)|\'|\")'), # Pattern to match old school node labels
        'quoted_label': re.compile(r"'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\""), # Pattern to match quoted node labels, with doubled quotes as escapes
        'branch_length': re.compile(rf'(\:)*({_NUMBER})'), # Pattern to match branch lengths without comments
    }
    ## patterns are matched in place (pattern.match(data, pos)) rather than on slices of data,
    ## so long tip names or annotations are never truncated by a fixed-size window
    match_at = lambda name, pos: patterns[name].match(data, pos) if pos >= 0 else None

    if not isinstance(data,str): ## tree string is not an instance of string (could be unicode) - convert
        data=str(data)
    data=data.rstrip() ## whitespace after the final semicolon (e.g. "; " or a Windows "\r") carries no data

    # Add in some checks that the data are in correct format
    if not data.endswith(";"):
        raise _parse_error(data, max(len(data) - 1, 0), 'Improperly formatted string: must end in semicolon')
    if data.count("(") != data.count(")"): ## recount outside quotes - quoted labels may contain parentheses
        opened = [] ## locate the first unmatched parenthesis to point at it
        quote = None
        for pos, char in enumerate(data):
            if quote is not None:
                if char == quote:
                    quote = None
            elif char in ('"', "'"):
                quote = char
            elif char == '(':
                opened.append(pos)
            elif char == ')':
                if not opened:
                    raise _parse_error(data, pos, "Improperly formatted string: unmatched ')'")
                opened.pop()
        if opened:
            raise _parse_error(data, opened[-1],
                               f"Improperly formatted string: {data.count('(')} '(' but {data.count(')')} ')'")

    if tre is None: ## calling without providing a tree object - create one
        tre=Tree(treeType)
    i=0 ## is an adjustable index along the tree string, it is incremented to advance through the string
    storedI=None ## store the i at the end of the loop, to make sure we haven't gotten stuck somewhere in an infinite loop
    tiedFields=[] ## annotation fields with tied states ("A+B") reduced to their first state

    while i < len(data): ## while there's characters left in the tree string - loop away
        if storedI == i: logger.debug(f'{i} >{data[i]}<')

        if storedI == i: ## index has not advanced - the string cannot be parsed here
            raise _stalled_parse_error(data, i)
        storedI=i # Store the current index at the end of the loop to check for infinite loops

        if data[i] == '(': ## look for new nodes
            logger.debug(f'{i} adding node')
            tre.add_node(i) ## add node to current node in tree ll
            i+=1 ## advance in tree string by one character

        if data[i-1] in ('(', ',') and data[i] in ('"', "'"): ## quoted tip name, which may hold any character, e.g. 'A, b' or 'O''Brien'
            match=match_at('quoted_label',i)
            if match is None:
                raise _parse_error(data, i, 'Unterminated quoted tip name')
            logger.debug(f'{i} adding leaf (quoted) {_unquote_label(match)}')
            tre.add_leaf(i,_unquote_label(match))
            i=match.end()

        match=match_at('beast_tip',i-1) ## look for tips in BEAST format (integers).
        if match:
            logger.debug(f'{i} adding leaf (BEAST) {match.group(2)}')
            tre.add_leaf(i,match.group(2)) ## add tip
            i+=len(match.group(2)) ## advance in tree string by however many characters the tip is encoded

        match=match_at('non_beast_tip',i-1)  ## look for tips with unencoded names - if the tips have some unusual format you'll have to modify this
        if match:
            logger.debug(f'{i} adding leaf (non-BEAST) {match.group(3)}')
            tre.add_leaf(i,match.group(3).strip('"').strip("'"))  ## add tip
            i+=len(match.group(3))+match.group().count("'")+match.group().count('"') ## advance in tree string by however many characters the tip is encoded

        match=match_at('multitype',i-1) ## look for multitype tree singletons.
        if match:
            logger.debug(f'{i} adding multitype node {match.group(1)}')
            i+=len(match.group(1))

        match=match_at('reticulation_start',i-1) ## look for beginning of reticulate branch
        if match:
            logger.debug(f'{i} adding outgoing reticulation branch {match.group(1)}')
            tre.add_reticulation(match.group(1)) ## add reticulate branch

            destination=None
            for k in tre.Objects: ## iterate over branches parsed so far
                if 'label' in k.traits and k.traits['label']==match.group(1): ## if there's a branch with a matching id
                    if destination is None: ## not set destination before
                        destination=k ## destination is matching node
                    else: ## destination seen before - raise an error (indicates reticulate branch ids are not unique)
                        raise _parse_error(data, i, f'Reticulate branch not unique: {match.group(1)} seen elsewhere in the tree')
            if destination: ## identified destination of this branch
                logger.debug(f'identified {match.group(1)} destination')
                tre.curNode.target=destination ## set current node's target as the destination
                setattr(destination,"contribution",tre.curNode) ## add contributing edge to destination
            else:
                logger.debug(f'destination of {match.group(1)} not identified yet')
            i+=len(match.group())-1

        match=match_at('reticulation_end',i-1) ## look for landing point of reticulate branch
        if match:
            logger.debug(f'{i} adding incoming reticulation branch {match.group(1)}')
            tre.curNode.traits['label']=match.group(1) ## set node label

            origin=None ## branch is landing, check if its origin was seen previously
            for k in tre.Objects: ## iterate over currently existing branches
                if isinstance(k,Reticulation) and k.name==match.group(1): ## check if any reticulate branches match the origin
                    if origin is None: ## origin not identified yet
                        origin=k ## origin is reticulate branch with the correct name
                    else: ## origin has been identified - shouldn't happen, implies that multiple reticulate branches exist with the same name
                        raise _parse_error(data, i, f'Reticulate branch not unique: {match.group(1)} seen elsewhere in the tree')
            if origin: ## identified origin
                logger.debug(f'identified {match.group(1)} origin')
                origin.target=tre.curNode ## set origin's landing at this node
                setattr(tre.curNode,"contribution",origin) ## add contributing edge to this node
            else:
                logger.debug(f'origin of {match.group(1)} not identified yet')
            i+=len(match.group())-1

        match=match_at('annotation',i)## look for MCC comments
        if match:
            logger.debug(f'{i} comment: {match.group(2)}')
            comment=match.group(2)
            try:
                fields=list(_split_annotation_fields(comment[1:]))
            except ValueError as err:
                raise _parse_error(data, i, f'Malformed annotation: {err}') from None
            for field in fields:
                if not field: ## nothing to read, e.g. an empty [&] comment
                    continue
                try:
                    parsed = _read_annotation_field(field)
                except ValueError as err:
                    raise _parse_error(data, i, *err.args) from None
                if parsed is None: ## FigTree display field
                    continue
                tr, val, tied = parsed
                tre.curNode.traits[tr] = val
                if tied: ## DO NOT ALLOW EQUIPROBABLE DOUBLE ANNOTATIONS (which are in format "A+B") - just get the first one
                    tiedFields.append(field)

            i+=len(match.group()) ## advance in tree string by however many characters it took to encode labels

        if i > 0 and data[i-1] == ')' and data[i] in ('"', "'"): ## quoted node label, e.g. IQ-TREE )'95.5/100': or ASTRAL )'[pp1=0.95;q1=0.6]':
            match=match_at('quoted_label',i)
            if match is None:
                raise _parse_error(data, i, 'Unterminated quoted node label')
            label=_unquote_label(match)
            logger.debug(f'quoted node label found: {label}')
            tre.curNode.traits['label']=label ## stored as text, exactly as written - support values are not split or converted
            i=match.end()

        match=match_at('node_label',i)## look for old school node labels

        if match:
            logger.debug(f'old school comment found: {match.group(1)}')
            tre.curNode.traits['label']=match.group(1)

            i+=len(match.group(1))

            if data[i] in ('"', "'"): ## name followed by a quoted label, e.g. ASTRAL N5'[pp1=1;pp2=0;pp3=0]' - kept apart from the name
                match=match_at('quoted_label',i)
                if match is None:
                    raise _parse_error(data, i, 'Unterminated quoted node label')
                tre.curNode.traits['quoted_label']=_unquote_label(match)
                i=match.end()

        micromatch=match_at('branch_length',i) ## look for branch lengths without comments
        if micromatch is not None:
            end=micromatch.end()
            if end < len(data) and (data[end].isalnum() or data[end] in '._'): ## e.g. 1.2e - the rest would otherwise become a node label
                token=re.match(r'[^,;\(\)\[\]:]*', data[micromatch.start(2):]).group()
                raise _parse_error(data, micromatch.start(2), f'Malformed branch length {token!r}')
            logger.debug(f'adding branch length ({i}) {float(micromatch.group(2))}')
            tre.curNode.length=float(micromatch.group(2)) ## set branch length of current node
            i+=len(micromatch.group()) ## advance in tree string by however many characters it took to encode branch length

        if data[i] == ',' or data[i] == ')': ## look for bifurcations or clade ends
            i+=1 ## advance in tree string
            tre.curNode=tre.curNode.parent

        if data[i] == ';': ## look for string end
            if i != len(data) - 1: ## a semicolon before the end would silently drop the rest of the tree
                raise _parse_error(data, i, 'Unexpected ";" before the end of the tree string')
            if tiedFields:
                logger.warning(f'{len(tiedFields)} annotation values held tied states (e.g. {tiedFields[0]}); kept the first state of each.')
            if sum(tre.get_parameter_list('length')) == 0.0:
                logger.warning(f"All branch lengths 0.0, setting branch lengths to make tree ultrametric with a height of 1.0")

                maxLvl = 0
                for k in tre.get_external():
                    k.height = 1.0

                    pathLen = len(k.get_path_to_root()) - 1
                    if pathLen > maxLvl:
                        maxLvl = pathLen

                heightless = tre.get_internal(lambda k: k.height is None)

                while len(heightless) > 0:
                    for node in heightless:
                        if len(node.children) == len([ch for ch in node.children if ch.height is not None]):
                            node.height = min((ch.height - (1.0 / maxLvl)) for ch in node.children)
                    heightless = tre.get_internal(lambda k: k.height is None)

                for k in tre.Objects:
                    k.length = (k.height - k.parent.height) if k.parent else 0.0

                tre.traverse_tree()

            return tre

def make_tree_JSON(jsonNode, jsonTranslationDict, treeType, tre=None,):
    """
    Build a :class:`baltic.tree.Tree` from an Auspice-style JSON node hierarchy.

    **Parameters**

    jsonNode : dict
        Root JSON node to parse.

    jsonTranslationDict : dict
        Mapping from ``baltic`` attribute names to JSON keys.

    treeType : {'divergence', 'time'}
        Interpretation of branch lengths in the resulting tree.

    tre : :class:`baltic.tree.Tree`, optional
        Existing tree object to populate recursively. If omitted, a new tree
        is created.

    **Returns**

    :class:`baltic.tree.Tree`
        Parsed tree object, already traversed.

    **Raises**

    KeyError
        If a JSON leaf lacks the name field named by *jsonTranslationDict*.

    **Notes**

    Internal nodes without a name in the JSON are assigned a generated
    ``NODE_<n>`` identifier, matching Auspice's own convention. This is why the
    second example below reports ``'NODE_0000001'`` for an unnamed root.

    **Examples**

    >>> import baltic as bt
    >>> json_tree = {
    ...     "name": "root",
    ...     "node_attrs": {"div": 0.0},
    ...     "children": [
    ...         {"name": "A", "node_attrs": {"div": 1.0}},
    ...         {"name": "B", "node_attrs": {"div": 1.2}},
    ...     ],
    ... }
    >>> translation = {"name": "name", "height": "div"}
    >>> ll = bt.make_tree_JSON(json_tree, translation, treeType="divergence")
    >>> sorted(tip.name for tip in ll.get_external())
    ['A', 'B']
    >>> unnamed_json_tree = {
    ...     "node_attrs": {"div": 0.0},
    ...     "children": [
    ...         {"name": "A", "node_attrs": {"div": 1.0}},
    ...         {"name": "B", "node_attrs": {"div": 1.2}},
    ...     ],
    ... }
    >>> ll = bt.make_tree_JSON(unnamed_json_tree, translation, treeType="divergence")
    >>> ll.root.index
    'NODE_0000001'
    """
    node_name_key = jsonTranslationDict['name']

    if 'children' in jsonNode: ## only nodes have children
        newNode=Node()
    else:
        if node_name_key not in jsonNode:
            raise KeyError(f"JSON leaf is missing required name field '{node_name_key}'.")
        newNode=Leaf(name = jsonNode[node_name_key])

    if tre is None:
        tre = Tree(treeType)
        tre.root = newNode

    if 'attr' in jsonNode:
        attr = jsonNode.pop('attr')
        jsonNode.update(attr)

    if node_name_key in jsonNode:
        node_name = jsonNode[node_name_key]
    elif 'children' in jsonNode:
        node_name = f"NODE_{len(tre.Objects) + 1:07d}"
        jsonNode[node_name_key] = node_name

    newNode.parent = tre.curNode ## set parent-child relationships
    tre.curNode.children.append(newNode)
    newNode.index = node_name ## indexing is based on name
    newNode.traits = {n: jsonNode[n] for n in list(jsonNode.keys()) if n != 'children'} ## set traits to non-children attributes
    tre.Objects.append(newNode)
    tre.curNode = newNode

    if 'children' in jsonNode:
        for child in jsonNode['children']:
            make_tree_JSON(jsonNode=child, jsonTranslationDict=jsonTranslationDict, treeType=treeType, tre=tre)
            tre.curNode = tre.curNode.parent
    return tre

if __name__ == '__main__':
    pass
