import re
import copy
import logging
import datetime as dt
import calendar
import warnings
import math
import numpy as np
from itertools import permutations
from scipy.stats import linregress
import matplotlib as mpl
from matplotlib.collections import LineCollection

logger = logging.getLogger("baltic.bt_utils")

def calendar_to_decimal_date(date, fmt="%Y-%m-%d", variable=False):
    """
    Convert a calendar date of a specified format into a decimal number.

    Parameters
    ----------
    date : str
        Date string to be converted.

    fmt : str, default="%Y-%m-%d"
        String encoding the format of the input date. Must be parsable
        by Python's `datetime <https://docs.python.org/3/library/datetime.html#>`__
        module.

    variable : bool, default=False
        Set to ``True`` when dates may be of variable lengths (e.g. when
        looping over ``["2025-01-01", "2025-02"]``). Will use highest
        precision available.

    Returns
    -------
    str

    Example
    -------
    >>> calendar_to_decimal_date("1253-07-06")
    1253.509589041096
    """

    if not fmt:
        return date

    # Detect delimiter in the format (e.g. '-', '/', etc.)
    delimiter = re.search(r"[^0-9A-Za-z%]", fmt)
    delimit = delimiter.group() if delimiter is not None else None

    # ------------------------------------------------------------------
    # 1. Determine how much of the date is present when variable=True
    # ------------------------------------------------------------------
    if variable:
        if delimit:
            dateL = len(date.split(delimit))
        else:
            # If no delimiter exists, treat the date as a single field (e.g. "2020")
            dateL = 1
        # Reduce fmt to available precision
        if dateL == 2: # "YYYY-MM"
            fmt = delimit.join(fmt.split(delimit)[:-1])
        elif dateL == 1: # "YYYY"
            fmt = delimit.join(fmt.split(delimit)[:-2])
        # If dateL == 3, fmt stays as full (YYYY-MM-DD)

    adatetime = dt.datetime.strptime(date, fmt)
    year = adatetime.year

    def to_decimal(d: dt.datetime): # Convert an exact datetime into a decimal year
        boy = dt.datetime(d.year, 1, 1)
        eoy = dt.datetime(d.year + 1, 1, 1)
        return d.year + (d - boy).total_seconds() / (eoy - boy).total_seconds()

    if not variable:
        return to_decimal(adatetime)

    if delimit:
        dateL = len(date.split(delimit))
    else:
        dateL = 1

    # CASE 1 — Full date available (YYYY-MM-DD)
    if dateL == 3:
        dec = to_decimal(adatetime)
        return dec, (dec, dec)

    # CASE 2 — Only year and month (YYYY-MM)
    if dateL == 2:
        year = adatetime.year
        month = adatetime.month

        # Min = first day of that month
        dmin = dt.datetime(year, month, 1)

        # Max = last day of that month
        last_day = calendar.monthrange(year, month)[1]
        dmax = dt.datetime(year, month, last_day)

        dec_min = to_decimal(dmin)
        dec_max = to_decimal(dmax)
        dec_mid = 0.5 * (dec_min + dec_max)

        return dec_mid, (dec_min, dec_max)

    # CASE 3 — Only year given (YYYY)
    if dateL == 1:
        # Min = Jan 1 of that year
        dmin = dt.datetime(year, 1, 1)
        # Max = Dec 31 of that year
        dmax = dt.datetime(year, 12, 31)

        dec_min = to_decimal(dmin)
        dec_max = to_decimal(dmax)
        dec_mid = 0.5 * (dec_min + dec_max)

        return dec_mid, (dec_min, dec_max)

    # Fallback (shouldn't happen)
    dec = to_decimal(adatetime)
    return dec, (dec, dec)


# def calendar_to_decimal_date(date,fmt="%Y-%m-%d",variable=False):
#     if not fmt:
#         return date
#     delimiter=re.search('[^0-9A-Za-z%]',fmt) ## search for non-alphanumeric symbols in fmt (should be field delimiter)
#     delimit=None
#     if delimiter is not None:
#         delimit=delimiter.group()

#     if variable: ## if date is variable - extract what is available
#         if delimit is not None:
#             dateL=len(date.split(delimit)) ## split date based on symbol
#         else:
#             dateL=1 ## no non-alphanumeric characters in date, assume dealing with an imprecise date (something like just year)

#         if dateL==2:
#             fmt=delimit.join(fmt.split(delimit)[:-1]) ## reduce fmt down to what's available
#         elif dateL==1:
#             fmt=delimit.join(fmt.split(delimit)[:-2])

#     adatetime=dt.datetime.strptime(date,fmt) ## convert to datetime object
#     year = adatetime.year ## get year
#     boy = dt.datetime(year, 1, 1) ## get beginning of the year
#     eoy = dt.datetime(year + 1, 1, 1) ## get beginning of next year
#     return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds())) ## return fractional year

def decimal_to_calendar_date(timepoint,fmt='%Y-%m-%d'):
    """
    Convert a decimal year value to a formatted calendar date string.

    Parameters
    ----------
    timepoint : float
        Decimal year to convert.

    fmt : str, optional
        Output date format passed to :func:`datetime.datetime.strftime`.

    Returns
    -------
    str
        Formatted calendar date.
    """
    year = int(timepoint)
    rem = timepoint - year

    base = dt.datetime(year, 1, 1)
    result = base + dt.timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    return dt.datetime.strftime(result,fmt)

def convert_date_format(dateString,startFormat,endFormat):
    """
    Reformat a date string from one ``datetime`` format to another.

    Parameters
    ----------
    dateString : str
        Input date string.

    startFormat : str
        Format used to parse *dateString*.

    endFormat : str
        Format used to render the output string.

    Returns
    -------
    str
        Reformatted date string.
    """
    return dt.datetime.strftime(dt.datetime.strptime(dateString,startFormat),endFormat)
    try: #TODO deal with stuff that comes after the return statement
        date_obj = dt.datetime.strptime(dateString, startFormat)
        return dt.datetime.strftime(date_obj, endFormat)
    except ValueError as e:
        raise ValueError('Error converting date "%s" from format "%s" to "%s": "%s"'%(dateString, startFormat, endFormat, e))

def state_collapse_tree(tree, switchFxn, keepLast=True, adjustEarlyHeights=False):
    """
    Return a deepcopied and reduced version of the tree provided where subtrees are labelled identically when branches evaluate switchFxn to False.
    Also known by the name of Phylotype maps/trees.
    """
    import copy

    localTree = copy.deepcopy(tree)

    localTree._partition_tree(partitionFxn = switchFxn) ## run tree labelling

    labels = set(localTree.get_parameter_list('partition', useTraitsDict=True)) ## get all unique labels in tree
    sortingFxn = lambda k: k.absoluteTime if localTree.treeType == 'time' else k.height ## determine whether to sort by height or absoluteTime

    labelToBranches = {label: [k for k in localTree.get_external() if k.traits['partition'] == label] for label in labels} ## map each unique label to the furthest tip with that label
    labelCounts = {label: len(labelToBranches[label]) for label in labels}

    labelToLastBranch = {label: sorted(labelToBranches[label], key=sortingFxn)[-1] for label in labels if labelCounts[label] > 0}

    if keepLast: ## keeping last tip
        localTree = localTree.reduce_tree(labelToLastBranch.values()) ## keep a single last tip for each unique label
        if keepLast and adjustEarlyHeights:
            logger.warning(f"Keeping last descendant of partition, adjustEarlyHeights parameter will be ignored, it is only used when keepLast == False.")
    else: ## keeping earliest tip
        labelToFirstBranch = {label: sorted(labelToBranches[label], key=sortingFxn)[0] for label in labels if labelCounts[label] > 0} ## get earliest branch

        localTree = localTree.reduce_tree(labelToFirstBranch.values()) ## keep a single earliest tip for each unique label

    
    for k in localTree.get_external(): ## iterate over remaining tips
        partition = k.traits['partition'] ## get partition label

        if keepLast == False and adjustEarlyHeights: ## kept earliest branch but want it adjusted
            lastBranch = labelToLastBranch[partition]

            k.length = lastBranch.height - k.parent.height ## new branch length is last height - current tip's parent height

            if localTree.treeType == 'time' and k.absoluteTime: ## absoluteTime set already, adjust too
                k.absoluteTime += lastBranch.absoluteTime - k.absoluteTime

        k.traits['size'] = labelCounts[partition] ## assign tip counts of this label to representative branch
        k.traits['members'] = labelToBranches[partition] ## remember descendants with label that may no longer be present

    localTree.traverse_tree()

    return localTree

def generate_calendar_timeline(startDateStr,endDateStr,spacing='monthly',dateFmt='%Y-%m-%d',roundDates=True):
    """
    Generate a list of calendar breakpoints between two dates.

    Parameters
    ----------
    startDateStr : str
        Start date of the interval.

    endDateStr : str
        End date of the interval.

    spacing : {'yearly', 'monthly', 'weekly'} or int, optional
        Calendar spacing to use. Integer values are interpreted as numbers of
        days.

    dateFmt : str, optional
        Date format used to parse inputs and format outputs.

    roundDates : bool, optional
        Whether to align the timeline to calendar boundaries when possible.

    Returns
    -------
    list[str]
        Sequence of formatted date strings.
    """
    assert spacing in ['yearly', 'monthly', 'weekly'] or isinstance(spacing, int), f"Invalid spacing {spacing}, must be int (for days) or str ('yearly', 'monthly' or 'weekly')"

    timeline = []
    startTime = dt.datetime.strptime(startDateStr, dateFmt)
    endTime = dt.datetime.strptime(endDateStr, dateFmt)

    if roundDates: ## rounding dates - additional breaks will be generated in the timeline to correspond with beginnings of months or years
        currentTime = dt.datetime(startTime.year, 1, 1) ## start from beginning of the year
        if spacing == 'yearly':
            timeline.append(dt.datetime.strftime(currentTime, dateFmt))
        elif isinstance(spacing,int):
            warnings.warn(f"Calendar timeline spacing defined as int (set to {spacing}) so roundDates (True by default) parameter ignored.")
    else: ## no rounding - timeline starts at the specified start date and is incremented at specified intervals
        currentTime = startTime
        dateStr = dt.datetime.strftime(currentTime, dateFmt)
        # timeline.append(dateStr)

    while currentTime < endTime:

        if startTime <= currentTime:
            dateStr = dt.datetime.strftime(currentTime, dateFmt)
            if dateStr not in timeline:
                timeline.append(dateStr)

        if isinstance(spacing,int):
            skip = dt.timedelta(days = spacing)

        elif spacing == 'weekly':
            skip = dt.timedelta(days = 7)
            if currentTime.month != (currentTime + skip).month: ## next month starts in a week
                daysInMonth = calendar.monthrange(currentTime.year, currentTime.month)[-1]
                lastDateOfMonth = dt.datetime(currentTime.year, currentTime.month, daysInMonth)
                dateStr = dt.datetime.strftime(lastDateOfMonth,dateFmt)
                if startTime < currentTime and currentTime != lastDateOfMonth and roundDates == True:
                    timeline.append(dateStr)

        if spacing == 'monthly':
            daysInMonth = calendar.monthrange(currentTime.year, currentTime.month)[-1]
            skip = dt.timedelta(days = daysInMonth)

        if spacing == 'yearly':
            skip = dt.timedelta(365 + calendar.isleap(currentTime.year))

        currentTime += skip

    dateStr = dt.datetime.strftime(currentTime,dateFmt)
    timeline.append(dateStr)

    return timeline




def plot_scale_bar(ax, xy, L = None, tree = None, alnL = None, textXY = None, unitText = None,
    style = 'simple', orientation = 'horizontal', ySpan = None, fancyWidth = 0.1, lineKwargs = None, textKwargs = None):
    """
    Plot a scale bar on the given axes.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes on which the scale bar will be plotted.

    xy : tuple[float, float]
        Coordinates of the starting point of the scale bar.

    L : float, optional
        Length of the scale bar. If not provided, it will be inferred from the tree or default to ``0.001``.

    tree : :class:`.Tree`, optional
        A ``baltic`` tree object used to infer the scale bar length and units if *L* is not provided.

    alnL : int, optional
        Alignment length used to convert branch lengths (in substitutions per site) to mutation counts.

    textXY : tuple[float, float], optional
        Coordinates for the scale bar label. If not provided, defaults to a position near the scale bar.

    unitText : str, optional
        Text describing the units of the scale bar. If not provided, defaults to "subs/site" for divergence trees
        or "years" for time trees.

    style : {'simple', 'fancy'}, optional
        Style of the scale bar. Defaults to ``'simple'``.

    orientation : {'horizontal', 'vertical'}, optional
        Orientation of the scale bar. Defaults to ``'horizontal'``.

    ySpan : float, optional
        Vertical span of the scale bar, used to calculate default label positions. If not provided, it will be
        inferred from the tree.

    lineKwargs : dict, optional
        Additional keyword arguments passed to the ``matplotlib`` line plotting function for the scale bar.

    textKwargs : dict, optional
        Additional keyword arguments passed to the ``matplotlib`` text plotting function for the scale bar label.

    Returns
    -------
    :obj:`matplotlib.axes.Axes`
        The modified matplotlib Axes object.

    Notes
    -----
    - If both *L* and *tree* are provided, *L* takes precedence.
    - If *alnL* is provided for a divergence tree, the scale bar will be labeled in mutation counts instead of substitutions per site.
    - The *style* parameter determines whether the scale bar has simple or fancy end markers.
    - The *orientation* parameter determines whether the scale bar is drawn horizontally or vertically.

    Raises
    ------
    ValueError
        If an invalid *style* or *orientation* is provided.

    Warnings
    --------
    - If neither *L* nor *tree* is provided, the scale bar defaults to a length of ``0.001`` with units of "subs/site".
    - If both *tree* and *ySpan* are provided, *ySpan* will be ignored in favor of the tree's *ySpan*.
    """

    assert style in ['simple', 'fancy'], f"Scale bar style {style} not recognised. Must be 'simple' or 'fancy'."
    assert orientation in ['horizontal', 'vertical'], f"Scale bar orientation {orientation} not recognised. Must be 'horizontal' or 'vertical'."

    if tree is None and ySpan is None:
        ySpan = 2e2
    elif tree is not None:
        if ySpan is not None:
            warnings.warn("Both tree and ySpan provided; using tree.ySpan.")
        ySpan = tree.ySpan

    localLineKwargs = dict(lineKwargs) if lineKwargs else {}
    localTextKwargs = dict(textKwargs) if textKwargs else {}

    if ('linewidth' not in localLineKwargs or 'lw' not in localLineKwargs): localLineKwargs['lw'] = 2
    if 'color' not in localLineKwargs: localLineKwargs['color'] = 'k'

    if orientation == 'vertical':
        if 'ha' not in localTextKwargs: localTextKwargs['ha'] = 'left'
        if 'va' not in localTextKwargs: localTextKwargs['va'] = 'center'
        localTextKwargs['rotation_mode'] = 'anchor'
    elif orientation == 'horizontal':
        if 'ha' not in localTextKwargs: localTextKwargs['ha'] = 'center'
        if 'va' not in localTextKwargs: localTextKwargs['va'] = 'top'

    n_mutations = None  # will hold equivalent mutation count if alnL is known

    if L is None:
        # No explicit L given – we choose it
        if tree is not None and tree.treeType == "divergence" and alnL is not None:
            # Divergence tree + alignment length known:
            # choose a NICE number of mutations, then convert to subs/site

            total_mut = tree.treeHeight * alnL           # total divergence in mutations
            target_mut = 0.05 * total_mut                # ~5% of total span

            if target_mut <= 0:
                target_mut = 1.0

            exponent = np.floor(np.log10(target_mut))
            fraction = target_mut / 10**exponent

            for nf in (1, 2, 5, 10):
                if fraction <= nf:
                    n_mutations = nf * 10**exponent
                    break
            else:
                n_mutations = 10 * 10**exponent

            # convert mutations -> subs/site
            L = n_mutations / alnL

        elif tree is not None:
            # No alnL: fall back to original subs/site logic
            proposed = tree.treeHeight
            exponent = np.floor(np.log10(proposed))
            fraction = (0.05 * proposed) / 10**exponent ## 5% of tree height

            for nf in (1, 2, 5, 10):  # nice numbers
                if fraction <= nf:
                    L = nf * 10**exponent
                    break
            else:
                L = 10 * 10**exponent

        else:
            # No L, no tree
            if alnL is not None:
                raise Exception(
                    "Neither scale bar length nor tree provided, "
                    "impossible to rescale scale bar to alignment length."
                )
            warnings.warn("Neither L nor tree provided; defaulting to 0.001 (subs/site).")
            L = 0.001

    else:
        # L was explicitly given
        if tree is not None:
            warnings.warn("Both L and tree provided; L takes precedence.")
        if alnL is not None and tree is not None and tree.treeType == "divergence":
            # Just record equivalent mutation count for labeling
            n_mutations = L * alnL

    if tree is not None and tree.treeType == 'time' and alnL is not None:
        warnings.warn("The tree provided has branch length units of time, the scale bar rescaling parameter alnL cannot be used and will be ignored.")

    if unitText is None:
        if tree is None:
            if alnL is None:
                warnings.warn("No units provided; assuming substitutions per site.")
                unitText = "subs/site"
            else:
                # No tree, alnL given: interpret L as subs/site and show mutations
                n_mutations = L * alnL if n_mutations is None else n_mutations
                unitText = "mutation" if np.isclose(n_mutations, 1.0) else "mutations"
        else:
            if tree.treeType == "divergence":
                if alnL is not None:
                    # Prefer to label in mutations when alnL is known
                    n_mutations = L * alnL if n_mutations is None else n_mutations
                    unitText = "mutation" if np.isclose(n_mutations, 1.0) else "mutations"
                else:
                    unitText = "subs/site"
            else:
                warnings.warn("No units provided; assuming branch lengths are years.")
                unitText = "years"

    x, y = xy

    xs = [x, x + L]
    ys = [y, y]

    if orientation == 'vertical':
        xs = [x, x]
        ys = [y, y + L]

    ax.plot(xs, ys, **localLineKwargs)

    if style == 'fancy':
        width = L * fancyWidth ## fancyWidth is a scalar of scale bar length
        left_xs, left_ys = [x, x], [y - width, y + width]
        right_xs, right_ys = [x + L, x + L], [y - width, y + width]

        if orientation == 'vertical':
            left_xs, left_ys = [x - width, x + width], [y, y]
            right_xs, right_ys = [x - width, x + width], [y + L, y + L]

        ax.plot(left_xs, left_ys, **localLineKwargs)
        ax.plot(right_xs, right_ys, **localLineKwargs)

    if textXY is None: ## set text position defaults / might be worth turning these coordinates as fractions in relation to the bar itself
        textX, textY = (x + L/2, y - ySpan * 0.02)
        if orientation == 'vertical':
            textX, textY = (x + ySpan * 0.002, y + L/2)
    else:
        textX, textY = textXY

    if alnL is not None and tree is not None and tree.treeType == "divergence":
        # Label in mutations if we know alnL
        if n_mutations is None:
            n_mutations = L * alnL
        # nice formatting
        if n_mutations >= 1:
            n_disp = int(round(n_mutations))
        else:
            n_disp = np.round(n_mutations, 2)
        scaleBarText = f'{n_disp}\n{unitText}'
    else:
        # Original behaviour: label in axis units
        scaleBarText = f'{L}\n{unitText}'

    ax.text(textX, textY, scaleBarText, **localTextKwargs)

    return ax

def _process_trait_prob_set(node, traitName):
    """
    Combine discrete trait states and probabilities into a lookup dictionary.

    Parameters
    ----------
    node : :class:`baltic.branchLike.BranchLike`
        Branch carrying ``traitName.set`` and ``traitName.set.prob`` entries.

    traitName : str
        Trait prefix used to locate the state and probability arrays.

    Returns
    -------
    dict
        Mapping from state label to posterior probability.
    """
    assert f"{traitName}.set" and f"{traitName}.set.prob" in node.traits, f"{traitName}.set or {traitName}.set.prob not found in node traits dict."

    stateSet = node.traits[f"{traitName}.set"]
    stateProbs = node.traits[f"{traitName}.set.prob"]

    stateSetDict = {key: value for key, value in zip(stateSet, stateProbs)}

    return stateSetDict

def branch_to_json(curNode, treeType, traits, mostRecentDate, treeDict=None):
    """
    Formats a baltic branchLike object into a dict that will be converted into auspice JSON
    """

    if treeDict is None:
        treeDict = {}
    
    nodeDict = {}
    nodeDict['node_attrs'] = {}

    #####
    if treeType == 'divergence':
        nodeDict['node_attrs']['div'] = curNode.height
    elif treeType == 'time':
        nodeDict['node_attrs']['num_date'] = {'value': curNode.absoluteTime}

    ##### rename nodes and leaves, recurse through children
    if curNode.is_node():
        if 'node_idx' in curNode.traits:
            nodeDict['name'] = curNode.traits['node_idx']
        else:
            logger.warning(f"Node does not have a name derived from pre-order traversal (str expected under trait key 'node_idx'). Importing this JSON into auspice.us will be fine but re-importing in baltic will have issues.")

        nodeDict['children'] = []

        for childNode in curNode.children:
            nodeDict['children'].append(branch_to_json(curNode=childNode, treeType=treeType, mostRecentDate=mostRecentDate, traits=traits, treeDict=treeDict))
    
    elif curNode.is_leaf():
        nodeDict['name'] = curNode.name
    
    else:
        logger.error(f"Attempted to convert baltic branchLike object of type {type(curNode)} which is not leaf or node which is currently not supported.")
        raise TypeError()
    #####
    
    for trait in traits: ## format traits
        if trait in curNode.traits:
            nodeDict['node_attrs'][trait] = {'value': curNode.traits[trait]}

            ### assign uncertainties
            if f"{trait}.set.prob" in curNode.traits: ## discrete traits
                traitProbs = _process_trait_prob_set(curNode, trait)
                nodeDict['node_attrs'][trait]['confidence'] = traitProbs
            
            elif f"{trait}_95%_HPD" in curNode.traits: ## continuous traits
                nodeDict['node_attrs'][trait]['confidence'] = curNode.traits[f"{trait}_95%_HPD"]

                if trait == 'height' and treeType == 'time' and mostRecentDate: ## special case - working with time tree
                    nodeDict['node_attrs']['num_date']['confidence'] = [mostRecentDate - value for value in curNode.traits["height_95%_HPD"]]
            
            elif trait in curNode.traits:
                nodeDict['node_attrs'][trait] = {'value': curNode.traits[trait]}
    ########
    return nodeDict

def plot_node_bar(ax, node, traitName, traitColourDict, xyFxn = None, height = 10, width = 0.2, otherThres = 0.0, connectNode = True, connectingCorner = 'lower middle', orientation = 'vertical', **kwargs):
    """
    Plot a stacked bar summarizing discrete trait probabilities for a node.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes on which the bar should be drawn.

    node : :class:`baltic.branchLike.BranchLike`
        Branch whose trait probabilities should be displayed.

    traitName : str
        Trait prefix used to locate ``.set`` and ``.set.prob`` values.

    traitColourDict : dict
        Mapping from trait state to display colour.

    xyFxn : callable, optional
        Function returning the anchor coordinates for the bar.

    height : float, optional
        Total span of the stacked bar.

    width : float, optional
        Width of the bar orthogonal to *height*.

    otherThres : float, optional
        Probability threshold below which states are grouped into ``other``.

    connectNode : bool, optional
        If ``True``, draw a dashed connector back to the node location.

    connectingCorner : str, optional
        Corner of the bar used as the connector origin.

    orientation : {'vertical', 'horizontal'}, optional
        Orientation of the stacked bar.

    ``**kwargs``
        Additional keyword arguments forwarded to
        :class:`matplotlib.patches.Rectangle`.
    """
    from matplotlib.patches import Rectangle

    assert f"{traitName}.set" and f"{traitName}.set.prob" in node.traits, f"{traitName}.set or {traitName}.set.prob not found in node traits dict."
    assert otherThres < 1.0, f"Threshold for assigning state to 'other' category ({otherThres}) should be <1.0."

    if xyFxn is None: xyFxn = lambda k: (k.x, k.y)

    stateSetDict = _process_trait_prob_set(node, traitName)

    stateOrder = sorted(stateSetDict.keys(), key = lambda state: -stateSetDict[state])

    otherCategory = [state for state in stateOrder if stateSetDict[state] <= otherThres]

    localKwargs = dict(kwargs)

    if len(otherCategory) > 0:
        stateSetDict['other'] = sum([stateSetDict[state] for state in otherCategory])
        for state in otherCategory:
            stateSetDict.pop(state)
            stateOrder.remove(state)
        stateOrder.append('other')

    if 'other' not in traitColourDict: traitColourDict['other'] = '#f5f5f5'

    probs = [stateSetDict[state] for state in stateOrder]
    cs = [traitColourDict[state] for state in stateOrder]

    for i, (state, prob) in enumerate(zip(stateOrder, probs)):
        x, y = xyFxn(node)

        cumulativeProb = sum(probs[:i])
        rectHeight = height * prob

        if 'zorder' not in localKwargs: localKwargs['zorder'] = 1

        if orientation == 'horizontal':
            rectBottom = x + height * cumulativeProb
            rect = Rectangle((rectBottom, y), width = rectHeight, height = width, fc = traitColourDict[state], **localKwargs)

        elif orientation == 'vertical':
            rectBottom = y + height * cumulativeProb
            rect = Rectangle((x, rectBottom), width = width, height = rectHeight, fc = traitColourDict[state], **localKwargs)

        ax.add_patch(rect)

    if connectNode:
        vaCorner, haCorner = connectingCorner.split(' ')
        assert vaCorner in ['upper', 'lower'], f"Vertical corner connection parameter {vaCorner} not recognised. Must be 'lower' or 'upper'"
        assert haCorner in ['left', 'middle', 'right'], f"Horizontal corner connection parameter {haCorner} not recognised. Must be 'left', 'middle' or 'right'"

        x, y = xyFxn(node)

        if vaCorner == 'upper':
            y += height if orientation == 'vertical' else width

        if haCorner == 'right':
            x += width if orientation == 'vertical' else height
        elif haCorner == 'middle':
            x += width/2 if orientation == 'vertical' else height/2

        xs, ys = zip(*[(x, y), (x, node.y), (node.x, node.y)])

        ax.plot(xs, ys, ls = '--', color = 'dimgray', zorder = 0)


def plot_node_treemap(ax, node, traitName, traitColourDict, height, width, centerFxn = None, area = 1.0, other_thres = 0.0, **kwargs):
    """
    Plot a treemap summarizing discrete trait probabilities for a node.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes on which the treemap should be drawn.

    node : :class:`baltic.branchLike.BranchLike`
        Branch whose trait probabilities should be displayed.

    traitName : str
        Trait prefix used to locate ``.set`` and ``.set.prob`` values.

    traitColourDict : dict
        Mapping from trait state to display colour.

    height, width : float
        Size of the treemap rectangle.

    centerFxn : callable, optional
        Function returning the rectangle center.

    area : float, optional
        Included for API compatibility with related plotting helpers.

    other_thres : float, optional
        Probability threshold below which states are grouped into ``other``.

    ``**kwargs``
        Additional keyword arguments forwarded to
        :class:`matplotlib.patches.Rectangle`.
    """

    assert other_thres < 1.0, f"Threshold for assigning state to 'other' category ({other_thres}) should be <1.0."

    import squarify
    from matplotlib.patches import Rectangle

    if centerFxn is None: centerFxn = lambda k: (k.x, k.y)

    stateSetDict = _process_trait_prob_set(node, traitName)

    stateOrder = sorted(stateSetDict.keys(), key = lambda state: -stateSetDict[state])

    otherCategory = [state for state in stateOrder if stateSetDict[state] <= other_thres]

    if len(otherCategory) > 0:
        stateSetDict['other'] = sum([stateSetDict[state] for state in otherCategory])
        for state in otherCategory:
            stateSetDict.pop(state)
            stateOrder.remove(state)
        stateOrder.append('other')

    if 'other' not in traitColourDict: traitColourDict['other'] = '#f5f5f5'

    probs = [stateSetDict[state] for state in stateOrder]
    cs = [traitColourDict[state] for state in stateOrder]

    x, y = centerFxn(node)

    x -= width/2
    y -= height/2

    values = squarify.normalize_sizes(probs, width, height)
    rects = squarify.squarify(values, x, y, width, height)

    for i,rect in enumerate(rects):
        rectPatch = Rectangle((rect['x'], rect['y']), rect['dx'], rect['dy'], fc = cs[i], **kwargs)
        ax.add_patch(rectPatch)

    ax.autoscale()
    return ax

def plot_node_piechart(ax, node, traitName, traitColourDict, centerFxn = None, radius = 0.5, other_thres = 0.0, **kwargs):
    """
    Plot a pie chart summarizing discrete trait probabilities for a node.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes on which the pie chart should be drawn.

    node : :class:`baltic.branchLike.BranchLike`
        Branch whose trait probabilities should be displayed.

    traitName : str
        Trait prefix used to locate ``.set`` and ``.set.prob`` values.

    traitColourDict : dict
        Mapping from trait state to display colour.

    centerFxn : callable, optional
        Function returning the chart center.

    radius : float, optional
        Pie chart radius.

    other_thres : float, optional
        Probability threshold below which states are grouped into ``other``.

    ``**kwargs``
        Additional keyword arguments forwarded to
        :meth:`matplotlib.axes.Axes.pie`.
    """
    assert f"{traitName}.set" and f"{traitName}.set.prob" in node.traits, f"{traitName}.set or {traitName}.set.prob not found in node traits dict."
    assert other_thres < 1.0, f"Threshold for assigning state to 'other' category ({other_thres}) should be <1.0."

    if centerFxn is None: centerFxn = lambda k: (k.x, k.y)

    stateSetDict = _process_trait_prob_set(node, traitName)

    stateOrder = sorted(stateSetDict.keys(), key = lambda state: -stateSetDict[state])

    otherCategory = [state for state in stateOrder if stateSetDict[state] <= other_thres]

    if len(otherCategory) > 0:
        stateSetDict['other'] = sum([stateSetDict[state] for state in otherCategory])
        for state in otherCategory:
            stateSetDict.pop(state)
            stateOrder.remove(state)
        stateOrder.append('other')
    
    if 'other' not in traitColourDict: traitColourDict['other'] = '#f5f5f5'

    probs = [stateSetDict[state] for state in stateOrder]
    cs = [traitColourDict[state] for state in stateOrder]

    ax.pie(x = probs, colors = cs, radius = radius, center = centerFxn(node), **kwargs)

    ax.autoscale()
    return ax

def plot_tmrca_posterior(ax, tmrcaFile, tmrcaName = 'age(root)', burnin = None, yCoord = None, fullViolin = True,
                         hpdLvl = 0.95, precision = 100, kdeWidth = 3, normalise=True, orientation = 'horizontal', connectNode = False, node = None, violinKwargs = {}, outlineKwargs = {}, connectionLineKwargs = {}):
    """
    Plot a KDE-based posterior density for a TMRCA statistic from a log file.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes on which the posterior should be drawn.

    tmrcaFile : str
        Path to the tab-delimited log file containing posterior samples.

    tmrcaName : str, optional
        Column name to extract from the log file.

    burnin : int, optional
        Minimum state value to retain from the log.

    yCoord : float, optional
        Anchor coordinate for plotting the density.

    fullViolin : bool, optional
        If ``True``, draw the full violin; otherwise draw a half violin.

    hpdLvl : float, optional
        Highest posterior density mass to report.

    precision : int, optional
        Number of x positions used to evaluate the KDE.

    kdeWidth : float, optional
        Width scaling applied to the KDE curve.

    normalise : bool, optional
        If ``True``, normalise the KDE height before scaling by *kdeWidth*.

    orientation : {'horizontal', 'vertical'}, optional
        Orientation of the violin plot.

    connectNode : bool, optional
        If ``True``, connect the posterior summary back to *node*.

    node : :class:`baltic.branchLike.BranchLike`, optional
        Branch to connect to when *connectNode* is enabled.

    violinKwargs, outlineKwargs, connectionLineKwargs : dict, optional
        Keyword arguments forwarded to the violin fill, outline, and connector
        line artists.

    Returns
    -------
    :obj:`matplotlib.axes.Axes`
        The modified matplotlib Axes object.
    """

    ### certain tree orientations will require xCoord too
    ### connect mean of HPD to node with dotted line (accommodate elbow in case KDE is on the x-axis)
    import csv
    from scipy.stats import gaussian_kde
    # from baltic.bt_utils import hpd, decimal_to_calendar_date


    func_logger = logging.getLogger("baltic.bt_utils.plot_tmrca_posterior")
    func_logger.setLevel(logging.INFO)
    func_logger.propagate = False

    if not func_logger.handlers:
        import sys
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[plot_tmrca_posterior] %(message)s"))
        handler.setLevel(logging.INFO)
        func_logger.addHandler(handler)

    tmrcaPosterior = []

    if burnin == None:
        burnin = 10e6
        logger.warning('No burnin set, defaulting to 10M states.')

    if node: assert (node.is_leaflike() or node.is_node()), f"Provided node object is {type(node)}, not a baltic branchLike object."
    handle = open(tmrcaFile, 'r')

    for l in csv.DictReader((line for line in handle if line.startswith('#') == False), delimiter = '\t'):
        state = int(l['state'])
        if state >= burnin:
            assert tmrcaName in l, f"'{tmrcaName}' not found in log file provided. Log file header looks like this: {', '.join(list(l.keys()))}"
            tmrcaPosterior.append(float(l[tmrcaName])) ## grab column with tmrca stat

    assert len(tmrcaPosterior) > 0, f"No TMRCA values loaded. Is burnin of {burnin} too high?"
    handle.close()

    if yCoord is None and node is None:
        yCoord = 0.0
    elif yCoord is None and node:
        yCoord = node.y
    elif yCoord is not None and node is not None:
        logger.warning(f"Both yCoord and node were provided, KDE will be positioned at yCoord value.")

    localViolinKwargs = dict(violinKwargs)
    if 'fc' in localViolinKwargs:
        localViolinKwargs['facecolor'] = localViolinKwargs['fc']
        localViolinKwargs.pop('fc')
    if 'ec' in localViolinKwargs:
        localViolinKwargs['edgecolor'] = localViolinKwargs['ec']
        localViolinKwargs.pop('ec')
    if 'facecolor' not in localViolinKwargs and 'fc' not in localViolinKwargs: localViolinKwargs['facecolor'] = 'gray'
    if 'edgecolor' not in localViolinKwargs and 'ec' not in localViolinKwargs: localViolinKwargs['edgecolor'] = 'none'
    if 'alpha' not in localViolinKwargs: localViolinKwargs['alpha'] = 0.1
    if 'zorder' not in localViolinKwargs: localViolinKwargs['zorder'] = 1

    localOutlineKwargs = dict(outlineKwargs)
    if 'color' not in localOutlineKwargs: localOutlineKwargs['color'] = 'gray'
    if 'linewidth' not in localOutlineKwargs and 'lw' not in localOutlineKwargs: localOutlineKwargs['linewidth'] = 2
    if 'zorder' not in localOutlineKwargs: localOutlineKwargs['zorder'] = 2

    localConnectionLineKwargs = dict(connectionLineKwargs)
    if 'color' not in localConnectionLineKwargs: localConnectionLineKwargs['color'] = 'dimgray'
    if 'linewidth' not in localConnectionLineKwargs and 'lw' not in localConnectionLineKwargs: localConnectionLineKwargs['linewidth'] = 2
    if 'zorder' not in localConnectionLineKwargs: localConnectionLineKwargs['zorder'] = 1
    if 'linestyle' not in localConnectionLineKwargs and 'ls' not in localConnectionLineKwargs: localConnectionLineKwargs['linestyle'] = '--'

    kde = gaussian_kde(tmrcaPosterior)
    hpdLo, hpdHi = hpd(tmrcaPosterior, hpdLvl)
    x_grid = np.linspace(hpdLo, hpdHi, precision)

    meanTmrca = np.mean(tmrcaPosterior)
    medianTmrca = np.median(tmrcaPosterior)

    func_logger.info(f"Node {tmrcaName} mean TMRCA: {meanTmrca:.3f} ({decimal_to_calendar_date(meanTmrca)}) median TMRCA: {medianTmrca:.3f} ({decimal_to_calendar_date(medianTmrca)})")
    func_logger.info(f"Node {tmrcaName} TMRCA 95% HPD: {hpdLo:.3f} - {hpdHi:.3f} / {decimal_to_calendar_date(hpdLo)} - {decimal_to_calendar_date(hpdHi)}")

    y_grid = kde(x_grid)
    y_max = y_grid.max()
    if normalise:
        y_grid /= y_max ## normalise KDE to peak at 1.0

    y_grid *= kdeWidth ## rescale

    upper_ys = [yCoord + y for y in y_grid]
    lower_ys = [yCoord - y for y in y_grid]
    constant_ys = [yCoord for _ in y_grid]

    if orientation == 'horizontal':
        plotKDE = ax.fill_between
    elif orientation == 'vertical':
        plotKDE = ax.fill_betweenx

    xs, uys, lys = x_grid, upper_ys, lower_ys

    if fullViolin:
        plotKDE(xs, uys, lys, **localViolinKwargs)
    else:
        plotKDE(xs, uys, constant_ys, **localViolinKwargs)

    if fullViolin:
        if orientation == 'vertical':
            ax.plot(uys, xs,**localOutlineKwargs)
            ax.plot(lys, xs,**localOutlineKwargs)
        elif orientation == 'horizontal':
            ax.plot(xs, uys,**localOutlineKwargs)
            ax.plot(xs, lys,**localOutlineKwargs)
    else:
        if orientation == 'vertical':
            ax.plot(uys, xs,**localOutlineKwargs)
        elif orientation == 'horizontal':
            ax.plot(xs, uys,**localOutlineKwargs)

    if connectNode and node:
        x, y = meanTmrca, yCoord

        mean_xs = [meanTmrca, meanTmrca]
        if fullViolin:
            if normalise:
                mean_ys = [yCoord - kde(meanTmrca) / y_max * kdeWidth, yCoord + kde(meanTmrca) / y_max * kdeWidth]
            else:
                mean_ys = [yCoord - kde(meanTmrca) * kdeWidth, yCoord + kde(meanTmrca) * kdeWidth]
        else:
            if normalise:
                mean_ys = [yCoord, yCoord + kde(meanTmrca).item() / y_max * kdeWidth]
            else:
                mean_ys = [yCoord, yCoord + kde(meanTmrca).item() * kdeWidth]    

            mean_ys = np.array(mean_ys, dtype=float).tolist()

        elbow_xs = [meanTmrca, meanTmrca, node.absoluteTime]
        elbow_ys = [yCoord, node.y, node.y]

        if orientation == 'vertical':
            x, y = y, x
            mean_xs, mean_ys = mean_ys, mean_xs
            elbow_xs, elbow_ys = elbow_ys, elbow_xs

        ax.scatter(x, y, s = 40, fc = localViolinKwargs['facecolor'], ec = 'none', zorder = 10)
        ax.scatter(x, y, s = 80, fc = localViolinKwargs['edgecolor'], ec = 'none', zorder = 9)

        ax.plot(mean_xs, mean_ys, **localConnectionLineKwargs)
        ax.plot(elbow_xs, elbow_ys, **localConnectionLineKwargs)
    elif connectNode and node is None:
        logger.warning(f"Need to specify node to connect to.")

    return ax

def plot_time_grid(ax, timeline, dateFmt='%Y-%m-%d', colourFxn=None, colour=None, edgeColourFxn=None, edgeColour=None, axis='x',**kwargs):
    """
    Shade alternating time intervals on an axis.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes on which the spans should be drawn.

    timeline : list[str] or range
        Time boundaries as calendar strings or numeric values.

    dateFmt : str, optional
        Date format used when *timeline* contains strings.

    colourFxn, edgeColourFxn : callable, optional
        Functions that map interval indices to face and edge colours.

    colour, edgeColour : color, optional
        Constant face and edge colours used when the corresponding functions
        are not provided.

    axis : {'x', 'y'}, optional
        Axis along which spans should be drawn.

    ``**kwargs``
        Additional keyword arguments forwarded to ``ax.axvspan`` or
        ``ax.axhspan``.

    Returns
    -------
    :obj:`matplotlib.axes.Axes`
        The modified matplotlib Axes object.
    """

    if colour is not None and colourFxn is not None:
        raise ValueError(
            "Cannot specify both colour and colourFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if colour is None and colourFxn is None:
        colourFxn = lambda k: "k"
    elif colourFxn is None:
        colourFxn = lambda k: colour

    if edgeColour is not None and edgeColourFxn is not None:
        raise ValueError(
            "Cannot specify both edgeColour and edgeColourFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if edgeColour is None and edgeColourFxn is None:
        edgeColourFxn = lambda k: "none"
    elif edgeColourFxn is None:
        edgeColourFxn = lambda k: edgeColour

    localKwargs = dict(kwargs)
    if 'alpha' not in localKwargs: localKwargs['alpha'] = 0.08

    if isinstance(timeline,list):
        try:
            timeline = [calendar_to_decimal_date(t,fmt=dateFmt) for t in timeline] ## convert timeline to
        except:
            warnings.warn(f"List of timeline dates are not recognised. Expected date format: {dateFmt}, first entry in list: {timeline[0]}.")
    else:
        assert isinstance(timeline,range), f"timeline is neither a list nor a range."

    if axis == 'x':
        [ax.axvspan(timeline[t], timeline[t+1], fc=colourFxn(t), ec=edgeColourFxn(t), **localKwargs) for t in range(0,len(timeline)-1,2)]
    elif axis == 'y':
        [ax.axhspan(timeline[t], timeline[t+1], fc=colourFxn(t), ec=edgeColourFxn(t), **localKwargs) for t in range(0,len(timeline)-1,2)]
    return ax

def format_time_grid(ax, timeline, inputDateFmt='%Y-%m-%d', outputFmtFxn=None, labelPosition='mid', axis='x', **kwargs):
    """
    Format tick positions and labels for a time-grid axis.

    Parameters
    ----------
    ax : :obj:`matplotlib.axes.Axes`
        Axes whose ticks should be updated.

    timeline : list[str]
        Ordered list of calendar dates defining grid boundaries.

    inputDateFmt : str, optional
        Date format used to parse entries in *timeline*.

    outputFmtFxn : callable, optional
        Function used to convert each timeline entry into a label string.

    labelPosition : {'left', 'mid'}, optional
        Whether labels should be placed on boundaries or interval midpoints.

    axis : {'x', 'y'}, optional
        Axis whose ticks should be updated.

    ``**kwargs``
        Additional keyword arguments forwarded to the tick label setters.

    Returns
    -------
    :obj:`matplotlib.axes.Axes`
        The modified matplotlib Axes object.
    """
    assert labelPosition in ['left', 'mid'], f"labelPosition {labelPosition} invalid. Must be 'left' or 'mid'"
    assert axis in ['x', 'y'], f"axis {axis} invalid. Must be 'x' or 'y'"

    if outputFmtFxn is None:
        outputFmtFxn = lambda date: convert_date_format(date, '%Y-%m-%d', '%b\n%Y') if convert_date_format(date, '%Y-%m-%d', '%m') == '01' else convert_date_format(date, '%Y-%m-%d', '%b')

    localKwargs = dict(kwargs)

    if axis == 'x':
        if labelPosition == 'left':
            ax.set_xticks([calendar_to_decimal_date(date, inputDateFmt) for date in timeline])
            ax.set_xticklabels([outputFmtFxn(date) for date in timeline],**localKwargs)
        elif labelPosition == 'mid':
            ax.set_xticks([np.mean([calendar_to_decimal_date(timeline[t], inputDateFmt), calendar_to_decimal_date(timeline[t+1], inputDateFmt)]) for t in range(len(timeline)-1)])
            ax.set_xticklabels([outputFmtFxn(date) for date in timeline[:-1]],**localKwargs)
    elif axis == 'y':
        if labelPosition == 'left':
            ax.set_yticks([calendar_to_decimal_date(date, inputDateFmt) for date in timeline])
            ax.set_yticklabels([outputFmtFxn(date) for date in timeline],**localKwargs)
        elif labelPosition == 'mid':
            ax.set_yticks([np.mean([calendar_to_decimal_date(timeline[t], inputDateFmt), calendar_to_decimal_date(timeline[t+1], inputDateFmt)]) for t in range(len(timeline)-1)])
            ax.set_yticklabels([outputFmtFxn(date) for date in timeline[:-1]],**localKwargs)

    ax.tick_params(axis=axis, size=0)
    return ax

def clean_axes(ax, hideSpines = ['left', 'top', 'right', 'bottom'], removeTickLabels = 'both'):
    """
    Remove selected spines, suppress ticks and ticklabels on x, y or both axes.
    """
    validSpines = ['left', 'top', 'right', 'bottom']
    assert set(validSpines) >= set(hideSpines), f"Spine {[val for val in hideSpines if val not in validSpines]} not recognised. Must belong to the set {validSpines}."

    validRemoveTickLabels = ['x', 'y', 'both', 'none']
    assert removeTickLabels in validRemoveTickLabels, f"removeTickLabels value {removeTickLabels} not recognised. Must be one of {', '.join(validRemoveTickLabels)}"

    if removeTickLabels in ['x', 'both']:
        ax.set_xticks([])
        ax.set_xticklabels([])
    if removeTickLabels in ['y', 'both']:
        ax.set_yticks([])
        ax.set_yticklabels([])

    [ax.spines[loc].set_visible(False) for loc in hideSpines]
    return ax



def untangle(tree, reference, min_shared=2, maxPolytomy=9):
    """
    Reorder internal node children across multiple trees to reduce tip crossing.

    Parameters
    ----------
    trees : list[:class:`.Tree`]
        Trees to reorder in sequence.

    costFxn : callable, optional
        Function used to score the discrepancy between the y-positions of a
        shared tip in adjacent trees. The function must accept a pair of
        numeric positions and return a scalar cost. If omitted, squared
        distance is used.

    iterations : int, optional
        Number of untangling passes to perform. Defaults to ``3``.

    Returns
    -------
    list[:class:`.Tree`]
        The input trees with child orderings updated in place.

    Raises
    ------
    RuntimeWarning
        If a node has ten or more children, making exhaustive permutation
        search impractical.
    """
    from itertools import permutations

    # Ensure coordinates exist
    reference._assign_tree_coordinates()
    tree._assign_tree_coordinates()

    y_ref = {k.name: k.y for k in reference.get_external()}
    y_tr  = {k.name: k.y for k in tree.get_external()}

    # Bottom-up so child decisions propagate correctly
    for node in sorted(tree.get_internal(), key=lambda n: -n.height):

        k = len(node.children)
        if k < 2:
            continue
        if k > maxPolytomy:
            continue  # avoid factorial explosion

        # Precompute descendant sets
        child_sets = [{ch.name} if ch.is_leaf() else ch.leaves for ch in node.children]

        # Reference means for each child-set
        ref_means = []
        for S in child_sets:
            shared = S & y_ref.keys()
            if len(shared) < min_shared:
                ref_means.append(None)
            else:
                ref_means.append(np.mean([y_ref[t] for t in shared if t in y_ref]))

        if any(m is None for m in ref_means):
            continue

        best_perm = None
        best_score = None

        for perm in permutations(range(k)):
            perm_sets = [child_sets[i] for i in perm]

            tree_means = []
            for S in perm_sets:
                shared = S & y_tr.keys()
                if len(shared) < min_shared:
                    break
                tree_means.append(np.mean([y_tr[t] for t in shared if t in y_tr]))
            else:
                # Score = L1 distance between ordered means
                score = sum(
                    abs(tree_means[i] - ref_means[i])
                    for i in range(k)
                )

                if best_score is None or score < best_score:
                    best_score = score
                    best_perm = perm

        if best_perm is not None:
            node.children = [node.children[i] for i in best_perm]

            # Commit geometry change
            tree._assign_tree_coordinates()
            y_tr = {k.name: k.y for k in tree.get_external()}

    return tree


def untangle_trees(
    trees,
    iterations=10,
    maxPolytomy=8,
    bidirectional=True
):
    """
    Untangle a list of trees for tanglegram visualisation.

    Parameters
    ----------
    trees : list[Tree]
        Trees ordered as they will appear in the tanglegram.
        Trees are modified in place.
    iterations : int
        Number of global passes along the chain.
    costFxn : callable
        Function mapping (y_ref, y_tree) -> cost.
    maxPoly : int
        Maximum polytomy size to brute-force.
    bidirectional : bool
        Whether to do backward passes as well as forward passes.

    Returns
    -------
    list[Tree]
        The same list, untangled.
    """

    logger.warning(f"Untangled trees can be misleading! Tangle lines can be perfectly parallel without much topological congruence between trees. Use at your own peril.")
    
    if len(trees) < 2:
        raise Exception(f"List of trees contains {len(trees)} trees, need at least 2.")

    # Ensure all trees start with valid coordinates
    for t in trees:
        t._assign_tree_coordinates()

    for _ in range(iterations):

        # ---------- forward pass ----------
        for i in range(1, len(trees)):
            untangle(
                tree=trees[i],
                reference=trees[i - 1],
                maxPolytomy=maxPolytomy
            )

        # ---------- backward pass ----------
        if bidirectional:
            for i in range(len(trees) - 2, -1, -1):
                untangle(
                    tree=trees[i],
                    reference=trees[i + 1],
                    maxPolytomy=maxPolytomy
                )

    return trees


def unnest(nodeList, towardsRoot = True):
    """
    Remove nested nodes from a selection until descendant sets no longer overlap.

    Parameters
    ----------
    nodeList : iterable[:class:`baltic.branchLike.BranchLike`]
        Nodes or leaf-like branches to filter.

    towardsRoot : bool, optional
        If ``True``, preferentially keep deeper nodes; otherwise keep more
        tip-proximal entries.

    Returns
    -------
    list
        Filtered list with nested overlaps removed.
    """
    nodeList = list(nodeList) ## make sure removal of nested nodes is not in-place
    assert all([(k.is_node() or k.is_leaflike()) for k in nodeList]), f"nodeList contains objects that are not baltic branch objects (node or leaflike): {', '.join([k for k in nodeList if k.is_node() == False and k.is_leaflike() == False])}"

    while any([A.leaves.isdisjoint(B.leaves) == False for A in nodeList for B in nodeList if A != B]): ## continue looping for as long as any pair of nodes are nested
        remove = set() ## store nodes for removal
        for A in nodeList: ## iterate over nodes once (A variable)
            for B in nodeList: ## iterate over nodes twice (B variable)
                if A != B and A.leaves.isdisjoint(B.leaves) == False: ## if descendant tips overlap between the two nodes
                    if towardsRoot:
                        remove.add(B if B.leaves.issubset(A.leaves) else A) ## keep nodes deeper in the tree - remove B if node B is subset of node A (remove nodes closer to tips)
                    else:
                        remove.add(A if B.leaves.issubset(A.leaves) else B) ## keep nodes closer to tips - remove A if node B is subset of node A (removing nodes deeper in the tree)

        for r in remove:
            nodeList.remove(r) ## remove designated nodes from list

    return nodeList ## when done return list of remaining nodes

def _root_to_tip(rootCandidate, tipDates, tipHeights, res, stat='r^2', forcePositive=True, frac=None):
    """
    Evaluate a root-to-tip regression for a candidate root position.

    Parameters
    ----------
    rootCandidate : :class:`baltic.branchLike.BranchLike`
        Candidate root branch or node.

    tipDates : sequence[float]
        Sampling dates of the tips.

    tipHeights : sequence[float]
        Root-to-tip distances corresponding to *tipDates*.

    res : dict
        Mutable result dictionary storing the current best regression summary.

    stat : {'r^2', 'correlation', 'sum of squares'}, optional
        Optimisation criterion used to compare candidate roots.

    forcePositive : bool, optional
        If ``True``, invalidate regressions with negative slopes.

    frac : float, optional
        Optional position along a branch associated with *rootCandidate*.

    Returns
    -------
    dict
        Updated regression summary.
    """
    slope, intercept, rval, _, _ = linregress(tipDates, tipHeights) ## run linear regression
    corr = np.corrcoef((tipDates, tipHeights))[0,1] ## correlation coefficient
    ssq = sum([(y - (slope * x + intercept))**2 for x, y in zip(tipDates, tipHeights)]) ## sum of squares

    if stat=='correlation': ## set stat to optimise
        localStat = corr
    elif stat=='sum of squares':
        localStat = ssq
    elif stat=='r^2':
        localStat = math.pow(rval, 2)
    else:
        raise ValueError(f'Unknown stat {stat} to optimise')

    optStat = res[stat] if stat in res else (np.inf if stat in ['sum of squares'] else -np.inf)

    invalidateRegression = True if forcePositive and slope < 0 else False

    if (localStat < optStat if stat in ['sum of squares'] else localStat > optStat): ## minimise sum of squares or maximise correlation/r^2
        optStat = localStat ## better root found
        res['correlation'] = -np.inf if invalidateRegression else corr
        res['sum of squares'] = np.inf if invalidateRegression else ssq
        res['slope'] = slope
        res['intercept'] = intercept
        res['r^2'] = -np.inf if invalidateRegression else math.pow(rval, 2)
        res['root'] = rootCandidate
        if frac is not None: res['frac'] = frac

    return res

def _rtt_worker(args):
    """
    Evaluate one candidate root in a Monte Carlo root-to-tip analysis.

    Parameters
    ----------
    args : tuple
        Tuple containing:

        - ``tree``: a tree object to copy and reroot.
        - ``root_index``: object index for the candidate root.
        - ``tipDates``: mapping of tip name to fixed sampling date.
        - ``uncertainDateRanges``: mapping of tip name to ``(min, max)`` date
          interval for uncertain tips.
        - ``n_mc``: number of Monte Carlo iterations.
        - ``stat``: optimization statistic, one of ``'r^2'``,
          ``'correlation'``, or ``'sum of squares'``.
        - ``forcePositive``: whether negative regression slopes should be
          rejected.

    Returns
    -------
    dict
        Best regression result found for the candidate root, with the root
        stored by index rather than by object reference.
    """

    (tree,
     root_index,
     tipDates,
     uncertainDateRanges,
     max_iters,
     stat,
     forcePositive) = args

    # deep-copy tree and reroot on candidate
    tree_copy = copy.deepcopy(tree)
    candidate = next(obj for obj in tree_copy.Objects if obj.index == root_index)
    cll = tree_copy.reroot(candidate)

    tips = cll.get_external()
    uncertainTips = [k for k in tips if k.name in uncertainDateRanges]

    # root-to-tip distances are fixed for a given root
    rootDistances = {k.name: k.height for k in tips}

    # start from a local copy of tipDates
    # (each worker gets its own copy due to pickling, but we clone anyway to be explicit)
    currentDates = dict(tipDates)

    # storage for the best result across iterations
    best_res = None
    best_score = None
    best_dates = None

    # tolerance for date convergence
    tol = 1e-6

    # main refinement loop
    for _ in range(max(1, int(max_iters))):

        # build regression inputs with current dates
        xs = [currentDates[k.name] for k in tips]               # dates
        ys = [rootDistances[k.name] for k in tips]               # heights

        res_local = {}
        res_local = _root_to_tip(
            candidate, xs, ys, res_local,
            stat=stat, forcePositive=forcePositive
        )

        # if regression is invalid (e.g., negative slope with forcePositive=True), bail out
        if not res_local:
            break

        # compute scalar score (same logic as before)
        if stat == "sum of squares":
            score = -res_local["sum of squares"]  # smaller SSE is better
        else:
            score = res_local[stat]               # larger r^2 / correlation is better


        # track best iteration so far
        if best_score is None or score > best_score:
            best_score = score
            best_res = res_local
            # store only the uncertain ones (like before with monte_carlo_dates)
            best_dates = {name: currentDates[name] for name in uncertainDateRanges}

        # if there are no uncertain dates, we are done after one regression
        if not uncertainDateRanges:
            break

        # update uncertain dates by projecting onto the regression line
        slope = res_local["slope"]
        intercept = res_local["intercept"]

        # if slope is zero, projection is undefined – nothing to refine
        if abs(slope) < 1e-12:
            break

        newDates = _adjust_tip_dates_by_regression(uncertainTips, slope, intercept)
        delta = 0.0
        for tipName in newDates:
            delta += abs(currentDates[tipName] - newDates[tipName])
            currentDates[tipName] = newDates[tipName]

        # stop if updates are tiny (converged)
        if delta < tol:
            break


    # if we never got a valid regression, return an empty-ish record
    if best_res is None:
        best_res = {}
        best_score = None
        best_dates = None

    # attach metadata expected by root_by_regression
    best_res["root_index"] = root_index
    best_res["score"] = best_score
    best_res["assigned_uncertain_dates"] = best_dates

    # IMPORTANT: replace 'root' node reference with index so we don't leak deep-copied nodes
    best_res["root"] = root_index

    return best_res

def _adjust_tip_dates_by_regression(
    uncertainTips,
    slope: float,
    intercept: float):
    """
    Receives a list of tips with uncertain dates and regression parameters, adjusts .absoluteTime (within uncertainty range .absoluteTimeRange) to be as close to regression line as possible.
    """
    adjustedDates = {}

    for k in uncertainTips:
        regressionDate = (k.height - intercept) / slope ## interpolate date of tip based on subs/site height
        minD, maxD = k.absoluteTimeRange ## get clamp

        adjustedDates[k.name] = min(max(regressionDate, minD), maxD) ## clamp tip date

    return adjustedDates


def project_to_polar(x, y, yRange, circleStart=0.0, circleFraction=1.0):
    """
    Convert rectangular tree coordinates to Cartesian coordinates on a circle.

    Parameters
    ----------
    x : float
        Radial distance from the origin.

    y : float
        Position along the non-informative tree axis.

    yRange : float
        Total span of the non-informative axis.

    circleStart : float, optional
        Starting angular offset as a fraction of a full turn.

    circleFraction : float, optional
        Fraction of the circle used by the layout.

    Returns
    -------
    tuple[float, float]
        Projected Cartesian coordinates.
    """

    # circle_start_radians = 2*math.pi * circleStart ## convert starting point to radians
    # circle_fraction_radians = 2*math.pi * circleFraction ## convert arc width to radians

    # rads = circle_start_radians + (circle_fraction_radians * y / yRange) ## compute position along circle
    rads = (circleStart + (circleFraction * y / yRange)) * 2*math.pi ## compute position along circle

    tx = math.sin(rads) * x ## convert to polar x coordinate, adjust radius by rectangular tree x coordinate
    ty = math.cos(rads) * x

    return (tx,ty)

def project_polar_vector(x,y,radians,length):
    """
    Translate a point by a vector specified in polar coordinates.

    Parameters
    ----------
    x, y : float
        Starting point.

    radians : float
        Direction of the vector in radians.

    length : float
        Vector length.

    Returns
    -------
    tuple[float, float]
        Endpoint coordinates.
    """

    new_x = x + length * math.cos(radians)
    new_y = y + length * math.sin(radians)

    return (new_x,new_y)


def desaturate(colour, desat=0.65, out="auto"):
    """
    Desaturate a colour by scaling its HSV saturation component.

    Parameters
    ----------
    colour : str or tuple
        Input colour specification.

    desat : float, optional
        Saturation multiplier in the interval ``[0, 1]``.

    out : {'auto', 'hex', 'rgb', 'rgba'}, optional
        Output format for the desaturated colour.

    Returns
    -------
    str or tuple
        Desaturated colour in the requested format.
    """
    if not (0 <= desat <= 1):
        raise ValueError(f"invalid desat value: {desat}, must be within interval [0, 1].")

    # Detect input kind for round-tripping
    in_is_str = isinstance(colour, str)
    in_is_tuple = isinstance(colour, (tuple, list)) and len(colour) in (3, 4)
    if not (in_is_str or in_is_tuple):
        raise TypeError(f"colour {colour} invalid, must be a name/hex string or RGB/RGBA tuple.")

    # Remember tuple arity and whether hex originally had alpha
    tuple_len = len(colour) if in_is_tuple else None
    hex_had_alpha = False
    if in_is_str and colour.strip().startswith("#"):
        cs = colour.strip()
        hex_had_alpha = len(cs) in (5, 9)  # #RGBA or #RRGGBBAA

    # Convert to RGBA
    rgba = np.array(mpl.colors.to_rgba(colour), dtype=float)

    # Desaturate in HSV (scale S)
    rgb = rgba[:3]
    hsv = mpl.colors.rgb_to_hsv(rgb.reshape(1, 1, 3))
    hsv[..., 1] *= desat
    new_rgb = mpl.colors.hsv_to_rgb(hsv).reshape(3)
    rgba[:3] = new_rgb

    # Decide output format
    if out == "auto":
        if in_is_tuple:
            return tuple(rgba[:tuple_len])  # keep RGB vs RGBA arity
        else:
            # strings round-trip to hex; include alpha if originally present or alpha<1
            keep_alpha = hex_had_alpha or (rgba[3] < 1.0)
            return mpl.colors.to_hex(rgba, keep_alpha=keep_alpha)
    elif out == "hex":
        # include alpha if <1
        return mpl.colors.to_hex(rgba, keep_alpha=(rgba[3] < 1.0))
    elif out == "rgb":
        return tuple(rgba[:3])
    elif out == "rgba":
        return tuple(rgba)
    else:
        raise ValueError(f"out {out} invalid, must be one of {'auto','hex','rgb','rgba'}.")

def make_cmap(colours, position=None, name="custom_cmap"):
    """
    Create a colormap from mixed color formats:
        - RGB float tuples (0–1)
        - RGB int tuples (0–255)
        - Hex strings "#RRGGBB" or "RRGGBB"
        - HTML/CSS names ("red", "steelblue")
        - Matplotlib shorthand ("r", "C0")
    """

    from matplotlib.colors import LinearSegmentedColormap, to_rgb

    normalized = []
    for c in colours:

        # tuple/list: could be float or int RGB
        if isinstance(c, (tuple, list)) and len(c) == 3:
            if all(isinstance(v, float) for v in c) and all(0 <= v <= 1 for v in c):
                normalized.append(tuple(c))
            elif all(isinstance(v, int) for v in c) and all(0 <= v <= 255 for v in c):
                normalized.append(tuple(v / 255.0 for v in c))
            else:
                raise ValueError(f"Invalid RGB tuple: {c}")

        # string: hex, HTML name, or mpl shorthand
        elif isinstance(c, str):
            # allow "ffaa00" without "#"
            if len(c) == 6 and all(ch in "0123456789abcdefABCDEF" for ch in c):
                c = "#" + c
            try:
                normalized.append(to_rgb(c))
            except ValueError:
                raise ValueError(f"Unrecognized color string: {c}")

        else:
            raise TypeError(f"Unsupported color format: {c}")

    if position is None:
        position = np.linspace(0, 1, len(colours))
    else:
        position = np.asarray(position, float)
        if len(position) != len(colours):
            raise ValueError("position must be same length as colors")
        if position[0] != 0 or position[-1] != 1:
            raise ValueError("position must start at 0 and end at 1")

    cdict = {"red": [], "green": [], "blue": []}
    for p, (r, g, b) in zip(position, normalized):
        cdict["red"].append((p, r, r))
        cdict["green"].append((p, g, g))
        cdict["blue"].append((p, b, b))

    return LinearSegmentedColormap(name, cdict)

def desaturate_cmap(cmap, desat=0.65):
    """
    Create a desaturated version of a matplotlib colormap.

    Parameters
    ----------
    cmap : :obj:`matplotlib.colors.Colormap`
        Colormap to desaturate.

    desat : float, optional
        Saturation multiplier applied to sampled colours.

    Returns
    -------
    :obj:`matplotlib.colors.ListedColormap`
        Desaturated colormap.
    """
    from matplotlib.colors import ListedColormap

    assert isinstance(cmap, mpl.colors.LinearSegmentedColormap) or isinstance(cmap, mpl.colors.ListedColormap), f"cmap type {type(cmap)} invalid, must be mpl.colors.LinearSegmentedColormap or mpl.colors.ListedColormap."
    base_colours = [tuple(c) for c in cmap(np.linspace(0, 1, 256))]
    desat_colours = [desaturate(c, desat = desat) for c in base_colours]

    return ListedColormap(desat_colours)

def hpd(data, level=0.95):
    """
    Compute the highest posterior density interval for a sample.

    Parameters
    ----------
    data : sequence[float]
        Posterior samples.

    level : float, optional
        Posterior mass to include in the interval. Defaults to ``0.95``.

    Returns
    -------
    tuple[float, float] or None
        Lower and upper bounds of the highest posterior density interval, or
        ``None`` if there are too few samples to estimate the interval.

    Notes
    -----
    Original implementation copyright (C) 2010 Joseph Heled.
    """
    d = list(data)
    d.sort()

    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2 :
        return None
    #raise RuntimeError("Not enough data. N data: %s"%(len(data)))

    i = 0
    r = d[i+nIn-1] - d[i]
    for k in range(len(d) - (nIn - 1)) :
        rk = d[k+nIn-1] - d[k]
        if rk < r :
            r = rk
            i = k

    assert 0 <= i <= i+nIn-1 < len(d)

    return (d[i], d[i+nIn-1])

def five_point_bezier(points, precision=50):
    """
    Quartic Bézier curve (5 control points).
    Returns arrays of x and y coordinates.
    """
    p0, p1, p2, p3, p4 = map(lambda x: np.array(x, dtype=float), points)

    t = np.linspace(0, 1, precision)

    Bx = ((1-t)**4 * p0[0] +
          4*(1-t)**3 * t * p1[0] +
          6*(1-t)**2 * t**2 * p2[0] +
          4*(1-t) * t**3 * p3[0] +
          t**4 * p4[0])

    By = ((1-t)**4 * p0[1] +
          4*(1-t)**3 * t * p1[1] +
          6*(1-t)**2 * t**2 * p2[1] +
          4*(1-t) * t**3 * p3[1] +
          t**4 * p4[1])

    return Bx, By

def draw_gradient_polygon(
    ax,
    polygonXY,
    extent,
    colour,
    minAlpha=0.0,
    maxAlpha=1.0,
    n=256,
    axis='y',
    origin='lower',
    reverse=False,
    zorder=0,
    interpolation='bicubic',
    addPatch=True,
    patchKwargs=None,
    imshowKwargs=None,
    ):
    """
    Draw an RGBA gradient image and clip it to a polygon.
    """
    import numpy as np
    import matplotlib as mpl
    from matplotlib.patches import Polygon

    assert 0.0 <= minAlpha <= 1.0 and 0.0 <= maxAlpha <= 1.0, f"minAlpha ({minAlpha}) and maxAlpha ({maxAlpha}) must be in range [0.0, 1.0]."
    assert n >= 2, f"Insufficient number of colours: {n}, must be >=2."
    assert axis in ['x', 'y'], f"Unrecognised axis {axis}. Must be 'x' or 'y'."

    polygonXY = np.asarray(polygonXY, dtype=float)
    if polygonXY.ndim != 2 or polygonXY.shape[1] != 2:
        raise ValueError("polygonXY must be an array-like of shape (N, 2)")

    xmin, xmax, ymin, ymax = extent

    # Build RGBA image. We keep the non-ramp dimension at 2 pixels because colour/alpha don't vary there.
    rgb = mpl.colors.colorConverter.to_rgb(colour)

    if axis == 'y':
        img = np.zeros((n, 2, 4), dtype=float)  # rows vary in y
    elif axis == 'x':
        img = np.zeros((2, n, 4), dtype=float)  # cols vary in x

    img[..., 0] = rgb[0]
    img[..., 1] = rgb[1]
    img[..., 2] = rgb[2]

    ramp = np.linspace(minAlpha, maxAlpha, n)
    if reverse: ramp = ramp[::-1]
    
    if axis == 'y':
        img[..., 3] = ramp[:, None]
    elif axis == 'x':
        img[..., 3] = ramp[None, :]
    
    localImshowKwargs = dict(imshowKwargs) if imshowKwargs else {}
    im = ax.imshow(
        img,
        extent=[xmin, xmax, ymin, ymax],
        origin=origin,
        aspect='auto',
        interpolation=interpolation,
        zorder=zorder,
        **localImshowKwargs,
    )

    localPatchKwargs = dict(patchKwargs) if patchKwargs else {}
    # Default: invisible patch used only for clipping
    patch = Polygon(
        polygonXY,
        facecolor='none',
        edgecolor='none',
        closed=localPatchKwargs.pop('closed', True),
        zorder=localPatchKwargs.pop('zorder', zorder),
        **localPatchKwargs,
    )
    if addPatch:
        ax.add_patch(patch)

    im.set_clip_path(patch)
    return im, patch

def get_path_effects(mainColour='k', outlineColour='w', mainWeight=0.5, outlineWeight=4):
    """
    Construct a simple stroked text/line path effect stack.

    Parameters
    ----------
    mainColour : color, optional
        Foreground colour for the inner stroke.

    outlineColour : color, optional
        Colour for the outer stroke.

    mainWeight : float, optional
        Line width of the inner stroke.

    outlineWeight : float, optional
        Line width of the outer stroke.

    Returns
    -------
    list
        Matplotlib path-effect objects.
    """

    import matplotlib.patheffects as path_effects

    effects = [path_effects.Stroke(linewidth=outlineWeight, foreground=outlineColour), 
               path_effects.Stroke(linewidth=mainWeight, foreground=mainColour)]

    return effects
