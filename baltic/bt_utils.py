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
    Convert calendar date into decimal year.
    If variable=True, return (midpoint, (min, max)) where min/max reflect
    the full possible range implied by the precision of the input date.
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
    year = int(timepoint)
    rem = timepoint - year

    base = dt.datetime(year, 1, 1)
    result = base + dt.timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    return dt.datetime.strftime(result,fmt)

def convert_date_format(dateString,startFormat,endFormat):
    return dt.datetime.strftime(dt.datetime.strptime(dateString,startFormat),endFormat)
    try: #TODO deal with stuff that comes after the return statement
        date_obj = dt.datetime.strptime(dateString, startFormat)
        return dt.datetime.strftime(date_obj, endFormat)
    except ValueError as e:
        raise ValueError('Error converting date "%s" from format "%s" to "%s": "%s"'%(dateString, startFormat, endFormat, e))

def state_collapse_tree(tree, switchFxn):
    """
    Return a deepcopied and reduced version of the tree provided where subtrees are labelled identically when branches evaluate switchFxn to False.
    Also known by the name of Phylotype maps/trees.
    """
    import copy
    
    local_tree = copy.deepcopy(tree)
    
    _partition_tree(local_tree, partitionFxn = switchFxn) ## run tree labelling
    
    labels = set(local_tree.get_parameter_list('partition',useTraitsDict=True)) ## get all unique labels in tree
    sortingFxn = lambda k: k.absoluteTime if local_tree.treeType == 'time' else k.height ## determine whether to sort by height or absoluteTime
    
    label_to_branches = {label: [k for k in local_tree.get_external() if k.traits['partition'] == label] for label in labels} ## map each unique label to the furthest tip with that label
    label_counts = {label: len(label_to_branches[label]) for label in labels}

    label_to_last_branch = {label: sorted(label_to_branches[label], key = sortingFxn)[-1] for label in labels if label_counts[label] > 0}
    
    local_tree = local_tree.reduce_tree(label_to_last_branch.values()) ## keep a single tip for each unique label

    for k in local_tree.get_external():
        k.traits['size'] = label_counts[k.traits['partition']] ## assign tip counts of this label to representative branch
        k.traits['members'] = label_to_branches[k.traits['partition']] ## remember descendants with label that may no longer be present
    
    return local_tree

def generate_calendar_timeline(startDateStr,endDateStr,spacing='monthly',dateFmt='%Y-%m-%d',roundDates=True):
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

def plot_tangled_chain(ax, treeList, colourMap=None, padding=None, treeSpaceFxn=None, treeSpace=None, treeKwargs={}, pointKwargs={}, **kwargs):
    localKwargs = dict(kwargs)
    localTreeKwargs = dict(treeKwargs)
    localPointKwargs = dict(pointKwargs)

    if treeSpace is not None and treeSpaceFxn is not None: ## treeSpace is how much space is left between consecutive trees, takes current tree if specified as a function
        raise ValueError(
            "Cannot specify both treeSpace and treeSpaceFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if treeSpace is None and treeSpaceFxn is None:
        treeSpaceFxn = lambda k: treeList[0].treeHeight * 0.20 ## 20% of first tree height is space between all trees
    elif treeSpaceFxn is None:
        treeSpaceFxn = lambda k: treeSpace

    if padding is None: ## padding is proportion of treeSpace protrudes beyond previous tree and before next tree (0 == line finishes at last tip of current tree and goes to root of next, 0.5 == line goes to middle between consecutive trees and switches abruptly)
        padding = 0.1
    else:
        assert 0.0 <= padding <= 0.5, f"Padding (given as {padding}) should be a float between 0 and 0.5."

    if colourMap is None: ## colourMap is dict that assigns colours to tips according to their y-axis order in first tree
        colourMap = {}

        cmap = mpl.cm.Spectral
        firstTreeTips = treeList[0].get_external()

        for i,k in enumerate(sorted(firstTreeTips, key = lambda q: q.y)):
            colourMap[k.name] = cmap(i/(len(firstTreeTips)-1))

    cumulativeX = 0 ## tracks x coordinate as we plot consecutive trees
    if 'coordinateFxn' in localTreeKwargs: warnings.warn(f"Custom x coordinate function for tree was specified but will be overriden for tangled chain visualisation.")
    if 'xCoordinateFxn' not in localTreeKwargs: localTreeKwargs['xCoordinateFxn'] = lambda k: k.x + cumulativeX

    if len(localPointKwargs)>0: ## tip points are required - override xCoordinateFxn, assign default colours if nothing specified
        if 'xCoordinateFxn' in localPointKwargs: warnings.warn(f"Custom x coordinate function for points was specified but will be overriden for tangled chain visualisation.")
        localPointKwargs['xCoordinateFxn'] = lambda k: k.x + cumulativeX
        if 'colour' not in localPointKwargs and 'colourFxn' not in localPointKwargs:
            warnings.warn(f"Point colours were not specified, defaulting to tangled chain colour defaults. This may cause issues if targetFxn is not set to identify tips.")
            localPointKwargs['colourFxn'] = lambda k: colourMap[k.name]

    connectionCoordinates = []
    connectionColours = []

    for curTree, nexTree in zip(treeList,treeList[1:]): ## iterate over pairs of consecutive trees
        curTree.plot_tree(ax,**localTreeKwargs) ## plot current tree
        if len(localPointKwargs)>0:
            curTree.plot_points(ax,**localPointKwargs) ## add points if specified

        spaceUnit = treeSpaceFxn(curTree)

        for curTip in curTree.get_external(): ## iterate over tips in current tree
            c = colourMap[curTip.name] if curTip.name in colourMap else 'lightgray'

            nexTip = nexTree.get_external(filterFxn = lambda k: k.name == curTip.name) ## identify matching tip
            if len(nexTip)>0:
                nexTip = nexTip[0]
                curX = localTreeKwargs['xCoordinateFxn'](curTip)
                curY = curTip.y

                lineAfterX = cumulativeX + curTree.treeHeight + spaceUnit*padding
                lineBeforeX = cumulativeX + curTree.treeHeight + spaceUnit*(1-padding)

                nexX = localTreeKwargs['xCoordinateFxn'](nexTip) + curTree.treeHeight + spaceUnit
                nexY = nexTip.y

                connectionCoordinates.append([(curX, curY),
                                              (lineAfterX, curY),
                                              (lineBeforeX, nexY),
                                              (nexX, nexY)]) ## coordinates of tangled line
                connectionColours.append(c) ## colour of tangled line

        cumulativeX += curTree.treeHeight + spaceUnit ## increment x-axis

    nexTree.plot_tree(ax,**localTreeKwargs) ## plot last tree
    if len(localPointKwargs)>0:
        nexTree.plot_points(ax,**localPointKwargs) ## plot its points

    if 'zorder' not in localKwargs: localKwargs['zorder'] = 0
    ax.add_collection(LineCollection(connectionCoordinates,color=connectionColours,**localKwargs)) ## add tangled lines

    return ax


def plot_scale_bar(ax, xy, L = None, tree = None, alnL = None, textXY = None, unitText = None, style = 'simple', orientation = 'horizontal', ySpan = None, lineKwargs = None, textKwargs = None):
    """
    Plots a scale bar at given coordinates, can plot a scale bar of required length L or inferred automatically from the tree or just defaults to 0.001.
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
        xs, ys = ys, xs
    
    ax.plot(xs, ys, **localLineKwargs)
    
    if style == 'fancy':
        width = 0.005 * ySpan
        left_xs, left_ys = [x, x], [y - width, y + width]
        right_xs, right_ys = [x + L, x + L], [y - width, y + width]

        if orientation == 'vertical':
            left_xs, left_ys = left_ys, left_xs
            right_xs, right_ys = right_ys, right_xs
        
        ax.plot(left_xs, left_ys, **localLineKwargs)
        ax.plot(right_xs, right_ys, **localLineKwargs)

    if textXY is None: ## set text position defaults / might be worth turning these coordinates as fractions in relation to the bar itself
        textX, textY = (x + L/2, y - ySpan * 0.02)
    else:
        textX, textY = textXY
    
    if orientation == 'vertical':
        textX, textY = textY, textX
    
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
    assert f"{traitName}.set" and f"{traitName}.set.prob" in node.traits, f"{traitName}.set or {traitName}.set.prob not found in node traits dict."

    stateSet = node.traits[f"{traitName}.set"]
    stateProbs = node.traits[f"{traitName}.set.prob"]
    
    stateSetDict = {key: value for key, value in zip(stateSet, stateProbs)}

    return stateSetDict

def plot_node_treemap(ax, node, traitName, traitColourDict, height, width, centerFxn = None, area = 1.0, other_thres = 0.0, **kwargs):
    
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
    
    probs = [stateSetDict[state] for state in stateOrder]
    cs = [traitColourDict[state] for state in stateOrder]
    
    ax.pie(x = probs, colors = cs, radius = radius, center = centerFxn(node), **kwargs)

    ax.autoscale()
    return ax

def plot_tmrca_posterior(ax, tmrcaFile, tmrcaName = 'age(root)', burnin = None, yCoord = None, fullViolin = True, 
                         hpdLvl = 0.95, precision = 100, kdeWidth = 3, orientation = 'horizontal', connectNode = False, node = None, violinKwargs = {}, outlineKwargs = {}):
    
    ### certain tree orientations will require xCoord too
    ### connect mean of HPD to node with dotted line (accommodate elbow in case KDE is on the x-axis)
    import csv
    from scipy.stats import gaussian_kde
    from baltic.bt_utils import hpd, decimal_to_calendar_date


    func_logger = logging.getLogger("baltic.tree.plot_tmrca_posterior")
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
        warnings.warn('No burnin set, defaulting to 10M states.')
    
    handle = open(tmrcaFile,'r')

    for l in csv.DictReader((line for line in handle if line.startswith('#') == False), delimiter = '\t'):
        state = int(l['state'])
        if state >= burnin:
            tmrcaPosterior.append(float(l[tmrcaName])) ## grab column with tmrca stat

    handle.close()
    
    if yCoord is None and node is None:
        yCoord = 0.0
    elif yCoord is None and node:
        yCoord = node.y
    elif yCoord is not None and node is not None:
        warnings.warn(f"Both yCoord and node were provided, KDE will be positioned at yCoord value.")
    
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
    if 'linewidth' not in localOutlineKwargs: localOutlineKwargs['linewidth'] = 2
    if 'zorder' not in localOutlineKwargs: localOutlineKwargs['zorder'] = 2
    
    kde = gaussian_kde(tmrcaPosterior)
    hpdLo, hpdHi = hpd(tmrcaPosterior, hpdLvl)
    x_grid = np.linspace(hpdLo, hpdHi, precision)

    meanTmrca = np.mean(tmrcaPosterior)
    medianTmrca = np.median(tmrcaPosterior)

    func_logger.info(f"Node {tmrcaName} mean TMRCA: {meanTmrca:.3f} ({decimal_to_calendar_date(meanTmrca)}) median TMRCA: {medianTmrca:.3f} ({decimal_to_calendar_date(medianTmrca)})")
    func_logger.info(f"Node {tmrcaName} TMRCA 95% HPD: {hpdLo:.3f} - {hpdHi:.3f} / {decimal_to_calendar_date(hpdLo)} - {decimal_to_calendar_date(hpdHi)}")
    
    y_grid = kde(x_grid)
    y_max = y_grid.max()
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
            mean_ys = [yCoord - kde(meanTmrca) / y_max * kdeWidth, yCoord + kde(meanTmrca) / y_max * kdeWidth]
        else:
            mean_ys = [yCoord, yCoord + kde(meanTmrca).item() / y_max * kdeWidth]
            mean_ys = np.array(mean_ys, dtype=float).tolist()
        
        elbow_xs = [meanTmrca, meanTmrca, node.absoluteTime]
        elbow_ys = [yCoord, node.y, node.y]
        
        if orientation == 'vertical':
            x, y = y, x
            mean_xs, mean_ys = mean_ys, mean_xs
            elbow_xs, elbow_ys = elbow_ys, elbow_xs
        
        ax.scatter(x, y, s = 40, fc = localViolinKwargs['facecolor'], ec = 'none', zorder = 10)
        ax.scatter(x, y, s = 80, fc = localViolinKwargs['edgecolor'], ec = 'none', zorder = 9)

        ax.plot(mean_xs, mean_ys, color = 'dimgray', ls = '-', zorder = 2)
        ax.plot(elbow_xs, elbow_ys, color = 'dimgray', ls = '--', zorder = 8)
    elif node is None:
        warnings.warn(f"Need to specify node to connect to.")
    
    return ax

def plot_time_grid(ax, timeline, dateFmt='%Y-%m-%d', colourFxn=None, colour=None, edgeColourFxn=None, edgeColour=None, axis='x',**kwargs):

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

    validRemoveTickLabels = ['x', 'y', 'both']
    assert removeTickLabels in validRemoveTickLabels, f"removeTickLabels value {removeTickLabels} not recognised. Must be one of {', '.join(validRemoveTickLabels)}"

    if removeTickLabels in ['x', 'both']:
        ax.set_xticks([])
        ax.set_xticklabels([])
    if removeTickLabels in ['y', 'both']:
        ax.set_yticks([])
        ax.set_yticklabels([])
    
    [ax.spines[loc].set_visible(False) for loc in hideSpines]
    return ax

def untangle(trees,costFxn=None,iterations=None):
    if iterations is None: iterations=3
    if costFxn is None: costFxn=lambda pair: math.pow(abs(pair[0]-pair[1]),2)

    y_positions={T: {k.name: k.y for k in T.getExternal()} for T in trees} ## get y positions of all the tips in every tree

    for iteration in range(iterations):
        logger.debug(f'Untangling iteration {iteration+1}')
        first_trees=list(range(len(trees)-1))+[-1] ## trees up to next-to-last + last
        next_trees=list(range(1,len(trees)))+[0] ## trees from second + first
        for cur,nex in zip(first_trees,next_trees): ## adjacent pairs
            tree1=trees[cur] ## fetch current tree
            tree2=trees[nex] ## fetch next tree
            logger.debug(f'{cur} vs {nex}')
            for k in sorted(tree2.getInternal(),key=lambda branch: branch.height): ## iterate through nodes of next tree by height (start from root)
                clade_y_positions=sorted([y_positions[tree2][tip] for tip in k.leaves]) ## sorted list of available y coordinates for node
                costs={} ## will store cost of all children permutations
                if len(k.children)>=10: raise RuntimeWarning('Node is too polytomic and untangling will take an astronomically long time')
                logger.debug(len(k.children))
                for permutation in permutations(k.children): ## iterate over permutations of node's children
                    clade_order=sum([[child.name] if child.is_leaf() else list(child.leaves) for child in permutation],[]) ## flat list of tip names as they would appear in permutation order
                    new_y_positions={clade_order[i]: clade_y_positions[i] for i in range(len(clade_y_positions))} ## assign available y positions in order

                    tip_costs=list(map(costFxn,[(y_positions[tree1][tip],new_y_positions[tip]) for tip in clade_order if tip in y_positions[tree1]]))
                    costs[permutation]=sum(tip_costs)/len(tip_costs) ## compute cost of this permutation in relation to next tree

                best=sorted(costs.keys(),key=lambda w: -costs[w])[0] ## get tree with smallest cost
                k.children=list(best) ## reorder children according to minimised cost

            tree2.drawTree() ## compute new y coordinates for nodes
            for k in tree2.getExternal(): ## iterate over tips
                y_positions[tree2][k.name]=k.y ## remember new coordinates

    return trees

def unnest(nodeList, towardsRoot = True):
    assert all([(k.is_node() or k.is_leaflike()) for k in nodeList]), f'nodeList contains objects that are not baltic branch objects (node or leaflike): {', '.join([k for k in nodeList if k.is_node() == False and k.is_leaflike() == False])}'
    
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
    slope, intercept, rval, _, _ = linregress(tipDates,tipHeights) ## run linear regression
    corr = np.corrcoef((tipDates, tipHeights))[0,1] ## correlation coefficient
    ssq = sum([(y - (slope * x + intercept))**2 for x, y in zip(tipDates, tipHeights)]) ## sum of squares

    if stat=='correlation': ## set stat to optimise
        localStat = corr
    elif stat=='sum of squares':
        localStat = ssq
    elif stat=='r^2':
        localStat = rval
    else:
        raise ValueError(f'Unknown stat {stat} to optimise')

    optStat = res[stat] if stat in res else (np.inf if stat in ['sum of squares'] else -np.inf)

    invalidateRegression = True if forcePositive and slope < 0 else False

    if (localStat < optStat if stat in ['sum of squares'] else localStat > optStat): ## minimise sum of squares or maximise correlation/r^2
        # if forcePositive and slope < 0: ## force positive True and slope<0
        #     pass ## if forcing positive root, slope must be positive; do nothing
        #     res['correlation'] = -np.inf
        #     res['sum of squares'] = np.inf
        #     res['slope'] = slope
        #     res['intercept'] = intercept
        #     res['r^2'] = -np.inf
        #     res['root'] = rootCandidate
        #     if frac is not None: res['frac'] = frac
        # else: ## force_positive is False or true and slope>0
        optStat = localStat ## better root found
#             new_root = [w for w in self.Objects if w.index==k.index][0]
        res['correlation'] = -np.inf if invalidateRegression else corr
        res['sum of squares'] = np.inf if invalidateRegression else ssq
        res['slope'] = slope
        res['intercept'] = intercept
        res['r^2'] = -np.inf if invalidateRegression else rval
        res['root'] = rootCandidate
        if frac is not None: res['frac'] = frac

    return res

def _rtt_worker(args):
    """
    Worker for root_by_regression.

    args = (
        tree,                # a picklable copy of the tree
        root_index,          # candidate root node index
        fixed_times,         # dict: tip.index -> fixed date
        uncertain_ranges,    # dict: tip.index -> (min_date, max_date)
        n_mc,                # number of Monte Carlo iterations
        stat,                # 'r^2', 'correlation', or 'sum of squares'
        forcePositive,       # as in your original function
    )
    """
    (tree,
     root_index,
     tipDates,
     uncertainDateRanges,
     n_mc,
     stat,
     forcePositive) = args

    tree_copy = copy.deepcopy(tree)

    # find candidate root node inside worker copy
    candidate = next(obj for obj in tree_copy.Objects if obj.index == root_index)

    # Reroot on this candidate for this sample
    cll = tree_copy.reroot(candidate)
    tips = cll.get_external()

    rootDistances = {k.name: k.height for k in tips} ## these don't change

    best_res = None
    best_score = None
    best_dates = None
    best_branch_frac = None

    rng = np.random.default_rng()

    for _ in range(n_mc): ## Monte Carlo iterations

        for k in uncertainDateRanges:
            tipDates[k] = rng.uniform(*uncertainDateRanges[k]) ## for uncertain tips sample from possible date range randomly

        xs = [tipDates[k.name] for k in tips]
        ys = [rootDistances[k.name] for k in tips]

        # Run root-to-tip regression
        res_local = {}
        res_local = _root_to_tip(
            candidate, xs, ys, res_local,
            stat=stat, forcePositive=forcePositive
        )

        if res_local: ## valid regression (when forcePositive == True but slope is negative res_local is None)
            # Compute a scalar score, taking into account direction of optimisation
            if stat == "sum of squares":
                # Smaller is better
                score = -res_local["sum of squares"]
            else:
                # Larger is better for r^2 or correlation
                score = res_local[stat]

            if best_score is None or score > best_score:
                best_res = res_local
                best_score = score
                best_dates = {k: tipDates[k] for k in uncertainDateRanges}

    if best_res is None:
        best_res = {}

    best_res["root_index"] = root_index
    best_res["score"] = best_score
    best_res["monte_carlo_dates"] = best_dates

    # IMPORTANT: we replace 'root' with index so we don't leak deep-copied nodes
    best_res["root"] = root_index

    return best_res


def project_to_polar(x,y,yRange,circleStart=0.0,circleFraction=1.0):

    # circle_start_radians = 2*math.pi * circleStart ## convert starting point to radians
    # circle_fraction_radians = 2*math.pi * circleFraction ## convert arc width to radians
    
    # rads = circle_start_radians + (circle_fraction_radians * y / yRange) ## compute position along circle
    rads = (circleStart + (circleFraction * y / yRange)) * 2*math.pi ## compute position along circle

    tx = math.sin(rads) * x ## convert to polar x coordinate, adjust radius by rectangular tree x coordinate
    ty = math.cos(rads) * x

    return (tx,ty)

def project_polar_vector(x,y,radians,length):

    new_x = x + length * math.cos(radians)
    new_y = y + length * math.sin(radians)

    return (new_x,new_y)


def desaturate(colour, desat = 0.65, out = "auto"):
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

def desaturate_cmap(cmap, desat = 0.65):
    from matplotlib.colors import ListedColormap

    assert isinstance(cmap, mpl.colors.LinearSegmentedColormap) or isinstance(cmap, mpl.colors.ListedColormap), f"cmap type {type(cmap)} invalid, must be mpl.colors.LinearSegmentedColormap or mpl.colors.ListedColormap."
    base_colours = [tuple(c) for c in cmap(np.linspace(0, 1, 256))]
    desat_colours = [desaturate(c, desat = desat) for c in base_colours]

    return ListedColormap(desat_colours)

def hpd(data, level = 0.95):
    """
    Return highest posterior density interval from a list,
    given the posterior density interval required.
    Copyright (C) 2010 Joseph Heled
    Author: Joseph Heled <jheled@gmail.com>
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