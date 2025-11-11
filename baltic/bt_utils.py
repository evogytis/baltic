import re
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

def calendar_to_decimal_date(date,fmt="%Y-%m-%d",variable=False):
    if not fmt:
        return date
    delimiter=re.search('[^0-9A-Za-z%]',fmt) ## search for non-alphanumeric symbols in fmt (should be field delimiter)
    delimit=None
    if delimiter is not None:
        delimit=delimiter.group()

    if variable: ## if date is variable - extract what is available
        if delimit is not None:
            dateL=len(date.split(delimit)) ## split date based on symbol
        else:
            dateL=1 ## no non-alphanumeric characters in date, assume dealing with an imprecise date (something like just year)

        if dateL==2:
            fmt=delimit.join(fmt.split(delimit)[:-1]) ## reduce fmt down to what's available
        elif dateL==1:
            fmt=delimit.join(fmt.split(delimit)[:-2])

    adatetime=dt.datetime.strptime(date,fmt) ## convert to datetime object
    year = adatetime.year ## get year
    boy = dt.datetime(year, 1, 1) ## get beginning of the year
    eoy = dt.datetime(year + 1, 1, 1) ## get beginning of next year
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds())) ## return fractional year

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
        outputFmtFxn = lambda date: bt_utils.convert_date_format(date, '%Y-%m-%d', '%b\n%Y') if bt_utils.convert_date_format(date, '%Y-%m-%d', '%m') == '01' else bt_utils.convert_date_format(date, '%Y-%m-%d', '%b')

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
    slope,intercept,rval,_,_ = linregress(tipDates,tipHeights) ## run linear regression
    corr = np.corrcoef((tipDates,tipHeights))[0,1] ## correlation coefficient
    ssq = sum([(y-(slope*x+intercept))**2 for x,y in zip(tipDates,tipHeights)]) ## sum of squares

    if stat=='correlation': ## set stat to optimise
        localStat = corr
    elif stat=='sum of squares':
        localStat = ssq
    elif stat=='r^2':
        localStat = rval
    else:
        raise ValueError(f'Unknown stat {stat} to optimise')

    optStat = res[stat] if stat in res else (np.inf if stat in ['sum of squares'] else -np.inf)

    if (localStat < optStat if stat in ['sum of squares'] else localStat > optStat): ## minimise sum of squares or maximise correlation/r^2
        if forcePositive and slope<0: ## force positive True and slope<0
            pass ## if forcing positive root slope must be positive, do nothing
        else: ## force_positive is False or true and slope>0
            optStat = localStat ## better root found
#             new_root = [w for w in self.Objects if w.index==k.index][0]
            res['correlation']=corr
            res['sum of squares']=ssq
            res['slope']=slope
            res['intercept']=intercept
            res['r^2']=rval
            res['root']=rootCandidate
            if frac is None: res['frac']=frac

    return res

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

def desaturate_cmap(cmap, desat = 0.65):
    from matplotlib.colors import ListedColormap

    assert isinstance(cmap, mpl.colors.LinearSegmentedColormap) or isinstance(cmap, mpl.colors.ListedColormap), f"cmap type {type(cmap)} invalid, must be mpl.colors.LinearSegmentedColormap or mpl.colors.ListedColormap."
    base_colours = [tuple(c) for c in cmap(np.linspace(0, 1, 256))]
    desat_colours = [desaturate(c, desat = desat) for c in base_colours]

    return ListedColormap(desat_colours)