import re
import logging
import datetime as dt
import math
import numpy as np
from itertools import permutations
from scipy.stats import linregress

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

    circle_start_radians = 2*math.pi * circleStart ## convert starting point to radians
    circle_fraction_radians = 2*math.pi * circleFraction ## convert arc width to radians
    
    rads = circle_start_radians + (circle_fraction_radians * y / yRange) ## compute position along circle

    tx = math.sin(rads) * x ## convert to polar x coordinate, adjust radius by rectangular tree x coordinate
    ty = math.cos(rads) * x

    return (tx,ty)

def project_polar_vector(x,y,radians,length):

    new_x = x + length * math.cos(radians)
    new_y = y + length * math.sin(radians)

    return (new_x,new_y)
