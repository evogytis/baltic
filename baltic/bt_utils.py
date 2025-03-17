import re
import datetime as dt

def decimalDate(date,fmt="%Y-%m-%d",variable=False):
    """
    Converts calendar dates in specified format to decimal date. 
    
    A decimal date represents the fraction of the year that has passed by the given date.
    
    Parameters:
    date (str): The date to be converted.
    fmt (str): The format of the input date string. Default is "%Y-%m-%d".
    variable (bool): If True, allows for variable date precision. Default is False.
    
    Returns:
    float: The decimal representation of the date.
    
    Notes:
    - If `variable` is True, the function adjusts the format to match the available precision in the date.
    - For example, a date like "2023" will be interpreted as 2023 Jan 01, while "2023-05" will be interpreted as 2023 May 01.
    
    Examples:
    >>> decimalDate("2023-05-23")
    2023.3890410958904
    >>> decimalDate("2023", fmt="%Y", variable=True)
    2023.0

    Docstring generated with ChatGPT 4o.
    """
    if fmt == "":
        return date
    delimiter=re.search('[^0-9A-Za-z%]',fmt) ## search for non-alphanumeric symbols in fmt (should be field delimiter)
    delimit=None
    if delimiter is not None:
        delimit=delimiter.group()

    if variable==True: ## if date is variable - extract what is available
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

def calendarDate(timepoint,fmt='%Y-%m-%d'):
    """
    Converts decimal dates to a specified calendar date format.
    
    A decimal date represents the fraction of the year that has passed by the given timepoint.
    This function converts it back to a calendar date in the given format.
    
    Parameters:
    timepoint (float): The decimal representation of the date.
    fmt (str): The desired format of the output date string. Default is '%Y-%m-%d'.
    
    Returns:
    str: The date in the specified calendar format.
    
    Examples:
    >>> calendarDate(2023.3923497267758)
    '2023-05-24'
    >>> calendarDate(2023.0, fmt='%Y')
    '2023'
    
    Docstring generated with ChatGPT 4o.
    """
    year = int(timepoint)
    rem = timepoint - year

    base = dt.datetime(year, 1, 1)
    result = base + dt.timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    return dt.datetime.strftime(result,fmt)

def convertDate(date_string,start,end):
    """
    Converts calendar dates between given formats.
    
    Parameters:
    x (str): The date string to be converted.
    start (str): The format of the input date string.
    end (str): The desired format of the output date string.
    
    Returns:
    str: The date converted to the new format.
    
    Examples:
    >>> convertDate('23-05-2023', '%d-%m-%Y', '%Y/%m/%d')
    '2023/05/23'
    >>> convertDate('2023/05/23', '%Y/%m/%d', '%B %d, %Y')
    'May 23, 2023'
    
    Docstring generated with ChatGPT 4o.
    """
    return dt.datetime.strftime(dt.datetime.strptime(date_string,start),end)
    try:
        date_obj = dt.datetime.strptime(date_string, start)
        return dt.datetime.strftime(date_obj, end)
    except ValueError as e:
        raise ValueError('Error converting date "%s" from format "%s" to "%s": "%s"'%(date_string, start, end, e))


def untangle(trees,cost_function=None,iterations=None,verbose=False):
    """
    Minimise y-axis discrepancies between tips of trees in a list.
    Only the tangling of adjacent trees in the list is minimised, so the order of trees matters.
    Trees do not need to have the same number of tips but tip names should match.
    
    Parameters:
    trees (list): A list of tree objects to untangle.
    cost_function (function or None): A function to calculate the cost of y-axis discrepancies between tips.
                                      Default is None, which uses the squared difference between y axis positions.
    iterations (int or None): The number of iterations to perform. Default is None, which sets the iterations to 3.
    verbose (bool): If True, prints verbose output during the process. Default is False.
    
    Returns:
    list: The list of untangled tree objects.
    
    Example:
    >>> untangled_trees = untangle(list_of_trees, iterations=5, verbose=True)
    
    Docstring generated with ChatGPT 4o.
    """
    from itertools import permutations

    if iterations==None: iterations=3
    if cost_function==None: cost_function=lambda pair: math.pow(abs(pair[0]-pair[1]),2)

    y_positions={T: {k.name: k.y for k in T.getExternal()} for T in trees} ## get y positions of all the tips in every tree

    for iteration in range(iterations):
        if verbose: print('Untangling iteration %d'%(iteration+1))
        first_trees=list(range(len(trees)-1))+[-1] ## trees up to next-to-last + last
        next_trees=list(range(1,len(trees)))+[0] ## trees from second + first
        for cur,nex in zip(first_trees,next_trees): ## adjacent pairs
            tree1=trees[cur] ## fetch current tree
            tree2=trees[nex] ## fetch next tree
            if verbose: print('%d vs %d'%(cur,nex))
            for k in sorted(tree2.getInternal(),key=lambda branch: branch.height): ## iterate through nodes of next tree by height (start from root)
                clade_y_positions=sorted([y_positions[tree2][tip] for tip in k.leaves]) ## sorted list of available y coordinates for node
                costs={} ## will store cost of all children permutations
                if len(k.children)>=10: raise RuntimeWarning('Node is too polytomic and untangling will take an astronomically long time')
                if verbose==True: print(len(k.children))
                for permutation in permutations(k.children): ## iterate over permutations of node's children
                    clade_order=sum([[child.name] if child.is_leaf() else list(child.leaves) for child in permutation],[]) ## flat list of tip names as they would appear in permutation order
                    new_y_positions={clade_order[i]: clade_y_positions[i] for i in range(len(clade_y_positions))} ## assign available y positions in order

                    tip_costs=list(map(cost_function,[(y_positions[tree1][tip],new_y_positions[tip]) for tip in clade_order if tip in y_positions[tree1]]))
                    costs[permutation]=sum(tip_costs)/len(tip_costs) ## compute cost of this permutation in relation to next tree

                best=sorted(costs.keys(),key=lambda w: -costs[w])[0] ## get tree with smallest cost
                k.children=list(best) ## reorder children according to minimised cost

            tree2.drawTree() ## compute new y coordinates for nodes
            for k in tree2.getExternal(): ## iterate over tips
                y_positions[tree2][k.name]=k.y ## remember new coordinates

    return trees