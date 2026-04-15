import re
import logging
import numpy as np
from baltic.baltic import make_tree
from baltic.bt_utils import calendar_to_decimal_date

logger = logging.getLogger("baltic.samogitia")

def posterior_tree_iterator(treesPath, burnin, mostRecentDate, tipRegex, dateFmt, treestringRegex):
    """
    Iterate over tree strings in a posterior ``.trees`` file after burn-in.

    Parameters
    ----------
    treesPath : str
        Path to the posterior tree file.

    burnin : int
        Number of initial sampled trees to skip.

    outputPath : str
        Output path associated with the downstream processing workflow.

    mostRecentDate : float or None
        Explicit most recent sampling date. Currently unused during parsing.

    tipRegex : str
        Regular expression used to extract dates from tip names.

    dateFmt : str
        Calendar date format used to parse matched tip dates.

    treestringRegex : str
        Regular expression used to identify tree lines and extract state
        numbers.

    Yields
    ------
    tuple
        Tuple ``(i, state, treeString, tips, maxDate)`` for each sampled tree.
    """
    tip_flag = False
    tips = {}
    tip_num = 0
    treeCounter = 0
    ############
    handle = open(treesPath, 'r')
    
    for line in handle:
        l = line.strip('\n')
        
        ############
        match = re.search('Dimensions ntax=([0-9]+);',l)
        if match:
            tip_num = int(match.group(1))
            logger.debug(f'File should contain {tip_num} taxa')

        ###############
        match = re.search(treestringRegex,l)
        if match:
            state = int(match.group(1))
            
            treeString_start = l.index('(')

            if treeCounter == 0: ## only need to get tip names on the first tree
                ####
                # do an assertion here that tipNames exist in tree
                ####
                tree = make_tree(l[treeString_start:], 'time')

                if len(tips) > 0:
                    tree.rename_tips(tips)
                
                tip_dates = []
                tip_names = []
                for k in tree.get_external():
                    tip_names.append(k.name)
                    match = re.search(tipRegex, k.name)
                    if match:
                        tip_dates.append(calendar_to_decimal_date(match.group(1), fmt = dateFmt, variable = False))

                if mostRecentDate:
                    maxDate = mostRecentDate
                    logger.warning(f"Using provided mostRecentDate ({mostRecentDate})")
                else:
                    assert len(tip_dates) > 0, f"Could not identify any tip dates from tip names."
                    maxDate = max(tip_dates)
                
            if state < burnin:
                continue

            yield treeCounter, state, l[treeString_start:], tips, maxDate
            treeCounter += 1
            
            logger.debug('Identified tree string')
        ##############
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
    
    handle.close()

def process_posterior_trees(treesPath, processFxn, workers = 4, burnin = None, outputPath = 'processed-output.log.txt', mostRecentDate = None, tipRegex=r'\|([0-9]+\-[0-9]+\-[0-9]+)', dateFmt='%Y-%m-%d', treestringRegex=r'tree [A-Za-z\_]+([0-9]+)', **kwargs):
    """
    Process a posterior tree set with a parallel worker function.

    Parameters
    ----------
    treesPath : str
        Path to the posterior ``.trees`` file.

    processFxn : callable
        Worker function applied to each sampled tree.

    workers : int, optional
        Maximum number of worker processes to use. Defaults to ``4``.

    burnin : int, optional
        Number of initial sampled trees to skip. Defaults to ``0``.

    outputPath : str, optional
        Path to the tab-delimited output file.

    mostRecentDate : float, optional
        Explicit most recent sampling date to pass to workers.

    tipRegex : str, optional
        Regular expression used to extract dates from tip names.

    dateFmt : str, optional
        Date format used to parse extracted tip dates.

    treestringRegex : str, optional
        Regular expression used to identify tree lines and their state values.

    **kwargs
        Additional keyword arguments forwarded to *processFxn*.
    """
    if burnin is None:
        burnin = 0
    
    import heapq
    from concurrent.futures import ProcessPoolExecutor, as_completed

    pq = []
    next_to_write = 0   # index of next tree result to write

    with ProcessPoolExecutor(max_workers = workers) as ex, open(outputPath, "w") as out:
        futures = []

        # submit each tree as we read it
        for i, state, treeString, tipRenameDict, maxDate in posterior_tree_iterator(treesPath, burnin = burnin, mostRecentDate = mostRecentDate, tipRegex = tipRegex, dateFmt = dateFmt, treestringRegex = treestringRegex): ## posterior tree iterator returns the index of the tree, its MCMC state number, treestring itself, tip renaming dict and most recent tip date
            
            if i == 0: ## at first tree
                _, _, header_val = processFxn(-1, state, treeString, tipRenameDict, maxDate, **kwargs, headerMode = True) ## call worker in header mode
                out.write(f"state\t{'\t'.join(header_val)}\n") ## write header to output file
                
            fut = ex.submit(processFxn, i, state, treeString, tipRenameDict, maxDate, **kwargs, headerMode = False) ## process tree as usual - calling worker in headerMode = False
            futures.append(fut)

        # process completed tasks
        for fut in as_completed(futures):
            try:
                i, state, value = fut.result()
            except Exception as e:
                print("Worker crashed:", e)
                raise
            
            # i, state, value = fut.result()
            heapq.heappush(pq, (i, state, value))

            # write anything we can write now
            while pq and pq[0][0] == next_to_write:
                _, st, val = heapq.heappop(pq)
                # print(state, '\t'.join(val))
                out.write(f"{st}\t{'\t'.join(map(str, val))}\n") ## output to file
                next_to_write += 1


def label_timepoints(times, edges, labels):
    """
    Assign string labels to times based on interval start times.

    Parameters
    ----------
    times : array-like
        Time points to label (e.g. from np.linspace).
    edges : array-like
        Start times of intervals (unordered OK).
        Label i applies on [edges[i], edges[i+1]),
        up to the last edge.
    labels : array-like
        Label for each interval start (same length as edges).

    Returns
    -------
    np.ndarray of object
        Labels for each time; None if time is < min(edges)
        or > max(edges).
    """
    times = np.asarray(times)
    edges = np.asarray(edges)
    labels = np.asarray(labels, dtype=object)

    # sort edges and reorder labels
    order = np.argsort(edges)
    edges = edges[order]
    labels = labels[order]

    # find index of last edge <= time
    idx = np.searchsorted(edges, times, side="right") - 1

    # default to None
    out = np.full(times.shape, None, dtype=object)

    # valid times: within [edges[0], edges[-1]]
    mask = (idx >= 0) & (times <= edges[-1])

    # assign labels for valid times
    out[mask] = labels[idx[mask]]

    return out

def trace_lineage_trait_worker(i, state, treeString, tipRenameDict, maxDate, tipNames, traitName, timeline, headerMode = False):
    """
    Track lineage trait states backward in time for one or more focal tips.

    Parameters
    ----------
    i : int
        Tree index within the posterior sample.

    state : int
        MCMC state associated with the sampled tree.

    treeString : str
        Serialized tree string for the sampled tree.

    tipRenameDict : dict
        Mapping from numeric tip identifiers to display names.

    maxDate : float
        Most recent sampling date used to assign absolute times.

    tipNames : str or list[str]
        Tip name or names whose lineage trait history should be extracted.

    traitName : str
        Trait to read along the path from each focal tip to the root.

    timeline : array-like
        Time points at which trait states should be reported.

    headerMode : bool, optional
        If ``True``, return output column names instead of values.

    Returns
    -------
    tuple
        Tuple ``(i, state, values)`` suitable for the posterior-processing
        pipeline.
    """

    tree = make_tree(treeString, 'time')
    if len(tipRenameDict) > 0:
        tree.rename_tips(tipRenameDict)
    tree.traverse_tree()
    tree.set_absolute_time(maxDate)
    
    if isinstance(tipNames, str): ## was given a single tip
        focalTips = tree.get_external(lambda k: k.is_leaf() and k.name == tipNames) ## get leaf object
        assert len(focalTips) == 1, f'Found {len(focalTips)} tips named {tipNames} in tree.'
    elif isinstance(tipNames, list): ## was given a list of tips
        focalTips = tree.get_external(lambda k: k.is_leaf() and k.name in tipNames) ## get leaf objects
        assert len(focalTips) == len(tipNames), f'Found {len(focalTips)} tips but {len(tipNames)} were expected in tree.\nFound in tree: {', '.join([k.name for k in focalTips if k.name in tipNames])}\nMissing from tree: {', '.join([k.name for k in focalTips if k.name not in tipNames])}'

    output_string = []
    
    for focalTip in focalTips: ## iterate over leaf objects

        if headerMode == False: ## processing tree and generating actual output
            path = focalTip.get_path_to_root() ## ignore last value (parent of root)
        
            branchTimes = np.array([k.absoluteTime for k in path][::-1])
            branchTraits = np.array([k.traits[traitName] for k in path][::-1])
            
            traitSwitchIdx = np.flatnonzero(branchTraits[1:] != branchTraits[:-1]) + 1 ## identify indices where trait switches
            traitSwitchIdx = np.concatenate(([0], traitSwitchIdx, [len(branchTimes)-1])) ## add first and last value to indices
            traitSwitchStates = branchTraits[traitSwitchIdx] ## grab (new) traits after switch

            interpolatedLabels = label_timepoints(timeline, branchTimes[traitSwitchIdx], traitSwitchStates)

            output_string += [traitState if traitState else '' for traitState in interpolatedLabels]
        else: ## doing header
            output_string += [f'{focalTip.name}__{t}' for t in timeline] ## only outputs output file header

    return i, state, output_string

def tmrca_worker(i, state, treeString, tipRenameDict, maxDate, tipNames, strictMRCA = False, headerMode = False):
    """
    Compute TMRCA values for one or more named tip sets in a sampled tree.

    Parameters
    ----------
    i : int
        Tree index within the posterior sample.

    state : int
        MCMC state associated with the sampled tree.

    treeString : str
        Serialized tree string for the sampled tree.

    tipRenameDict : dict
        Mapping from numeric tip identifiers to display names.

    maxDate : float
        Most recent sampling date used to assign absolute times.

    tipNames : list[str] or dict[str, list[str]]
        Tip set definitions for which TMRCA values should be reported.

    strictMRCA : bool, optional
        If ``True``, only report TMRCAs for nodes whose descendant set exactly
        matches the requested tip set.

    headerMode : bool, optional
        If ``True``, return output column names instead of values.

    Returns
    -------
    tuple
        Tuple ``(i, state, values)`` suitable for the posterior-processing
        pipeline.
    """
    tree = make_tree(treeString, 'time')
    if len(tipRenameDict) > 0:
        tree.rename_tips(tipRenameDict)
    tree.traverse_tree()
    tree.set_absolute_time(maxDate)
    
    if isinstance(tipNames, list):
        tipNames = {'tmrca': tipNames} ## convert to dict format for later processing
    
    descendants = {}
    for tmrca in tipNames:
        assert len(tipNames[tmrca]) >= 2, f'Need at least 2 tip names to identify TMRCA, {len(tipNames[tmrca])} provided.'
        descendants[tmrca] = tree.get_external(lambda k: k.name in tipNames[tmrca]) ## if dict - key is tmrca name and value is a list of str corresponding to tip names
        assert len(descendants[tmrca]) == len(tipNames[tmrca]), f'Not all tips could be found: {', '.join([name for name in tipNames if name not in [k.name for k in descendants[tmrca]]])}.'

    out = []
    for tmrca in descendants: ## iterate over required tips
        if headerMode == False:
            commonAncestor = tree.find_MRCA(descendants[tmrca])
            if (strictMRCA and commonAncestor.leaves == set(tipNames[tmrca])) or strictMRCA == False: ## want strict mrca and node matches OR don't need strict mrca
                assert commonAncestor.absoluteTime is not None, f'Node {commonAncestor.index} does not have a set absolute time. Function received maxDate = {maxDate}'
                out += [str(commonAncestor.absoluteTime)]
            elif strictMRCA: ## need strictMRCA which wasn't found, write nan
                out += ['nan']
        else: ## header mode - just output tmrca name
            out += [tmrca]
    
    return i, state, out

def tree_length_worker(i, state, treeString, tipRenameDict, maxDate, headerMode = False):
    # i, state, treeString, tipRenameDict, maxDate
    tree = make_tree(treeString, 'time')
    tree.traverse_tree()

    outputLine = []
    if headerMode == False:
        outputLine.append(sum(tree.get_parameter_list('length')))
    else:
        outputLine.append('treeLength')

    return i, state, outputLine