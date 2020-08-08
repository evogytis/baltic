import argparse
import re
import datetime as dt
import baltic as bt
import sys
import collections

def overlap(a,b):
    """
    Return the elements shared by two lists in the following format:
    [overlap],[list 1 remainder],[list 2 remainder]
    """
    a_multiset = collections.Counter(a)
    b_multiset = collections.Counter(b)

    overlap = list((a_multiset & b_multiset).elements())
    a_remainder = list((a_multiset - b_multiset).elements())
    b_remainder = list((b_multiset - a_multiset).elements())

    return overlap, a_remainder, b_remainder

############ arguments
samogitia = argparse.ArgumentParser(description="samogitia.py analyses trees drawn from the posterior distribution by BEAST.\n")

samogitia.add_argument('-b','--burnin', default=0, type=int, help="Number of states to remove as burnin (default 0).\n")
samogitia.add_argument('-nc','--nocalibration', default=True, action='store_false', help="Use flag to prevent calibration of trees into absolute time (default True). Should be used if tip names do not contain information about *when* each sequence was collected.\n")
samogitia.add_argument('-t','--treefile', type=open, required=True, help="File with trees sampled from the posterior distribution (usually with suffix .trees).\n")
samogitia.add_argument('-a','--analyses', type=str, required=True, nargs='+', help="Analysis to be performed, can be a list separated by spaces.\n")
samogitia.add_argument('-o','--output', type=argparse.FileType('w'), default='samogitia.out.txt', help="Output file name (default samogitia.out.txt).\n")
samogitia.add_argument('-s','--states', type=str, default='0-inf', help="Define range of states for analysis.\n")
samogitia.add_argument('-df','--date_format', type=str, default='%Y-%m-%d', help="Define date format encoded in tips (default \'%%Y-%%m-%%d\').\n")
samogitia.add_argument('-tf','--tip_format', type=str, default='\|([0-9]+)\-*([0-9]+)*\-*([0-9]+)*$', help="Define regex for capturing dates encoded in tips (default \'\|([0-9]+)\-*([0-9]+)*\-*([0-9]+)*$\'.\n")

args = vars(samogitia.parse_args())
burnin, treefile, analyses, outfile, calibration, states, dformat, tformat = args['burnin'], args['treefile'], args['analyses'], args['output'], args['nocalibration'], args['states'], args['date_format'], args['tip_format']

lower,upper=states.split('-')
lower=int(lower)
if upper=='inf':
    upper=float('inf')
else:
    upper=int(upper)

try:
    for line in open('../docs/banner_samogitia.txt','r'):
        sys.stderr.write('%s'%(line))
except:
    pass

## keeps track of things in the tree file at the beginning
plate=True
taxonlist=False
treecount=0 ## counts tree
##

tips={} ## remembers tip encodings

############################# progress bar stuff
Ntrees=10000 ## assume 10 000 trees in posterior sample
barLength=30
progress_update=Ntrees/barLength ## update progress bar every time a tick has to be added
threshold=progress_update ## threshold at which to update progress bar
processingRate=[] ## remember how quickly script processes trees
#############################

available_analyses=['treeLength','RC','Sharp','tmrcas','transitions','subtrees'] ## analysis names that are possible

assert analyses,'No analyses were selected.'
for queued_analysis in analyses: ## for each queued analysis check if samogitia can do anything about them (i.e. whether they're known analysis types)
    assert queued_analysis in available_analyses,'%s is not a known analysis type\n\nAvailable analysis types are: \n* %s\n'%(queued_analysis,'\n* '.join(available_analyses))

begin=dt.datetime.now() ## start timer

for line in treefile: ## iterate through each line
    ###################################################################################
    if plate==True and 'state' not in line.lower():
        cerberus=re.search('Dimensions ntax\=([0-9]+)\;',line) ## Extract useful information from the bits preceding the actual trees.
        if cerberus is not None:
            tipNum=int(cerberus.group(1))

        if 'Translate' in line:
            taxonlist=True ## taxon list to follow

        if taxonlist==True and ';' not in line and 'Translate' not in line: ## remember tip encodings
            cerberus=re.search('([0-9]+) ([\'\"A-Za-z0-9\?\|\-\_\.\/]+)',line)
            tips[cerberus.group(1)]=cerberus.group(2).strip("'")

    if 'tree STATE_' in line and plate==True: ## starting actual analysis
        plate=False
        assert (tipNum == len(tips)),'Expected number of tips: %s\nNumber of tips found: %s'%(tipNum,len(tips)) ## check that correct numbers of tips have been parsed
    ################################################################################### start analysing trees
    cerberus=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line) ## search for crud at the beginning of the line that's not a tree string

    if cerberus is not None: ## tree identified
        ################################################################# at state 0 - create the header for the output file and read the tree (in case the output log file requires information encoded in the tree)
        if treecount==0: ## At tree state 0 insert header into output file
            ll=bt.tree() ## empty tree object
            start=len(cerberus.group()) ## index of where tree string starts in the line
            treestring=str(line[start:]) ## grab tree string
            bt.make_tree(treestring,ll) ## read tree string
            if lower==0 and upper==float('inf'): ## only add a header if not doing a chunk
                outfile.write('state') ## begin the output log file
                ########################################### add header to output log file
                if 'treeLength' in analyses:
                    outfile.write('\ttreeLength')
                ###########################################
                if 'RC' in analyses:
                    outfile.write('\tN\tS\tuN\tuS\tdNdS')
                ###########################################
                if 'tmrcas' in analyses:
                    tmrcas={'A':[],'B':[],'C':[]} ## dict of clade names
                    ll.renameTips(tips)
                    for k in ll.Objects: ## iterate over branches
                        if isinstance(k,bt.leaf): ## only interested in tips
                            if 'A' in k.name: ## if name of tip satisfies condition
                                tmrcas['A'].append(k.numName) ## add the tip's numName to tmrca list - tips will be used to ID the common ancestor
                            elif k.name in Conakry_tips:
                                tmrcas['B'].append(k.numName)
                            tmrcas['C'].append(k.numName)

                    outfile.write('\t%s'%('\t'.join(sorted(tmrcas.keys()))))
                ###########################################
                if 'transitions' in analyses:
                    outfile.write('state\ttotalChangeCount\tcompleteHistory_1')
                ###########################################
                if 'subtrees' in analyses:
                    import copy
                ###########################################
                ## your custom header making code goes here
                ## if 'custom' in analyses:
                ##     trait='yourTrait'
                ##     trait_vals=[]
                ##     for k in ll.Objects:
                ##         if k.traits.has_key(trait):
                ##             trait_vals.append(k.traits[trait])
                ##     available_trait_values=sorted(bt.unique(trait_vals))
                ##     for tr in available_trait_values:
                ##         outfile.write('\t%s.time'%(tr))
                ###########################################
                treecount+1
                outfile.write('\n') ## newline for first tree
        #################################################################
        if int(cerberus.group(1)) >= burnin and lower <= int(cerberus.group(1)) < upper: ## After burnin start processing
            ll=bt.tree() ## ll is the tree object
            start=len(cerberus.group()) ## find start of tree string in line
            treestring=str(line[start:]) ## get tree string
            bt.make_tree(treestring,ll) ## pass it to make_tree function
            ll.traverse_tree() ## Traverse the tree - sets the height of each object in the tree
            #### renaming tips
            if len(tips)>0:
                ll.renameTips(tips) ## Rename tips so their name refers to sequence name
            else:
                for k in ll.Objects:
                    if isinstance(k,leaf):
                        k.name=k.numName ## otherwise every tip gets a name that's the same as tree string names
            #### calibration
            dateCerberus=re.compile(tformat) ## search pattern + brackets on actual calendar date
            if calibration==True: ## Calibrate tree so everything has a known position in actual time
                tipDatesRaw=[dateCerberus.search(x).group(1) for x in tips.values()]
                tipDates=[]
                for tip_date in tipDatesRaw:
                    tipDates.append(bt.decimalDate(tip_date,fmt=dformat,variable=True))
                maxDate=max(tipDates) ## identify most recent tip
                ll.setAbsoluteTime(maxDate)
            outfile.write('%s'%cerberus.group(1)) ## write MCMC state number to output log file
            ################################################################################
            if 'treeLength' in analyses:
                treeL=sum([k.length for k in ll.Objects]) ## do analysis
                outfile.write('\t%s'%(treeL)) ## output to file
            ###################################################
            if 'RC' in analyses: ## 'RC' was queued as an analysis
                Ns=[] ## empty list
                Ss=[]
                uNs=[]
                uSs=[]
                for k in ll.Objects: ## iterate over branch objects in the tree
                    if k.traits.has_key('N'): ## if branch has a trait labelled "N"...
                        Ns.append(k.traits['N']) ## add it to empty list
                        Ss.append(k.traits['S']) ## likewise for every other trait
                        uNs.append(k.traits['b_u_N'])
                        uSs.append(k.traits['b_u_S'])
                tNs=sum(Ns) ## sum of numbers in list
                tSs=sum(Ss)
                tuNs=sum(uNs)
                tuSs=sum(uSs)
                dNdS=(tNs/tSs)/(tuNs/tuSs) ## calculate dNdS
                outfile.write('\t%s\t%s\t%s\t%s\t%s'%(tNs,tSs,tuNs,tuSs,dNdS)) ## output to file, separated by tabs
            ###################################################
            if 'tmrcas' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                nodes={x:None for x in tmrcas.keys()} ## each TMRCA will correspond to a single object
                score={x:len(ll.Objects)+1 for x in tmrcas.keys()} ## this will be used to score the common ancestor candidate
                for required in tmrcas.keys(): ## iterate over TMRCAs
                    searchNodes=sorted([nd for nd in ll.Objects if nd.branchType=='node' and len(nd.leaves)>=len(tmrcas[required])],key=lambda n:len(n.leaves)) ## common ancestor candidates must have at least as many descendants as the list of search tips
                    for k in searchNodes: ## iterate over candidates
                        common,queryLeft,targetLeft=overlap(k.leaves,tmrcas[required]) ## find how many query tips exist as descendants of candidate nodes
                        if len(targetLeft)==0 and len(queryLeft)<=score[required]: ## all of query tips must be descended from common ancestor, every extra descendant of common ancestor not in query contributes to a score
                            nodes[required]=k ## if score improved - assign new common ancestor
                            score[required]=len(queryLeft) ## score is extra descendants not in the list of known tips

                outTMRCA=['%.6f'%(nodes[n].absoluteTime) for n in sorted(nodes.keys())] ## fetch absoluteTime of each identified common ancestor
                outfile.write('\t%s'%('\t'.join(outTMRCA)))
            ###################################################
            if 'Sharp' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                assert len(analyses)==1,'More that one analysis queued in addition to Sharp, which is inadvisable'
                outSharp=[]
                for k in ll.Objects:
                    if k.traits.has_key('N'):
                        N=k.traits['N']
                        S=k.traits['S']
                        halfBranch=k.length*0.5
                        if isinstance(k,node):
                            all_leaves=[tips[lf] for lf in k.leaves]
                            t=min(map(bt.decimalDate,[dateCerberus.search(x).group(1) for x in all_leaves]))-k.absoluteTime+halfBranch
                        else:
                            t=halfBranch

                        outSharp.append('(%d,%d,%.4f)'%(N,S,t))
                outfile.write('\t%s'%('\t'.join(outSharp)))
            ###################################################
            if 'transitions' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                assert len(analyses)==1,'More that one analysis queued in addition to transitions, which is inadvisable'
                outTransitions=[]
                for k in ll.Objects:
                    if k.traits.has_key('location.states') and k.parent.traits.has_key('location.states'):
                        cur_value=k.traits['location.states']
                        par_value=k.parent.traits['location.states']
                        if cur_value!=par_value:
                            outTransitions.append('{1,%s,%s,%s}'%(ll.treeHeight-k.height-0.5*k.length,par_value,cur_value))
                outfile.write('\t%d\t%s'%(len(outTransitions),'\t'.join(outTransitions)))
            ###################################################
            if 'subtrees' in analyses:
                traitName='location.states'
                assert [k.traits.has_key(traitName) for k in ll.Objects].count(True)>0,'No branches have the trait "%s"'%(traitName)
                for k in ll.Objects:
                    if k.traits.has_key(traitName) and k.parent.traits.has_key(traitName) and k.parent.index!='Root' and k.traits[traitName]!=k.parent.traits[traitName] and k.traits[traitName]=='human':
                        proceed=False ## assume we still can't proceed forward
                        kloc=k.traits[traitName]
                        if isinstance(k,bt.leaf): ## if dealing with a leaf - proceed
                            N_children=1
                            proceed=True
                        else:
                            N_children=len(k.leaves)
                            if [ch.traits[traitName] for ch in k.children].count(kloc)>0:
                                proceed=True
                        #print k.index,k.parent.index,k.traits,k.parent.traits,k.traits[traitName]

                        if proceed==True: ## if at least one valid tip and no hanging nodes
                            subtree=copy.deepcopy(ll.traverseWithinTrait(k,traitName))
                            subtree_leaves=[x.name for x in subtree if isinstance(x,bt.leaf)]

                            if len(subtree_leaves)>0:
                                mostRecentTip=max([bt.decimalDate(x.strip("'").split('|')[-1]) for x in subtree_leaves])
                                while sum([len(nd.children)-sum([1 if ch in subtree else 0 for ch in nd.children]) for nd in subtree if isinstance(nd,bt.node) and nd.index!='Root'])>0: ## keep removing nodes as long as there are nodes with children that are not entirely within subtree
                                    for nd in sorted([q for q in subtree if isinstance(q,bt.node)],key=lambda x:(sum([1 if ch in subtree else 0 for ch in x.children]),x.height)): ## iterate over nodes in subtree, starting with ones that have fewest valid children and are more recent

                                        child_status=[1 if ch in subtree else 0 for ch in nd.children] ## check how many children of current node are under the right trait value

                                        if sum(child_status)<2 and nd.index!='Root': ## if less than 2 children in subtree (i.e. not all children are under the same trait state)
                                            #print 'removing: %d, children in: %s'%(nd.index,[location_to_country[ch.traits[traitName]] for ch in nd.children])
                                            grand_parent=nd.parent ## fetch grandparent of node to be removed
                                            grand_parent.children.remove(nd) ## remove node from its parent's children

                                            if sum(child_status)==0: ## node has no valid children - current grandparent will be removed on next iteration, since it will only have one child
                                                pass
                                            else: ## at least one child is still valid - reconnect the one valid child to grandparent
                                                child=nd.children[child_status.index(1)] ## identify the valid child
                                                child.parent=grand_parent ## child's parent is now its grandparent
                                                grand_parent.children.append(child) ## child is now child of grandparent
                                                child.length+=nd.length ## child's length now includes it's former parent's length
                                            subtree.remove(nd) ## remove node from subtree
                                outfile.write('\t{%s,%s,%s,%s,%d}'%(k.absoluteTime,mostRecentTip,k.parent.traits[traitName],k.traits[traitName],len(subtree_leaves)))
                                sys.stderr.write('\t{%s,%s,%s,%s,%d}'%(k.absoluteTime,mostRecentTip,k.parent.traits[traitName],k.traits[traitName],len(subtree_leaves)))
                                ##########
                                ## Comment out to output stats rather than trees
                                ##########
#                                 if len(subtree)>0: ## only proceed if there's at least one tip in the subtree
#                                     local_tree=bt.tree() ## create a new tree object where the subtree will be
#                                     local_tree.Objects=subtree ## assign branches to new tree object
#                                     local_tree.root.children.append(subtree[0]) ## connect tree object's root with subtree
#                                     subtree[0].parent=local_tree.root ## subtree's root's parent is tree object's root
#                                     #local_tree.root.absoluteTime=subtree[0].absoluteTime-subtree[0].length ## root's absolute time is subtree's root time
#                                     local_tree.sortBranches() ## sort branches, draw small tree
#                                     subtreeString=local_tree.toString()
#                                     outfile.write('\t%s'%(subtreeString))
            ###################################################
            ## your analysis and output code goes here, e.g.
            ## if 'custom' in analyses:
            ##     out={x:0.0 for x in available_trait_values}
            ##     for k in ll.Objects:
            ##         if k.traits.has_key(trait):
            ##             out[k.traits[trait]]+=k.length
            ##     for tr in available_trait_values:
            ##         outfile.write('\t%s'%(out[tr]))
            ###################################################
            outfile.write('\n') ## newline for post-burnin tree
        treecount+=1 ## increment tree counter
        ################################################################################
        if treecount==threshold: ## tree passed progress bar threshold
            timeTakenSoFar=dt.datetime.now()-begin ## time elapsed
            timeElapsed=float(divmod(timeTakenSoFar.total_seconds(),60)[0]+(divmod(timeTakenSoFar.total_seconds(),60)[1])/float(60))
            timeRate=float(divmod(timeTakenSoFar.total_seconds(),60)[0]*60+divmod(timeTakenSoFar.total_seconds(),60)[1])/float(treecount+1) ## rate at which trees have been processed
            processingRate.append(timeRate) ## remember rate
            ETA=(sum(processingRate)/float(len(processingRate))*(Ntrees-treecount))/float(60)/float(60) ## estimate how long it'll take, given mean processing rate

            excessiveTrees=treecount
            if treecount>=10000:
                excessiveTrees=10000
            if timeElapsed>60.0: ## took over 60 minutes
                reportElapsed=timeElapsed/60.0 ## switch to hours
                reportUnit='h' ## unit is hours
            else:
                reportElapsed=timeElapsed ## keep minutes
                reportUnit='m'

            sys.stderr.write('\r') ## output progress bar
            sys.stderr.write("[%-30s] %4d%%  trees: %5d  elapsed: %5.2f%1s  ETA: %5.2fh (%6.1e s/tree)" % ('='*(excessiveTrees/progress_update),treecount/float(Ntrees)*100.0,treecount,reportElapsed,reportUnit,ETA,processingRate[-1]))
            sys.stderr.flush()

            threshold+=progress_update ## increment to next threshold
        ################################################################################
        if 'End;' in line:
            pass
outfile.close()
sys.stderr.write('\nDone!\n') ## done!
