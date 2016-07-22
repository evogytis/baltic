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
austechia = argparse.ArgumentParser(description="samogitia.py analyses trees drawn from the posterior distribution by BEAST.\n")

austechia.add_argument('-b','--burnin', default=0, type=int, help="Number of states to remove as burnin (default 10 million).\n")
austechia.add_argument('-c','--calibrate', default=True, type=bool, help="Set the tree to have absolute time information (default True). Should only be used if tip names contain information about *when* each sequence was collected.\n")
austechia.add_argument('-t','--treefile', type=open, help="File with trees sampled from the posterior distribution (usually with suffix .trees).\n")
austechia.add_argument('-a','--analyses', type=str, nargs='+', help="Analysis to be performed, can be a list separated by spaces.\n")
austechia.add_argument('-o','--output', type=argparse.FileType('w'), default='samogitia.out.txt', help="Output file name (default samogitia.out.txt).\n")

args = vars(austechia.parse_args())
burnin, treefile, analyses, outfile, calibration = args['burnin'], args['treefile'], args['analyses'], args['output'], args['calibrate']

for line in open('banner_samogitia.txt','r'):
    sys.stderr.write('%s'%(line))

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

available_analyses=['treeLength','RC','Sharp','tmrcas'] ## analysis names that are possible

assert analyses,'No analyses were selected.'
for queued_analysis in analyses: ## for each queued analysis check if austechia can do anything about them (i.e. whether they're known analysis types)
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
            outfile.write('state') ## begin the output log file
            ########################################### add header to output log file
            if 'treeLength' in analyses:
                outfile.write('\ttreeLength')
			###########################################
            if 'RC' in analyses:
                outfile.write('\tN\tS\tuN\tuS\tdNdS')
            ###########################################
            if 'tmrcas' in analyses:
                tmrcas={'Kailahun':[],'Root':[],'Coyah':[]}
                Coyah_tips=['EBOV|KG87||GIN|Boke|2015-06-19','EBOV|GUI_CTS_2015_0050||GIN|Boke|2015-06-20','EBOV|KG168||GIN|Boke|2015-06-29','EBOV|KG90||GIN|Boke|2015-06-19','EBOV|KG35||GIN|Boke|2015-06-08','EBOV|KG167||GIN|Boke|2015-06-29','EBOV|KG91||GIN|Boke|2015-06-20','EBOV|KG88||GIN|Boke|2015-06-19','EBOV|GUI_CTS_2015_0051||GIN|Boke|2015-06-21','EBOV|GUI_CTS_2015_0052||GIN|Boke|2015-06-25','EBOV|KG80||GIN|Boke|2015-06-18','EBOV|EM_COY_2015_017574||GIN|Boke|2015-06-10','EBOV|KG45||GIN|Boke|2015-06-09','EBOV|IPDPFHGINSP_GUI_2015_6899||GIN|Boke|2015-05-14','EBOV|EM_COY_2015_016483||GIN|Boke|2015-05-13','EBOV|EM_COY_2015_016743||GIN|Boke|2015-05-19','EBOV|KG12||GIN|Boke|2015-05-27','EBOV|EM_COY_2015_017091||GIN|Boke|2015-05-28','EBOV|IPDPFHGINSP_GUI_2015_5339||GIN|Conakry|2015-04-08','EBOV|EM_COY_2015_015802||GIN|Coyah|2015-04-14','EBOV|EM_COY_2015_014261||GIN|Forecariah|2015-04-02','EBOV|IPDPFHGINSP_GUI_2015_4909||GIN|Conakry|2015-03-29','EBOV|IPDPFHGINSP_GUI_2015_5117||GIN|Dubreka|2015-04-03','EBOV|EM_COY_2015_014100||GIN|Conakry|2015-03-26','EBOV|IPDPFHGINSP_GUI_2015_4786||GIN|Conakry|2015-03-26','EBOV|EM_COY_2015_014102||GIN|Conakry|2015-03-26','EBOV|EM_COY_2015_014098||GIN|Conakry|2015-03-26','EBOV|EM_COY_2015_013731||GIN|Coyah|2015-03-14','EBOV|EM_COY_2015_016642||GIN|Forecariah|2015-05-17','EBOV|IPDPFHGINSP_GUI_2015_7070||GIN|Forecariah|2015-05-18','EBOV|EM_FORE_2015_875||GIN|Forecariah|2015-06-17','EBOV|EM_FORE_2015_816||GIN|Forecariah|2015-06-15','EBOV|EM_COY_2015_015982||GIN|Forecariah|2015-04-20','EBOV|EM_COY_2015_016236||GIN|Forecariah|2015-05-03','EBOV|EM_COY_2015_016293||GIN|Forecariah|2015-05-06','EBOV|EM_COY_2015_016263||GIN|Forecariah|2015-05-04','EBOV|IPDPFHGINSP_GUI_2015_6505||GIN|Forecariah|2015-05-03','EBOV|CON12930||GIN|Conakry|2015-10-13','EBOV|EM_COY_2015_017018||GIN|Forecariah|2015-05-26','EBOV|EM_COY_2015_016449||GIN|Forecariah|2015-05-11','EBOV|PL4432|KU296481|SLE|Kambia|2015-03-07','EBOV|PL5249|KU296744|SLE|Kambia|2015-03-24','EBOV|PL5294|KU296764|SLE|Kambia|2015-03-24','EBOV|KT5317|KU296352|SLE|WesternRural|2015-02-24','EBOV|EM_COY_2015_013576||GIN|Coyah|2015-03-09','EBOV|EM_COY_2015_017021||GIN|Fria|2015-05-25','EBOV|1327|KR534591|GIN|Coyah|2014-10-07','EBOV|1320|KR534590|GIN|Coyah|2014-10-07','EBOV|1551|KR534562|GIN|Conakry|2014-10-18','EBOV|1249|KR534541|GIN|Conakry|2014-10-04','EBOV|1374|KR534556|GIN|Coyah|2014-10-10','EBOV|1652|KR534571|GIN|Coyah|2014-10-24','EBOV|1316|KR534546|GIN|Coyah|2014-10-07','EBOV|1331|KR534548|GIN|Kerouane|2014-10-07','EBOV|1298|KR534545|GIN|Conakry|2014-10-06','EBOV|1278|KR534544|GIN|Coyah|2014-10-04','EBOV|1341||GIN|Coyah|2014-10-08','EBOV|1686|KR534572|GIN|Coyah|2014-10-24','EBOV|1689|KR534573|GIN|Coyah|2014-10-25','EBOV|1436|KR534558|GIN|Coyah|2014-10-13','EBOV|1213|KR534583|GIN|Conakry|2014-10-02','EBOV|1193|KR534537|GIN|Conakry|2014-09-28','EBOV|1648|KR534569|GIN|Kindia|2014-10-23','EBOV|1047|KR534528|GIN|Kindia|2014-09-20','EBOV|1690|KR534574|GIN|Coyah|2014-10-25','EBOV|1277|KR534543|GIN|Coyah|2014-10-04','EBOV|1333|KR534549|GIN|Coyah|2014-10-07','EBOV|1129|KR534581|GIN|Conakry|2014-09-25','EBOV|EM_COY_2015_013857||GIN|Forecariah|2015-03-18']
                
                ll.renameTips(tips)
                for k in ll.Objects:
                    if isinstance(k,bt.leaf):
                        if 'Kailahun' in k.name:
                            tmrcas['Kailahun'].append(k.numName)
                        elif k.name in Coyah_tips:
                            tmrcas['Coyah'].append(k.numName)
                        tmrcas['Root'].append(k.numName)
                outfile.write('\t%s'%('\t'.join(sorted(tmrcas.keys()))))
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
            outfile.write('\n') ## newline
        
        #################################################################
        if int(cerberus.group(1)) >= burnin: ## After burnin start processing
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
            dateCerberus=re.compile('\|([0-9\-]+)$') ## search pattern + brackets on actual calendar date
            if calibration==True: ## Calibrate tree so everything has a known position in actual time
                tipDatesRaw=[dateCerberus.search(x).group(1) for x in tips.values()]
                tipDates=map(bt.decimalDate,tipDatesRaw)
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
                nodes={x:None for x in tmrcas.keys()}
                score={x:len(ll.Objects)+1 for x in tmrcas.keys()}
                for k in ll.Objects:
                    for required in tmrcas.keys():
                        if isinstance(k,bt.node) and len(k.leaves)>=len(tmrcas[required]):
                            common,queryLeft,targetLeft=overlap(k.leaves,tmrcas[required])
                            if len(targetLeft)==0 and len(queryLeft)<=score[required]:
                                nodes[required]=k
                                score[required]=len(queryLeft)
                                
                outTMRCA=['%.6f'%(nodes[n].absoluteTime) for n in sorted(nodes.keys())]
                outfile.write('\t%s'%('\t'.join(outTMRCA)))
            ###################################################
            if 'Sharp' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                outSharp=[]
                for k in ll.Objects:
                    if k.traits.has_key('N'):
                        N=k.traits['N']
                        S=k.traits['S']
                        halfBranch=k.length*0.5
                        if isinstance(k,node):
                            all_leaves=[tips[lf] for lf in k.leaves]
                            t=min(map(decimalDate,[dateCerberus.search(x).group(1) for x in all_leaves]))-k.absoluteTime+halfBranch
                        else:
                            t=halfBranch
                            
                        outSharp.append('(%d,%d,%.4f)'%(N,S,t))
                outfile.write('\t%s'%('\t'.join(outSharp)))
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
            outfile.write('\n') ## newline
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