## baltic: the Backronymed Adaptable Lightweight Tree Import Code

baltic was initially developed to extract various statistics from molecular phylogenies derived from [BEAST](https://github.com/beast-dev/beast-mcmc) in a customised way. My [influenza B virus reassortment paper](https://dx.doi.org/10.1093/molbev/msu287) used early versions of baltic’s code to look at how the human influenza B virus segment diversity is structured according to genomic background. I’ve since split up the various bits of code into three parts:

--------------------

## baltic
[`baltic.py`](baltic.py) is the tree parser itself. It uses three main classes - node, leaf and tree to import, manipulate and plot BEAST trees with their rich variety of comments. Node and leaf classes have references to the usual set of parameters you would find in a phylogeny - index of character in string designating the branch (a unique identifier within the tree), length, height, position in time (absoluteTime), X and Y coordinates, a dictionary encoding BEAST branch comments (traits) and a reference to their parent (None for root) and a string designating their type (branchType). The node class additionally have a children attribute, which is a list of the node’s children, another list called leaves that contains tip names that descend from that node, a numChildren attribute, which is the length of the leaves list and a childHeight attribute which tracks when the last tip descending from the node existed. The leaf class has two extra attributes called name and numName. Trees drawn from the posterior distribution will usually encode tips as numbers to save space and require a translation map to convert back into actual interpretable tip names. In baltic numName will be the exact name for the tip that was used in the tree string, with functions to allow the translation of numName into name.

baltic.py evolved from a [short linked list script on StackOverflow](http://stackoverflow.com/questions/280243/python-linked-list/280286#280286) and underwent a major overhaul in order to correct [an article](https://dx.doi.org/10.1126/science.aaa5646) that was wrong. The code should be fairly legible (and commented) and highly adaptable to suit anyone’s needs.

--------------------
### Usage basics

By convention baltic is imported as bt:

`import baltic as bt`

When called with a tree string the `make_tree()` function return a baltic tree object:
```python
treeString='((A:1.0,B:2.0):1.0,C:3.0);'

myTree = bt.make_tree(treeString)
```
Note that if you're not using newick, nexus or nextstrain JSON trees (`loadNewick`, `loadNexus`, and `loadJSON` functions respectively) you'll have to write some code to parse out the tree string. baltic will warn the user if it can't parse something. If this happens you should check if your tip names or annotations contain characters that should **never** be found outside of functional tree string bits, such as commas, parentheses or semicolons. Alternatively, it may be that the regexes that are used to parse out tip names or annotations don't cover some special character you use to define your taxa and will require some editing of baltic.py to alleviate the problem. Feel free to raise an issue if this happens.

`make_tree()` is a function that parses the tree string and builds the data structure that is the phylogenetic tree. It works exactly like all other tree parsers:

- Every time an opening parenthesis (`(`) is encountered in the tree string a new instance of `node` class is created. The new class' `.index` attribute is set to the index along the tree string where it was encountered, giving that particular class a unique identifier within the tree string. The `.parent` attribute is set to whatever the previous object encountered was, similarly, since the last encountered object could only be another node, the current node is added to its parents list of children. Finally we set our new node as the 'current' node of the tree and append the node to the list of objects (`.Objects`, which are branches) contained in the tree.

- Every time a string is encountered which may or may not be surrounded by quotation marks (`'` or `"`) or have the beginning of an annotation block (`[`) we create a new `leaf` class. It also receives an `.index` identifier, like the `node` class. Unlike the `node` class, however, the `.numName` attribute is also set as the string that defined the tip. In BEAST trees it will be the number that identifies the tip, but it could also be a regular string.

- Next baltic looks for annotations, which are the blocks in the format `[&parameter1=1.0,parameter2=0.0]`. These are transformed into the `.traits` dictionary for the branch. In this example the branch being parsed would receive a dictionary with two keys: `cur_branch.traits = {'parameter1' : 1.0, 'parameter2' : 0.0}`.

- Annotations should be followed by branch lengths preceded by a colon (`:`). The branch length is assigned to the current branch's `.length` attribute.

- Forks in the tree string are defined as commas (`,`) and ends of clades are defined by closing parentheses (`)`) and both mean that whatever comes next is in relation to the parent branch of whatever branch we were dealing with earlier.

- Finally tree strings are finished with a semi colon (`;`).


Before you can run any analysis you will usually have to traverse the tree such that branch lengths which are available in the tree string are transformed into branch heights:

`myTree.traverse_tree()`

This takes the `.length` attributes of each branch and sums or subtracts them, as appropriate during a tree traversal and the `.height` attribute of each branch (`node` or `leaf` object) is set, where the root of the tree has `.height = 0.0` and the most recent tip is the highest object in the tree. The tree traversal will also set the tree's `.treeHeight` attribute.

If your tree happens to have branch lengths in units of time you can use the `.height` attribute to modify the `.absoluteTime` attribute of each branch, such that the entire tree is calibrated and position correctly in time. This involves finding the most recent tip in the tree (in absolute time), subtracting the `.treeHeight` and adding the `.height` of the each branch.


Most analytic operations will involve looking at each branch individually without referring back to the tree structure much, beyond immediate children or parents of a particular branch. This is done by iterating over the tree's `.Objects` list, which contains all the branches in the tree. If you want to print out the height of each internal branch whose parent had a different trait value you would do it as:

```python
for k in myTree.Objects:
   if isinstance(k,bt.node) ## (or, alternatively if k.branchType=='node')
       if k.traits[myTrait] != k.parent.traits[myTrait]:
           print k.height
```

--------------------

## samogitia

<img src="docs/figures/coa_samogitia.jpg" width=200px>

[`samogitia.py`](samogitia.py) is the heavy-lifting, tree file-wrangling script in the collection. It’s main role is to parse BEAST tree files, use baltic to create tree data structures, which samogitia then manipulates to create BEAST-like log files that can usually be imported into [Tracer](http://tree.bio.ed.ac.uk/software/tracer/) or used in another program.

--------------------

## austechia

<img src="docs/figures/coa_austechia.png" width=200px>

[`austechia.ipynb`](austechia.ipynb) is the fancy Jupyter notebook that takes tree files, usually MCC trees from BEAST, and plots them. It is meant to be part teaching tool to get people to think about how trees are plotted, to allow for highly customisable representations of trees (e.g. Fig 6 in my [MERS-CoV paper](http://dx.doi.org/10.1093/ve/vev023)) and to improve the aesthetics situation in phylogenetics.

--------------------

## galindia

<img src="docs/figures/coa_galindia.jpg" width=200px>

[`galindia.ipynb`](galindia.ipynb) is a notebook that uses baltic to plot JSONs from [nextstrain.org](nextstrain.org) in order to allow customisable, static, publication-ready figures for phylogenies coming from nextstrain's augur pipeline.

--------------------

## curonia

<img src="docs/figures/coa_curonia.jpg" width=200px>

[`curonia.ipynb`](curonia.ipynb) generalises the [notebook](https://github.com/ebov/space-time/blob/master/Scripts/notebooks/EBOV_phylogeography_animation.ipynb) used to animate the [spread of Ebola virus in West Africa](https://www.youtube.com/watch?v=j4Ut4krp8GQ). This notebook should require minimal manual editing to produce similarly styled animation of other study systems.

--------------------
## baltic was used in the following publications:
- Munnink BBO, Nieuwenhuijse DF, Stein M, O'Toole A, Haverkarte M, ..., Meijer A, Rambaut A, van Dissel J, Sikkema R, Timen A, Koopmans M, 2020. _Rapid SARS-CoV-2 whole genome sequencing for informed public health decision making in the Netherlands_, __bioRxiv__: 2020.04.21.050633; [doi](https://doi.org/10.1101/2020.04.21.050633).
- Venkatesh D, Brouwer A, Ellis R, Goujgoulova G, Seekings J, Brown IH, Lewis NS, 2020. _Regional transmission and reassortment of 2.3.4.4b highly pathogenic avian influenza (HPAI) viruses in Bulgarian poultry 2017/18_, __bioRxiv__: 2020.04.14.040592; [doi](https://doi.org/10.1101/2020.04.14.040592).
- Suárez-Esquivel M, Hernández-Mora G, Ruiz-Villalobos N, ..., Thomson NR, Moreno E, Guzmán-Verri C, 2020. _Persistence of Brucella abortus lineages revealed by genomic characterization and phylodynamic analysis_, __PLOS Neglected Tropical Diseases__ 14(4): e0008235; [doi](https://doi.org/10.1371/journal.pntd.0008235).
- Moncla LH, Bedford T, Dussart P, ..., Guan Y, Friedrich TC, Horwood PF, 2020. _Quantifying within-host evolution of H5N1 influenza in humans and poultry in Cambodia_, __PLOS Pathogens__ 16(1): e1008191; [doi](https://doi.org/10.1371/journal.ppat.1008191).
- Dudas G, Bedford T, 2019. _The ability of single genes vs full genomes to resolve time and space in outbreak analysis_. __BMC Evolutionary Biology__ 19: 232; [doi](https://doi.org/10.1186/s12862-019-1567-0).
- Müller NF, Stolz U, Dudas G, Stadler T, Vaughan TG, 2019. _Bayesian inference of reassortment networks reveals fitness benefits of reassortment in human influenza viruses_, __bioRxiv__: 726042; [doi](https://doi.org/10.1101/726042).
- Wiley R, Fakoli L, Letizia AG, Welch SR, Ladner JT, Albariño CG, Fallah M, Palacios G, 2019. _Lassa virus circulating in Liberia: a retrospective genomic characterisation_, __Lancet Infect Dis__
- Theys K, Lemey P, Vandamme AM, Baele G, 2019. _Advances in Visualization Tools for Phylogenomic and Phylodynamic Studies of Viral Diseases_, __Front. Public Health__, 7: 208; [doi](https://doi.org/10.3389/fpubh.2019.00208).
- Müller NF, Dudas G, Stadler T, 2019. _Inferring time-dependent migration and coalescence patterns from genetic sequence and predictor data in structured populations_, __Virus Evolution__, 5(2): vez030; [doi](https://doi.org/10.1093/ve/vez030).
- Mbala-Kingebeni P, Aziza A, Di Paola N, Wiley MR, ..., Peeters M, Palacios G, Ahuka-Mundeke S, 2019. _Medical countermeasures during the 2018 Ebola virus disease outbreak in the North Kivu and Ituri Provinces of the Democratic Republic of the Congo: a rapid genomic assessment_, __Lancet Infect Dis__ 19(6): 648-657; [doi](https://doi.org/10.1016/S1473-3099%2819%2930118-5).
- Work TM, Dagenais J, Stacy BA, Ladner JT, ..., Rameyer RA, Taylor DR, Waltzek TB, 2019. _A novel host-adapted strain of Salmonella Typhimurium causes renal disease in olive ridley turtles (Lepidochelys olivacea) in the Pacific_, __Scientific Reports__ 9: 9313; [doi](https://doi.org/10.1038/s41598-019-45752-5).
- Poen MJ, Venkatesh D, Bestebroer TM, ..., Brown IH, Fouchier RAM, Lewis NS, 2019. _Co-circulation of genetically distinct highly pathogenic avian influenza A clade 2.3.4.4 (H5N6) viruses in wild waterfowl and poultry in Europe and East Asia, 2017–18_, __Virus Evolution__ 5(1): vez004; [doi](https://doi.org/10.1093/ve/vez004).
- van Vuren JP, Ladner JT, Grobbelaar AA, Wiley MR, Lovett S, Allam M, Ismail A, le Roux C, Weyer J, Moolla N, Storm N, Kgaladi J, Sanchez-Lockhart M, Conteh O, Palacios G, Paweska JT, 2019. _Phylodynamic Analysis of Ebola Virus Disease Transmission in Sierra Leone_. __Viruses__, 11(1), 71; [doi](https://doi.org/10.3390/v11010071).
- Bell SM, Katzelnick L, Bedford T, 2018. _Dengue antigenic relationships predict evolutionary dynamics_, __eLife__ 8: e42496; [doi](https://doi.org/10.7554/eLife.42496).
- Wille M, Latorre-Margalef N, ..., Raghwani J, Pybus OG, Olsen B, Waldenström J, 2018. _Where do all the subtypes go? Temporal dynamics of H8–H12 influenza A viruses in waterfowl_, __Virus Evolution__, 4(2): vey025; [doi](https://doi.org/10.1093/ve/vey025).
- Dokubo EK, Wendland A, Mate SE, Ladner JT, ..., Palacios G, Fallah MP, 2018. _Persistence of Ebola virus after the end of widespread transmission in Liberia: an outbreak report_, __Lancet Infect Dis__ 18: 1015–1024; [doi](https://doi.org/10.1016/S1473-3099(18)30417-1).
- Lee JM, Huddleston J, Doud MB, Hooper KA, Wu NC, Bedford T, Bloom JD, 2018. _Deep mutational scanning of hemagglutinin helps predict evolutionary fates of human H3N2 influenza variants_, __PNAS__ 115(35): 8276-8285; [doi](https://doi.org/10.1073/pnas.1806133115).
- Venkatesh D, Poen MJ, Bestebroer TM, ..., Brown IH, Fouchier RAM, Lewis NS, 2018. _Avian influenza viruses in wild birds: virus evolution in a multi-host ecosystem_, __J Virol__ 92: e00433-18; [doi](https://doi.org/10.1128/JVI.00433-18).
- Chu DKW, Hui Kenrie PY, Perera RAPM, Miguel E, Niemeyer D, Zhao J, Channappanavar R, Dudas G, Oladipo JO, Traoré A, Fassi-Fihri O, Ali A, Demissie GF, Muth D, Chan MCW, Nicholls JM, Meyerholz DK, Kuranga SA, Mamo G, Zhou Z, So RTY, Hemida MG, Webby RJ, Roger F, Rambaut A, Poon LLM, Perlman S, Drosten C, Chevalier V, Peiris M, 2018. _MERS coronaviruses from camels in Africa exhibit region-dependent genetic diversity_. __PNAS__ 115(12): 3144-3149; [doi](https://doi.org/10.1073/pnas.1718769115).
- Whitmer SLM, Ladner JT, Wiley MR, Patel K, Dudas G, Rambaut A, Sahr F, Prieto K, Shepard SS, Carmody E, Knust B, Naidoo D, Deen G, Formenty P, Nichol ST, Palacios G, Ströher U, 2018. _Active Ebola Virus Replication and Heterogeneous Evolutionary Rates in EVD Survivors_. __Cell Reports__ 22(5): 1159-1168; [doi](https://doi.org/10.1016/j.celrep.2018.01.008).
- Dudas G, Carvalho L, Rambaut A, Bedford T. _MERS-CoV spillover at the camel-human interface_, 2017. __eLife__ 7: e31257; [doi](https://doi.org/10.7554/eLife.31257).
- Langat P, Raghwani J, Dudas G, ..., Russell C, Pybus OG, McCauley J, Kellam P, Watson SJ. _Genome-wide evolutionary dynamics of influenza B viruses on a global scale_, 2017. __PLOS Pathogens__ 13(12): e1006749; [doi](https://doi.org/10.1371/journal.ppat.1006749).
- Grubaugh ND, Ladner JT, Kraemer MUG, Dudas G, Tan AL, Gangavarapu K, Wiley MR, White S, Thézé J, ..., Bedford T, Pybus OG, Isern S, Palacios G, Andersen KG. _Multiple introductions of Zika virus into the United States revealed through genomic epidemiology_, 2017. __Nature__ 546: 401–405; [doi](https://doi.org/10.1038/nature22400).
- Dudas G, Carvalho LM, Bedford T, Tatem AJ, ..., Suchard M, Lemey P, Rambaut A. _Virus genomes reveal the factors that spread and sustained the West African Ebola epidemic_, 2017. __Nature__ 544(7650): 309-315; [doi](https://doi.org/10.1038/nature22040).
- Bell SM, Bedford T. _Modern-Day SIV viral diversity generated by extensive recombination and cross-species transmission_, 2017. __PLOS Pathogens__ 13(7): e1006466; [doi](https://doi.org/10.1371/journal.ppat.1006466).
- Holmes EC, Dudas G, Rambaut A, Andersen KG. _The Evolution of Ebola virus: Insights from the 2013-2016 Epidemic_, 2016. __Nature__ 538(7624): 193-200; [doi](https://doi.org/10.1038/nature19790).
- Whitmer SLM, Albariño C, Shepard SS, Dudas G, ..., Nichol ST, Ströher U. _Preliminary Evaluation of the Effect of Investigational Ebola Virus Disease Treatments on Viral Genome Sequences_, 2016. __Journal of Infectious Diseases__: jiw177; [doi](https://doi.org/10.1093/infdis/jiw177).
- Rambaut A, Dudas G, Carvalho LM, Park DJ, Yozwiak NL, Holmes EC, Andersen KG. _Comment on “Mutation rate and genotype variation of Ebola virus from Mali case sequences”_, 2016. __Science__ 353(6300): 658-658; [doi](https://doi.org/10.1126/science.aaf3823).
- Lewis NS, Russell CA, Langat P, ..., Dudas G, ..., Watson SJ, Brown IH, Vincent AL. _The global antigenic diversity of swine influenza A viruses_, 2016. __eLife__ 5: e12217; [doi](https://doi.org/10.7554/eLife.12217).
- Quick J, Loman NJ, Duraffour S, Simpson JT, Severi E, Cowley L ..., Dudas G, ..., Günther S, Carroll MW. _Real-time, portable genome sequencing for Ebola surveillance_, 2016. __Nature__ 530(7589): 228-232; [doi](https://doi.org/10.1038/nature16996).
- Dudas G, Rambaut A, _MERS-CoV recombination: implications about the reservoir and potential for adaptation_, 2016. __Virus Evolution__ 2(1): vev023; [doi](https://doi.org/10.1093/ve/vev023).
- Ladner JT, Wiley MR, Mate S, Dudas G, ... Palacios G. _Evolution and Spread of Ebola Virus in Liberia, 2014-2015_, 2015. __Cell Host & Microbe__ 18(6): 659-669; [doi](https://doi.org/10.1016/j.chom.2015.11.008).
- Park DJ, Dudas G, Wohl S, Goba A, Whitmer SLM, ..., Sabeti PC. _Ebola Virus Epidemiology, Transmission, and Evolution during Seven Months in Sierra Leone_, 2015. __Cell__ 161(7): 1516-1526; [doi](https://doi.org/10.1016/j.cell.2015.06.007).
- Carroll MW, Matthews DA, Hiscox JA, ... Dudas G, ... Günther S. _Temporal and spatial analysis of the 2014-2015 Ebola virus outbreak in West Africa_, 2015. __Nature__ 524(7563): 97-101; [doi](https://doi.org/10.1038/nature14594).
- Dudas G, Bedford T, Lycett S, Rambaut A. _Reassortment between Influenza B Lineages and the Emergence of a Coadapted PB1–PB2–HA Gene Complex_, 2015. __Molecular Biology and Evolution__ 32(1): 162-172; [doi](https://doi.org/10.1093/molbev/msu287).
- Gire SK, Goba A, Andersen KG, ... Dudas G, ... Sabeti PC. _Genomic surveillance elucidates Ebola virus origin and transmission during the 2014 outbreak_, 2014. __Science__ 345(6202): 1369-1372; [doi](https://doi.org/10.1126/science.1259657).

--------------------

Copyright 2016 [Gytis Dudas](https://twitter.com/evogytis). Licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).
