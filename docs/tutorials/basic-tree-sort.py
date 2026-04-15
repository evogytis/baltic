import baltic as bt
from baltic import bt_utils

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.use("Agg")

treestring = '((("A":1,"B":1):1,"C":2):1,"D":3);'
ll = bt.make_tree(treestring, 'divergence')
ll.treeStats()

fig = plt.figure(figsize=(5,5),facecolor='w')
gs = GridSpec(1,1)

ax = plt.subplot(gs[0])

def sort_function(node):
    """
    Equivalent to anonymous function lambda k: -1 if k.is_leaf() and k.name == 'C' else 1
    """
    if node.is_leaf() and node.name == 'C':
        return -1
    else:
        return 1

ll.sort_branches(sortFxn = sort_function)

ll.plot_tree(ax)
ll.plot_text(ax, size = 26)
ll.plot_points(ax)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.clean_axes(ax)

plt.savefig('basic-tree-sort.png', bbox_inches='tight')