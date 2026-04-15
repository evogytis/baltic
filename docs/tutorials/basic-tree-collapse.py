import baltic as bt
from baltic import bt_utils

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.use("Agg")

treestring = '((("A":1,"B":1):1,"C":2):1,"D":3);'
ll = bt.make_tree(treestring, 'divergence')
ll.treeStats()

fig = plt.figure(figsize=(10,5),facecolor='w')
gs = GridSpec(1,2)

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

def remove_AB_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_node() and k.leaves == set(['A', 'B'])
    """
    return node.is_node() and node.leaves == set(['A', 'B'])

collapse_ll = ll.collapse_branches(collapseIfFxn = remove_AB_function)
collapse_ll.sort_branches(descending = False)

ll.plot_tree(ax)
collapse_ll.plot_tree(ax2)

def target_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_node()
    """
    return node.is_node()

ll.plot_text(ax, size = 26)
ll.plot_points(ax, targetFxn = target_function)

collapse_ll.plot_text(ax2, size = 26)
collapse_ll.plot_points(ax2, targetFxn = target_function)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.plot_scale_bar(ax2, xy = (0.5, 0.1), tree = collapse_ll, style = 'fancy', textKwargs={'fontsize': 18})

bt_utils.clean_axes(ax)
bt_utils.clean_axes(ax2)

plt.savefig('basic-tree-collapse.png', bbox_inches='tight')