import baltic as bt
from baltic import bt_utils

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

treestring = '((("A":1,"B":1):1,"C":2):1,"D":3);'
ll = bt.make_tree(treestring, 'divergence')
ll.treeStats()

fig = plt.figure(figsize=(10,5),facecolor='w')
gs = GridSpec(1,2)

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

descendants = []
for k in ll.get_branches():
    if k.is_leaf() and k.name in ['A', 'C']:
        descendants.append(k)

commonAncestor = ll.find_MRCA(descendants)

subtree_ll = ll.subtree(commonAncestor)

ll.plot_tree(ax)
subtree_ll.plot_tree(ax2)

def target_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_node()
    """
    return node.is_node()

ll.plot_text(ax, size = 26)
ll.plot_points(ax, targetFxn = target_function)

subtree_ll.plot_text(ax2, size = 26)
subtree_ll.plot_points(ax2, targetFxn = target_function)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.plot_scale_bar(ax2, xy = (0.5, 0.1), tree = subtree_ll, style = 'fancy', textKwargs={'fontsize': 18})

bt_utils.clean_axes(ax)
bt_utils.clean_axes(ax2)

plt.savefig('basic-tree-subtree.png', bbox_inches='tight')