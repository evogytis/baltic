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

keep = []
for k in ll.get_external():
    if k.name in ['A', 'C', 'D']:
        keep.append(k)

reduced_ll = ll.reduce_tree(keep)
reduced_ll.sort_branches(descending=False)

ll.plot_tree(ax)
reduced_ll.plot_tree(ax2)

def target_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_node()
    """
    return node.is_node()

ll.plot_text(ax, size = 26)
ll.plot_points(ax, targetFxn = target_function)

reduced_ll.plot_text(ax2, size = 26)
reduced_ll.plot_points(ax2, targetFxn = target_function)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.plot_scale_bar(ax2, xy = (0.5, 0.1), tree = reduced_ll, style = 'fancy', textKwargs={'fontsize': 18})

bt_utils.clean_axes(ax)
bt_utils.clean_axes(ax2)

plt.savefig('basic-tree-reduced.png', bbox_inches='tight')