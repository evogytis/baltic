import baltic as bt
from baltic import bt_utils

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

treestring = '((("A":1,"B":1):1,"C":2):1,"D":3);'
ll = bt.make_tree(treestring, 'divergence')
ll.treeStats()

fig = plt.figure(figsize=(5,5),facecolor='w')
gs = GridSpec(1,1)

ax = plt.subplot(gs[0])

descendants = []
for k in ll.get_branches():
    if k.is_leaf() and k.name in ['A', 'B']:
        descendants.append(k)

commonAncestor = ll.find_MRCA(descendants)
ll.collapse_subtree_to_clade(commonAncestor, 'AB')

ll.plot_tree(ax)

def target_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_leaflike()
    """
    return node.is_leaflike()

def x_coordinate_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_leaflike()
    """
    if isinstance(node, bt.Clade):
        return node.lastHeight
    else:
        return node.height

ll.plot_text(ax, targetFxn = target_function, xCoordinateFxn = x_coordinate_function, size = 26)
ll.plot_points(ax)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.clean_axes(ax)

plt.show()