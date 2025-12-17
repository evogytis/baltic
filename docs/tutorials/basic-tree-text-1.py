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

ll.plot_tree(ax)

def text_function(node):
    """
    Equivalent to anynomous function lambda k: f"node at height {k.height}"
    """
    return f"node at height {node.height}"

def target_function(node):
    """
    Equivalent to anonymous function lambda k: k.is_node()
    """
    return node.is_node()

def x_coordinate_function(node):
    """
    Equivalent to anonymous function lambda k: k.height + 0.1
    """
    return node.height + 0.1

ll.plot_text(ax, xCoordinateFxn = x_coordinate_function, size = 26) ## by default only targets leaf objects
ll.plot_text(ax, xCoordinateFxn = x_coordinate_function, targetFxn = target_function, textContentFxn = text_function, size = 12)
ll.plot_points(ax, targetFxn = target_function)

bt_utils.clean_axes(ax)

plt.show()