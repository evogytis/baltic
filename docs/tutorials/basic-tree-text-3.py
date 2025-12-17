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
    Equivalent to anynomous function lambda k: f"this branch leads to tip {k.name}\n(and has branch length > 1.0)"
    """
    return f"this branch leads to tip {node.name}\n(and has branch length > 1.0)"

def target_function(node):
    """
    Equivalent to anonymous function lambda k: k > 1.0
    """
    return node.length > 1.0

def x_coordinate_function(node):
    """
    Equivalent to anonymous function lambda k: k.height - node.length * 0.5
    """
    return node.height - node.length * 0.5

ll.plot_text(ax, size = 20)
ll.plot_text(ax, xCoordinateFxn = x_coordinate_function, targetFxn = target_function, textContentFxn = text_function, 
             horizontalalignment = 'center', verticalalignment = 'bottom', size = 10)

ll.plot_points(ax, pointSize = 10, xCoordinateFxn = x_coordinate_function, targetFxn = target_function)
bt_utils.clean_axes(ax)

plt.show()