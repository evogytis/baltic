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

inwardSpace = -5.0
ll.plot_tree(ax, inwardSpace = inwardSpace, treeType = 'circular')
ll.plot_text(ax, inwardSpace = inwardSpace, treeType = 'circular', size = 26)
ll.plot_points(ax, inwardSpace = inwardSpace, treeType = 'circular')

bt_utils.clean_axes(ax)

plt.savefig('basic-tree-circular-inward.png', bbox_inches='tight')