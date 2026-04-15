import baltic as bt
from baltic import bt_utils

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.use("Agg")

treestring = '((("A":1,"B":1):1,"C":1):1,"D":1);'
ll = bt.make_tree(treestring, 'divergence')
ll.treeStats()

fig = plt.figure(figsize=(5,5),facecolor='w')
gs = GridSpec(1,1)

ax = plt.subplot(gs[0])

ll.plot_tree(ax)
ll.plot_points(ax)

ll.plot_aligned_tip_labels(ax, xSpace = 0.5, size = 30)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.clean_axes(ax)

plt.savefig('basic-tree-aligned-tips.png', bbox_inches='tight')