import baltic as bt
from baltic import bt_utils

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.use("Agg")

treestring = '((("A":1,"B":1):1,"C":1.5):1,"D":0.5);'
ll = bt.make_tree(treestring, 'divergence')
ll.treeStats()

fig = plt.figure(figsize=(10,5),facecolor='w')
gs = GridSpec(1,2)

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

ll.plot_tree(ax) ## plot tree as is
ll.plot_text(ax, size = 26)
ll.plot_points(ax)

ll.midpoint_root() ## reroot

ll.plot_tree(ax2) ## plot again after midpoint rooting
ll.plot_text(ax2, size = 26)
ll.plot_points(ax2)

bt_utils.plot_scale_bar(ax, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})
bt_utils.plot_scale_bar(ax2, xy = (0.5, 0.1), tree = ll, style = 'fancy', textKwargs={'fontsize': 18})

bt_utils.clean_axes(ax)
bt_utils.clean_axes(ax2)

plt.savefig('basic-tree-reroot.png', bbox_inches='tight')