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

ll.plot_tree(ax, treeType = 'circular')
ll.plot_text(ax, treeType = 'circular', size = 26)
ll.plot_points(ax, treeType = 'circular')

scale_y = ll.ySpan * 0.05
scale_start = ll.treeHeight * 0.25
scale_end = scale_start + 1.0
scale_bar = [
    ll.project_circular_point(x, scale_y)
    for x in (scale_start, scale_end)
]
scale_xs, scale_ys = zip(*scale_bar)
ax.plot(scale_xs, scale_ys, color='k', lw=3)

bt_utils.clean_axes(ax)

plt.savefig('basic-tree-circular.png', bbox_inches='tight')