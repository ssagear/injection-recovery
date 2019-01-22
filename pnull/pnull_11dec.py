

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import csv

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

a = np.arange(26)
ybins = [x / 20. for x in range(8)]
count = 0

natural_val = np.linspace(0, 1.415, num=8, endpoint=True)
natural_val = natural_val.round(decimals=2)
per_segments = [10**x for x in natural_val]
per_segments = [i.round(decimals=2) for i in per_segments]
rprs_segments = [x / 20 for x in range(8)]

block_percent = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]


with open('/Users/sheilasagear/Dropbox/K2/surfacetest/pnull/block_percent/pnull_blockpercent.csv') as csvfile:
    reader = csv.reader(csvfile)
    block_lst = list(reader)
    for i in range(len(block_lst)):
        for j in range(len(block_lst[0])):
            block_lst[i][j] = float(block_lst[i][j])

    for i in range(len(block_percent)):
        for j in range(len(block_percent[0])):
            block_percent[i][j] = block_percent[i][j]+block_lst[i][j]


block_percent = np.asarray(block_percent)
block_percent = np.flip(block_percent, axis=0)

eta = np.linspace(0, 1, 100000)
n_stars = 500
p_target = .05


p = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]
p_null = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]
p_close = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]
p_null_lst = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]
p_ind = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]
eta_p = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]


for i in range(len(block_percent)):
    for j in range(len(block_percent[0])):
        p[i][j] = eta*block_percent[i][j]
        p_null[i][j] = (1-p[i][j])**n_stars
        p_close[i][j] = find_nearest(p_null[i][j], p_target)
        p_null_lst[i][j] = p_null[i][j].tolist()
        p_ind[i][j] = p_null_lst[i][j].index(p_close[i][j])
        eta_p[i][j] = eta[p_ind[i][j]]
        if eta_p[i][j] == 0.:
            eta_p[i][j] = 1.

eta_p = np.asarray(eta_p)
print(np.mean(eta_p))
print(eta_p)


matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['xtick.minor.width'] = 0

fig, ax = plt.subplots()
heatmap = plt.pcolor(per_segments, ybins, eta_p, cmap='Blues_r')
ax.set_xscale('log')

ax.xaxis.set_ticks(per_segments)
ax.xaxis.set_ticklabels(per_segments)

ax.yaxis.set_ticks(ybins)
ax.yaxis.set_ticklabels(ybins)

ax.set_title('Upper Limit of $\eta$, $P_{null}$ = 0.05')
ax.set_xlabel('Period (days)')
ax.set_ylabel('Radius (stellar radii)')


cbar = plt.colorbar(heatmap)
heatmap.set_clim(0.3, 1)
cbar.set_label('$\eta$', rotation=270, labelpad=13)

#this is modified from https://stackoverflow.com/questions/25071968/heatmap-with-text-in-each-cell-with-matplotlibs-pyplot
def show_values(fig, ax, pc, fmt="%.2f", **kw):
    from itertools import izip
    pc.update_scalarmappable()
    ax1 = ax.axes
    for p, color, value in izip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
        x, y = p.vertices[:-2, :].mean(0)
        if np.all(color[:3] > 0.5):
            color = (0.0, 0.0, 0.0)
        else:
            color = (1.0, 1.0, 1.0)
        ax1.text(x, y, fmt % value, ha="center", va="center", color=color, **kw)

show_values(fig, ax, heatmap)

plt.savefig('/Users/sheilasagear/Dropbox/ssagear_k2/plots/pnull_5percent.png', dpi = 400)

plt.show()
