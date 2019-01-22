
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim


a = np.arange(26)
ybins = [x / 20. for x in range(8)]
count = 0

natural_val = np.linspace(0, 1.415, num=8, endpoint=True)
natural_val = natural_val.round(decimals=2)
per_segments = [10**x for x in natural_val]
per_segments = [i.round(decimals=2) for i in per_segments]
rprs_segments = [x / 20 for x in range(8)]

block_percent = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]


with open('/Users/sheilasagear/Dropbox/K2/heatmap_inclination_csv/heatmaps_inprogress/inclination/heatmaps.txt') as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    for x in content:
        count+=1

    for f in content:#for each csv file
        with open(str(f)) as csvfile:
            reader = csv.reader(csvfile)
            block_lst = list(reader)
            for i in range(len(block_lst)):
                for j in range(len(block_lst[0])):
                    block_lst[i][j] = float(block_lst[i][j])

            block_lst = np.asarray(block_lst)

            for i in range(len(block_percent)-1):
                for j in range(len(block_percent[0])-1):
                    block_percent[i][j] = block_percent[i][j]+block_lst[i][j]

for i in range(len(block_percent)):
    for j in range(len(block_percent[0])):
        block_percent[i][j] = block_percent[i][j]/count

block_percent = np.asarray(block_percent)
print(np.mean(block_percent))
print(block_percent)

np.savetxt('/Users/sheilasagear/Dropbox/K2/surfacetest/pnull/block_percent/pnull_blockpercent.csv', block_percent, fmt='%f', delimiter=',')


block_percent = np.flip(block_percent, axis=0)

matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['xtick.minor.width'] = 0

fig, ax = plt.subplots()
heatmap = ax.pcolormesh(per_segments, ybins, block_percent, cmap='Blues_r')
ax.set_xscale('log')

ax.xaxis.set_ticks(per_segments)
ax.xaxis.set_ticklabels(per_segments)

ax.yaxis.set_ticks(ybins)
ax.yaxis.set_ticklabels(ybins)

ax.set_title('Fraction of Transits Recovered')
ax.set_xlabel('Period (days)')
ax.set_ylabel('Radius (stellar radii)')


cbar = plt.colorbar(heatmap)
heatmap.set_clim(0.0, 0.03)
cbar.set_label('Fraction recovered', rotation=270, labelpad=13)


#this is modified from https://stackoverflow.com/questions/25071968/heatmap-with-text-in-each-cell-with-matplotlibs-pyplot
def show_values(fig, ax, pc, fmt="%.2f", **kw):
    from itertools import izip
    pc.update_scalarmappable()
    ax1 = ax.axes
    for p, color, value in izip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
        x, y = p.vertices[:-1, :].mean(0)
        if np.all(color[:3] > 0.5):
            color = (0.0, 0.0, 0.0)
        else:
            color = (1.0, 1.0, 1.0)
        ax1.text(x, y, fmt % value, ha="center", va="center", color=color, **kw)

show_values(fig, ax, heatmap)

plt.savefig('/Users/sheilasagear/Dropbox/ssagear_k2/plots/merged_heatmap_1_with_vals.png', dpi = 400)


plt.show()
