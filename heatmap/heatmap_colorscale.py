import matplotlib.pyplot as plt
import matplotlib
import numpy as np


countsT = np.loadtxt(open("/Users/sheilasagear/Dropbox/K2/heatmap_inclination_csv/heatmaps_inprogress/inclination/heatmap229227257_fullplanets.csv", "rb"), delimiter=",", skiprows=1)

matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['xtick.minor.width'] = 0

fig, ax = plt.subplots()
heatmap = ax.pcolormesh(per_segments, ybins, countsT, cmap='Blues_r')
ax.set_xscale('log')

ax.xaxis.set_ticks(per_segments)
ax.xaxis.set_ticklabels(per_segments)

ax.yaxis.set_ticks(ybins)
ax.yaxis.set_ticklabels(ybins)

ax.set_title('Fraction of Transits Recovered (EPIC ' + str(targetname) + ')')
ax.set_xlabel('Period (days)')
ax.set_ylabel('Radius (stellar radii)')


cbar = plt.colorbar(heatmap)
heatmap.set_clim(0.0, 1.0)
cbar.set_label('Fraction recovered', rotation=270, labelpad=13)

plt.show()
