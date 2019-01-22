
import csv
import numpy as np
import matplotlib.pyplot as plt
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim



per_segments = np.logspace(-1, 1.415, num=8, endpoint=True, base=10)
rprs_segments = [x / 20 for x in range(8)]

block_percent = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]

count = 0

with open('/Users/sheilasagear/Dropbox/K2/surfacetest/heatmap/heatmaps.txt') as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    for x in content:
        count+=1

    print(content[0])

    for f in content:#for each csv file
        with open('/Users/sheilasagear/Dropbox/K2/heatmap_inclination_csv/' + str(f)) as csvfile:
            reader = csv.reader(csvfile)
            block_lst = list(reader)
            for i in range(len(block_lst)):
                for j in range(len(block_lst[0])):
                    block_lst[i][j] = float(block_lst[i][j])

            print(block_lst)

            for i in range(len(block_percent)):
                for j in range(len(block_percent[0])):
                    block_percent[i][j] = block_percent[i][j]+block_lst[i][j]

for i in range(len(block_percent)):
    for j in range(len(block_percent[0])):
        block_percent[i][j] = block_percent[i][j]/count


print(count)
print(block_percent)

np.savetxt('Pnull_blockpercent_5apr18.csv', block_percent, fmt='%f', delimiter=',')


xlabels = ['0.3', '0.57', '1.1', '2.0', '3.8', '7.3', '13.7', '26.0', ' ', ' ']
ylabels = ['0.5', '1', '1.5', '2', '2.5', '3', '3.5']


plt.imshow(block_percent, cmap = 'Blues_r', extent = [0, 25, 0, .4], aspect='auto')
plt.ylim(.05, .35)
plt.xlim(1, 25)
plt.xticks(np.arange(0, 28.125, 3.125), xlabels)
plt.yticks(np.arange(.05, .4, .05), ylabels)
plt.xlabel('Period (days)')
plt.ylabel('~$\mathregular{R_e}$')
plt.title('Fraction of Transits Recovered')
plt.clim(0,1)
plt.colorbar()

#plt.ticks(a)
#plt.xaxis.ticklabels(labels)

plt.show()
