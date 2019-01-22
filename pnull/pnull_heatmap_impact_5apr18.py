

import matplotlib.pyplot as plt
import numpy as np
import csv

#229227260
#213244700
#229227268

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]



per_segments = np.logspace(-1, 1.415, num=8, endpoint=True, base=10)
rprs_segments = [x / 20 for x in range(8)]

block_percent = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]


with open('/Users/sheilasagear/Dropbox/K2/surfacetest/Pnull_blockpercent_5apr18.csv') as csvfile:
    reader = csv.reader(csvfile)
    block_lst = list(reader)
    for i in range(len(block_lst)):
        for j in range(len(block_lst[0])):
            block_lst[i][j] = float(block_lst[i][j])

    for i in range(len(block_percent)):
        for j in range(len(block_percent[0])):
            block_percent[i][j] = block_percent[i][j]+block_lst[i][j]




eta = np.arange(0, 1, .001)
t_p = .035
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
        p[i][j] = eta*t_p*block_percent[i][j]
        p_null[i][j] = (1-p[i][j])**n_stars
        p_close[i][j] = find_nearest(p_null[i][j], p_target)
        p_null_lst[i][j] = p_null[i][j].tolist()
        p_ind[i][j] = p_null_lst[i][j].index(p_close[i][j])
        eta_p[i][j] = eta[p_ind[i][j]]


np.savetxt('Pnull_eta_5apr18.csv', eta_p, fmt='%f', delimiter=',')


xlabels = ['0.3', '0.57', '1.1', '2.0', '3.8', '7.3', '13.7', '26.0', ' ', ' ']
ylabels = ['0.5', '1', '1.5', '2', '2.5', '3', '3.5']

print(np.mean(eta_p))

plt.imshow(eta_p, cmap = 'Blues_r', extent = [0, 25, 0, .4], aspect='auto')
plt.ylim(.05, .35)
plt.xlim(1, 25)
plt.xticks(np.arange(0, 28.125, 3.125), xlabels)
plt.yticks(np.arange(.05, .4, .05), ylabels)
plt.xlabel('Period (days)')
plt.ylabel('~$\mathregular{R_e}$')
plt.title('Upper Limit of $\eta$; $\mathregular{P_{null}}$ = ' + str(p_target))
plt.clim(0,1)
plt.colorbar()

plt.tight_layout()

plt.show()








"""
p_1 = eta*t_p*deff_1
print(p_1)

p_2 = eta*t_p*deff_2
#print(p_2)

p_3 = eta*t_p*deff_3

p_null = (1-p_1)**n_stars
#print(p_null)

p=.01

p_close = find_nearest(p_null, p)

p_null_lst = p_null.tolist()

p_ind = p_null_lst.index(p_close)

eta_p = eta[p_ind]


fig, ax = plt.subplots()

ax.plot(eta, p_null)
ax.set_ylim(0, 1)
ax.set_xlabel('$\eta$')
ax.set_ylabel('$\mathregular{P_{null}}$')
ax.set_title('Period 1.1 to 2.4 days; $\mathregular{R_e}$ 3 to 3.5; ' + str(n_stars) + ' stars')

ax.scatter(eta_p, p_close)

ax.text(.77, .15, '$\mathregular{P_{null}}$ = ' + str(p) + '\n$\eta$ = ' + str(eta_p))

plt.show()

"""
