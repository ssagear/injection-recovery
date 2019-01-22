

import matplotlib.pyplot as plt
import numpy as np

import operator
import functools

#229227260
#213244700
#229227268

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def multiply(numbers):
    numbers = list(numbers)
    total = 1
    for x in numbers:
        total *= x  
    return total  

deff_1 = 1
deff_2 = 1
deff_3 = 1
eta = np.arange(0, 1, .01)
t_p = float(.035)
n_stars = 500

p_lst = []

filecount = sum(1 for line in open('all_txt_files.txt'))

n_stars = filecount

with open('all_txt_files.txt') as f:
    content = f.readlines()

    index = 0
    while (index < filecount):
        file = str(content[index])
        file = file.strip('\n')
        file = file.strip('\r')

        targetname = file[25:]
        targetname = targetname[:-35]
        #title('EPIC ' + targetname + ' Light Curve')
    
        #print("TARGET: EPIC " + str(targetname))
    
        campaign_no = file[37:]
        campaign_no = campaign_no[:-31]
        #print("CAMPAIGN " + str(campaign_no))

        #open csv(file)
        import csv
        row = 1
        col = 3
        with open('/Users/sheilasagear/Dropbox/K2/heatmap_csv/' + '229227254' + '.csv', 'rb') as f:
            csv = csv.reader(f)
            csv_lst = list(csv)
            efficiency = csv_lst[row][col]
            efficiency = float(efficiency)

            print(efficiency)

            p_n = eta*t_p*efficiency

            p_nnull = 1-p_n

            p_lst.append(list(p_nnull)) # a list of lists of p_nnull as a function of eta
            #p_lst.append(list(eta*t_p*efficiency)

        index += 1
print(p_lst[1])

"""
for i in p_lst:
    for j in range(len(p_lst)):
        p_lst[i][j]*p_lst[i+1][j]
"""
#something like this

#[functools.reduce(operator.mul, i) for i in p_lst]

print(composite_prob)

#print(p_lst)



"""        
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

#ax.scatter(eta_p, p_close)

#ax.text(.77, .15, '$\mathregular{P_{null}}$ = ' + str(p) + '\n$\eta$ = ' + str(eta_p))

plt.show()
"""
