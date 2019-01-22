

import matplotlib.pyplot as plt
import numpy as np

#229227260
#213244700
#229227268

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

deff_1 = 1
deff_2 = 1
deff_3 = 1
eta = np.arange(0, 1, .001)
t_p = .035
n_stars = 500

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
