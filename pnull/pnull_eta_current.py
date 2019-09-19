#this makes the pnull line plot

import matplotlib.pyplot as plt
import numpy as np

#229227260
#213244700
#229227268

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


deff_1 = 0.044248
#deff_1 = .90
#deff_2 = 1
#deff_3 = 1
eta = np.arange(0, 1, .001)
t_p = 1
n_stars = 374
impact = 1
#i don't know what impact is, so I set it to 1

#the probablility that a planet transits: takes into account
#1. eta, occurrence rate of planets (ie probability of there being a planet)
#2. t_p, transit probability. How likely is it that that planet transits? In our case,
#we already simulated non-transiting planets, so we set this to 1.
p_1 = eta*t_p*impact*deff_1

#p_2 = eta*t_p*deff_2
#print(p_2)

#p_3 = eta*t_p*deff_3

p_null = (1-p_1)**n_stars
#print(p_null)

#####
p=.01
#####

p_close = find_nearest(p_null, p)

p_null_lst = p_null.tolist()

p_ind = p_null_lst.index(p_close)

eta_p = eta[p_ind]


fig, ax = plt.subplots()

ax.plot(eta, p_null)
ax.set_ylim(0, 1)
ax.set_xlabel('$\eta$')
ax.set_ylabel('$\mathregular{P_{null}}$')
ax.set_title('$P_{null}$ vs $\eta$: detection efficiency 5%')

ax.scatter(eta_p, p_close)

ax.text(.60, .88, '$\mathregular{P_{null}}$ = ' + str(p) + '\n$\eta$ = ' + str(round(eta_p, 4)))
ax.text(.60, .83, 'Planets from .9 to 1 $R_{p}/R_{s}$')
ax.text(.60, .78, 'and 1 to 1.58 day periods')

fig.savefig('pnull_eta_1percent_5percentdetectionefficiency.png')
plt.show()
