import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import math


per_lst = []
rad_lst = []
imp_lst = []

def randomInc(n):
    randNums = np.random.uniform(low=0.5, high=1, size=n)
    incs = np.arccos(2*randNums - 1)
    return incs

def loguniform(low=float(1), high=float(26)):
    return 10**(np.random.uniform(np.log10(low), np.log10(high)))

def createT():#time array, transit mid-time, impact parameter

    #random period and radius: uniformly weighted
    #you can inject a specific period if you want
    period = loguniform(low=float(1), high=float(26))
    rprs = np.random.uniform(low=.01, high=1)
    rprs = float(rprs)

    Mc = 2e29#kilograms - mass of L dwarf
    Rs = 7.1492e7#meters - radius of L dwarf (1 jupiter radius)
    periodsec = 86400*period
    a1 = sc.G*Mc*periodsec**2
    a2 = (4*sc.pi)**2
    a = np.cbrt(a1/a2)#semimajor axis, meters

    i = randomInc(1)
    impact = (a*math.cos(i))/Rs

    return period, rprs, impact


for i in range(10000):
    incs = randomInc(1)

    period, rprs, impact = createT()
    per_lst.append(period)
    rad_lst.append(rprs)
    imp_lst.append(impact)


perbins = list(np.arange(1, 26.5, 0.5))
fig, ax = plt.subplots()
ax.hist(per_lst, bins=perbins)
ax.set_xlabel('Periods (days)')
ax.set_ylabel('Frequency')
ax.set_title('Frequency of Period - 10,000 Injections')
ax.set_xlim(1, 26)
#plt.show()

radbins = list(np.arange(0, 1.02, 0.02))
fig, ax = plt.subplots()
ax.hist(rad_lst, bins=radbins)
ax.set_xlabel('Planet Radius (stellar radii)')
ax.set_ylabel('Frequency')
ax.set_title('Frequency of Radius - 10,000 Injections')
ax.set_xlim(0.1, 1)
ax.set_ylim(0, 300)
plt.show()

impbins = list(np.arange(0, 2.03, 0.02))
fig, ax = plt.subplots()
ax.hist(imp_lst, bins=impbins)
ax.set_xlabel('Impact parameter $b$')
ax.set_ylabel('Frequency')
ax.set_title('Frequency of Impact Parameter - 10,000 Injections')
ax.set_xlim(0, 2)
ax.set_ylim(0, 18)
plt.axvline(1.5, color='r')
ax.axvspan(1.5, 2, alpha=0.5, color='r')
plt.text(0.69, 16.8, 'Mini-Neptunes out of transit', color='r')
plt.show()
