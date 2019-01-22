import ktransit
import matplotlib.pyplot as plt
import numpy as np

import bls

import pylab

M = ktransit.LCModel()
M.add_star(
        rho=1.5, # mean stellar density in cgs units
        ld1=0.2, # ld1--4 are limb darkening coefficients 
        ld2=0.4, # if only ld1 and ld2 are non-zero then a quadratic limb darkening law is used
        ld3=0.0, # if all four parameters are non-zero we use non-linear flavour limb darkening
        ld4=0.0, 
        dil=0.0, # a dilution factor: 0.0 -> transit not diluted, 0.5 -> transit 50% diluted
        zpt=0.0  # a photometric zeropoint, incase the normalisation was wonky
        )
M.add_planet(
        T0=1.0,     # a transit mid-time  
        period=1.0, # an orbital period in days
        impact=0.1, # an impact parameter
        rprs=0.1,   # planet stellar radius ratio  
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
        occ=0.0)    # a secondary eclipse depth in ppm

M.add_planet() # you can add as many planets as you like (up to 10)

M.add_data(
        time=np.arange(0,10,0.0188),                                 # timestamps to evaluate the model on
        itime=np.zeros_like(np.arange(0,10,0.0188))+0.0188 )      # integration time of each timestamp

tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt
plt.plot(M.time,tmod)


u = [0.0]*len(M.time)
v = [0.0]*len(M.time)

u = np.array(u)
v = np.array(v)

nf = 1000.0
fmin = 0.2
df = 0.001
nb = 200
qmi = 0.001
qma = 0.3
    
results = bls.eebls(M.time, tmod, u, v, nf, fmin, df, nb, qmi, qma)
            
#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

print('BLS period: ' + str(results[1]))
cd ..
SR_array = results[0]
max_SR = max(SR_array)
avg_SR = np.mean(SR_array)
sd_SR = np.std(SR_array)

#normalize SR_array between 0 and 1
SR_array = [i/max_SR for i in SR_array]

    #print(max_SR, avg_SR, sd_SR)

#np.arange(freq start, freq stop (must be calculated), df (freq step))
    #freq = np.arange(.2, 1.2, .001, dtype=None)

freq = fmin + np.arange(nf) * df

plt.clf()

plt.plot(freq, SR_array)
plt.show()
