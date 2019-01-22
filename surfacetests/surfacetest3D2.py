#!/usr/bin/python

from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

import bls, ktransit, math, numpy, pylab, os
import matplotlib.pyplot as plt
#axes = plt.gca()

from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#axe = Axes3D(plt.gcf())

import numpy as np

from matplotlib import cm

#Depending on step size and limits:
#detect_x and undetect_x are a maximum of 95*199 = 18905

detect_count = 0
detect_SDE = []*18905
detect_per = []*18905
detect_rad = []*18905

undetect_count = 0
undetect_SDE = []*18905
undetect_per = []*18905
undetect_rad = []*18905

detect_diff = []*18905
undetect_diff = []*18905

diff = [[0 for x in range(96)] for y in range(200)]

rindex = -1
pindex = -1


rad_diff_values = [[] for y in range(96)]
per_diff_values = [[] for y in range(200)]


rad_range = numpy.arange(.3, .4, .1)
per_range = numpy.arange(5, 6, 1)


##########################
#BEGIN LOOP
##########################


for p in per_range:
    pindex += 1
    rindex = 0
    for r in rad_range:
        rindex += 1

        
        #if not os.path.exists("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r)):
            #os.makedirs("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r))        


##########################
#Import time, flux
##########################


        with open(r'/Users/sheilasagear/Dropbox/K2/surfacetest/hlsp_k2sff_k2_lightcurve_229227248-c06_kepler_v1_llc-default-aper.txt') as f:
            data1 = f.read()
    
#converts txt file to string (and removes heading: 30 characters) 
#in order to isolate each (x,y) coordinate as an element in list 'data'
        datastr = str(data1)
        datastr = datastr[30:]
        data1 = datastr.split('\n')
    
#removes comma after each (x,y) coordinate;
#isolates x and y values as indicies of list 'data'
        index = -1
        while (index < len(data1)):
            tempstring = str(data1[index])
            data1[index] = tempstring.rstrip(',')
            data1[index] = tempstring.split(', ')
            index+=1
            data1.pop
    
        data2 = sum(data1, [])
    
        index = 0
        while (index < len(data2)):
            if index % 2 == 1:
                data2[index] = data2[index].rstrip(',')
            index+=1
            data2.pop()
    
#converts str data points to float
        data_final = [float(i) for i in data2]
    
#defines x and y values by index of 'data'
        time = data_final[0::2]
        flux = data_final[1::2]
    
    
#normalizes flux values to 1.0 (at avg of flux values)
        flux_count = len(flux)
        sumflux = sum(flux)
        flux_avg = sumflux/flux_count
    
        flux = [i/flux_avg for i in flux]
        flux = [i-1 for i in flux]
    
        #xlabel('Time')
        #ylabel('Corrected Flux (normalized to 0)')  

        #print(len(time))
        #print(len(flux))
        
        time.pop()
        
        #plot (time,flux)
        #show()

##########################
#Create ktransit Data: CODE FROM GITHUB
##########################

#M.time = ktransit time stamps. We need to create an array the same size as the K2 data.  

#How many data points do we have for K2?
#The period is p

        num_time = len(time)

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
        period=p, # an orbital period in days
        impact=0.1, # an impact parameter
        rprs=r,   # planet stellar radius ratio  
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
        occ=0.0)    # a secondary eclipse depth in ppm


        M.add_data(time=numpy.array(time[:])),
#the ktransit time array doesn't matter - this is just to create a flux array (tmod) the same size as K2 timestamps

        tmod = M.transitmodel# the out of transit data will be 0.0 unless you specify zpt

        pylab.cla()
        #xlabel('Time')
        #ylabel('Corrected Flux (normalized to 0)')    
        title('KTransit Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
        #graph = plt.plot(M.time,tmod)
        #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_1_2/Period" + str(p) + "Radius" + str(r) + "/ktransitLCPer" + str(p) + "Rad" + str(r) + '.png')

        #pylab.show(graph)

        #print(M.time)
        #print(tmod)


##########################
#Inject ktransit LC into K2 data
##########################


        merged_flux = tmod + flux

        pylab.cla()
        #xlabel('Time')        
        #ylabel('Merged Flux (normalized to 0)')
        #title('EPIC 212820594 Merged Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
        #merged_LC = plt.scatter(time, merged_flux)
        #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_1_2/Period" + str(p) + "Radius" + str(r) + "/merged_fluxPer" + str(p) + "Rad" + str(r) + '.png')

        #plt.show(merged_LC)


##########################
#BLS routine
##########################

        u = [0.0]*len(time)
        v = [0.0]*len(time)

        u = numpy.array(u)
        v = numpy.array(v)

        results = bls.eebls(time, merged_flux, u, v, 1000.0, .3, .001, 200.0, .01, .1)

        #print(results)
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

        print('Planet radius/star radius:' + str(r))
        print('Period from BLS is: ' + str(results[1]))
        print('Set period is: ' + str(p))

##########################
#SR plot
##########################

        SR_array = results[0]

        max_SR = max(SR_array)
        avg_SR = numpy.mean(SR_array)
        sd_SR = numpy.std(SR_array)

#normalize SR_array between 0 and 1
        SR_array = [i/max_SR for i in SR_array]

        #print(max_SR, avg_SR, sd_SR)

        freq = numpy.arange(.3, 1.3, .001, dtype=None)
#numpy.arange(frew start, freq stop (must be calculated), df (freq step)

        pylab.cla()
        #xlabel('Frequency')
        #ylabel('Signal Residue')    
        #title('EPIC 212820594 Merged Signal Residue: Period' + str(p) + ' Radius Ratio ' + str(r))
        #SR_freq_plot = plt.plot(freq, SR_array)
        #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_1_2/Period" + str(p) + "Radius" + str(r) + "/SR_freq_plotPer" + str(p) + "Rad" + str(r) + '.png')

        #pylab.show(SR_freq_plot)

#Best period: 

##########################
#Phase Folding LC
##########################

        start = time[0]
        end = time[len(time)-1]

        time_folded = [i for i in time if i < results[1]]

        phase = [math.fmod(((i-start)/results[1]),1) for i in time]

        pylab.cla()
        #xlabel('Phase')
        #ylabel('Corrected Flux (normalized to 0)')
        #title('EPIC 212820594 Merged Phase Folded Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
        #folded_plt = plt.scatter(phase, merged_flux)
        #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_1_2/Period" + str(p) + "Radius" + str(r) + "/folded_pltPer" + str(p) + "Rad" + str(r) + '.png')
        
        #pylab.show(folded_plt)


##########################
#Calculate SDE (SNR)
########################## 

        SDE = (max_SR-avg_SR)/sd_SR
        
        print('Signal Detection Efficiency: ' + str(SDE))
        if SDE >= 6:
            print('SDE above 6: Transit Detectable')
            detect_count += 1
            detect_per.append(p)
            detect_rad.append(r)
            detect_diff.append(results[1]-p)
            
        else:
            print('SDE below 6: Transit Undetectable')
            undetect_count += 1
            undetect_SDE.append(SDE)
            undetect_per.append(p)
            undetect_rad.append(r)
            undetect_diff.append(results[1]-p)


        diff[pindex][rindex] = (results[1]-p)
        print('Difference in period is '+ str(diff[pindex][rindex]))
        
        rad_diff_values[rindex].append(diff[pindex][rindex])
        
        per_diff_values[pindex].append(diff[pindex][rindex])
            


##########################
#END LOOP
##########################


for i in range(len(detect_SDE)):
    if detect_SDE[i] == 0:
        del detect_SDE[i]

for i in range(len(detect_per)):
    if detect_per[i] == 0:
        del detect_per[i]

for i in range(len(detect_rad)):
    if detect_rad[i] == 0:
        del detect_rad[i]

for i in range(len(undetect_SDE)):
    if undetect_SDE[i] == 0:
        del undetect_SDE[i]

for i in range(len(undetect_per)):
    if undetect_per[i] == 0:
        del undetect_per[i]

for i in range(len(undetect_rad)):
    if undetect_rad[i] == 0:
        del undetect_rad[i]

#print(detect_count, detect_SDE, detect_per, detect_rad)
#print(undetect_count, undetect_SDE, undetect_per, undetect_rad)

"""
pylab.cla()
plt.scatter(detect_per, detect_rad)
axes.set_ylim([0,max(rad_range)])
axes.set_xlim([0,max(per_range)])
plt.xlabel('Orbital Period (days)')
plt.ylabel('Planet Radius/Star Radius')
plt.title('Detectable Transits (SDE > 6)')
plt.show()
"""

"""
pylab.cla()

for x, y in zip(per_range, per_diff_values):
    if x in detect_per:
        plt.scatter([x] * len(y), y, c='r')
    else:
        plt.scatter([x] * len(y), y, c='b')
    #axes.set_ylim([0,5])
    plt.xlabel('Orbital Period (days)')
    plt.ylabel('BLS Period - ktransit Period')
    plt.title('EPIC 229227248 & ktransit Period')


plt.show()
    
pylab.cla()
"""
"""
for x, y in zip(rad_range, per_diff_values):
    plt.scatter([x] * len(y), y, c='b')
    if x in undetect_rad:
        plt.scatter([x] * len(y), y, c='b')
    else:
        plt.scatter([x] * len(y), y, c='r')
    #axes.set_ylim([0,5])
    plt.xlabel('Planet Radius/Star Radius')
    plt.ylabel('BLS Period - ktransit Period')
    plt.title('EPIC 229227248 & ktransit Radius')

plt.show()
"""
"""
    axes.set_ylim([min(per_diff_values),max(per_diff_values)])
    axes.set_xlim([0,max(per_range)])
    plt.xlabel('Orbital Period (days)')
    plt.ylabel('abs(BLS Period - Actual Period)')
    plt.title('BLS Period - Actual Period vs Actual Period')
"""
"""
axe.scatter(detect_per, detect_rad, detect_diff, zdir='z', s=20, c=detect_SDE)
axe.set_xlabel('Orbital Period')
axe.set_ylabel('Planet Radius/Star Radius')
axe.set_zlabel('BLS Period - ktransit Period')
axe.set_title('Detectable Transits')
"""



detect_per = np.asarray(detect_per)
detect_rad = np.asarray(detect_rad)
detect_diff = np.asarray(detect_diff)
detect_SDE = np.asarray(detect_SDE)

# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=plt.figaspect(0.5))

#===============
#  First subplot
#===============
# set up the axes for the first plot
ax = fig.add_subplot(111, projection='3d')

ax.scatter(detect_per, detect_rad, detect_diff, c=detect_SDE)

plt.show()
