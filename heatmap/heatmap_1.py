#!/usr/bin/python

from __future__ import division

import bls, ktransit, math, pylab, os, csv, untrendy, PyAstronomy

import numpy as np

from scipy.stats import gamma

from astropy.io import fits

import matplotlib.pyplot as plt

from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

from ktransit import FitTransit

from scipy import signal

from numpy import mean, sqrt, square, arange

from PyAstronomy import pyasl

#import plotly.plotly as py

#import plotly.graph_objs as go

tested_rprs = []
tested_period = []
recovered_rprs = []
recovered_period = []
recoveredarr = []
unrecoveredarr = []
unrecovered_period = []
unrecovered_rprs = []

flux = []
merged_flux = []
time = []
campaign_no = []
targetname = []

tmod = []
mergedfluxDetrend = []
SDE = []

BLSper = []
RMSarr = []

z = [0, 0, 0, 0]

z[0] = -3.22185
z[1] = 1.02655
z[2] = 1.6613
z[3] = .225603

p = np.poly1d(z)


fig, ax = plt.subplots(1, 1, figsize=[15,10])
ax.set_xlabel('$\mathregular{R_p}$/$\mathregular{R_s}$', fontsize=14)
ax.set_xlim(0, 0.4)
ax.set_ylim(0, 26)

ax.set_ylabel('Orbital Period (days)', fontsize=14)
ax.set_title('EPIC 229227260')
ax.hold(True)


iterations = 10


filecount = 1 #sum(1 for line in open('all_txt_files.txt'))

for i in range(iterations):
    print(i)

    
###########################################################################

    
#opens list of file names, creates list of file names
    with open('all_txt_files1.txt') as f:
        content = f.readlines()

#converts each file name to string, removes \n from end, opens each object's txt file (loop)
#and stores (x,y) in list 'data'
    tempindex = 0
    while (tempindex < filecount):
        file = str(content[tempindex])
        file = file.strip('\n')
        file = file.strip('\r')

        with open(file) as f:
            data1 = f.read()

#converts txt file to string (and removes heading: 30 characters) 
#in order to isolate each (x,y) coordinate as an element in list 'data'
        datastr = str(data1)
        datastr = datastr[30:]
        data1 = datastr.split('\n')


#removes comma after each (x,y) coordinate;
#islates x and y values as indicies of list 'data'
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
    
#convers str data points to float
        data_final = [float(a) for a in data2]

#defines x and y values by index of 'data'
        time = data_final[0::2]
        flux = data_final[1::2]
        
#normalizes flux values to 1.0 (at avg of flux values)
        flux_count = len(flux)
        sumflux = sum(flux)
        flux_avg = sumflux/flux_count
        
        flux = [i/flux_avg for i in flux]

#targets titled by EPIC names
        targetname = file[25:]
        targetname = targetname[:-35]
    #title('EPIC ' + targetname + ' Light Curve')
    
        print("TARGET: EPIC " + str(targetname))
    
        campaign_no = file[37:]
        campaign_no = campaign_no[:-31]
        print("CAMPAIGN " + str(campaign_no))
    
#normalize to 0.0
        flux = [i-1 for i in flux]


#SIGMA CLIPPING
        flux = sigma_clip(flux, sigma=2, iters=1)


        flux = flux.filled(fill_value=0)
    
#uncomment if extra
        
        #plt.plot(time, flux)
        #show()
        
        period = np.random.uniform(low=1, high=26)
#period = 8.0
        tested_period.append(period)
        print('Inj period = ' + str(period))
        
        #rprs = np.random.choice(np.arange(.1, .9, .1), 1)
        rprs = np.random.uniform(low=.01, high=.4)
#rprs = .6
        rprs = float(rprs)
        tested_rprs.append(rprs)
        print('Inj rprs = ' + str(rprs))

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
        T0=2425,     # a transit mid-time  
        period=period, # an orbital period in days
        impact=0.0, # an impact parameter
        rprs=rprs,   # planet stellar radius ratio  
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
        occ=0.0)    # a secondary eclipse depth in ppm

        M.add_data(time=np.array(time[:]))      # integration time of each timestamp

        tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt
#plt.plot(M.time,tmod)



        RMS = sqrt(mean(square(flux)))

        print('RMS: ' + str(RMS))
                


        merged_flux = flux + tmod
        
        #plt.plot(time, merged_flux)
        #show()

        trend = untrendy.median(time, merged_flux)

        mergedfluxDetrend = np.zeros(len(time))

        for i in range(len(time)):
            mergedfluxDetrend[i] = merged_flux[i]-trend[i]

        #plt.plot(time, mergedfluxDetrend)
        #show()

###########################################################################

        u = [0.0]*len(time)
        v = [0.0]*len(time)
        u = np.array(u)
        v = np.array(v)

    #time, flux, u, v, number of freq bins (nf), min freq to test (fmin), freq spacing (df), number of bins (nb), min transit dur (qmi), max transit dur (qma)

        nf = 1000.0
        fmin = .035
        df = 0.001
        nbins = 300
        qmi = 0.001
        qma = 0.3

        results = bls.eebls(time, mergedfluxDetrend, u, v, nf, fmin, df, nbins, qmi, qma)
            
#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6


        print('DEPTH: ' + str(results[3]))
    
        print('BLS period: ' + str(results[1]))


###################################################################################

        SR_array = results[0]
        max_SR = max(SR_array)
        avg_SR = np.mean(SR_array)
        sd_SR = np.std(SR_array)
#normaze SR_array between 0 and 1
        SR_array = [i/max_SR for i in SR_array]

#np.arange(freq start, freq stop (must be calculated), df (freq step))
    #freq = np.arange(.2, 1.2, .001, dtype=None)

        freq = fmin + np.arange(nf) * df

        #plt.plot(freq, SR_array)
        #show()


###################################################################################

        
        SDE = (max_SR-avg_SR)/sd_SR

        print('Signal Detection Efficiency: ' + str(SDE))
        if SDE >= 6:
            print('SDE above 6: Transit Detectable')
        else:
            print('SDE below 6: Transit Undetectable')




        high = results[3]*results[4]
        low = high - results[3]
        
        fit = np.zeros(nbins) + high # H
        fit[results[5]:results[6]+1] = low # L

        depth = high - low

        print('Guess rprs: ' + str(p(depth)))


        fitT = FitTransit()
        fitT.add_guess_star(rho=1.5)    
        fitT.add_guess_planet(
        period=results[1], impact=0.0, 
        T0=2425, rprs=p(depth))
        fitT.add_data(time=time, flux=mergedfluxDetrend)

        vary_star = []      # free stellar parameters
        vary_planet = (['period',       # free planetary parameters
            'rprs'])                # free planet parameters are the same for every planet you model

        fitT.free_parameters(vary_star, vary_planet)
        fitT.do_fit()                   # run the fitting

        bestFplanet = fitT.fitresultplanets.items()
        bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
        bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs
        
        fitT.print_results()    


        #fig = ktransit.plot_results(time,mergedfluxDetrend,fitT.transitmodel)
        #show()

        print('\n')
        print('\n')



        
        if abs(bestFperiod-period) < .03*period and abs(bestFrprs-rprs) < .03*rprs:
            print('YES')
            recovered = True
            recoveredarr.append(True)
            recovered_period.append(period)
            recovered_rprs.append(rprs)
            #ax.scatter(rprs, period, c='b', s=8)

        else:
            recovered = False
            recoveredarr.append(False)
            unrecovered_period.append(period)
            unrecovered_rprs.append(rprs)
            #ax.scatter(rprs, RMS, c='b', s=8)





        
        tempindex+=1

period_segments = []
rprs_segments = []

period_segments = [0, 5]#, 10, 15, 20, 25]
rprs_segments = [0, 0.1]#, 0.2, 0.3, 0.4]

block_count = [[False for i in period_segments] for j in rprs_segments]
count = []
block_percent = [[0 for i in period_segments] for j in rprs_segments]

num = 0
den = 0



count = []


for i in range(iterations):
    for p in range(len(period_segments) - 1):
        for r in range(len(rprs_segments) - 1):
            total = 0
            if tested_period[i] > period_segments[p] and tested_period[i] < period_segments[p+1]:
                if tested_rprs[i] > rprs_segments[r] and tested_rprs[i] < rprs_segments[r+1]:
                    total += 1
                    count = []
                    if tested_rprs[i] in recovered_rprs and tested_period[i] in recovered_period:
                        count.append(True)
                    else:
                        count.append(False)
            num = 0
            for b in range(len(count)):
                if count[b] == True:
                    num += 1

            if total != 0:
                block_percent[p][r] = num/total
            else:
                block_percent[p][r] = 0

            print(total)


print(block_percent)
                    
                

"""
plt.savefig('/Users/sheilasagear/Desktop/EPIC229227260-1-heatmap.png')

show()
"""
