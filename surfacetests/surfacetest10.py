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




tested_rprs = []
tested_period = []
recovered_rprs = []
recovered_period = []
recoveredarr = []

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


T0 = 3.0
impact = 0.0
rho = 1.5


z = [0, 0, 0, 0]

z[0] = -3.22185
z[1] = 1.02655
z[2] = 1.6613
z[3] = .225603

p = np.poly1d(z)



fig, ax = plt.subplots(1, 1, figsize=[15,10])
ax.set_xlabel('Period')
ax.set_ylabel('RMS')
ax.set_title('Period vs RMS (blue=recovered)')
ax.hold(True)



filecount = sum(1 for line in open('all_txt_files.txt'))

for i in range(1):
    
############################################################################################

    
#opens list of file names, creates list of file names
    with open('all_txt_files.txt') as f:
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
        
        while len(flux) > len(time):
            flux.pop()
        while len(time) > len(flux):
            time.pop()

############################################################################################

    
        period = np.random.uniform(low=1, high=26) #min 3 transits
        tested_period.append(period)
        print('Inj period = ' + str(period))
        
        rprs = np.random.uniform(low=.05, high=.9)
        tested_rprs.append(rprs)
        print('Inj rprs = ' + str(rprs))

        
        M = ktransit.LCModel()
        M.add_star(
            rho=rho, # mean stellar density in cgs units
            ld1=0.2, # ld1--4 are limb darkening coefficients 
            ld2=0.4, # if only ld1 and ld2 are non-zero then a quadratic limb darkening law is used
            ld3=0.0, # if all four parameters are non-zero we use non-linear flavour limb darkening
            ld4=0.0, 
            dil=0.0, # a dilution factor: 0.0 -> transit not diluted, 0.5 -> transit 50% diluted
            zpt=0.0  # a photometric zeropoint, incase the normalisation was wonky
        )
        M.add_planet(
            T0=T0,     # a transit mid-time  
            period=period, # an orbital period in days
            impact=impact, # an impact parameter
            rprs=rprs,   # planet stellar radius ratio  
            ecosw=0.0,  # eccentricity vector
            esinw=0.0,
            occ=0.0)    # a secondary eclipse depth in ppm

        M.add_data(time=np.array(time[:]))      # integration time of each timestamp

        tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt

        #fig.clf()
        #fig, ax = plt.subplots(1, 1, figsize=[15,10])
        #ax.plot(time, tmod, c='k')
        #ax.set_xlabel('BJD')
        #ax.set_ylabel('Injected flux')
        #ax.set_title('Period: ' + str(period) + ' | RPRS: ' + str(rprs))
        #show()

        ###########################################################################


        tmod = np.asarray(tmod)
        time = np.asarray(time)
        flux = np.asarray(flux)

        merged_flux = np.zeros(len(time))
        
        for x in range(len(time)):
            merged_flux[x] = flux[x] + tmod[x]
        
        #merged_flux = tmod + flux


        RMS = sqrt(mean(square(flux)))
        RMSarr.append(RMS)
        print('RMS = ' + str(RMS))

        #fig.clf()
        #fig, ax = plt.subplots(1, 1, figsize=[15,10])
        #ax.scatter(time, merged_flux, c='k', s=2)
        #ax.set_xlabel('BJD')
        #ax.set_ylabel('Merged flux')
        #ax.set_title('Merged Light Curve')
        #show()
        

        ###########################################################################
        


        #trend = untrendy.median(time, merged_flux)

        mergedfluxDetrend = np.zeros(len(time))
        
        #for x in range(len(time)):
        #    mergedfluxDetrend[x] = merged_flux[x]-trend[x]

        #mergedfluxDetrend = sigma_clip(mergedfluxDetrend, sigma=4, iters=1)

        #mergedfluxDetrend = mergedfluxDetrend.filled(fill_value=0)


        for i in range(1, len(time)):
            mergedfluxDetrend[i] = merged_flux[i]-merged_flux[i-1]

        #plt.scatter(time, mergedfluxDetrend)
        #plt.show()


        
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


        print('BLS period: ' + str(results[1]))

        

###################################################################################

        SR_array = results[0]
        max_SR = max(SR_array)
        avg_SR = np.mean(SR_array)
        sd_SR = np.std(SR_array)
#normalize SR_array between 0 and 1
        SR_array = [i/max_SR for i in SR_array]

#np.arange(freq start, freq stop (must be calculated), df (freq step))
    #freq = np.arange(.2, 1.2, .001, dtype=None)

        freq = fmin + np.arange(nf) * df


        plt.plot(freq, SR_array)
        plt.title('EPIC ' + str(targetname))
        plt.show()


###################################################################################

        
        SDE = (max_SR-avg_SR)/sd_SR

        print('Signal Detection Efficiency: ' + str(SDE))
        if SDE >= 6:
            print('SDE above 6: Transit Detectable')
        else:
            print('SDE below 6: Transit Undetectable')



###################################################################################

        phases = PyAstronomy.pyasl.foldAt(time, results[1], getEpoch=False)

        sortPhase = np.argsort(phases)
        phases = phases[sortPhase]
        fluxFolded = mergedfluxDetrend[sortPhase]


###################################################################################

        
        fitT = FitTransit()
        
        fitT.add_guess_star(rho=rho)
        
        fitT.add_guess_planet(
                period=results[1], impact=impact, 
                T0=T0, rprs=.5)#need a guess rprs
        fitT.add_data(time=time, flux=mergedfluxDetrend)
                    

        vary_star = []      # free stellar parameters
        vary_planet = (['period',       # free planetary parameter, 
        'rprs'])  
                
        fitT.free_parameters(vary_star, vary_planet)
        fitT.do_fit()

        #print(fitT.fitresultstellar.items())
        #print(fitT.fitresultplanets.items())

        #bestFstellar = fitT.fitresultstellar.items()
        #bestFrho = bestFstellar[0][1]#Best Fit Rho
                    
        bestFplanet = fitT.fitresultplanets.items()
        bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
        bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs
                    

        fitT.print_results()

        fig.clf()
        fig = ktransit.plot_results(time,merged_flux,fitT.transitmodel)
        fig.show()

        #exit(0)


###################################################################################

        
        if abs(bestFperiod-period) < .03*period:
            recovered = True
            recoveredarr.append(True)
            ax.scatter(period, RMS, c='b', s=4)

        else:
            recovered = False
            recoveredarr.append(False)
            ax.scatter(period, RMS, c='r', s=4)
        
        tempindex += 1



###################################################################################
#END LOOP
###################################################################################

        
        
fig.savefig(r'/Users/sheilasagear/Desktop/per_vs_rms5.png')
