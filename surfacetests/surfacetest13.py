#!/usr/bin/python

from __future__ import division

#from MCcubed import rednoise

import bls, ktransit, math, pylab, os, csv, untrendy

import numpy as np

from scipy.stats import gamma

from astropy.io import fits

import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = plt.gca()

import pylab
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

from ktransit import FitTransit

from scipy import signal

from numpy import mean, sqrt, square, arange

import PyAstronomy

from PyAstronomy import pyasl

T0 = 0.0
impact = 0.0
rho = 1.5

RMSarr = []

filecount = sum(1 for line in open('all_txt_files.txt'))

for i in range(1):

    
###########################################################################

    
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


        p = np.random.randint(low=1, high=26)
        print('period: ' + str(p))
        r = np.random.choice(np.arange(.1, .9, .1), 1)
        print('rprs: ' + str(r))

        

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

        M.add_data(time=np.array(time[:]))      # integration time of each timestamp

        tmod = M.transitmodel



        merged_flux = tmod + flux

        
        RMS = sqrt(mean(square(flux)))
        RMSarr.append(RMS)
        
        print('RMS = ' + str(RMS))

        
        """
        trend = untrendy.median(time, flux)

        mergedfluxDetrend = np.zeros(len(time))
        length = len(merged_flux)
        
        for x in range(length):
            mergedfluxDetrend[x] = merged_flux[x]-trend[x]

        plt.plot(time, mergedfluxDetrend)
        show()
        """

        #mergedfluxDetrend = np.zeros(len(time))
        
        #for i in range(len(time)):
        #    mergedfluxDetrend[i] - merged_flux[i]-merged_flux[i-1]

        trend = untrendy.median(time, merged_flux)

        mergedfluxDetrend = merged_flux-trend

        mergedfluxDetrend = sigma_clip(mergedfluxDetrend, sigma=3, iters=1)

        mergedfluxDetrend = mergedfluxDetrend.filled(fill_value=0)


        plt.scatter(time, mergedfluxDetrend)
        show()
        
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

        #plt.clf()

    #plt.plot(freq, SR_array)
    #plt.title('EPIC ' + str(targetname))
    #plt.show()


        SDE = (max_SR-avg_SR)/sd_SR

        print('Signal Detection Efficiency: ' + str(SDE))
        if SDE >= 6:
            print('SDE above 6: Transit Detectable')
        else:
            print('SDE below 6: Transit Undetectable')


        print('\n')
        print('\n')

        tempindex+=1
