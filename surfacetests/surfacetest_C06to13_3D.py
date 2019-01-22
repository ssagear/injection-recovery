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

import os

from mpl_toolkits.mplot3d import axes3d, Axes3D


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

z = [0, 0, 0, 0]

z[0] = -3.22185
z[1] = 1.02655
z[2] = 1.6613
z[3] = .225603

p = np.poly1d(z)


fig = plt.figure()
ax = Axes3D(fig)

ax.set_xlabel('Period', fontsize=14)
ax.set_ylabel('~Re', fontsize=14)

ax.set_zlabel('Injected SNR', fontsize=14)
ax.set_zlim(-.0002, .0002)
ax.hold(True)


with open('/Users/sheilasagear/Dropbox/K2/K2_targets/C06to13fits/C06to13lc.txt') as f:
     names = f.readlines()

     

filecount = sum(1 for line in open('/Users/sheilasagear/Dropbox/K2/K2_targets/C06to13fits/C06to13lc.txt'))

for i in range(1):
    print(i)

    
###########################################################################

    
#opens list of file names, creates list of file names
    
    for n in names: #for each LC
    #pulls time, flux from each C14 LC

        time = []
        flux = []
    
        n = n.strip('\n')
        print(n)
        if os.path.isfile(n):
            f = fits.open(n, ignore_missing_end=True)
            

            bestaper = f[1].data

            print(len(bestaper))
            
        
            for i in range(3300):
                time.append(bestaper[i][0])
                flux.append(bestaper[i][1])

    

    #normalize flux values to 1.0
            flux_count = len(flux)
            sumflux = sum(flux)
            flux_avg = sumflux/flux_count
    
            flux = [i/flux_avg for i in flux]

#targets titled by EPIC names

            targetname = n[79:]
            targetname = targetname[:-23]
        #title('EPIC ' + targetname + ' Light Curve')
            print("TARGET: EPIC " + str(targetname))
        
    #campaign_no = n[36:]
    #campaign_no = campaign_no[:-31]
            campaign_no = n[90:]
            campaign_no = campaign_no[:-19]
            print("CAMPAIGN " + str(campaign_no))
        
#normalize to 0.0
            flux = [i-1 for i in flux]

    #print('flux = ' + str(type(flux)))

            flux = np.asarray(flux)
            time = np.asarray(time)

    #print(type(flux))


#SIGMA CLIPPING
            flux = sigma_clip(flux, sigma=4, iters=1)
            flux = flux.filled(fill_value=0)

#uncomment if extra time stamp
    #time.pop()


            """
            fig, ax = plt.subplots(1, 1, figsize=[15,10])
            ax.scatter(time, flux, c='k', s=2)
            ax.set_xlabel('BJD')
            ax.set_ylabel('Flux')
            ax.set_title('EPIC ' + str(targetname))
            """
    
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




##########################
#IMPORT CDPP FROM FITS
##########################
            f = fits.open('/Users/sheilasagear/Dropbox/K2/K2_targets/C06to13fits/hlsp_k2sff_k2_lightcurve_' + str(targetname) + '-c' + str(campaign_no) + '_kepler_v1_llc.fits')

            #print(f.info())
            
            #bestaper = f.data
            #besthead = f.header

            CDPP = f[1].header['QCDPP6']

            SNR = (rprs**2)*np.sqrt(1/period)*(1/CDPP)
            print('SNR = ' + str(SNR))
        
            #QCDPP = besthead['QCDPP6']
            #print('QCDPP = ' + str(QCDPP))
            #QCDPParr.append(QCDPP)



            T0=time[0]+((time[-1]-time[0])/2)
            print('T0 is ' + str(T0))
            

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
        T0=T0,     # a transit mid-time  
        period=period, # an orbital period in days
        impact=0.0, # an impact parameter
        rprs=rprs,   # planet stellar radius ratio  
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
            occ=0.0)    # a secondary eclipse depth in ppm

            M.add_data(time=np.array(time[:]))      # integration time of each timestamp

            tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt
#plt.plot(M.time,tmod)
                


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
            T0=T0, rprs=p(depth))
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


            print('INJECTED PERIOD: ' + str(period))
            print('RECOVERED PERIOD: ' + str(abs(bestFperiod-period)))
            print('INJECTED RPRS: ' + str(rprs))
            print('RECOVERED RPRS: ' + str(abs(bestFrprs-rprs)))
            

        
            if abs(bestFperiod-period) < .1*period and abs(bestFrprs-rprs) < .1*rprs:
                print('YES')
                recovered = True
                recoveredarr.append(True)
                if SNR < 1:
                    ax.scatter(period, rprs, SNR, c='r', s=8)

            else:
                recovered = False
                recoveredarr.append(False)
                if SNR < 1:
                    ax.scatter(period, rprs, SNR, c='b', s=8)





        



plt.savefig('/Users/sheilasagear/Desktop/scatter_test.png')

show()
