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

#from sklearn.metrics import mean_squared_error


bestFperarr = []
bestFrprsarr = []
SDEarr = []
diffarr = []
targetnamearr = []

tested_rprs = []
tested_period = []
recovered_rprs = []
recovered_period = []

rprsarr = []
lowarr = []

recoveredarr = []
SNR = []

time = []
flux = []

SNRforLC = []

QCDPParr = []

RMSarr = []

deptharr = []


flux = []
merged_flux = []
time = []
campaign_no = []
targetname = []

tmod = []
mergedfluxDetrend = []
SDE = []

BLSper = []


T0 = 1.0
impact = 0.0 #np.random.uniform(low=0.0, high=1.0)
rho = 1.5



#fig, ax = plt.subplots(1, 1, figsize=[15,10])

#ax.hold(True)

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
        print(i)
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



        #fig, ax = plt.subplots(1, 1, figsize=[15,10])
        #ax.scatter(time, flux, c='k', s=2)
        #ax.set_xlabel('BJD')
        #ax.set_ylabel('Normalized flux')
        #ax.set_title('EPIC ' + str(targetname))
        #show()
        

###########################################################################


##########################
#IMPORT CDPP FROM FITS
##########################
        #f = fits.open('/Users/sheilasagear/OneDrive/K2_Research/Cycle2_FITS/CYCLE2FITS/hlsp_k2sff_k2_lightcurve_' + str(targetname) + '-c0' + str(campaign_no) + '_kepler_v1_llc.fits')

        #bestaper = f.data
        #besthead = f.header
        
        #QCDPP = besthead['QCDPP6']
        #print('QCDPP = ' + str(QCDPP))
        #QCDPParr.append(QCDPP)



    
###########################################################################


        #np.random.uniform(low=.01, high=1)
    
        period = np.random.uniform(low=1, high=26)
        #period = 12.0
        tested_period.append(period)

        print('period = ' + str(period))
        
        rprs = np.random.uniform(low=.05, high=1.0)
        #rprs = 0.5
        tested_rprs.append(rprs)

        
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

        
        if len(tmod) != len(flux):
            flux = flux[0:len(flux)-1]
        
        merged_flux = tmod + flux


        RMS = sqrt(mean(square(flux)))
        #RMS = sqrt(mean_squared_error(merged_flux, tmod))
        RMSarr.append(RMS)
        
        print('RMS = ' + str(RMS))

        #fig.clf()
        #fig, ax = plt.subplots(1, 1, figsize=[15,10])
        #ax.scatter(time, merged_flux, c='k', s=2)
        #ax.set_xlabel('BJD')
        #ax.set_ylabel('Merged flux')
        #ax.set_title('Merged Light Curve')
        #show()
        

        #exit(0)

        
        
        """
        mergedfluxDetrend = list()
                
        for i in range(len(trend)):
            val = merged_flux[i]-merged_flux[i-1]
            mergedfluxDetrend.append(val)
        """

        trend = untrendy.median(time, flux)

        mergedfluxDetrend = np.zeros(len(time))
        length = len(merged_flux)
        
        for x in range(length):
            mergedfluxDetrend[x] = merged_flux[x]-trend[x]

        mergedfluxDetrend = sigma_clip(mergedfluxDetrend, sigma=4, iters=1)

        mergedfluxDetrend = mergedfluxDetrend.filled(fill_value=0)



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

        BLSper.append(results[1])


        

###########################################################################


        
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


###########################################################################


##########################
#Folding and Binning
##########################

    #pylab.cla()
        """
        f0 = 1.0/results[1]
        n = len(time)
        ibi = np.zeros(nbins)
        y = np.zeros(nbins)
        phase = np.linspace(0.0, 1.0, nbins)
        for i in range(n):
            ph = u[i]*f0
            ph = ph-int(ph)
            j = int(nbins*ph)
            ibi[j] = ibi[j] + 1.0
            y[j] = y[j] + v[i]
    #binned and folded plot
    #pylab.cla()
    #plt.scatter(phase, y/ibi, s=3)
    #plt.title("EPIC " + str(targetname) + " folded and binned LC")
    #pylab.show()
        """
    

###########################################################################

            
###########################
#BLS Overlay
###########################

        
        high = results[3]*results[4]
        low = high - results[3]
        
        fit = np.zeros(nbins) + high # H
        fit[results[5]:results[6]+1] = low # L

        depth = high - low
        deptharr.append(depth)


    #plt.plot(phase, fit)
    #plt.xlabel(r"Phase")
    #plt.ylabel(r"Mean value of flux")
    #plt.title("SDE " + str(SDE) + "; BLS period " + str(results[1]))
    #plt.ylim(-.1, .1)
    #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r) + "/folded_pltPer" + str(p) + "Rad" + str(r) + 'BLSoverlay.png')



###########################################################################

        tempindex += 1
        
###########################################################################


#z = np.polyfit(deptharr, tested_rprs, 3)
#p = np.poly1d(z)


z = [0, 0, 0, 0]

z[0] = -3.22185
z[1] = 1.02655
z[2] = 1.6613
z[3] = .225603

p = np.poly1d(z)

#fig.clf()

fig, ax = plt.subplots(1, 1, figsize=[15,10])
ax.set_xlabel('Period')
ax.set_ylabel('RMS')
ax.set_title('Period vs RMS (blue=recovered)')
ax.hold(True)


"""
deptharr = sigma_clip(deptharr, sigma=4, iters=1)

x = np.linspace(-.4, .7, 100)

bbox_props = dict(boxstyle="square,pad=0.3", facecolor='none', edgecolor='black')


ax.scatter(deptharr, tested_rprs, color='k', s=2)
ax.plot(x, p(x), color = 'red')
ax.set_ylim(-.1, 1.1)
ax.set_xlim(-0.5, 0.8)
ax.set_ylabel('RPRS')
ax.set_xlabel('BLS depth')
ax.set_title('p(x) = ' + str(z[0]) + '$x^3$  + ' + str(z[1]) + '$x^2$ +' + str(z[2]) + '$x$ + ' + str(z[3]))

fig.savefig(r'/Users/sheilasagear/Desktop/rprs_vs_depth2.png')


exit(0)
"""

#print(tested_period)
#print(len(tested_period))

#exit(0)
        
###########################################################################


#fig.clf()
#fig, ax = plt.subplots(1, 1, figsize=[20,15])
#ax.set_ylabel('RMS')
#ax.set_xlabel('Period')
#ax.set_title('Per vs RMS (blue=recovered)')
#ax.hold(True)



depthi = 0
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

        #time = [float(x) for x in time]
        #flux = [float(x) for x in flux]


#SIGMA CLIPPING
        flux = sigma_clip(flux, sigma=2, iters=1)
        
    #uncomment if extra time stamp
        while len(flux) > len(time):
            flux.pop()
        while len(time) > len(flux):
            time.pop()


        time = np.asarray(time)
        flux = np.asarray(flux)

        #print(len(flux))
        #print(len(time))
        #exit(0)





###########################################################################


        T0 = 1.0 #np.random.uniform(low=.01, high=1)
    
        period = tested_period[depthi]
        print('Set period = ' + str(period))
        
        impact = 0.0 #np.random.uniform(low=0.0, high=1.0)
        
        rprs = tested_rprs[depthi]
        print('Set rprs = ' + str(rprs))
        
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


        

        
###########################################################################


        merged_flux = tmod + flux

    
        #trend = untrendy.median(time, merged_flux)

        trend = untrendy.median(time, flux)

        mergedfluxDetrend = np.zeros(len(time))
        length = len(merged_flux)
        
        for x in range(length):
            mergedfluxDetrend[x] = merged_flux[x]-trend[x]

        mergedfluxDetrend = sigma_clip(mergedfluxDetrend, sigma=4, iters=1)

        mergedfluxDetrend = mergedfluxDetrend.filled(fill_value=0)


        #even out len
        #while len(mergedfluxDetrend) > len(time):
        #    mergedfluxDetrend.pop()
        #while len(time) > len(mergedfluxDetrend):
        #    time.pop()

###

        phases = PyAstronomy.pyasl.foldAt(time, results[1], getEpoch=False)

        sortPhase = np.argsort(phases)
        phases = phases[sortPhase]
        fluxFolded = merged_flux[sortPhase]

###
        
        
##########################
#Detrending Flux and Levenberg-Marquardt Fitting:
#if SDE>6 
##########################

        print('Set rprs ' + str(tested_rprs[depthi]))
        print('Loop rprs ' + str(rprs))
        print('rprs guess ' + str(p(deptharr[depthi])))

        print('Set period ' + str(tested_period[depthi]))
        print('Loop period ' + str(period))
        print('period guess ' + str(BLSper[depthi]))
            
        fitT = FitTransit()
        fitT.add_guess_star(rho=rho)    
        fitT.add_guess_planet(
            period=BLSper[depthi], impact=impact, 
            T0=T0, rprs=(.5))#need a guess rprs
        fitT.add_data(time=time, flux=mergedfluxDetrend)
                
                    
        vary_star = ['zpt']      # not sure how to avoid free stellar parameters? ideally would not vary star at all
        vary_planet = (['period', 'rprs'])
                    
        fitT.free_parameters(vary_star, vary_planet)
        fitT.do_fit()

    #print(fitT.fitresultstellar.items())
    #print(fitT.fitresultplanets.items())
                    
    #bestFstellar = fitT.fitresultstellar.items()
    #bestFrho = bestFstellar[0][1]#Best Fit Rho
            
        bestFplanet = fitT.fitresultplanets.items()
        bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
        bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs
    #bestFimpact = bestFplanet[0][1]['impact']
    #bestFT0 = bestFplanet[0][1]['T0']
                    
        fitT.print_results()
        #fig = ktransit.plot_results(time,merged_flux,fitT.transitmodel)
        #fig.show()
    



        q = untrendy.median(time, flux)

        """
        
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(11, 5))
#fig.subplots_adjust(hspace=.01)

        bbox_props = dict(boxstyle="square,pad=0.4", facecolor='none', edgecolor='black')

        ax1.scatter(time, flux, color='k', s=2)
        ax1.plot(time, q, color='mediumaquamarine')
        ax1.set_xlabel('Time (days)')
        ax1.set_ylabel('Normalized Flux')
        ax1.set_title('EPIC ' + str(targetname))


        ax2.plot(freq, SR_array, color='black')
        ax2.set_ylabel('Power')
        ax2.set_xlabel('Frequency')
        ax2.set_title('Box Least-Squares')
        ax2.set_ylim(0, 1.3)
        ax2.text(0.2, 1.05, 'SDE: ' + str(SDE) + '\n' + 'BLS Period: ' + str(results[1]), bbox=bbox_props, fontsize=7)
#fitting?????


        ax3.scatter(phases, fluxFolded, color='k', s=2)
        ax3.set_xlabel('Phase')
        ax3.set_ylabel('Normalized Flux')
        ax3.set_title('Phase Folded: Period ' + str(results[1]))

        
        #ax4.scatter(time, mergedfluxDetrend, color='k', s=2)
        #ax4.plot(time, fitT.transitmodel, color='mediumaquamarine')
        #ax4.set_ylim(-0.25,0.25)

        #ax4.set_xlabel('Time (days)')
        #ax4.set_ylabel('Detrended Flux')
        #ax4.set_title('Levenberg-Marquardt')
        #ax4.annotate('Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), xy = (time[1], 0.04), bbox = bbox_props)
        #ax4.text(time[1], 0.12, 'Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), bbox=bbox_props, fontsize=7)
        """


        #plt.tight_layout()
        #show()

        #f.savefig(r'/Users/sheilasagear/Dropbox/ssagear_k2/cycle2/injections2/EPIC' + str(targetname) + 'period' + str(period) + 'rprs' + str(rprs) + 'png')


        #exit(0)

        
###########################################################################
        
        if abs(bestFperiod-tested_period[depthi]) < .03*tested_period[depthi]:
            recovered = True
            recoveredarr.append(True)
            ax.scatter(tested_period[depthi], RMSarr[depthi], c='b', s=4)

        else:
            recovered = False
            recoveredarr.append(False)
            ax.scatter(tested_period[depthi], RMSarr[depthi], c='r', s=4)

        print(str(depthi))
        
        depthi += 1
        tempindex += 1



###########################################################################


#plt.show()

fig.savefig(r'/Users/sheilasagear/Desktop/per_vs_rms1.png')

#pylab.clf()

#plt.scatter(RMSarr, QCDPParr, c='k')
#plt.show()


