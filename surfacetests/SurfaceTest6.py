
from __future__ import division

import bls, ktransit, math, pylab, os, batman

from scipy.stats import gamma

import untrendy

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
fig = plt.figure()
axes = plt.gca()

from mpl_toolkits.mplot3d import Axes3D
#axe = Axes3D(plt.gcf())
#ax = fig.add_subplot(111, projection='3d')

from matplotlib import pyplot
import pylab
from pylab import xlabel, ylabel, title, plot, show, ylim, xlim

from astropy.stats import sigma_clip

#fige = pylab.figure()
#ax = Axes3D(fig)

from ktransit import FitTransit

from scipy import signal

import obspy
from obspy.signal.detrend import polynomial

perarr = []
radarr = []
bestFperarr = []
bestFrprsarr = []
SDEarr = []
diffarr = []
targetnamearr = []

tested_rprs = []
tested_period = []


rprsarr = []
lowarr = []




##########################
#INITIALIZATION
##########################

#Depending on step size and limits:
#detect_x and undetect_x are a maximum of 95*199 = 18905
#Radius ratio from .05 to 1 with a step size of .01
#Period from .1 to 20 days with a step size of .1

detect_count = 0
detect_SDE = []*18905
detect_per = []*18905
detect_rad = []*18905
detect_BLS_per = []*18905

undetect_count = 0
undetect_SDE = []*18905
undetect_per = []*18905
undetect_rad = []*18905

detect_diff = []*18905
undetect_diff = []*18905

diff = [[0 for x in range(96)] for y in range(200)]

maskindex = []


rad_diff_values = [[] for y in range(96)]
per_diff_values = [[] for y in range(200)]

###########################################
#These limits and step size can be changed!
rad_range = np.arange(.05, 1, .08)
per_range = np.arange(.5, 25, .3)
###########################################


total_rad = len(rad_range)
total_per = len(per_range)
total = total_rad*total_per

recovered = []
SNR = []

##########################
#BEGIN LOOP
##########################

#LATER: loop through C6 and C7 light curves; open and inject range of transits into each 

#count number of files being graphed
filecount = sum(1 for line in open('all_txt_files.txt'))

#opens list of file names, creates list of file names
with open('all_txt_files.txt') as f:
    content = f.readlines()

    SNRforLC = []
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
    data_final = [float(i) for i in data2]

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
        
#uncomment if extra time stamp
    time.pop()
                    
##########################
#IMPORT CDPP FROM FITS
##########################
    #f = fits.open('/Users/sheilasagear/OneDrive/K2_Research/Cycle2_FITS/CYCLE2FITS/hlsp_k2sff_k2_lightcurve_' + str(targetname) + '-c0' + str(campaign_no) + '_kepler_v1_llc.fits')

    #bestaper = f[1].data
    #besthead = f[1].header
        
    #CDPP = besthead['QCDPP6']
    #print(CDPP)

    rindex = -1
    pindex = -1


        
##########################
#Create ktransit Data: CODE FROM GITHUB
##########################

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

    T0 = 1.0 #np.random.uniform(low=.01, high=1)
    period = np.random.uniform(low=.05, high=100)
    tested_period.append(period)
    impact = 0.0 #np.random.uniform(low=0.0, high=1.0)
    rprs = np.random.uniform(low=.05, high=1.0)
    tested_rprs.append(rprs)

    print('Set period is ' + str(period))
    print('Set rprs is ' + str(rprs))

    rprsarr.append(rprs)
            
            
    M.add_planet(
                T0=T0,     # a transit mid-time  
            period=period, # an orbital period in days
            impact=impact, # an impact parameter
            rprs=rprs,   # planet stellar radius ratio  
            ecosw=0.0,  # eccentricity vector
            esinw=0.0,
    occ=0.0)    # a secondary eclipse depth in ppm
            

    M.add_data(time=np.array(time[:]))

    tmod = M.transitmodel# the out of transit data will be 0.0 unless you specify zpt
    #pylab.cla()

    #plot(time, tmod)
    #title('KTRANSIT injected LC')
    #show()
        
######################################
#UNCOMMENT TO SAVE KTRANSIT LIGHT CURVE
######################################

    """
            xlabel('Time')
                ylabel('Corrected Flux (normalized to 0)')    
                title('KTransit Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
                graph = plt.plot(M.time,tmod)
                pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r) + "/ktransitLCPer" + str(p) + "Rad" + str(r) + '.png')
    """

    #pylab.show(graph)


##########################
#Inject ktransit LC into K2 data
##########################

    if len(tmod) != len(flux):
        flux = flux[0:len(flux)-1]
    merged_flux = tmod# + flux
        #plus flux

    #pylab.cla()
                
######################################
#UNCOMMENT TO SAVE MERGED LIGHT CURVE
######################################

    """
            xlabel('Time')        
            ylabel('Merged Flux (normalized to 0)')
            title('EPIC 212820594 Merged Light Curve: Period' + str(p) + ' Radius Ratio ' + str(r))
            merged_LC = plt.scatter(time, merged_flux, s=3)
            pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r) + "/merged_fluxPer" + st(p) + "Rad" + str(r) + '.png')
    """
    #plt.show(merged_LC)


##########################
#BLS routine
##########################

    mergedfluxDetrend = tmod
    
    ###mergedfluxDetrend, ferr = untrendy.untrend(time, merged_flux)

    #trend, ferr = untrendy.untrend(time, flux)

    #mergedfluxDetrend = list()
                
    #for i in range(len(trend)):
    #    val = merged_flux[i]-merged_flux[i-1]
    #    mergedfluxDetrend.append(val)

    u = [0.0]*len(time)
    v = [0.0]*len(time)

    u = np.array(u)
    v = np.array(v)

    #time, flux, u, v, number of freq bins (nf), min freq to test (fmin), freq spacing (df), number of bins (nb), min transit dur (qmi), max transit dur (qma)

    nf = 1000.0
    fmin = 0.2
    df = 0.001
    nbins = 200
    qmi = 0.001
    qma = 0.3
    
    results = bls.eebls(time, flux, u, v, nf, fmin, df, nbins, qmi, qma)
            
#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

    print('BLS period: ' + str(results[1]))
    
##########################
#SR plot
##########################

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

    #pylab.cla()

#UNCOMMENT TO SAVE SIGNAL RESIDUE PLOT

    """
            xlabel('Frequency')
            ylabel('Signal Residue')    
            title('EPIC 212820594 Merged Signal Residue: Period' + str(p) + ' Radius Ratio ' + str(r))
            SR_freq_plot = plt.plot(freq, SR_array)
            pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/LC_5/Period" + str(p) + "Radius" + str(r) + "/SR_freq_plotPer" + str(p) + "Rad" + str(r) + '.png')
    """
        
    #pylab.show(SR_freq_plot)

        
##########################
#Calculate SDE
########################## 

    SDE = (max_SR-avg_SR)/sd_SR
    
    print('Signal Detection Efficiency: ' + str(SDE))
    if SDE >= 6:
        print('SDE above 6: Transit Detectable')
                #detect_count += 1
                #detect_SDE.append(SDE)
                #detect_per.append(p)
                #detect_rad.append(r)
                #detect_BLS_per.append(results[1])
        SDEarr.append(True)
                
                #detect_diff.append(results[1]-p)
    else:
        print('SDE below 6: Transit Undetectable')
                #undetect_count += 1
                #undetect_SDE.append(SDE)
                #undetect_per.append(p)
                #undetect_rad.append(r)
                #undetect_diff.append(results[1]-p)
        SDEarr.append(False)

            
            #rad_diff_values[rindex].append(diff[pindex][rindex])
        
            #per_diff_values[pindex].append(diff[pindex][rindex])
        

########
#data centering?
########



    

##########################
#Folding and Binning
##########################

    #pylab.cla()

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
    

###########################
#BLS Overlay
###########################

        
    high = results[3]*results[4]
    low = high - results[3]
        
    fit = np.zeros(nbins) + high # H
    fit[results[5]:results[6]+1] = low # L
    
    #plt.plot(phase, fit)
    #plt.xlabel(r"Phase")
    #plt.ylabel(r"Mean value of flux")
    #plt.title("SDE " + str(SDE) + "; BLS period " + str(results[1]))
    #plt.ylim(-.1, .1)
    #pylab.savefig("/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/Pipeline/EPIC229227254-2sigclip/Period" + str(p) + "Radius" + str(r) + "/folded_pltPer" + str(p) + "Rad" + str(r) + 'BLSoverlay.png')


    #pylab.show()


    print('low = ' + str(low))
    lowarr.append(low)


##########################
#Detrending Flux and Levenberg-Marquardt Fitting:
#if SDE>6 
##########################

#start values are correct values

    if True: # SDE >= 6:

        #polynomial detrending
        #mergedfluxDetrend = polynomial(merged_flux, order=5)

    #untrendy
        """
        mergedfluxDetrend, ferr = untrendy.untrend(time, merged_flux)

        trend, ferr = untrendy.untrend(time, flux)

        mergedfluxDetrend = list()
                
        for i in range(len(trend)):
            val = merged_flux[i]-merged_flux[i-1]
            mergedfluxDetrend.append(val)
        """


                #plot(time, mergedfluxDetrend)
                #show()
        

        fitT = FitTransit()
        fitT.add_guess_star(rho=1.5)    
        fitT.add_guess_planet(
            period=results[1], impact=0.0, 
            T0=1.0, rprs=0.5)#need a guess rprs
        fitT.add_data(time=time, flux=mergedfluxDetrend)
                
                    
        vary_star = ['rho']      # not sure how to avoid free stellar parameters? ideally would not vary star at all
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
                    
                #save figure
                #fig = ktransit.plot_results(time,mergedfluxDetrend,fitT.transitmodel)
                #fig.savefig('/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/ktransitfits_untrendy/' + str(targetname) + 'fitPer' + str(p) + 'Rprs' + str(r) + '.png')



        #if not os.path.exists("/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/ktransitfits_untrendy1/EPIC" + str(targetname)):
        #    os.makedirs("/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/ktransitfits_untrendy1/EPIC" + str(targetname))

        """
        pylab.cla()
        fig, ax = plt.subplots(1, 1, figsize=[15,10])
        ax.scatter(time, mergedfluxDetrend, color='k', s=2)
        ax.plot(time, fitT.transitmodel, color='mediumaquamarine')
        ax.set_ylim(-0.3,0.3)
        ax.set_xlabel('Time (BJD - 2450000)')
        ax.set_title('EPIC ' + str(targetname) + ' || Injected Period ' + str(period) + ' RPRS ' + str(rprs))
        bbox_props = dict(boxstyle="square,pad=0.3", facecolor='none', edgecolor='black')
        ax.annotate('Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs) + '\n' + 'Impact: ' + str(bestFimpact) + '\n' + 'T0: ' + str(bestFT0), xy = (time[5], 0.15), bbox = bbox_props)
        #fig.savefig('/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/ktransitfits_untrendy1/EPIC' + str(targetname) + '/fitper' + str(period) + 'rprs' + str(rprs) + '.png')
        show()

        """
                    
##########################
#CALCULATE SNR
##########################

    #SNR.append((rprs**2)*(((time[-1]-time[0])/period)**(1/2))*(1/CDPP))
                
            
##########################
#This needs to be changed
#based on what we consider
#to be a recovered transit:
#if best-fit period & rprs
#falls within a certain error?
    """
    if True:#period
        if True:#radius
            recovered += 1
    """
    if abs(bestFperiod-period) < .3:
        if abs(bestFrprs-rprs) < .1:
            recovered.append(True)
        else:
            recovered.append(False)
    else:
        recovered.append(False)


    #RECOVERED IS BOOLEAN!!! IT IS EITHER RECOVERED OR IT IS NOT

    
    """
    perarr.append(p)
    radarr.append(r)
    bestFperarr.append(bestFperiod)
    bestFrprsarr.append(bestFrprs)
    SDEarr.append(SDE)
    diffarr.append(diff[pindex][rindex])
    targetnamearr.append(targetname)
    """

    q = untrendy.median(time, flux)

    #SNR.append(sum(SNRforLC)/len(SNRforLC))

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 6))
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
    ax2.text(0.05, 1.05, 'SDE: ' + str(SDE) + '\n' + 'BLS Period: ' + str(bestFrprs), bbox=bbox_props, fontsize=7)
#fitting?????

    
    ax3.scatter(phase, y/ibi, color='k', s=2)
    ax3.set_xlabel('Phase')
    ax3.set_ylabel('Normalized Flux')
    ax3.set_title('Phase Folded')
    

    ax4.scatter(time, mergedfluxDetrend, color='k', s=2)
    ax4.plot(time, fitT.transitmodel, color='mediumaquamarine')
    ax4.set_ylim(-0.25,0.25)

    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Detrended Flux')
    ax4.set_title('Levenberg-Marquardt')
    #ax4.annotate('Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), xy = (time[1], 0.04), bbox = bbox_props)
    ax4.text(time[1], 0.12, 'Best fitting planet parameters: ' + '\n' + 'Period: ' + str(bestFperiod) + '\n' + 'RPRS: ' + str(bestFrprs), bbox=bbox_props, fontsize=7)

    plt.tight_layout()
    show()

    
        
    tempindex+=1

##########################
#END LOOP
##########################





#print(frac_recovered)



#print(frac_recovered)

            
######################################
#UNCOMMMENT TO SAVE DIFF, SDE, PER, RAD TEXT OUTPUT
#CHANGE PATH TO APPROPRIATE BLANK TEXT FILE
######################################


#print(len(perarr), len(radarr), len(bestFperarr), len(bestFrprsarr), len(SDEarr), len(diffarr))


#np.savetxt('/Users/sheilasagear/OneDrive/K2_Research/bls-ktransit/SDE_text/data213244700.txt', np.column_stack((perarr, radarr, bestFperarr, bestFrprsarr, SDEarr, diffarr)), fmt='%.5f')


#deletes extras from arrays

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

##########################
#INITIAL LIGHT CURVE
##########################
"""
pylab.cla()

xlabel('Time')
ylabel('Corrected Flux (normalized to 0)')
title('EPIC ' + targetname + ' Light Curve')
        

plot (time,flux)
show()
"""

##########################
#2D SCATTER ALL LIGHT CURVES SDE>6
##########################

"""
pylab.cla()
plt.scatter(detect_per, detect_rad, s=3)
axes.set_ylim([0,max(rad_range)])
axes.set_xlim([0,max(per_range)])
plt.xlabel('Orbital Period (days)')
plt.ylabel('Planet Radius/Star Radius')
plt.title('Detectable Transits (SDE > 6)')
plt.show()
"""

##########################
#3D SCATTER: ERROR, PERIOD, RADIUS FOR SDE>6
##########################

#print(len(detect_per))
#print(len(detect_rad))
#print(len(detect_diff))

"""
pylab.cla()
ax.scatter(detect_per, detect_rad, detect_diff, zdir='z', s=20, c=detect_SDE)
ax.set_xlabel('Orbital Period')
ax.set_ylabel('Planet Radius/Star Radius')
ax.set_zlabel('BLS Period - ktransit Period')
ax.set_title('Detectable Transits')

show()
"""

##########################
#3D SCATTER: ERROR, PERIOD, RADIUS FOR RECOVERED TRANSITS
##########################

"""
pylab.cla()
ax.scatter(per_range, rad_range, frac_recovered, zdir='z', s=20, c='blue')
ax.set_xlabel('Orbital Period')
ax.set_ylabel('Radius Ratio')
ax.set_zlabel('Fraction of Transits Recovered')
ax.set_title('Recovered Transits?')

show()
"""

fig, ax = plt.subplots(1, 1, figsize=[15,10])

print('SNR = ' + str(SNR))

print('Recovered = ' + str(recovered))

#ax.scatter(SNR, tested_period, s=3, c=recovered)
#show()


"""
#parameters???
gamma.cdf(SNR, 17.6, 1.00, 0.49)

#this doesn't even make sense
rv = gamma(frac_recovered)
ax.plot(SNR, rv.cdf(SNR), 'k-', lw=2, label='frozen pdf')
show()
"""


lowarr = sigma_clip(lowarr, sigma=2, iters=1)

#ax.scatter(rprsarr, lowarr, c=SDEarr)
#show()



