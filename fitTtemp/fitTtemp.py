import ktransit
from ktransit import FitTransit

import numpy as np

import bls, ktransit, math, pylab, os, batman

from scipy.stats import gamma

import untrendy

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = plt.gca()

import csv

#from mpl_toolkits.mplot3d import Axes3D
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




time = []
flux = []

file = 'hlsp_k2sff_k2_lightcurve_212826600-c06_kepler_v1_llc-default-aper.csv'

with open(r'/Users/sheilasagear/OneDrive/K2_Research/CorrectedLC_EPICID/Cycle2_4_CSV/hlsp_k2sff_k2_lightcurve_212826600-c06_kepler_v1_llc-default-aper.csv') as f:
    reader = csv.reader(f)

    for row in reader:
        time.append(row[0])
        flux.append(row[1])

    for i in range(len(time)):
        time[i] = float(time[i])
        flux[i] = float(flux[i])

    #print(time)
    #print('\n')
    #print(flux)

#normalize flux values to 1.0
    flux_count = len(flux)
    sumflux = sum(flux)
    flux_avg = sumflux/flux_count
    
    flux = [i/flux_avg for i in flux]

#targets titled by EPIC names

    targetname = file[25:]
    targetname = targetname[:-35]
    #title('EPIC ' + targetname + ' Light Curve')
    print("TARGET: EPIC " + str(targetname))
        
    campaign_no = file[36:]
    campaign_no = campaign_no[:-31]
    print("CAMPAIGN " + str(campaign_no))
        
#normalize to 0.0
    flux = [i-1 for i in flux]

    print('flux = ' + str(type(flux)))

    flux = np.asarray(flux)

    #print(type(flux))


#SIGMA CLIPPING
    flux = sigma_clip(flux, sigma=3, iters=1)

#uncomment if extra time stamp
    #time.pop()


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

    #T0 = np.random.uniform(low=.01, high=1)
    #period = np.random.uniform(low=.05, high=100)
    #impact = np.random.uniform(low=0.0, high=1.0)
    #rprs = np.random.uniform(low=.05, high=1.0)
    T0 = 2.0
    period = 5
    impact = 0.0
    rprs = .2
            
    M.add_planet(
                T0=T0,     # a transit mid-time  
            period=period, # an orbital period in days
            impact=impact, # an impact parameter
            rprs=rprs,   # planet stellar radius ratio  
            ecosw=0.0,  # eccentricity vector
            esinw=0.0,
    occ=0.0)    # a secondary eclipse depth in ppm
            

    M.add_data(time=np.array(time[:])),

    tmod = M.transitmodel# the out of transit data will be 0.0 unless you specify zpt
    #pylab.cla()

    plot(time, tmod)
    title('KTRANSIT injected LC')
    #show()
        

    if len(tmod) != len(flux):
        flux = flux[0:len(flux)-1]
        
    merged_flux = tmod + flux

    mergedfluxDetrend = tmod

    """
    mergedfluxDetrend, ferr = untrendy.untrend(time, merged_flux)

    trend, ferr = untrendy.untrend(time, flux)

    mergedfluxDetrend = list()
    
    for i in range(len(trend)):
        val = merged_flux[i]-merged_flux[i-1]
        mergedfluxDetrend.append(val)

    """
    

    u = [0.0]*len(time)
    v = [0.0]*len(time)

    u = np.array(u)
    v = np.array(v)

    nbins = 200
        
    results = bls.eebls(time, mergedfluxDetrend, u, v, 1000.0, .1, .001, nbins, .001, .3)
            
#RESULTS:
#power, best_period, best_power, depth, q, in1, in2
#0      1            2           3      4  5    6

    print('BLS period: ' + str(results[1]))


    fitT = FitTransit()
    fitT.add_guess_star(rho=7.0)    
    fitT.add_guess_planet(
        period=results[1], impact=0.0,
        T0=2.0, rprs=0.3)
    fitT.add_data(time=time, flux=mergedfluxDetrend)

    vary_star = ['rho']      # free stellar parameters
    vary_planet = (['period', 'T0',       # free planetary parameters,
                        'rprs'])                # free planet parameters are the same for every planet you model

    fitT.free_parameters(vary_star, vary_planet)
    fitT.do_fit()                   # run the fitting

    fitT.print_results()   

    fig = ktransit.plot_results(time,mergedfluxDetrend,fitT.transitmodel)

    fig.savefig('fitTemp.png')
