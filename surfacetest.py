from __future__ import division

import numpy as np

import bls, ktransit, pylab, os, untrendy

from astropy.io import fits

from astropy.stats import sigma_clip

import matplotlib.pyplot as plt

########################################
#time and flux arrays must be np arrays
########################################

def data_K2(fitspath):

    time = []
    flux = []

    fitspath = fitspath.strip('\n')
    #print(fitspath)
    if os.path.isfile(fitspath):
        f = fits.open(fitspath, ignore_missing_end=True)

        bestaper = f[1].data

        for i in range(3300):
            time.append(bestaper[i][0])
            flux.append(bestaper[i][1])

        flux_count = len(flux)
        sumflux = sum(flux)
        flux_avg = sumflux/flux_count

        flux = [i/flux_avg for i in flux]

        #targets titled by EPIC names
        targetname = fitspath[79:]
        targetname = targetname[:-23]

        campaign_no = fitspath[90:]
        campaign_no = campaign_no[:-19]

        #normalize to 0.0
        flux = [i-1 for i in flux]

        flux = np.asarray(flux)
        time = np.asarray(time)

        #SIGMA CLIPPING
        flux = sigma_clip(flux, sigma=4, iters=1)
        flux = flux.filled(fill_value=0)

        return time, flux, targetname, campaign_no

def CDPP(fitspath):

    f = fits.open(fitspath)
    QCDPP = f[1].header['QCDPP6']
    return QCDPP

def createT(time, T0):#time array, transit mid-time, impact parameter

    #random period and radius: uniformly weighted
    #you can inject a specific period if you want
    period = np.random.uniform(low=1, high=26)
    print('inj period: ' + str(period))
    rprs = np.random.uniform(low=.01, high=.4)
    rprs = float(rprs)
    print('inj rprs: ' + str(rprs))

    #this is Tom Barclay's ktransit package I use for injection and fitting (https://github.com/mrtommyb/ktransit)
    M = ktransit.LCModel()
    M.add_star(
    rho=1.5, # mean stellar density in cgs units
    ld1=0.2, # ld1--4 are limb darkening coefficients
    ld2=0.4, # assuming quadratic limb darkening
    ld3=0.0,
    ld4=0.0,
    dil=0.0, # a dilution factor: 0.0 -> transit not diluted, 0.5 -> transit 50% diluted
    zpt=0.0  # photometric zeropoint
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

    return tmod, period, rprs

#tmod = the noiseless LC with the transits you want to inject (flux array)


def addT(flux, tmod):#flux array, transit flux array you want to inject
    merged_flux = flux + tmod
    return merged_flux

#dt=x (2 in trappist ntbk) - untrendy median parameter
def detrend(time, merged_flux):
    trend = untrendy.median(time, merged_flux, dt=0.6)
    mergedfluxD = np.zeros(len(time))
    for i in range(len(time)):
        mergedfluxD[i] = merged_flux[i]-trend[i]

    return mergedfluxD


#dan foreman-mackey's python BLS from Kovacs' fortran subroutine
#https://github.com/dfm/python-bls
def BLS(time, mergedfluxD):#time array, detrended flux arr (with transits injected)

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

    results = bls.eebls(time, mergedfluxD, u, v, nf, fmin, df, nbins, qmi, qma)

    #RESULTS:
    #power, best_period, best_power, depth, q, in1, in2
    #0      1            2           3      4  5    6

    SR_array = results[0]
    max_SR = max(SR_array)
    avg_SR = np.mean(SR_array)
    sd_SR = np.std(SR_array)
    #normaze SR_array between 0 and 1
    SR_array = [i/max_SR for i in SR_array]

    freq = fmin + np.arange(nf)*df

    #Signal Detection Efficiency
    SDE = (max_SR-avg_SR)/sd_SR

    #depth
    high = results[3]*results[4]
    low = high - results[3]
    fit = np.zeros(nbins) + high # H
    fit[results[5]:results[6]+1] = low # L
    depth = high - low

    return results, SR_array, freq, SDE, depth
    #plotting freq vs. SR_array gives you a power spectrum
    #depth is the BLS transit depth



def fitT(time, mergedfluxD, guessper, guessrprs, T0):#time arr, merged (transits added) flux arr detrended, guess period and rp/rs from BLS, transit mid-time

    fitT = ktransit.FitTransit()
    fitT.add_guess_star(rho=1.5)
    fitT.add_guess_planet(
    period=guessper, impact=0.0,
    T0=T0, rprs=guessrprs)
    fitT.add_data(time=time, flux=mergedfluxD)

    vary_star = []      # free stellar parameters
    vary_planet = (['period',       # free planetary parameters
        'rprs'])                # free planet parameters are the same for every planet you model

    fitT.free_parameters(vary_star, vary_planet)
    fitT.do_fit()                   # run the fitting

    bestFplanet = fitT.fitresultplanets.items()
    bestFperiod = bestFplanet[0][1]['period']#Best Fit Period
    bestFrprs = bestFplanet[0][1]['rprs']#Best Fit rprs

    fitT.print_results()

    return bestFperiod, bestFrprs




def is_recovered(period, fitper, rprs, fitrprs):#true (injected) period, L-M fitting period, true rp/rs, L-M fitting rp/rs
    if abs(fitper-period) < .03*period and abs(fitrprs-rprs) < .03*rprs:
        return True
    else:
        return False


def main():

    #this is a polynomial fit of BLS transit depth vs true rp/rs - gives us a rp/rs to guess when starting transit fitting
    z = [0, 0, 0, 0]
    z[0] = -3.22185
    z[1] = 1.02655
    z[2] = 1.6613
    z[3] = .225603
    p = np.poly1d(z)

    path = '/Users/sheilasagear/Dropbox/K2/K2_targets/C06to13fits/hlsp_k2sff_k2_lightcurve_229228324-c07_kepler_v1_llc.fits'

    time, flux, name, cnum = data_K2(path)

    print('TARGET: ' + str(name))
    print('CAMPAIGN: ' + str(cnum))

    QCDPP = CDPP(path)

    print('QCDPP: ' + str(QCDPP))

    T0 = time[0]+((time[-1]-time[0])/2)
    tmod, period, rprs = createT(time, T0)

    merged_flux = addT(flux, tmod)
    plt.plot(time, merged_flux)
    plt.show()
    mergedfluxD = detrend(time, merged_flux)
    plt.plot(time, mergedfluxD)
    plt.show()

    results, power, freq, SDE, depth = BLS(time, mergedfluxD)

    guessper = results[1]
    guessrprs = p(depth)

    fitper, fitrprs = fitT(time, mergedfluxD, guessper, guessrprs, T0)

    isrec = is_recovered(period, fitper, rprs, fitrprs)
    print(isrec)


for i in range(1):#however many planets you want to try
