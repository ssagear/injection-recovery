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

from itertools import izip_longest

import scipy.constants as sc

import matplotlib


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

total_per = []
total_rprs = []

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

natural_val = np.linspace(-.523, 1.415, num=8, endpoint=True)
natural_val = natural_val.round(decimals=2)

per_segments = [10**x for x in natural_val]
rprs_segments = [x / 20 for x in range(8)]

p = np.poly1d(z)

r_list = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]

r_per_list = []
r_rprs_list = []

u_list = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]

u_per_list = []
u_rprs_list = []

def loguniform(low=float(0.3), high=float(26)):
    return 10**(np.random.uniform(np.log10(low), np.log10(high)))


def randomInc(n):

    randNums = np.random.uniform(low=0.5, high=1, size=n)
    incs = np.arccos(2*randNums - 1)

    return incs


#################################
iterations = 20
#################################

filecount = sum(1 for line in open('/Users/sheilasagear/Dropbox/K2/pipeline/lcs_txt/epicids_injection_10oct.txt'))

print(filecount)

with open('/Users/sheilasagear/Dropbox/K2/pipeline/lcs_txt/epicids_injection_10oct.txt') as f:
    content = f.readlines()

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

    plt.plot(time, flux)
    plt.title('229227256 without processing')
    plt.show()

#normalizes flux values to 1.0 (at avg of flux values)
    #flux_count = len(flux)
    #sumflux = sum(flux)
    #flux_avg = sumflux/flux_count

    #flux = [i/flux_avg for i in flux]

#targets titled by EPIC names
    targetname = file[82:]
    targetname = targetname[:-35]
    print('TARGET: ' + str(targetname))

    campaign_no = file[37:]
    campaign_no = campaign_no[:-31]

#normalize to 0.0
    flux = [i-1 for i in flux]


#SIGMA CLIPPING
    flux = sigma_clip(flux, sigma=2, iters=1)


    flux = flux.filled(fill_value=0)


    plt.plot(time, flux)
    plt.show()
    exit(0)

    for i in range(iterations):
        print(i)

        period = loguniform()
        tested_period.append(period)
        print('Inj period = ' + str(period))

        rprs = np.random.uniform(low=.01, high=.4)
        rprs = float(rprs)
        tested_rprs.append(rprs)
        print('Inj rprs = ' + str(rprs))

        Mc = 9e28#kilograms - mass of L dwarf
        Rs = 7.1492e7#meters - radius of L dwarf

        periodsec = 86400*period
        a1 = sc.G*Mc*periodsec**2
        a2 = (4*sc.pi)**2
        a = np.cbrt(a1/a2)#semimajor axis, meters

        i = randomInc(1)

        impact = (a*math.cos(i))/Rs
        impact = 0.0

        print('Inj impact parameter = ' + str(impact))
        print(filecount)
        print(tempindex)

        total_per.append(period)
        total_rprs.append(rprs)

        #if impact > (1 + rprs):
        #    continue

    #uncomment if extra



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
        impact=impact, # an impact parameter
        rprs=rprs,   # planet stellar radius ratio
        ecosw=0.0,  # eccentricity vector
        esinw=0.0,
        occ=0.0)    # a secondary eclipse depth in ppm

        M.add_data(time=np.array(time[:]))      # integration time of each timestamp

        tmod = M.transitmodel # the out of transit data will be 0.0 unless you specify zpt
    #plt.plot(M.time,tmod)



        RMS = sqrt(mean(square(flux)))




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

        nf = 10000.0
        fmin = .035
        df = 0.001
        nbins = 300
        qmi = 0.001
        qma = 0.3

        results = bls.eebls(time, mergedfluxDetrend, u, v, nf, fmin, df, nbins, qmi, qma)

    #RESULTS:
    #power, best_period, best_power, depth, q, in1, in2
    #0      1            2           3      4  5    6


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

        if abs(bestFperiod-period) < .10*period and abs(bestFrprs-rprs) < .10*rprs:
            print('RECOVERED')
            recovered = True
            recoveredarr.append(True)
            recovered_period.append(period)
            recovered_rprs.append(rprs)
            for p_i in range(len(per_segments)-1):
                for r_i in range(len(rprs_segments)-1):
                    if period > per_segments[p_i] and period < per_segments[p_i+1]:
                        if rprs > rprs_segments[r_i] and rprs < rprs_segments[r_i+1]:
                            r_list[r_i][p_i] += 1
        else:

            recovered = False
            recoveredarr.append(False)
            unrecovered_period.append(period)
            unrecovered_rprs.append(rprs)
            #ax.scatter(rprs, RMS, c='b', s=8)

            for p_i in range(len(per_segments)-1):
                for r_i in range(len(rprs_segments)-1):
                    if period > per_segments[p_i] and period < per_segments[p_i+1]:
                        if rprs > rprs_segments[r_i] and rprs < rprs_segments[r_i+1]:
                            u_list[r_i][p_i] += 1






    tempindex+=1


block_percent = [[0 for i in range(len(per_segments))] for j in range(len(rprs_segments))]

for r in range(len(rprs_segments)-1):
    for p in range(len(per_segments)-1):
        if u_list[r][p]+r_list[r][p] > 0:
            block_percent[r][p] = r_list[r][p]/(u_list[r][p]+r_list[r][p])

block_percent = block_percent[::-1]

per_segments = [i.round(decimals=2) for i in per_segments]

a = np.arange(26)

ybins = [x / 20 for x in range(8)]

print(recovered_period)
print(recovered_rprs)

counts, _, _ = np.histogram2d(recovered_period, recovered_rprs, bins=(per_segments, ybins))
counts_tot, _, _ = np.histogram2d(total_per, total_rprs, bins=(per_segments, ybins))

for i in range(len(counts.T)):
    for j in range(len(counts.T[i])):
        counts.T[i][j] = counts.T[i][j]/counts_tot.T[i][j]
        if np.isnan(counts.T[i][j]):
            counts.T[i][j] = 0

np.savetxt('/Users/sheilasagear/Dropbox/K2/heatmap_inclination_csv/' + str(targetname) + '-cutoff.csv', counts.T, fmt='%f', delimiter=',')

matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['xtick.minor.width'] = 0

fig, ax = plt.subplots()
heatmap = ax.pcolormesh(per_segments, ybins, counts.T, cmap='Blues_r')
ax.set_xscale('log')

ax.xaxis.set_ticks(per_segments)
ax.xaxis.set_ticklabels(per_segments)

ax.yaxis.set_ticks(ybins)
ax.yaxis.set_ticklabels(ybins)

ax.set_title('Fraction of Transits Recovered (EPIC 229227256)')
ax.set_xlabel('Period (days)')
ax.set_ylabel('Radius (stellar radii)')


cbar = plt.colorbar(heatmap)
heatmap.set_clim(0.0, 1.0)
cbar.set_label('Fraction recovered', rotation=270, labelpad=13)


plt.savefig('/Users/sheilasagear/Dropbox/ssagear_k2/plots/in_progress/heatmap_229227256_cutoff.png')

plt.show()
