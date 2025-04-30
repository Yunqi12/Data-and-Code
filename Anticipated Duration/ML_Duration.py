## adapted from Liu (2022)

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import obspy
from obspy.core.utcdatetime import UTCDateTime as UTC
import random
import datetime
import scipy
from scipy.optimize import minimize
from scipy import linalg
from scipy.special import factorial as fac
#from confidence import get_CI
from geopy.distance import geodesic
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import os

import emcee
import corner

font = {'family':'sans-serif','sans-serif':['Helvetica'],
        'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('font', **font)

np.array([0.98, 1.51])+0.24

## Read catalog using pandas dataframe and convert to numpy array
def get_CI(mean, cl, sd):
  loc = stats.norm.ppf(1 - cl/2)
  rng_val = stats.norm.cdf(loc - mean/sd)

  lwr_bnd = value - rng_val
  upr_bnd = value + rng_val

  return_val = (lwr_bnd, upr_bnd)
  return(return_val)

def read_cat(filename):
    pdcat = pd.read_csv(filename, dtype='str')
    c    = pdcat.to_numpy()
    time = c[:,0]
    lat  = c[:,1].astype('float')
    lon  = c[:,2].astype('float')
    dep  = c[:,3].astype('float')
    mag  = c[:,4].astype('float')
    evid = c[:,11]
    obdt = []
    pydt = []
    relt = []
    for i in range(len(time)):
        obdt.append(UTC(time[i]))
        pydt.append(UTC(time[i]).datetime)
        relt.append((UTC(time[i]) - UTC(time[0])) / 86400)
    return (evid, obdt, pydt, relt, lat, lon, dep, mag)


def epoch_Mc(mag, obspyDT, Nt=10, plot='no', title=''):
    Ndt = (obspyDT[-1]-obspyDT[0])/Nt

    epochs = []
    for i in range(int(Nt)):
        epochs.append(np.where(np.array(obspyDT)>=obspyDT[0]+Ndt*i)[0][0])

    Mcs = []
    for i in range(Nt):
        if i==Nt-1:
            sub_mag = mag[epochs[i]:-1]
        else:
            sub_mag = mag[epochs[i]:epochs[i+1]]
        Mcs.append(maxc_Mc(sub_mag, plot=plot, title=title))
    return epochs, Ndt, Mcs


def omorilaw(x, params):
    c, K, p = params
    return K/((x+c)**p)
    
def omorilaw_integrate(x, params):
    c, K, p = params
    return K/(1-p)*((x+c)**(1-p)-(x-1+c)**(1-p))


def omori_syn(c, p, tmin, tmax, N):
    """
    Adapted from tgeobel's GitHub: https://github.com/tgoebel/aftershocks
    Felzer et al. 2002, Triggering of 1999 Mw 7.1 Hector Mine earthquake
    - define create power-law distributed aftershock time vector between tmin and tmax
    - tmin can be = 0, Omori-K parameter is a fct. of all four parameter (tmax-tmin), p, and N
    INPUT:  c, p       - omori parameters describing time shift for complete recording and rate decay exponent
                       - in alphabetical order
           tmin, tmax  - time window for aftershock catalog
           N           - total number of aftershocks
    """
    vRand = np.random.random_sample( N)
    #===========================================================================
    #          case1:  p != 1
    #===========================================================================
    #if p != 1.0: #abs(p - 1) < 1e-6:
    p += 1e-4 # this will make it unlikely for p to be exactly 1
    a1 = (tmax + c)**(1-p)
    a2 = (tmin + c)**(1-p)
    a3 = vRand*a1 + (1-vRand)*a2#
    otimes = a3**(1/(1-p))-c
#     else: # p == 1
#         a1 = np.log( tmax + c)
#         a2 = np.log( tmin + c)
#         a3 = vRand*a1 + (1-vRand)*a2
#         otimes = np.exp( a3) - c
    otimes.sort()
    return otimes


def ogata_logL(otimes, params):
    c, K, p = params
    S, T = min(otimes), max(otimes)
    n = len(otimes)
    if abs(p - 1) < 1e-8:
        A = np.log(T+c) - np.log(S+c)
    else:
        A = ((T+c)**(1-p) - (S+c)**(1-p))/(1-p)
    L = -n*np.log(K) + p*np.sum(np.log(otimes+c)) + K*A
    return L


def hol_logL(otimes, params):
    c, p = params
    S, T = min(otimes), max(otimes)
    n = len(otimes)
    if abs(p - 1) < 1e-8:
        D = (1/c) / (np.log(1+T/c) - np.log(1+S/c))
    else:
        D = ((1-p)/c) / ((1+T/c)**(1-p) - (1+S/c)**(1-p))
    L = -n*np.log(D) + p*np.sum(np.log(1+otimes/c))
    return L



def hol_post(otimes, params):
    c, p = params
    S, T = min(otimes), max(otimes)
    n = len(otimes)
    if abs(p - 1) < 1e-8:
        D = (1/c) / (np.log(1+T/c) - np.log(1+S/c))
    else:
        D = ((1-p)/c) / ((1+T/c)**(1-p) - (1+S/c)**(1-p))
    prob = np.prod(D/((1+otimes/c)**p)) * 1/c * 1/p
    return prob



def hol_getK(otimes, params):
    c, p = params
    S, T = min(otimes), max(otimes)
    n = len(otimes)
    if abs(p - 1) < 1e-8:
        D = (1/c) / (np.log(1+T/c) - np.log(1+S/c))
    else:
        D = ((1-p)/c) / ((1+T/c)**(1-p) - (1+S/c)**(1-p))
    L = -n*np.log(D) + p*np.sum(np.log(1+otimes/c))
    K = D * n * c**p
    return L, K

# Read catalog
cat = read_cat('/Users/yunqhuang/Desktop/Code_Seismica/Results/Python/WP_50k_AS.txt')
evid, obdt, pydt, relt, lat, lon, dep, mag = np.array(cat)


# covered time period
t_start = UTC('2021091')
t_termi = UTC('20300101')




# MS info
ms_id = np.argmax(mag[obdt >= t_start])
ms = {'id':evid[ms_id], 'obdt':obdt[ms_id], 'pydt':pydt[ms_id], 'tAMS':relt[ms_id],
      'lat':lat[ms_id], 'lon':lon[ms_id], 'dep':dep[ms_id], 'mag':mag[ms_id]}

# save meta parameters
meta = {'t_start'   : t_start,
        't_termi'   : t_termi,
        'Mcut'      : 0.8,
        'rmax'      : 10**(0.25*ms['mag']-.22),  # max radius of influence (Gardner & Knopoff, 1967)
        'nbin'      : 100,
        'c_bound'   : [1e-4,   2],
        'K_bound'   : [   2, 1e3],
        'p_bound'   : [  .2,  2],
        'c0'        : .5,
        'K0'        : 50,
        'p0'        : 1.1,
        'ylim'      : [1e-3, 1e5],
        'xlim'      : [1e-3, 1e3],
        'syn_c'     : 0.6,
        'syn_p'     : 1.3,
        'syn_tStart': 1e-2,
        'syn_tEnd'  : 1e3,
        'syn_N'     : 4000,}


## get aftershocks within a radius
aR = []
for i in range(len(evid)):
    aR.append(geodesic((ms['lat'], ms['lon']),(lat[i], lon[i])).km)
aR = np.array(aR)
rd_id = aR <= meta['rmax']
meta['rmax'] = float('inf')

# selections
mc_id  = mag  >= meta['Mcut']
as_id  = obdt >  ms['obdt']
end_id = obdt <  t_termi
rd_id  = aR   <= meta['rmax']
select = mc_id * as_id * end_id * rd_id
evid, obdt, pydt, relt, lat, lon, dep, mag = np.array(cat).T[select].T

print(' Mainshock magnitude: %.2f \n' % ms['mag'],
      'Mainshock time: %s \n'         % ms['obdt'],
      'Minimum magnitude: %.2f \n'    % meta['Mcut'],
      'Maximum radius: %.2f \n'       % meta['rmax'],
      'Start time: %s \n'             % meta['t_start'],
      'End time: %s \n'               % meta['t_termi'],
      '# events selected: %d \n'      % select.sum())


relt = relt.astype('float')
lat = lat.astype('float')
lon = lon.astype('float')
dep = dep.astype('float')
mag = mag.astype('float')
relt = relt - ms['tAMS']


# Choose dataset:
data = '1'

if data == '0':
    # Synthetic dataset
    c = meta['syn_c']
    p = meta['syn_p']
    synt = omori_syn(c, p, meta['syn_tStart'], meta['syn_tEnd'], meta['syn_N'])
    otimes = np.array(sorted(synt))

elif data == '1':
    # Real dataset
    otimes = np.array(sorted(relt))



# Calc likelihood:
o_objFunc = lambda X: ogata_logL(otimes, X)
h_objFunc = lambda X: hol_logL(otimes, X)

disp = 0
method = 'SLSQP'


# Ogata 1989: MLE
o_Bounds  = np.array([meta['c_bound'], meta['K_bound'], meta['p_bound']])
o_Par0    = np.array([meta['c0'], meta['K0'], meta['p0']])
ogata_fit = scipy.optimize.minimize(o_objFunc, o_Par0, bounds=o_Bounds, tol = 1e-4, method=method, options={'disp': disp, 'maxiter':500})
print(ogata_fit)

# Holschneider et al., 2012: Bayesian
h_Bounds   = np.array([meta['c_bound'], meta['p_bound']])
h_Par0     = np.array([meta['c0'], meta['p0']])
holsch_fit = scipy.optimize.minimize(h_objFunc, h_Par0, bounds=h_Bounds, tol = 1e-4, method=method, options={'disp': disp, 'maxiter':500})
finalL, K  = hol_getK(otimes, holsch_fit['x'])
print('\n',holsch_fit)
print('       K:',K)


meta['Ogata_fit'] = list(ogata_fit['x'])
meta['Bayes_fit'] = [holsch_fit['x'][0], K, holsch_fit['x'][1]]

meta

Cs = np.linspace(meta['c_bound'][0], meta['c_bound'][1])
Ps = np.linspace(meta['p_bound'][0], meta['p_bound'][1])

L = np.zeros([len(Cs), len(Ps)])
for i in range(len(Cs)):
    for j in range(len(Ps)):
        L[i,j] = h_objFunc((Cs[i],Ps[j]))
L = L.T
        

plot_res = 'both'

bins = np.logspace(np.log10(otimes[0]), np.log10(otimes[-1]), meta['nbin'])
count, bine = np.histogram(otimes, bins=bins)
bin_loc = (bine[1:] + bine[:-1]) / 2
occ_dens = count/np.diff(bins)


plt.figure(figsize=[10,8])
plt.scatter(bin_loc, occ_dens, s=100, ec='k', fc='lightgrey')
plt.xscale('log')
plt.yscale('log')
plt.ylim(meta['ylim'])
plt.xlabel('Days since mainshock')
plt.ylabel('# events/day')


linspace_array = np.linspace(np.max(bin_loc), 100000,len(bin_loc))



lgdstr = 'Ogata {c,K,p}={%.2f, %.2f, %.2f}'
if plot_res=='o':
    plt.plot(bin_loc, omorilaw(bin_loc, meta['Ogata_fit']), c='r')
    plt.plot(linspace_array, omorilaw(linspace_array, meta['Ogata_fit']), c='black')
elif plot_res=='h':
    plt.plot(bin_loc, omorilaw(bin_loc, meta['Bayes_fit']), c='b')
    plt.plot(linspace_array, omorilaw(linspace_array, meta['Bayes_fit']), c='black')
elif plot_res=='both':
    plt.plot(bin_loc, omorilaw(bin_loc, meta['Bayes_fit']), c='black')
    plt.plot(linspace_array, omorilaw(linspace_array, meta['Bayes_fit']), c='black',label= 'Omori Law')

aa = np.concatenate((bin_loc, linspace_array))
slope = np.array([0.0204,0.0844,0.0510])

linspace_array = np.linspace(1, 100000, 100000)

bb = np.ones(len(aa)) * slope[0];
plt.plot(aa, bb, c= 'blue',label= 'Smoothed Leonard (2008)');
bb = np.ones(len(aa)) * slope[1];
plt.plot(aa, bb, c= 'red',label= 'This Study');
bb = np.ones(len(aa)) * slope[2];
plt.plot(aa, bb, linestyle='--',c= 'g',label= 'Combined Result');


## Calculate the Omori's Law Intergrate (events/day)
plt.legend(loc='upper right',fontsize=14)
plt.title('Mc = 0.80 (Ml)')
plt.show()