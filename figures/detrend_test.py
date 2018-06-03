import numpy as np
import matplotlib.pyplot as pl
from scope.scopemath import PSF, PLD
import scope
from tqdm import tqdm
import itertools
from everest.pool import Pool
from everest.missions.k2 import CDPP
from everest.config import EVEREST_SRC
import os
import os.path

# astroML format for consistent plotting style
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=10, usetex=True)

# Number of targets to simulate
niter = 5

# Magnitude and motion arrays
mags = np.arange(10., 16., .5)
m_mags = np.arange(0., 21., 1)

sK2 = scope.Target()
sK2.GenerateLightCurve(ncadences=1)

# Plot several different motion vectors
_, kepler_kp, kepler_cdpp6 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'kepler.cdpp'), unpack = True)
fig, ax = pl.subplots(1, figsize=(7.5,5))
ax.plot(kepler_kp, kepler_cdpp6, color='palegreen', marker='.', alpha = 0.01, zorder = -1)
ax.set_rasterization_zorder(-1)
bins = np.arange(7.5,18.5,0.5)
by = np.zeros_like(bins) * np.nan
for b, bin in enumerate(bins):
    i = np.where((kepler_cdpp6 > -np.inf) & (kepler_cdpp6 < np.inf) & (kepler_kp >= bin - 0.5) & (kepler_kp < bin + 0.5))[0]
    if len(i) > 10:
        by[b] = np.median(kepler_cdpp6[i])
ax.scatter(bins, by, color='palegreen', label = 'Kepler', edgecolors='k', s=50)
for m_mag, color in zip([1, 2, 5, 10], ['mediumblue', 'darkviolet', 'k', 'orange']):
    cdpp = [[] for mag in mags]
    for i, mag in enumerate(mags):
        for iter in range(niter):
            aperture = np.ones((13,13))
            '''
            aperture = np.zeros((13,13))
            if m_mag == 1:
                aperture[5:8,5:8] = 1
            elif m_mag == 2:
                aperture[4:9,4:9] = 1
            elif m_mag == 5:
                aperture[3:10,3:10] = 1
            elif m_mag == 10:
                aperture[2:11,2:11] = 1
            else:
                aperture = np.ones((13,13))
            '''
            fpix = np.load('/Users/nksaunders/Documents/Research/scope/scope/batch/plot_run7/%2dmag%.2fmotion%.2f.npz' % (iter, mag, m_mag))['fpix']
            # ferr = np.ones(len(fpix))
            ferr = np.load('/Users/nksaunders/Documents/Research/scope/scope/batch/error_run/mag12roll%i.npz' % m_mag)['ferr']
            det_flux, rawflux = PLD(fpix, ferr, [], sK2.t, aperture)
            cdpp[i].append(CDPP(det_flux))
    cdpp = np.nanmean(np.array(cdpp), axis = 1)
    ax.plot(mags, cdpp, '.', color = color, label = 'Synthetic (%dx motion)' % m_mag)
    ax.plot(mags, cdpp, '-', color = color)
ax.set_xlabel('Kepler Magnitude')
ax.set_ylabel('CDPP [ppm]')
ax.set_ylim(-30, 2500)
ax.set_xlim(8, 18)
ax.legend(loc = 'best')

pl.show()

import pdb; pdb.set_trace()

pl.savefig('detrend_test.pdf', format='pdf', bbox_inches='tight')
