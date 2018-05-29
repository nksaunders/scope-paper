import numpy as np
import matplotlib.pyplot as pl
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

_, kp, cdpp6r, _, _, _, _, _, _ = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c03_nPLD.cdpp'), unpack = True, skiprows = 2)
fig, ax = pl.subplots(1, figsize=(7.5,5))
ax.scatter(kp, cdpp6r, color='palegreen', alpha = 0.05, zorder = -1)
ax.set_rasterization_zorder(-1)
bins = np.arange(7.5,18.5,0.5)
by = np.zeros_like(bins) * np.nan
for b, bin in enumerate(bins):
    i = np.where((cdpp6r > -np.inf) & (cdpp6r < np.inf) & (kp >= bin - 0.5) & (kp < bin + 0.5))[0]
    if len(i) > 10:
        by[b] = np.median(cdpp6r[i])
ax.scatter(bins, by, label = 'Raw K2', color='palegreen', edgecolors='k', s=50)
for iter in range(niter):
    cdpp = np.zeros_like(mags)
    for i, mag in enumerate(mags):
        # flux = np.load('batch/plot_run7/%2dmag%.2fmotion%.2f.npz' % (iter, mag, 1.))['flux']
        # perform aperture masking
        fpix = np.load('/Users/nksaunders/Documents/Research/scope/scope/batch/plot_run7/%2dmag%.2fmotion%.2f.npz' % (iter, mag, 1.))['fpix']
        # crop out extra pixels
        crop = np.array([f[2:10,2:10] for f in fpix])
        # sum into flux
        flux = np.sum(crop.reshape((len(crop)), -1), axis=1)
        cdpp[i] = CDPP(flux)
    if iter == 0:
        ax.scatter(mags, cdpp, color='mediumblue', marker='.', label = 'Synthetic (1x motion)')
    else:
        ax.scatter(mags, cdpp, color='mediumblue', marker='.')
ax.set_xlabel('Kepler Magnitude')
ax.set_ylabel('CDPP [ppm]')
ax.set_ylim(-30, 1500)
ax.set_xlim(8, 18)
ax.legend(loc = 'best')

pl.savefig('k2_benchmark.pdf', format='pdf', bbox_inches='tight')
