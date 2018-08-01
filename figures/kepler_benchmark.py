import scope
import numpy as np
import matplotlib.pyplot as pl
from everest.missions.k2 import CDPP
from everest.config import EVEREST_SRC
import os
import os.path

def Generate():

    mags = np.arange(10., 16., .5)
    for iter in range(5):
        print("Iteration %i" % iter)

        for m in mags:
            print("Running mag = %.2f" % m)

            sK2 = scope.Target(variable=True)

            # check if lc exists
            if os.path.isfile('kepler_benchmark/%2dmag%.2fmotion0.0.npz' % (iter, m)):
                print("Mag = %.2f, m_mag = 0.0 already exists!" % m)

            # create missing lc
            else:
                fpix, flux, ferr = sK2.GenerateLightCurve(mag=m, roll=0, background_level=10, ncadences=500, apsize=5)
                np.savez('kepler_benchmark/%2dmag%.2fmotion0.0.npz' % (iter, m), fpix=fpix, flux=flux, ferr=ferr)

def Plot():

    # astroML format for consistent plotting style
    from astroML.plotting import setup_text_plots
    setup_text_plots(fontsize=10, usetex=True)
    pl.rcParams['agg.path.chunksize'] = 1000

    # Compare zero-motion synthetic data to original Kepler raw CDPP
    _, kepler_kp, kepler_cdpp6 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'kepler.cdpp'), unpack = True)
    fig, ax = pl.subplots(1)
    ax.scatter(kepler_kp, kepler_cdpp6, color='red', alpha = 0.01, zorder = -1)
    ax.set_rasterization_zorder(-1)
    bins = np.arange(7.5,18.5,0.5)
    by = np.zeros_like(bins) * np.nan
    for b, bin in enumerate(bins):
        i = np.where((kepler_cdpp6 > -np.inf) & (kepler_cdpp6 < np.inf) & (kepler_kp >= bin - 0.5) & (kepler_kp < bin + 0.5))[0]
        if len(i) > 10:
            by[b] = np.median(kepler_cdpp6[i])
    ax.plot(bins, by, 'ro', label = 'Kepler', markeredgecolor = 'k')

    mags = np.arange(10., 16., .5)

    for iter in range(5):
        cdpp = np.zeros_like(mags)
        for i, m in enumerate(mags):
            # flux = np.load('batch/plot_run7/%2dmag%.2fmotion%.2f.npz' % (iter, mag, 0.))['flux']

            # load in fpix
            fpix = np.load('kepler_benchmark/%2dmag%.2fmotion0.0.npz' % (iter, m))['fpix']
            # crop out extra pixels
            crop = np.array([f[2:3,2:3] for f in fpix])
            # sum into flux
            flux = np.sum(crop.reshape((len(crop)), -1), axis=1)

            # calculate CDPP
            cdpp[i] = CDPP(flux)
            # import pdb; pdb.set_trace()

        if iter == 0:
            ax.plot(mags, cdpp, 'b.', label = 'Synthetic (0x motion)')
        else:
            ax.plot(mags, cdpp, 'b.')

    ax.set_xlabel('Kepler Magnitude')
    ax.set_ylabel('CDPP [ppm]')
    ax.set_ylim(-10, 500)
    ax.set_xlim(8, 18)
    ax.legend(loc = 'best')

    # pl.show()
    pl.savefig('kepler_benchmark.pdf', format='pdf', bbox_inches='tight')

# Generate()
Plot()
