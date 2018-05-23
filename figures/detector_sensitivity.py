import scope;
import numpy as np
import matplotlib.pyplot as pl
from tqdm import tqdm

# astroML format for consistent plotting style
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=10, usetex=True)

star = scope.Target()
fpix, flux, ferr = star.GenerateLightCurve(ncadences=1)

sub_pixel = np.zeros((994,994))

star.DisplayDetector()
detector = 1+(star.detector/(np.max(np.abs(star.detector))))

countx = 0; county = 0
indx = 0; indy = 0

for i in tqdm(range(994)):
    for j in range(994):
        sub_pixel[i][j] = fpix[0][indx][indy] * (detector[i][j])

        county += 1
        if county >= 142:
            indy += 1
            county = 0

        if j == 993:
            indy = 0

    countx += 1
    if countx >= 142:
        indx += 1
        countx = 0

    if i == 993:
        indx = 0

sub_pixel += np.abs(np.min(sub_pixel))

pl.clf()
fig, (ax0, ax1, ax2) = pl.subplots(1, 3, figsize=(12,4))

fig0 = ax0.imshow(detector, cmap='gray', origin='lower', extent=[0,6,0,6])
ax0.set_xlabel('(a)')
fig1 = ax1.imshow(sub_pixel, vmin=0, vmax=50000, origin='lower', extent=[0,6,0,6])
ax1.set_xlabel('(b)')
fig2 = ax2.imshow(fpix[0], vmin=0, vmax=50000, origin='lower')
ax2.set_xlabel('(c)')

pl.colorbar(fig0, ax=ax0)
pl.colorbar(fig1, ax=ax1)
pl.colorbar(fig2, ax=ax2)

pl.savefig('detector_sensitivity', format='eps')

pl.show()
