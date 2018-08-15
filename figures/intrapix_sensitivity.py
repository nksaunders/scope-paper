import numpy as np
import matplotlib.pyplot as pl
import scope

from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=16, usetex=True)

def intra(x, y, cx, cy):

    xsens = 0; ysens = 0
    for i in range(len(cx)):
        xsens += cx[i] * x ** i
        ysens += cy[i] * y ** i

    '''
    sens = (cx[0] + cx[1] * x + cx[2] * x**2) \
         * (cy[0] + cy[1] * y + cy[2] * y**2)
    '''
    sens = xsens * ysens
    return sens


x = np.arange(-.5,.5,.01)
y = np.arange(-.5,.5,.01)

cx = [1, 0, -.01, 0, -.05, 0, 0, 0, 0]
cy = [1, 0, -.1, 0, .75, 0, -1, 0, -2.5]
#cy = [1, 0, -.01, 0, -.05, 0, 0, 0, 0]

fig = pl.figure(figsize=(10,5))

det = np.zeros((len(x), len(y)))

for i, vx in enumerate(x):
    for j, vy in enumerate(y):
        det[i][j] = intra(vx, vy, cx, cy)

det /= np.max(det)
det /= 100

pl.plot(x, np.sum(det, axis=0))
pl.plot(y, np.sum(det, axis=1))

pl.xlabel('x and y (pixels)')
pl.ylabel('Relative Sensitivity')
pl.savefig('intra_xy_other.pdf', format='pdf', bbox_inches='tight')

pl.close()

fig = pl.figure(figsize=(12.5, 10))

pl.axis('off')
pl.imshow(det * 100, cmap='viridis'); pl.colorbar()
pl.savefig('intra_image_other.pdf', format='pdf', bbox_inches='tight')
