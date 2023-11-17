from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt

template_file = '../data/k_nmf_derived.newdefault.fits'

# open files, template and coefficients
hdulist = fits.open(template_file)

# read in templates
lam	 		= hdulist[11].data
tspec_v0		= hdulist[1].data      # no smoothing
tspec_v0_nl		= hdulist[2].data      # no smoothing without lines
tspec_v0_nd		= hdulist[3].data      # no smoothing without dust
tspec_v0_nd_nl	= hdulist[4].data      # no smoothing without lines and dust
tspec_v300		= hdulist[5].data      # smoothing
tspec_v300_nl	= hdulist[6].data      # smoothing without lines
tspec_v300_nd	= hdulist[7].data      # smoothing without dust 
tspec_v300_nd_nl = hdulist[8].data     # smoothing without line without dust
lspec_v300		= hdulist[9].data      # add lines ?

def plot_template(X, Y, title, log):

	plt.plot(X, Y[0], label='1')
	plt.plot(X, Y[1], label='2')
	plt.plot(X, Y[2], label='3')
	plt.plot(X, Y[3], label='4')
	plt.plot(X, Y[4], label='5')
	if log:
		plt.xscale('log')
		plt.yscale('log')
	plt.xlim(1000, 20000)
	plt.legend(loc='lower right', fontsize=10)
	plt.xlabel('$\lambda$ (A)')
	plt.title(title)
	plt.show()


plt.figure(1, figsize=(12,10))
plt.suptitle('K correct templates')
plot_template(lam, tspec_v0, 'True spectrum', True)





