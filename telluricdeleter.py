import pyfits
import numpy as np
from scipy import signal
from scipy import ndimage
from scipy import interpolate
import os.path
import fileinput
from astropy.modeling import models, fitting

f1=  raw_input('Path to fits file:')
f2=  raw_input('Name of fits file:')


hdu = pyfits.open(f1+'/'+f2)
spec = hdu[0].data.copy()
hdr = hdu[0].header.copy()
hdu.close()

crval = float(hdr['CRVAL1'])
cdelt = float(hdr['CDELT1'])
nlam = float(hdr['NAXIS1'])
lam = np.arange(crval, crval + cdelt * nlam - 1e-4, cdelt)
waves=lam


template_spectrum = signal.savgol_filter(spec, 21, 9)
noise = np.abs(spec - template_spectrum)
noise = ndimage.filters.gaussian_filter1d(noise, 100.0)
not_telluric = np.ones(spec.shape, dtype=np.bool)



#telluricWaves = {'B': (6855, 6935), 'A': (7590, 7685)}
telluricWaves = [(5880,5965),(6830,6900)]
lamA= np.arange(telluricWaves[0][0],telluricWaves[0][1],cdelt)
lamB= np.arange(telluricWaves[1][0],telluricWaves[1][1],cdelt)
#lamC= np.arange(telluricWaves[2][0],telluricWaves[2][1],cdelt)


# For each telluric region
for wavereg in telluricWaves:
        in_telluric_region = np.logical_and(waves >= wavereg[0],
                                            waves <= wavereg[1])
        not_telluric = np.logical_and(not_telluric,
        np.logical_not(in_telluric_region))



#for regionA
region_A = np.logical_and(waves >= telluricWaves[0][0],
                                            waves <= telluricWaves[1][1])
#for regionB
region_B = np.logical_and(waves >= telluricWaves[1][0],
                                            waves <= telluricWaves[1][1])
#for regionC
#region_C = np.logical_and(waves >= telluricWaves[2][0],
#                                            waves <= telluricWaves[2][1])


# noise per region
#A
template_spectrumA = signal.savgol_filter(spec[region_A], 21,3)
noiseA = np.abs(spec[region_A] - template_spectrumA)
noiseA = ndimage.filters.gaussian_filter1d(noiseA, 100.0)
#B
template_spectrumB = signal.savgol_filter(spec[region_B], 21, 3)
noiseB = np.abs(spec[region_B] - template_spectrumB)
noiseB = ndimage.filters.gaussian_filter1d(noiseB, 100.0)
#C
#template_spectrumC = signal.savgol_filter(spec[region_C], 15, 3)
#noiseC = np.abs(spec[region_C] - template_spectrumC)
#noiseC = ndimage.filters.gaussian_filter1d(noiseC, 100.0)


sgspec= signal.savgol_filter(spec,51,7)
m = not_telluric.sum()
intpr = interpolate.splrep(waves[not_telluric], sgspec[not_telluric],w=1/noise[not_telluric], k=3,  s=2*m)



newreg=np.logical_and(waves >= 5850,
                                            waves <= 6862)

smoothedspec = interpolate.splev(waves, intpr)
smoothedspec[not_telluric]= spec[not_telluric]






noiseA = np.abs(spec[region_A] - template_spectrumA)
noiseB = np.abs(spec[region_B] - template_spectrumB)
#noiseC = np.abs(spec[region_C] - template_spectrumC)

smoothedspec[region_A]=smoothedspec[region_A]+noiseA
smoothedspec[region_B]=smoothedspec[region_B]+noiseB
#smoothedspec[region_C]=smoothedspec[region_C]+noiseC

smoothedspec[newreg]=spec[newreg]



hdu[0].data=smoothedspec

output = f1+'/'+'realfinal.fits'



hdu.writeto(output, clobber='True')
