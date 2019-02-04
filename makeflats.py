global iraf
from pyraf import iraf
import numpy as np
import pyfits
from glob import glob
import os
iraf.pysalt()
iraf.saltspec()
iraf.saltred()
iraf.set(clobber='YES')
iraf.noao()
iraf.twodspec()
iraf.longslit()


def tofits(filename, data, hdr=None, clobber=False):
    """simple pyfits wrapper to make saving fits files easier."""
    from pyfits import PrimaryHDU, HDUList
    hdu = PrimaryHDU(data)
    if hdr is not None:
        hdu.header = hdr
    hdulist = HDUList([hdu])
    hdulist.writeto(filename, clobber=clobber, output_verify='ignore')


def get_ims(fs, imtype):
    imtypekeys = {'sci': 'OBJECT', 'arc': 'ARC', 'flat': 'FLAT'}
    ims = []
    grangles = []
    for f in fs:
        if pyfits.getval(f, 'OBSTYPE') == imtypekeys[imtype]:
            ims.append(f)
            grangles.append(pyfits.getval(f, 'GR-ANGLE'))
    return np.array(ims), np.array(grangles)

def get_scis_and_arcs(fs):
    scifs, scigas = get_ims(fs, 'sci')
    arcfs, arcgas = get_ims(fs, 'arc')

    ims = np.append(scifs, arcfs)
    gas = np.append(scigas, arcgas)
    return ims, gas



fs = glob('mbxgp*.fits')
if len(fs) == 0:
        print "WARNING: No flat-fields to combine and normalize."
        # Fail gracefully by going up a directory
# make a flats directory
if not os.path.exists('flats'):
        os.mkdir('flats')

# Figure out which images are flats and which grating angles were used
allflats, grangles = get_ims(fs, 'flat')

for ga in np.unique(grangles):
        # grab the flats for this gr angle
        flats = allflats[grangles == ga]
        # run imcombine with average and sigclip, weighted by exposure time
        flatlist = ''
        for f in flats:
        	flatlist += '%s[%i],' % (f,1)
        	# Add the exptime keyword to each extension
       		pyfits.setval(f, 'EXPTIME', ext=1,
                  	value=pyfits.getval(f, 'EXPTIME'))

       		 # set the output combined file name
       	combineoutname = 'flats/flt%05.2fcom.fits' % (ga)
        if os.path.exists(combineoutname):
        	os.remove(combineoutname)
        # initialize the iraf command
        iraf.unlearn(iraf.imcombine)
        print(flatlist)
        # don't forget to remove the last comma in the filelist
        iraf.imcombine(input=flatlist[:-1], output=combineoutname,
               combine='average', reject='sigclip', lsigma=3.0,
               hsigma=3.0, weight='exposure', expname='EXPTIME')

        pyfits.setval(combineoutname, 'DISPAXIS', value=1)
        # We want to make an illumination correction file
        # before running response:
        illumoutname = 'flats/flt%05.2fill.fits' % (ga)
        iraf.unlearn(iraf.illumination)
        iraf.illumination(images=combineoutname,
                                  illuminations=illumoutname, interactive=False,
                                  naverage=-40, order=11, low_reject=3.0,
                                  high_reject=3.0, niterate=5, mode='hl')

        # Flag any pixels in the illumination correction< 0.1
        illumhdu = pyfits.open(illumoutname, mode='update')
        illumhdu[0].data[illumhdu[0].data <= 0.1] = 0.0
        illumhdu[0].data[np.isnan(illumhdu[0].data)]=0.0
        illumhdu.flush()

        # Get 40 pixels out of the middle of the image and
        # median them to run response
        combinehdu = pyfits.open(combineoutname)
        ny = combinehdu[0].data.shape[0]
	flat1dtop = combinehdu[0].data[ny / 2 - 21: ny / 2 + 20, :]
	flat1dbot = illumhdu[0].data[ny / 2 - 21: ny / 2 + 20, :]
	flat1dbotsig = np.std(flat1dbot - 1.0)
        flat1dbot[abs(flat1dbot - 1.0) > 5.0 * flat1dbotsig] = 1.0
        # divide out the illumination correction before running response
        flat1d = np.median(flat1dtop
                           / flat1dbot,
                           axis=0)
        # close the illumination file because we don't need it anymore
        illumhdu.close()

        # File stage m1d for median 1-D
        flat1dfname = 'flats/flt%05.2fm1d.fits' % (ga)
        tofits(flat1dfname, flat1d, hdr=combinehdu[0].header.copy())

        # run response
        # r1d = response1d
        resp1dfname = 'flats/flt%05.2fr1d.fits' % (ga)
        iraf.response(flat1dfname, flat1dfname, resp1dfname, order=31,
                     interactive=False, naverage=-5, low_reject=3.0,
                     high_reject=3.0, niterate=5, mode='hl')

        resp1dhdu = pyfits.open(resp1dfname)
        resp1d = resp1dhdu[0].data.copy()
        resp1dhdu.close()

        # After response divide out the response function
        # normalize the 1d resp to its median
        resp1d /= np.median(resp1d)

        # Chuck any outliers
        flatsig = np.std(resp1d - 1.0)
        resp1d[abs(resp1d - 1.0) > 5.0 * flatsig] = 1.0
        resp = flat1d / resp1d

        resp2dfname = 'flats/flt%05.2fres.fits' % (ga)
        resp2d = combinehdu[0].data.copy() / resp
        tofits(resp2dfname, resp2d, hdr=combinehdu[0].header.copy())
        combinehdu.close()
        # close the combined flat because we don't need it anymore
        combinehdu.close()

        pyfits.setval(resp2dfname, 'DISPAXIS', value=1)

        # Reset any pixels in the flat field correction< 0.1
        # We could flag bad pixels here if we want, but not right now
        flathdu = pyfits.open(resp2dfname, mode='update')
        flathdu[0].data[flathdu[0].data <= 0.1] = 0.0
        flathdu[0].data[np.isnan(flathdu[0].data)]=0.0
        flathdu.flush()
        flathdu.close()


