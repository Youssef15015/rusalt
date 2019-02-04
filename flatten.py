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
with np.errstate(divide='ignore'):
	np.float64(1.0)/0.0

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



fs = glob('mbxgpP*.fits')
if len(fs) == 0:
    print "WARNING: No images to flat-field."
    # Change directories to fail more gracefully
if not os.path.exists('flts'):
    os.mkdir('flts')
    # Make sure there are science images or arcs and what grating angles were
    # used
scifs, scigas = get_ims(fs, 'sci')
arcfs, arcgas = get_ims(fs, 'arc')

ims = np.append(scifs, arcfs)
gas = np.append(scigas, arcgas)
# For each science and arc image
for i, f in enumerate(ims):
    thishdu = pyfits.open(f)
    ga = gas[i]
    flatfile = 'flats/flt%05.2fres.fits' % (ga)
    if len(glob(flatfile)) == 0:
           if masterflatdir is None:
               print("No flat field image found for %s"% f)
               continue
           # Check for the master flat directory
           flatfile = masterflatdir+'/flt%05.2fres.fits' % (ga)
           if len(glob(flatfile)) == 0:
                   # Still can't find one? Abort!!
                   print("No flat field image found for %s"% f)
                   continue

    # open the corresponding response file
    resphdu = pyfits.open(flatfile)
    # divide out the illumination correction and the flatfield
    # make sure divzero = 0.0
    with np.errstate(divide='ignore'):
        thishdu[1].data /= resphdu[0].data.copy()
    # replace the infinities with 0.0
    thishdu[1].data[np.isinf(thishdu[1].data)] = 0.0
    resphdu.close()

# save the updated file
if f in scifs:
    typestr = 'sci'
else:
    typestr = 'arc'
# get the image number
# by salt naming convention, these should be the last 4 characters
# before the '.fits'
imnum = f[-9:-5]
outname = 'flts/' + typestr + '%05.2fflt%04i.fits' % (float(ga),
                                                             int(imnum))
thishdu.writeto(outname)
thishdu.close()
