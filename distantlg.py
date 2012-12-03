from __future__ import division

import numpy as np

NSA_VERSION = '0.1.2'  # used to find the download location/file name
NSAFILENAME = 'nsa_v{0}.fits'.format(NSA_VERSION.replace('.', '_'))

def get_nsa(fn=None):
    """
    Download the NASA Sloan Atlas if it hasn't been already, open it, and
    return the data.

    Parameters
    ----------
    fn : str or None
        The name of the file to load (or to save as if its not present).  If
        None, the convention from the NSA web site will be used.

    Returns
    -------
    nsadata
        The data as an `astropy.io.fits` record array.
    """
    import os
    import sys
    from urllib2 import urlopen

    #from astropy.io import fits
    import pyfits as fits


    if fn is None:
        fn = NSAFILENAME

    if os.path.exists(fn):
        print 'Loading NSA from local file', fn
    else:
        # download the file if it hasn't been already
        NSAurl = 'http://sdss.physics.nyu.edu/mblanton/v0/' + NSAFILENAME

        with open(fn, 'w') as fw:
            u = urlopen(NSAurl)
            try:
                l = int(u.info()['content-length'])  # bytes
                print 'Downloading NSA from', NSAurl, 'to', fn, '\nSize:', l * 1024 ** -2, 'MB'

                for i in range(100):
                    fw.write(u.read(int(l / 100.)))
                    sys.stdout.write('\r{0}%'.format(i + 1))
                    sys.stdout.flush()

                sys.stdout.write('\n')  # leave the 100% part
                sys.stdout.flush()
                #get the little bit that might be left over due to rounding of l/100
                end = u.read()
                if end != '':
                    fw.write(end)

            finally:
                u.close()

    # use pyfits from astropy to load the data
    return fits.getdata(fn, 1)