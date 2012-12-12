from __future__ import division

"""
These functions are for the "Distant local groups" project WIYN-related work.
"""

import numpy as np


def select_fops(host, faintlimit=14, brightlimit=12.5, randomize=True):
    """
    Selects candidate FOP stars from USNO-B

    Parameters
    ----------
    host : NSAHost
    faintlimit : number
    brightlimit : number
    randomize : bool
        Randomize the order of the catalog and the very end

    Returns
    -------
        cat : table
            The USNO-B catalog with the selection applied
    """

    cat = host.get_usnob_catalog()

    mag = cat['R2']

    magcuts = (brightlimit < mag) & (mag < faintlimit)

    # only take things with *both* R mags
    bothRs = (cat['R1'] != 0) & (cat['R2'] != 0)

    res = cat[bothRs & magcuts]
    if randomize:
        res = res[np.random.permutation(len(res))]

    return res


def select_sky_positions(host, nsky=500, sdsscat=None, usnocat=None, nearnesslimitarcsec=15):
    """
    Produces sky positions uniformly covering a circle centered at the host

    Parameters
    ----------
    host : NSAHost
    sdsscat
    usnocat
    nsky : int
        Number of sky positions to generate

    Returns
    -------
    ra : array
    dec : array
    """
    from scipy.spatial import cKDTree

    if sdsscat is None:
        sdsscat = host.get_sdss_catalog()
    if usnocat is None:
        usnocat = host.get_usnob_catalog()

    neardeg = nearnesslimitarcsec / 3600.

    skdt = cKDTree(np.array([sdsscat['ra'], sdsscat['dec']]).T)
    ukdt = cKDTree(np.array([usnocat['RA'], usnocat['DEC']]).T)

    raddeg = host.environsarcmin / 60.

    ras = np.array([])
    decs = np.array([])

    i = -1
    while len(ras) < nsky:
        i += 1

        rs = raddeg * 2 * np.arccos(np.random.rand(nsky)) / np.pi
        thetas = 2 * np.pi * np.random.rand(nsky)

        ra = host.ra + rs * np.sin(thetas)
        dec = host.dec + rs * np.cos(thetas)

        dsdss = skdt.query(np.array([ra, dec]).T)[0]
        dusno = ukdt.query(np.array([ra, dec]).T)[0]

        msk = (dsdss > neardeg) & (dusno > neardeg)

        ras = np.append(ras, ra[msk])
        decs = np.append(decs, dec[msk])

        if i > 100:
            raise ValueError('Could not produce {nsky} sky positions after {i} iterations!'.format(nsky=nsky, i=i))

    return ras[:nsky], decs[:nsky]


def construct_whydra_file(fnout, host, lst, texp=1.5, wl=7000, obsdatetime=None, objcat=None, fopcat=None, skyradec=None):
    import time

    from astropy.time import Time

    from targeting import usno_vs_sdss_offset, select_targets

    if obsdatetime is None:
        obsdatetime = Time(time.time(), format='unix', scale='utc')

    if objcat is None:
        objcat = select_targets(host)
        print len(objcat), 'objects'
    if fopcat is None:
        fopcat = select_fops(host)
        print len(fopcat), 'FOPS'
    if skyradec is None:
        skyradec = select_sky_positions(host)

    if len(objcat) > 2000:
        print('whydra cannot handle > 2000 objects, truncating')
        objcat = objcat[:1999]
    if len(fopcat) > 2000:
        print('whydra cannot handle > 2000 FOPS, truncating')
        fopcat = fopcat[:1999]
    if len(skyradec[0]) > 2000:
        print('whydra cannot handle > 2000 sky points, truncating')
        skyradec = skyradec[0][:1999], skyradec[1][:1999]

    #determine the SDSS vs. USNO offset
    dra, ddec = usno_vs_sdss_offset(host.get_sdss_catalog(), host.get_usnob_catalog())
    print 'USNO/SDSS offsets:', dra * 3600, ddec * 3600

    with open(fnout, 'w') as fw:
        #first do the header
        fw.write('FIELD NAME: NSA{0}\n'.format(host.nsaid))
        fw.write('INPUT EPOCH: 2000.00\n')
        fw.write('CURRENT EPOCH: {0:.1f}\n'.format(obsdatetime.jyear))
        fw.write('SIDEREAL TIME: {0:.2f}\n'.format(lst))
        fw.write('EXPOSURE LENGTH: {0:.2f}\n'.format(texp))
        fw.write('WAVELENGTH: {0:f}\n'.format(int(wl)))
        fw.write('CABLE: RED\n')
        fw.write('WEIGHTING: WEAK\n')

        #first add host as center, and as object
        fw.write(_whydra_file_line(0001, 'Center'.format(host.nsaid), host.ra, host.dec, 'C'))
        fw.write('\n')

        fw.write(_whydra_file_line(1000, 'NSA{0}'.format(host.nsaid), host.ra, host.dec, 'O'))
        fw.write('\n')

        i = 2000
        for obj in objcat:
            ln = _whydra_file_line(i, 'SDSS', obj['ra'], obj['dec'], 'O')
            i += 1
            fw.write(ln)
            fw.write('\n')

        i = 5000
        for fop in fopcat:
            ln = _whydra_file_line(i, 'USNO' + fop['id'], fop['RA'] - dra, fop['DEC'] - ddec, 'F')
            fw.write(ln)
            i += 1
            fw.write('\n')

        i = 8000
        j = 1
        for skyra, skydec in zip(*skyradec):
            ln = _whydra_file_line(i, 'sky{0}'.format(j), skyra, skydec, 'S')
            fw.write(ln)
            i += 1
            j += 1
            fw.write('\n')


def _whydra_file_line(i, name, radeg, decdeg, code):
    from warnings import warn
    from astropy.coordinates import Angle
    from astropy.units import degree, hour

    i = int(i)
    if i > 9999:
        raise ValueError('i too large: ' + str(i))

    if len(name) > 20:
        warn('Name {0} too long - trimming'.format(name))
        name = name[:20]

    if code not in 'COSFE':
        raise ValueError('invalid whydra line code ' + code)

    rastr = Angle(radeg, degree).format(hour, sep=':', pad=True, precision=3)
    decstr = Angle(decdeg, degree).format(degree, sep=':', pad=True, precision=2, alwayssign=True)

    if name == 'SDSS':
        name = 'J' + rastr[:-1].replace(':', '') + decstr[:-1].replace(':', '')

    return '{i:04} {name:20} {ra} {dec} {code}'.format(i=i, name=name, ra=rastr,
        dec=decstr, code=code)
