from __future__ import division

"""
These functions are for the "Distant local groups" project target selection.
"""
#important note: SDSS 'type' field: 3=galaxy, 6=star


import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u

NSA_VERSION = '0.1.2'  # used to find the download location/file name
NSAFILENAME = 'nsa_v{0}.fits'.format(NSA_VERSION.replace('.', '_'))

SDSS_SQL_URL = 'http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'
SDSS_IMAGE_LIST_URL = 'http://skyserver.sdss3.org/dr10/en/tools/chart/list.aspx'
SDSS_FINDCHART_URL = 'http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx'

USNOB_URL = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi'

GAMA_URL = 'http://www.gama-survey.org/dr1/data/GamaCoreDR1_v1.csv.gz'
#gama sted limit: r < 19.8

TELL_IF_USING_CACHED = False  # If True, show an informational message about where the NSA is coming from

_cachednsa={}
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
    from urllib2 import urlopen

    from astropy.io import fits
    from hosts import download_with_progress_updates

    if fn is None:
        fn = NSAFILENAME

    if fn in _cachednsa:
        if TELL_IF_USING_CACHED:
            print 'Using cached NSA for file', fn
        return _cachednsa[fn]

    if os.path.exists(fn):
        if TELL_IF_USING_CACHED:
            print 'Loading NSA from local file', fn
    else:
        # download the file if it hasn't been already
        NSAurl = 'http://sdss.physics.nyu.edu/mblanton/v0/' + NSAFILENAME

        with open(fn, 'w') as fw:
            msg = 'Downloading NSA from ' + NSAurl + ' to ' + fn
            u = urlopen(NSAurl)
            try:
                download_with_progress_updates(u, fw, msg=msg)
            finally:
                u.close()

    # use pyfits from astropy to load the data
    res = fits.getdata(fn, 1)
    _cachednsa[fn] = res

    return res


_cachedgama = {}
def get_gama(fn=None):
    """
    Download or load the GAMA survey data

    Parameters
    ----------
    fn : str or None
        A file to load from or None to use the default.
    Returns
    -------
    gamadata
        The data as an astropy `Table`.  Table will also have `decmax`,
        `decmin`, `ramax`, and `ramin`.
    """
    import os
    from urllib2 import urlopen

    from astropy.io import ascii
    from hosts import download_with_progress_updates

    if fn is None:
        fn = os.path.join('catalogs', os.path.split(GAMA_URL)[-1])

    if fn in _cachedgama:
        #print 'Using cached GAMA for file', fn
        return _cachedgama[fn]

    if not os.path.exists(fn):
        with open(fn, 'w') as fw:
            msg = 'Downloading GAMA from ' + GAMA_URL + ' to ' + fn
            u = urlopen(GAMA_URL)
            try:
                download_with_progress_updates(u, fw, msg=msg)
            finally:
                u.close()

    tab = _cachedgama[fn] = ascii.read(fn, delimiter=',', guess=False)

    tab.ramax = np.max(tab['RA_J2000'])
    tab.ramin = np.min(tab['RA_J2000'])
    tab.decmax = np.max(tab['DEC_J2000'])
    tab.decmin = np.min(tab['DEC_J2000'])

    return tab


def construct_sdss_query(ra, dec, radius=1, into=None, magcut=None):
    """
    Generates the query to send to the SDSS to get the full SDSS catalog around
    a target.

    Parameters
    ----------
    ra : float
        The center/host RA in degrees
    dec : float
        The center/host Dec in degrees
    radius : float
        The radius to search out to in degrees
    into : str or None
        The name of the table to construct in your `mydb` if you want to use
        this with CasJobs, or None to have no "into" in the SQL. This also
        adjust other parts of the query a little to work with CasJobs instead
        of the direct query.
    magcut : 2-tuple or None
        if not None, adds a magnitude cutoff.  Should be a 2-tuple
        ('magname', faintlimit). Ignored if None.

    Returns
    -------
    query : str
        The SQL query to send to the SDSS skyserver


    """
    from textwrap import dedent

    query_template = dedent("""
    SELECT  p.objId  as objID,
    p.ra, p.dec, p.type, p.flags, p.specObjID, dbo.fPhotoTypeN(p.type) as phot_sg,
    p.modelMag_u as u, p.modelMag_g as g, p.modelMag_r as r,p.modelMag_i as i,p.modelMag_z as z,
    p.modelMagErr_u as u_err, p.modelMagErr_g as g_err, p.modelMagErr_r as r_err,p.modelMagErr_i as i_err,p.modelMagErr_z as z_err,
    p.psfMag_u as psf_u, p.psfMag_g as psf_g, p.psfMag_r as psf_r, p.psfMag_i as psf_i, p.psfMag_z as psf_z,
    p.fibermag_r, p.fiber2mag_r,
    p.petroMag_r + 2.5*log10(2*PI()*p.petroR50_r*p.petroR50_r) as sb_petro_r,
    p.expMag_r, p.expMag_r + 2.5*log10(2*PI()*p.expRad_r*p.expRad_r + 1e-20) as sb_exp_r,
    p.deVMag_r, p.deVMag_r + 2.5*log10(2*PI()*p.deVRad_r*p.deVRad_r + 1e-20) as sb_deV_r,
    p.lnLExp_r, p.lnLDeV_r, p.lnLStar_r,
    p.extinction_u as Au, p.extinction_g as Ag, p.extinction_r as Ar, p.extinction_i as Ai, p.extinction_z as Az,
    ISNULL(s.z, -1) as spec_z, ISNULL(s.zErr, -1) as spec_z_err, ISNULL(s.zWarning, -1) as spec_z_warn, s.class as spec_class, s.subclass as spec_subclass


    {into}
    FROM {funcprefix}fGetNearbyObjEq({ra}, {dec}, {radarcmin}) n, PhotoPrimary p
    LEFT JOIN SpecObj s ON p.specObjID = s.specObjID
    WHERE n.objID = p.objID{magcutwhere}
    """)
    #if using casjobs, functions need 'dbo' in front of them for some reason
    if into is None:
        intostr = ''
        funcprefix = 'dbo.'
    else:
        intostr = 'INTO MyDB.' + into
        funcprefix = ''

    intostr = '' if into is None else ('INTO MyDB.' + into)

    if magcut is None:
        magcutwhere = ''
    else:
        magcutwhere = ' and p.{0} < {1}'.format(*magcut)

    return query_template.format(ra=float(ra), dec=float(dec),
        radarcmin=radius * 60., into=intostr, funcprefix=funcprefix,
        magcutwhere=magcutwhere)


def construct_usnob_query(ra, dec, radius=1, verbosity=1, votable=False, baseurl=USNOB_URL):
    """
    Generate a USNO-B query for the area around a target.

    Parameters
    ----------
    ra : float
        The center/host RA in degrees
    dec : float
        The center/host Dec in degrees
    radius : float
        The radius to search out to in degrees
    verbosity : int
        The USNO verbosity level
    votable : bool
        If True, query gets a VOTable, otherwise ASCII

    Returns
    -------
    url : str
        The url to query to get the catalog.
    """
    from urllib import urlencode

    parameters = urlencode([('CAT', 'USNO-B1'),
                            ('RA', ra),
                            ('DEC', dec),
                            ('SR', radius),
                            ('VERB', verbosity),
                            ('cftype', 'XML/VO' if votable else 'ASCII'),
                            ('slf', 'ddd.ddd/dd.ddd'),
                            ('skey', 'Mag')])

    return baseurl + '?' + parameters


def download_sdss_query(query, fn=None, sdssurl=SDSS_SQL_URL, format='csv',
                   dlmsg='Downloading...', inclheader=True):
    """
    Runs the provided query on the given SDSS `url`, and either returns the
    result or saves it as a file.

    Parameters
    ----------
    query : str
        The SQL query string.
    fn : str or None
        The filename to save the result to or None to return it from
        this function.
    sdssurl : str
        The URL to send the query to - defaults to whatever
        `SDSS_SQL_URL` is (defined at the top of this file)
    format : str
        The format to return the query.  As far as I know, SDSS only
        accepts 'csv', 'xml', and 'html'
    dlmsg : str or None
        A string to print when the download begins.  If None, there will
        also be no progress updates on the download.
    inclheader : bool or str
        Whether or not to include a header with information about the query
        in the resulting file.  If a string, that will be at the end of the
        header.

    Returns
    -------
    result : str, optional
        If `fn` is None, this will contain the result of the query.

    Raises
    ------
    ValueError
        If the SQL query results in an error or returns no rows

    Notes
    -----
    This way of querying the SDSS has time and # of row limits - if
    you exceed them you'll get an error and

    """
    import urllib2
    import datetime
    from urllib import urlencode
    from StringIO import StringIO

    from hosts import download_with_progress_updates

    parameterstr = urlencode([('cmd', query.strip()), ('format', format)])
    url = sdssurl + '?' + parameterstr

    #either open the requested file or a buffer to later return the values
    if fn is None:
        fw = StringIO()
    else:
        fw = open(fn, 'w')

    try:
        q = urllib2.urlopen(url)
        try:
            #first read the initial two lines to check for errors
            firstline = q.readline()
            secondline = q.readline()
            if 'error' in firstline.lower() or 'error' in secondline.lower():
                rest = q.read()
                raise ValueError('SQL query returned an error:\n' + firstline +
                                 secondline + rest)

            if 'No objects have been found' == firstline:
                raise ValueError('No objects were returned from the request!')

            if inclheader:
                dtstr = str(datetime.datetime.today())
                fw.write('#Retrieved on {0} from {1}\n'.format(dtstr, sdssurl))
                fw.write('#Query:\n#{0}\n'.format(query.strip().replace('\n', '\n#')))
                if isinstance(inclheader, basestring):
                    fw.write('#{0}\n'.format(inclheader.replace('\n', '\n#')))

            fw.write(firstline)
            fw.write(secondline)
            if dlmsg is None:
                fw.write(q.read())
            else:
                download_with_progress_updates(q, fw, msg=dlmsg)

        finally:
            q.close()

        if fn is None:
            # f should be a StringIO object, so we return its value
            return fw.getvalue()
    finally:
        fw.close()

# the color cuts specified in the BOSSANOVA proposal
bossanova_color_cuts = {'g-r': (None, 1.3), 'r-i': (None, 0.7)}

def select_targets(host, band='r', faintlimit=21, brightlimit=15,
    galvsallcutoff=20, inclspecqsos=False, removespecstars=True,
    removegalsathighz=True, removegama='now', photflags=True,
    outercutrad=250, innercutrad=20, colorcuts={},
    randomize=True, removeallsdss=False):
    """
    Selects targets from the SDSS catalog.

    Parameters
    ----------
    hostorcat : NSAHost
        The host object to select targets for
    band : str
        The name of the photometric band to use for magnitude cuts
    faintlimit : number
        The faint cutoff in `band` magnitudes
    brightlimit : number
        The bright cutoff in `band` magnitudes
    galvsallcutoff : number
        The cutoff below which both stars and galaxies are used
    inclspecqsos : bool
        Whether or not to include SDSS spectroscopic targets classified as QSOs
    removespecstars : bool
        Whether or not to ignore SDSS spectroscopic targets classified as stars
    removegalsathighz : bool
        Whether or not to ignore SDSS spectroscopic targets classified as
        galaxies but with redshfits > 3*zerr +z_host
    removegama : str
        If not empty string/False, don't select targets that are already in
        GAMA.  Can be 'all' or 'now' to remove targets that will *eventually* be
        in GAMA, or just those currently observed. Note that the tolerance is
        currently 1 arcsec, which seems fine based on a check of one field.
    photflags : bool
        Apply the extended object recommended photometry flags (see
        http://www.sdss3.org/dr9/tutorials/flags.php)
    outercutrad : Quantity or None
        A separation angle or distance beyond which to not select targets
        or None to not do this cut
    innercutrad : number or None
        A separation angle in *kpc* inside which to not select targets
        (or negative for arcmin), or None to not cut
    colorcuts: dict of 2-tuples
        A dictionary mapping SDSS colors to the range of colors to accept as
        ``(bluest, reddest)``. E.g., {'g-r': (-1, 2)}.  Default is to do no
        color cuts.
    randomize : bool
        Randomize the order of the catalog and the very end
    removeallsdss : bool
        If True, removes *all* SDSS spectra from the target selection

    Returns
    -------
        cat : table
            The SDSS catalog with the selection applied
    """
    from astropy.table import Column, MaskedColumn
    from math import cos, radians

    cat = host.get_sdss_catalog()

    mag = cat[band]

    #raw magnitude cuts
    magcuts = (brightlimit < mag) & (mag < faintlimit)

    #color cuts if any are present
    colorcutmsk = np.ones_like(magcuts)  # start by accepting everything
    if colorcuts:
        for k, v in colorcuts.iteritems():
            c1, c2 = k.split('-')
            color = cat[c1] - cat[c2]

            bluec, redc = v
            if bluec is None:
                bluec = -float('inf')
            if redc is None:
                redc = float('inf')

            assert bluec < redc, 'blue cut larger than red cut!: ' + str(v)

            #this now adds the cuts to the mask for this color
            colorcutmsk = colorcutmsk & (bluec < color) & (color < redc)


    #type==3 is an imaging-classified galaxy - but only do it if you're brighter than galvsallcutoff
    nonphotgal = (cat['type'] == 3) | (mag > galvsallcutoff)

    #base selection is based on the above
    msk = magcuts & colorcutmsk & nonphotgal

    cdec = cos(radians(host.dec))


    if 'rhost' not in cat.colnames:
        dra = cat['ra'] - host.ra
        ddec = cat['dec'] - host.dec
        rhost = ((dra * cdec) ** 2 + ddec ** 2) ** 0.5
        colcls = MaskedColumn if hasattr(rhost, 'mask') else Column
        cat.add_column(colcls(name='rhost', data=rhost))

    #negative for arcmin
    if outercutrad is not None:
        if hasattr(outercutrad, 'unit') and outercutrad.unit.is_equivalent(u.Mpc):
            outercutraddeg = np.degrees(outercutrad.to(u.Mpc).value / host.distmpc)
        elif hasattr(outercutrad, 'unit') and outercutrad.unit.is_equivalent(u.degree):
            outercutraddeg = outercutrad.to(u.degree).value
        elif isinstance(outercutrad, float) or isinstance(outercutrad, int):  # pre-Quantity
            if outercutrad < 0:  # arcmin
                outercutraddeg = -outercutrad / 60.
            else:  # kpc
                outercutraddeg = np.degrees(outercutrad / (1000 * host.distmpc))
        else:
            raise ValueError('Invalid outercutrad')
    else:
        outercutraddeg = cat['rhost'].max()

    outercutrad = cat['rhost'] < outercutraddeg

    msk = msk & outercutrad

    if innercutrad is not None:
        if innercutrad < 0:  # arcmin
            innercutraddeg = -innercutrad / 60.
        else:  # kpc
            innercutraddeg = np.degrees(innercutrad / (1000 * host.distmpc))

        innercutrad = cat['rhost'] > innercutraddeg

        msk = msk & innercutrad

    if photflags:
        flags = cat['flags']
        binned1 = (flags & 0x10000000) != 0  # BINNED1 detection
        nsaturated = (0x0000000000040000 & flags) == 0  # not saturated
        nbce = (0x0000010000000000 & flags) == 0  # not BAD_COUNTS_ERROR
        #photqual = (flags & 0x8100000c00a0) == 0  # not NOPROFILE, PEAKCENTER,
            # NOTCHECKED, PSF_FLUX_INTERP, SATURATED, or BAD_COUNTS_ERROR
        #deblendnopeak = ((flags & 0x400000000000) == 0)  # | (psfmagerr_g <= 0.2)  # DEBLEND_NOPEAK
        msk = msk & binned1 & nsaturated & nbce

    #below are "overrides" rather than selection categories:

    #include SDSS spectroscopy QSOs
    specqsos = cat['spec_class'] == 'QSO'
    msk[specqsos] = inclspecqsos
    if inclspecqsos:
        print 'Found', sum(specqsos), 'QSO candidates'

    if removespecstars:
        specstars = cat['spec_class'] == 'STAR'
        msk[specstars] = False

    if removegalsathighz:
        gals = cat['spec_class'] == 'GALAXY'
        #"high" z means more than 3sigma above the host's redshift
        highzgals = gals & ((cat['spec_z']) > (host.zspec + 3 * host.zdisterr))
        lowzgals = gals & ((cat['spec_z']) <= (host.zspec + 3 * host.zdisterr))
        msk[highzgals] = False

    if removeallsdss:
        sdssspecs = ~cat['spec_class'].mask
        print 'Removing ALL objects with SDSS spec: {0} of {1} objects'.format(sdssspecs.sum(), len(sdssspecs))
        msk[sdssspecs] = False

    if removegama:
        g = get_gama()
        if (host.dec + outercutraddeg > g.decmax or
            host.dec - outercutraddeg < g.decmin or
            host.ra + outercutraddeg > g.ramax or
            host.ra - outercutraddeg < g.ramin):
            pass#print 'Host not in GAMA area - not looking at GAMA'
        else:
            print 'Found host', host.name, 'in GAMA!'
            if removegama == 'all':
                future = True
            elif removegama == 'now':
                future = False
            else:
                raise ValueError('invalid removegama')

            gamamatchmsk = find_gama(cat, host, outercutraddeg, tol=1 / 3600.,
                matchfuture=future)[0]

            msk = msk & ~gamamatchmsk
            print 'Removing', np.sum(gamamatchmsk),'GAMA objects'



    res = cat[msk]
    if randomize:
        res = res[np.random.permutation(len(res))]

    return res


def find_gama(cat, host, raddeg, tol, matchfuture=True):
    """
    Find GAMA objects that match the given ra/decs within a tolerance

    Parameters
    ----------
    catra : record array
        The catalog of objects to match
    host : NSAHost
        The host
    raddeg : array
        The radius of the field in degrees
    tol : float
        The distance (in degrees) of a "close enough" match.
    matchfuture : bool
        If True, include things that are planned for future GAMA releases.
        Otherwise, only accepts things that are currently in GAMA

    Returns
    -------
    msk : bool array
        An array that's true if the catalog entry has a matching GAMA object
    gamacat : Table
        The corresponding GAMA entries
    ds : float array
        On sky-distances (w/o cos(dec)) between GAMA and SDSS
    """
    from scipy import spatial

    catra = cat['ra']
    catdec = cat['dec']
    ra0 = host.ra
    dec0 = host.dec

    g = get_gama()
    gamacoverage = ( (g['RA_J2000'] < (ra0 + raddeg)) &
                     (g['RA_J2000'] > (ra0 - raddeg)) &
                     (g['DEC_J2000'] < (dec0 + raddeg)) &
                     (g['DEC_J2000'] > (dec0 - raddeg)) )

    if matchfuture:
        gamaspec = g['Z_QUALITY'] > 2
    else:
        gamaspec = (g['Z_HELIO'] > -2) & (g['Z_QUALITY'] > 2)

    gm = g[gamacoverage & gamaspec]
    kdt = spatial.cKDTree(np.array([gm['RA_J2000'], gm['DEC_J2000']]).T)

    ds, idx = kdt.query(np.array([catra, catdec]).T)

    msk = ds < tol
    return msk, gm[idx[msk]], ds[msk]


def usno_vs_sdss_offset(sdsscat, usnocat, plots=False, raiseerror=0.5):
    """
    Determines the offset between a USNO-B1 catalog and a sdss catalog covering
    the same fields by matching nearby stars and computing the median.

    Parameters
    ----------
    sdsscat : astropy.table.Table
        Table from the SDSS catalog - must have 'ra' and 'dec'
    usnocat : astropy.table.Table
        Table from the SDSS catalog - must have 'RA' and 'DEC'
    plots : bool
        True to show diagnostic plots
    raiseerror : float or None
        The distance to raise an error if the resulting offset is more than the
        given number of arcmin.

    Returns
    -------
    dra : array
    ddec : array

    Raises
    ------
    ValueError
        If the separation is larger that `raiseerror`
    """
    from scipy.spatial import cKDTree
    from math import cos

    sra = sdsscat['ra']
    sdec = sdsscat['dec']
    ura = usnocat['RA']
    udec = usnocat['DEC']

    kdt = cKDTree(np.array([sra, sdec]).T)

    d, si = kdt.query(np.array([ura, udec]).T)

    dra = ura - sra[si]
    ddec = udec - sdec[si]

    newura = ura - np.median(dra)
    newudec = udec - np.median(ddec)

    d2, si2 = kdt.query(np.array([newura, newudec]).T)

    dra2 = ura - sra[si2]
    ddec2 = udec - sdec[si2]
    d2off = np.hypot(dra2, ddec2)

    cdec = cos(np.radians(np.median(udec)))

    if plots:
        plt.figure()
        bins = np.linspace(0, 3, 200)
        plt.hist(d * 3600, bins=bins, histtype='step')
        plt.hist(d2off * 3600, bins=bins, histtype='step')
        plt.xlabel('d [arcmin]')
        plt.figure()
        plt.plot(dra * 3600, ddec * 3600, '.b', ms=1, alpha=.5)
        plt.plot(dra2 * 3600, ddec2 * 3600, '.r', ms=1, alpha=.5)
        plt.scatter([0], [0], color='k', zorder=2)
        plt.scatter([np.median(dra) * 3600], [np.median(ddec) * 3600], color='r', alpha=.95, zorder=2)
        plt.scatter([np.median(dra2) * 3600], [np.median(ddec2) * 3600], color='g', alpha=.95, zorder=2)
        plt.xlim(-1 / cdec, 1 / cdec)
        plt.ylim(-1, 1)

    dres = np.hypot(np.median(dra2) / cdec, np.median(ddec2))
    if raiseerror is not None and ((dres * 3600) > raiseerror):
        raise ValueError('median separation from USNO to SDSS is {0} arcsec'.format(np.median(d) * 3600))

    return np.median(dra2), np.median(ddec2)


def sampled_imagelist(ras, decs, n=25, names=None, url=SDSS_IMAGE_LIST_URL,
    copytoclipboard=False, posttoimglist=3.):
    """
    Returns the text to be pasted into the sdss image list page.  Also opens
    the page (if `url` is not None) and copies the text to the clipboard if on
    a mac or linux.

    Parameters
    ----------
    ras : array-like
        RA of objects to be marked in degrees
    decs : array-like
        Dec of objects to be marked in degrees
    n : int
        Maximum number of objects (randomly sampled if this is greater than
        `ras` or `decs` length)
    url : str or None
        The URL to the SDSS image list page or None to not open in a web
        browser.
    copytoclipboard : bool
        If True, copies the list of images to the clipboard for use on the SDSS
        web site
    posttoimglist : bool or float
        If True, makes a form to post to the URL site. If a float, gives the
        number of seconds to wait until deleting the temporary file (to gives
        the browser time to load)

    Returns
    -------
    text : str
        The table to be pasted into the image list text box

    """
    import webbrowser
    import platform
    import tempfile
    import time
    import os

    if len(ras) != len(decs):
        raise ValueError('ras and decs not the same size!')

    ras = np.array(ras, copy=False)
    decs = np.array(decs, copy=False)

    idx = None
    if len(ras) > n:
        idx = np.random.permutation(len(ras))[:n]
        ras = ras[idx]
        decs = decs[idx]

    if names is None:
        names = [str(i) for i in range(len(ras))]
    elif idx is not None:
        names = np.array(names, copy=False)[idx]

    text = ['name ra dec']
    for nmi, rai, deci in zip(names, ras, decs):
        text.append('{0} {1} {2}'.format(nmi, rai, deci))
    text = '\n'.join(text)

    if copytoclipboard:
        if platform.system() == 'Darwin':
            clipproc = os.popen('pbcopy', 'w')
            clipproc.write(text)
            clipproc.close()
        elif platform.system() == 'Linux':
            clipproc = os.popen('xsel -i', 'w')
            clipproc.write(text)
            clipproc.close()
        else:
            print("Not on a mac or linux, so can't use clipboard. ")

    if url:
        if posttoimglist:
            page = _imglist_post_templ.format(url=url,text=text)
            tf = tempfile.NamedTemporaryFile(delete=False)
            tf.write(page)
            tf.flush()
            fiurl = 'file://' + os.path.abspath(tf.name)
            webbrowser.open(fiurl)
            if isinstance(posttoimglist, float):
                time.sleep(posttoimglist)
        else:
            webbrowser.open(url)

    return text
_imglist_post_templ = """
<html>
<head>
<title>SDSS sampled_imagelist form</title>
</head>
<body>
<h1> SDSS sampled_imagelist form </h1>
<p>Using URL {url}</p>

<form action="{url}"
method="post">
<TEXTAREA name="paste">
{text}
</TEXTAREA>
<br>
Scale: <input class="in" type="text" value="0.4" name="scale">
<br>
Opt: <input class="in" type="text" value="" name="opt">
<br>
<input type="submit">
</form>
</body>
</html>
"""
