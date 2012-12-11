from __future__ import division

"""
These functions are for the "Distant local groups" project target selection.
"""

import sys

import numpy as np
from matplotlib import pyplot as plt


NSA_VERSION = '0.1.2'  # used to find the download location/file name
NSAFILENAME = 'nsa_v{0}.fits'.format(NSA_VERSION.replace('.', '_'))

SDSS_SQL_URL = 'http://skyserver.sdss3.org/dr9/en/tools/search/x_sql.asp'
SDSS_FINDCHART_URL = 'http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx'
SDSS_IMAGE_LIST_URL = 'http://skyserver.sdss3.org/dr9/en/tools/chart/list.asp'
USNOB_URL = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi'

#SDSS 'type': 3=galaxy, 6=star


def download_with_progress_updates(u, fw, nreports=100, msg=None, outstream=sys.stdout):
    """
    Download a file and give progress updates on the download.

    Parameters
    ----------
    u : result of `urllib2.urlopen`
        The file-like object to read from
    fw : writeable file-like object
        The file object to fill with the content of the download.
    nreports : int
        The number of times to update the download percentage if size available
    msg : str or None
        A message to print when the download stars or None for no message
    outstream : file-like
        The stream to write the updates to
    """
    nreports = int(nreports)
    if 'content-length' in u.headers:
        l = int(u.headers['content-length'])  # bytes
    else:
        l = None
    if msg is not None:
        outstream.write(msg)
        if l is None:
            outstream.write('\nUnknown Size\n')
        else:
            outstream.write('\nSize: {0} kB\n'.format(l / 1024.))
        outstream.flush()
    else:
        outstream.write('\n')  # prepare for the percentage report

    if l is None:
        i = 0
        buf = 'notempty'
        while buf:
            buf = u.read(1024)
            fw.write(buf)
            outstream.write('\r{0} kB downloaded'.format(i + 1))
            outstream.flush()
            i += 1
    else:
        for i in range(nreports):
            fw.write(u.read(int(l / nreports)))
            outstream.write('\r{0}%'.format((i + 1) * 100 / nreports))
            outstream.flush()

    outstream.write('\n')  # leave the 100% part
    outstream.flush()
    #get the little bit that might be left over due to rounding
    end = u.read()
    if end != '':
        fw.write(end)


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

    from astropy.io import fits
    #can also do this if you don't have astropy:
    #import pyfits as fits

    if fn is None:
        fn = NSAFILENAME

    if os.path.exists(fn):
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
    return fits.getdata(fn, 1)


def construct_sdss_query(ra, dec, radius=1, into=None):
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

    Returns
    -------
    query : str
        The SQL query to send to the SDSS skyserver

    """
    from textwrap import dedent

    query_template = dedent("""
    SELECT  p.objId  as objID,
    p.ra, p.dec, p.type, p.flags, p.specObjID, p. fracDeV_r,
    p.modelMag_u as u, p.modelMag_g as g, p.modelMag_r as r,p.modelMag_i as i,p.modelMag_z as z,
    p.modelMagErr_u as u_err, p.modelMagErr_g as g_err, p.modelMagErr_r as r_err,p.modelMagErr_i as i_err,p.modelMagErr_z as z_err,
    p.psfMag_u as psf_u, p.psfMag_g as psf_g, p.psfMag_r as psf_r, p.psfMag_i as psf_i, p.psfMag_z as psf_z,
    p.extinction_u as Au, p.extinction_g as Ag, p.extinction_r as Ar, p.extinction_i as Ai, p.extinction_z as Az,
    ISNULL(s.z, -1) as spec_z, ISNULL(s.zErr, -1) as spec_z_err, ISNULL(s.zWarning, -1) as spec_z_warn, s.class as spec_class, s.subclass as spec_subclass

    {into}
    FROM {funcprefix}fGetNearbyObjEq({ra}, {dec}, {radarcmin}) n, PhotoPrimary p
    LEFT JOIN SpecObj s ON p.specObjID = s.specObjID
    WHERE n.objID = p.objID
    """)
    #if using casjobs, functions need 'dbo' in front of them for some reason
    if into is None:
        intostr = ''
        funcprefix = 'dbo.'
    else:
        intostr = 'INTO mydb.' + into
        funcprefix = ''

    intostr = '' if into is None else ('INTO mydb.' + into)

    return query_template.format(ra=float(ra), dec=float(dec),
        radarcmin=radius * 60., into=intostr, funcprefix=funcprefix)


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



class NSAHost(object):
    """
    A host for targetting extracted from the Nasa-Sloan Atlas.

    Parameters
    ----------
    nsaid : int
        The NSA ID# of this host
    environsradius : float
        The distance to consider as the edge of the "environs" - if positive,
        this will be taken as arcmin, otherwise -kpc
    fnsdss : str or None
        The filename for the SDSS data table.  If None, the default will
        be used.
    fnusnob : str or None
        The filename for the USNO-B data table.  If None, the default
        will be used.

    Attributes
    ----------
    nsaid : int
        The NSA ID# of this host
    nsaindx : int
        The index into the NSA array for this host - *not* the same as `nsaid`
    ra : float
        RA in degrees of the host
    dec : float
        Declination in degrees of the host
    zdist : float
        'ZDIST' field from NSA
    zdisterr : float
        'ZDIST_ERR' field from NSA
    zspec : float
        'Z' field from NSA (`z` is the photometric band)
    F,N,u,g,r,i,z : float
        magnitudes from the NSA

    Notes
    -----
    This assumes the NSA catalog is sorted by NSAID.  This is True as of
    this writing but could in theory change.  In that case the catalog should
    be pre-sorted or something
    """
    def __init__(self, nsaid, environsradius=-35, fnsdss=None, fnusnob=None):
        from os import path

        self.nsaid = nsaid  # The NSA ID #

        nsa = get_nsa()

        # find the object that's in the right order for the requested ID
        self.nsaindx = np.searchsorted(nsa['NSAID'], nsaid)
        obj = nsa[self.nsaindx]

        # make sure its actually the right object
        if obj['NSAID'] != nsaid:
            raise ValueError('NSAID #{0} not present in the catalog'.format(nsaid))

        #now populate various things
        self.ra = obj['RA']
        self.dec = obj['DEC']
        self.zdist = obj['ZDIST']
        self.zdisterr = obj['ZDIST_ERR']
        self.zspec = obj['Z']

        for i, band in enumerate('FNugriz'):
            setattr(self, band, obj['ABSMAG'][i])

        if environsradius > 0:
            self.environskpc = environsradius
        else:
            self.environsarcmin = -environsradius

        if fnsdss is None:
            self.fnsdss = path.join('target_catalogs',
                'NSA{0}_sdss.dat'.format(self.nsaid))
        else:
            self.fnsdss = fnsdss

        if fnusnob is None:
            self.fnusnob = path.join('target_catalogs',
                'NSA{0}_usnob.dat'.format(self.nsaid))
        else:
            self.fnusnob = fnusnob

    @property
    def distmpc(self):
        """
        Distance in Mpc (given WMAP7 cosmology/H0)
        """
        from astropy.cosmology import WMAP7

        return WMAP7.luminosity_distance(self.zdist)

    @property
    def disterrmpc(self):
        """
        Distance error in Mpc (given WMAP7 cosmology/H0)
        """
        from astropy.cosmology import WMAP7

        dist = WMAP7.luminosity_distance(self.zdist)
        dp = abs(dist - WMAP7.luminosity_distance(self.zdist + self.zdisterr))
        dm = abs(dist - WMAP7.luminosity_distance(self.zdist - self.zdisterr))

        return (dp + dm) / 2

    @property
    def environskpc(self):
        """
        The environs radius in kpc
        """
        from math import radians
        return radians(self._environsarcmin / 60.) * self.distmpc * 1000
    @environskpc.setter
    def environskpc(self, val):
        from math import degrees
        self._environsarcmin = degrees(val / (1000 * self.distmpc)) * 60.

    @property
    def environsarcmin(self):
        """
        The environs radius in arcminutes
        """
        return self._environsarcmin
    @environsarcmin.setter
    def environsarcmin(self, val):
        self._environsarcmin = val

    def sdss_environs_query(self, dl=False):
        """
        Constructs an SDSS query to get the SDSS objects around this
        host and possibly downloads the catalog.

        .. note ::
            If you do `usecas`=True with this, be sure to put the result
            in whatever file `fnsdss` points to.

        Parameters
        ----------
        usecas : bool
            If True, download the catalog to `fnsdss`.  Otherwise, just
            return the query.

        Returns
        -------
        query : str
            The SQL query if `dl` is False.  Otherwise, `None` is
            returned.

        Raises
        ------
        ValueError
            if the requested id is not in the catalog.
        """
        raddeg = self.environsarcmin / 60.

        query = construct_sdss_query(self.ra, self.dec, raddeg,
            into=None if dl else ('NSA{0}_environs'.format(self.nsaid)))

        if dl:
            msg = 'Downloading NSA ID{0} to {1}'.format(self.nsaid, self.fnsdss)
            download_sdss_query(query, fn=self.fnsdss, dlmsg=msg,
                inclheader='Environs of NSA Object {0}'.format(self.nsaid))
        else:
            return query

    def usnob_environs_query(self, dl=True):
        """
        Constructs a query to get USNO-B objects around this host, and
        possibly downloads the catalog.

        .. note::
            If `dl` is True, this also fiddles with the catalog a bit to
            make the header easier to read by `astropy.io.ascii`.

        Parameters
        ----------
        dl : bool
            If True, download the catalog

        """
        import urllib2

        raddeg = self.environsarcmin / 60.

        usnourl = construct_usnob_query(self.ra, self.dec, raddeg)

        if dl:
            u = urllib2.urlopen(usnourl)
            try:
                with open(self.fnusnob, 'w') as fw:
                    download_with_progress_updates(u, fw,
                        msg='Downloading USNO-B to ' + self.fnusnob)
            finally:
                u.close()
        else:
            return usnourl

    def get_usnob_catalog(self):
        """
        Loads and retrieves the data for the USNO-B catalog associated with this
        host.

        Returns
        -------
        cat : astropy.table.Table
            The USNO-B catalog
        """
        if getattr(self, '_cached_usnob', None) is None:
            from astropy.io import ascii

            with open(self.fnusnob) as f:
                for l in f:
                    if l.startswith('#1') and 'id' in l:
                        colnames = [nm.strip() for nm in l.replace('#1', '').split('|') if nm.strip() != '']
                        break
                else:
                    raise ValueError('USNO-B catalog does not have header - wrong format?')

            self._cached_usnob = ascii.read(self.fnusnob, names=colnames, guess=False)

        return self._cached_usnob

    def get_sdss_catalog(self):
        """
        Loads and retrieves the data for the SDSS environs catalog associated
        with this host.

        Returns
        -------
        cat : astropy.table.Table
            The SDSS catalog
        """
        if getattr(self, '_cached_sdss', None) is None:
            from astropy.io import ascii

            self._cached_sdss = ascii.read(self.fnsdss, delimiter=',')

        return self._cached_sdss



def mark_findingchart(ras, decs, fnout, arcsecperpix=0.396127, borderpix=15,
                      grid=True, baseurl=SDSS_FINDCHART_URL):
    """
    Downloads an SDSS finding chart marked with the requested locations. The
    edges of the image are determined from the outermost objects to be marked.                                                                                                                                                                                                                                               b

    Parameters
    ----------
    ras : array-like
        RA of objects to be marked in degrees
    decs : array-like
        Dec of objects to be marked in degrees
    fnout : str
        File name to save the image as.  will get '.jpg' added if not given.
    arcsecperpix : float
        Number of arcseconds per pixel for the output image
    borderpix : int
        Number of additional buffer pixels around the edge of the image.
    grid : bool
        If True, the image will have the scale grid/bar.
    baseurl : str
        The URL for the finding chart web site.

    Returns
    -------
    fn : str
        The saved output file name.

    """
    import urllib2
    from urllib import urlencode

    if len(ras) != len(decs):
        raise ValueError('ras and decs not the same size!')

    ras = np.array(ras, copy=False)
    decs = np.array(decs, copy=False)

    rmx = np.max(ras)
    rmn = np.min(ras)
    dmx = np.max(decs)
    dmn = np.min(decs)

    racen = (rmx + rmn) / 2.
    deccen = (dmx + dmn) / 2.

    rarngasec = (rmx - rmn) * 3600
    decrngasec = (dmx - dmn) * 3600
    w = int(np.ceil(arcsecperpix * rarngasec)) + int(2 * borderpix)
    h = int(np.ceil(arcsecperpix * decrngasec)) + int(2 * borderpix)

    qr = ['RA DEC']
    for radeci in zip(ras, decs):
        qr.append('{0} {1}'.format(*radeci))
    qstr = '\n'.join(qr)
    print qstr

    parameters = urlencode([('ra', racen), ('dec', deccen),
                            ('scale', arcsecperpix),
                            ('width', w), ('height', h),
                            ('opt', 'G' if grid else ''),
                            ('query', qstr)])
    url = baseurl + '?' + parameters

    fn = fnout if fnout.endswith('.jpg') else (fnout + '.jpg')

    print 'Getting', w, 'by', h, 'image'
    u = urllib2.urlopen(url)
    try:
        with open(fn, 'w') as fw:
            download_with_progress_updates(u, fw, msg='Downloading to ' + fn)
    finally:
        u.close()

    return fn


def sampled_imagelist(ras, decs, n=25, url=SDSS_IMAGE_LIST_URL):
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

    Returns
    -------
    text : str
        The table to be pasted into the image list text box

    """
    import webbrowser
    import platform
    import os

    if len(ras) != len(decs):
        raise ValueError('ras and decs not the same size!')

    ras = np.array(ras, copy=False)
    decs = np.array(decs, copy=False)

    if len(ras) > n:
        idx = np.random.permutation(len(ras))[:n]
        ras = ras[idx]
        decs = decs[idx]

    text = ['name ra dec']
    for i, (rai, deci) in enumerate(zip(ras, decs)):
        text.append('{0} {1} {2}'.format(i, rai, deci))
    text = '\n'.join(text)

    if url:
        if platform.system() == 'Darwin':
            clipproc = os.popen('pbcopy', 'w')
            clipproc.write(text)
            clipproc.close()
            webbrowser.open(url)
        elif platform.system() == 'Linux':
            clipproc = os.popen('xsel -i', 'w')
            clipproc.write(text)
            clipproc.close()
        else:
            print ("Not on a mac or linux, so can't use clipboard. "
                   " Instead, returning the query and you can do what "
                   "you want with it at", url)
    return text

def select_targets(host, band='r', faintlimit=21, brightlimit=15,
    galvsallcutoff=19, inclspecqsos=False, removespecstars=True,
    removegalsathighz=True, photflags=True, rcut=300):
    """
    Selects targets from the SDSS catalog.

    Parameters
    ----------
    host : NSAHost
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
    photflags : bool
        Apply the extended object recommended photometry flags (see
        http://www.sdss3.org/dr9/tutorials/flags.php)
    rcut : number or None
        A separation angle in kpc beyond which to not select targets
        (or negative for arcmin), or None to use the whole catalog

    Returns
    -------
        cat : table
            The SDSS catalog with the selection applied
    """
    cat = host.get_sdss_catalog()

    mag = cat[band]

    #raw magnitude cuts
    magcuts = (brightlimit < mag) & (mag < faintlimit)

    #type==3 is an imaging-classified galaxy - but only do it if you're brighter than galvsallcutoff
    nonphotgal = (cat['type'] == 3) | (mag > galvsallcutoff)

    #TODO: APPLY FLAGS!
    flagcuts = np.ones_like(nonphotgal)

    #base selection is based on the above
    msk = magcuts & nonphotgal & flagcuts

    if rcut is not None:
        dra = cat['ra'] - host.ra
        ddec = cat['dec'] - host.dec
        rhost = (dra ** 2 + ddec ** 2) ** 0.5

        if rcut < 0:  # arcmin
            rcutdeg = -rcut / 60.
        else:  # kpc
            rcutdeg = np.degrees(rcut / (1000 * host.distmpc))

        rcut = rhost < rcutdeg

        msk = msk & rcut

    if photflags:
        flags = cat['flags']
        binned1 = (flags & 0x10000000) != 0  # BINNED1 detection
        photqual = (flags & 0x8100000c00a0) == 0  # not NOPROFILE, PEAKCENTER,
            # NOTCHECKED, PSF_FLUX_INTERP, SATURATED, or BAD_COUNTS_ERROR
        deblendnopeak = ((flags & 0x400000000000) == 0)  # | (psfmagerr_g <= 0.2)  # DEBLEND_NOPEAK
        msk = msk & binned1 & photqual  # & deblendnopeak

    #below are "overrides" rather than selection categories:

    #include SDSS spectroscopy QSOs
    specqsos = cat['spec_class'] == 'QSO'
    msk[specqsos] = inclspecqsos

    if removespecstars:
        specstars = cat['spec_class'] == 'STAR'
        msk[specstars] = False

    if removegalsathighz:
        gals = cat['spec_class'] == 'GALAXY'
        #"high" z means more than 3sigma above the host's redshift
        highzgals = gals & ((cat['spec_z']) > (host.zspec + 3 * host.zdisterr))
        msk[highzgals] = False

    return cat[msk]


def select_fops(host, faintlimit=14, brightlimit=10):
    """
    Selects candidate FOP stars from USNO-B

    Parameters
    ----------
    host : NSAHost
    faintlimit : number
    brightlimit : number

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

    return cat[bothRs & magcuts]


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
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)

    dres = np.hypot(np.median(dra2), np.median(ddec2))
    if raiseerror is not None and ((dres * 3600) > raiseerror):
        raise ValueError('median separation from USNO to SDSS is {0} arcsec'.format(np.median(d) * 3600))

    return np.median(dra2), np.median(ddec2)


def select_sky_positions(host, nsky=100):
    """
    Produces sky positions uniformly covering a circle centered at the host

    Parameters
    ----------
    host : NSAHost
    nsky : int
        Number of sky positions to generate

    Returns
    -------
    ra : array
    dec : array
    """
    raddeg = host.environsarcmin / 60.

    rs = raddeg * 2 * np.arccos(np.random.rand(nsky)) / np.pi
    thetas = 2 * np.pi * np.random.rand(nsky)

    return (host.ra + rs * np.sin(thetas)), (host.dec + rs * np.cos(thetas))


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


def construct_whydra_file(fnout, host, lst, texp=1.5, wl=7000, obsdatetime=None, objcat=None, fopcat=None, skyradec=None):
    import time

    from astropy.time import Time

    if obsdatetime is None:
        obsdatetime = Time(time.time(), format='unix', scale='utc')

    if objcat is None:
        objcat = select_targets(host)
    if fopcat is None:
        fopcat = select_fops(host)
    if skyradec is None:
        skyradec = select_sky_positions(host)

    if len(objcat) > 2000:
        raise ValueError('whydra cannot handle > 2000 objects')
    if len(fopcat) > 2000:
        raise ValueError('whydra cannot handle > 2000 FOPS')
    if len(skyradec[0]) > 2000:
        raise ValueError('whydra cannot handle > 2000 sky points')

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
        fw.write('WAVELENGTH: {0}\n'.format(int(wl)))
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
