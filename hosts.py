"""
Defines the hosts for the distant local group project
"""
from __future__ import division

import sys

import numpy as np
from matplotlib import pyplot as plt


class NSAHost(object):
    """
    A host for targetting extracted from the Nasa-Sloan Atlas.

    Parameters
    ----------
    nsaid : int
        The NSA ID# of this host
    name : str or sequence of strings or None
        The name of this host (or None to go with "NSA###").  If a list,
        the first will be the `name` attribute, while all others will be
        in `altname`
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
    def __init__(self, nsaid, name=None, environsradius=-35, fnsdss=None, fnusnob=None):
        from targeting import get_nsa
        from os import path
        from astropy.coordinates import ICRSCoordinates

        self.nsaid = nsaid  # The NSA ID #

        if name is None:
            name = 'NSA{0}'.format(self.nsaid)

        if isinstance(name, basestring):
            self.name = name
            self.altnames = []
        else:
            self.name = name[0]
            self.altnames = name[1:]

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
        self.mstar = obj['MASS']

        galcoord = ICRSCoordinates(self.ra, self.dec, unit=('deg', 'deg')).galactic
        self.l = galcoord.l.degrees
        self.b = galcoord.b.degrees

        for i, band in enumerate('FNugriz'):
            setattr(self, band, obj['ABSMAG'][i])

        if environsradius > 0:
            self.environskpc = environsradius
        else:
            self.environsarcmin = -environsradius

        if fnsdss is None:
            self.fnsdss = path.join('catalogs',
                '{0}_sdss.dat'.format(self.name))
        else:
            self.fnsdss = fnsdss

        if fnusnob is None:
            self.fnusnob = path.join('catalogs',
                '{0}_usnob.dat'.format(self.name))
        else:
            self.fnusnob = fnusnob

        self.sdssquerymagcut = None

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
    def distmod(self):
        """
        Distance modulus (mags)
        """
        from math import log10

        return 5 * log10(self.distmpc * 100000)

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

    def physical_to_projected(self, distkpc):
        """
        Returns the angular distance (in arcmin) given a projected physical distance (in kpc)
        """
        return np.degrees(distkpc / (1000 * self.distmpc)) * 60

    def projected_to_physical(self, angle):
        """
        Returns the projected physical distance (in kpc) given an angular distance (in arcmin)
        """
        if hasattr(angle, 'degrees'):
            angle = angle.degrees

        return np.radians(angle) * 1000 * self.distmpc

    def sdss_environs_query(self, dl=False, usecas=False, magcut=None):
        """
        Constructs an SDSS query to get the SDSS objects around this
        host and possibly downloads the catalog.

        .. note ::
            If you do `usecas`=True with this, be sure to put the result
            in whatever file `fnsdss` points to.

        Parameters
        ----------
        dl : bool
            If True, download the catalog to `fnsdss`.  Otherwise, just
            return the query.
        usecas : bool
            If True, includes an `INTO` in the SQL for use with casjobs.
            Ignored if `dl` is True
        magcut : float or None
            `magcut` as accepted by `targeting.construct_sdss_query` or
            None to use `self.sdssquerymagcut`

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
        from os.path import exists
        from targeting import construct_sdss_query, download_sdss_query

        raddeg = self.environsarcmin / 60.
        usecas = False if dl else usecas
        magcut = self.sdssquerymagcut if magcut is None else magcut

        query = construct_sdss_query(self.ra, self.dec, raddeg,
            into=('{0}_environs'.format(self.name)) if usecas else None,
            magcut=magcut)

        if dl:
            if exists(self.fnsdss):
                print 'File', self.fnsdss, 'exists - not downloading anything.'
            else:
                msg = 'Downloading NSA ID{0} to {1}'.format(self.nsaid, self.fnsdss)
                download_sdss_query(query, fn=self.fnsdss, dlmsg=msg,
                    inclheader='Environs of NSA Object {0}'.format(self.nsaid))
        else:
            return query

    def usnob_environs_query(self, dl=False):
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
        from os.path import exists
        import urllib2
        from targeting import construct_usnob_query

        raddeg = self.environsarcmin / 60.

        usnourl = construct_usnob_query(self.ra, self.dec, raddeg, verbosity=1)

        if dl:
            if exists(self.fnusnob):
                print 'File', self.fnusnob, 'exists - not downloading anything.'
            else:
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
                        #now if there's a group of "S/G" columns, add the appropriate mag suffix
                        for i in range(len(colnames)):
                            if colnames[i] == 'S/G':
                                colnames[i] = colnames[i] + '_' + colnames[i - 1]
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
            from astropy.table import Column

            self._cached_sdss = tab = ascii.read(self.fnsdss, delimiter=',')

            #add UBVRI converted from SDSS mags
            U, B, V, R, I = sdss_to_UBVRI(*[tab[b] for b in 'ugriz'])
            pU, pB, pV, pR, pI = sdss_to_UBVRI(*[tab['psf_' + b] for b in 'ugriz'])

            for b in 'UBVRI':
                tab.add_column(Column(name=b, data=locals()[b]))
            for b in 'UBVRI':
                tab.add_column(Column(name='psf_' + b, data=locals()['p' + b]))


        return self._cached_sdss

    def open_on_nsasite(self):
        """
        Uses `webbrowser` to open the page for this host on the NSA web site.
        """
        import webbrowser

        urltempl = 'http://www.nsatlas.org/getAtlas.html?search=nsaid&nsaID={0}&submit_form=Submit'

        webbrowser.open(urltempl.format(self.nsaid))


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


def load_all_hosts(hostsfile='hosts.dat', existinghosts='globals', usedlgname=False, keyonname=False):
    """
    Loads all the hosts in the specified host file and resturns them
    as a dictionary.

    If `existinghosts` is 'globals', it will be taken from from the global
    variable named 'hosts' if it exists.  Otherwise it is a list with all
    the existing NSA hosts
    """
    d = {}

    if existinghosts == 'globals':
        existinghosts = globals().get('hosts', None)

    ids = []
    if existinghosts is not None:
        for h in existinghosts:
            if isinstance(h, NSAHost):
                ids.append(h.nsaid)

    i = 1
    with open(hostsfile) as f:
        f.readline()  # header
        for l in f:
            while i in ids:
                i += 1
            nsanum, ra, dec, z = l.split()
            nsanum = int(nsanum)
            hnm = 'h' + str(i)
            d[hnm] = NSAHost(nsanum, 'DLG' + str(i) if usedlgname else None)
            if keyonname:
                h = d[hnm]
                d[h.name] = h
                del d[hnm]
            i += 1

    return d

def sdss_to_UBVRI(u, g, r, i, z):
    """
    Converts SDSS ugriz to UBVRI - uses Jordi+ from
    http://www.sdss3.org/dr9/algorithms/sdssUBVRITransform.php

    Returns (U, B, V, R, I)
    """
    UmB   =     (0.52)*(u-g)    + (0.53)*(g-r) - (0.82)
    Bmg   =     (0.313)*(g-r)  + (0.219)
    Vmg   =     (-0.565)*(g-r) - (0.016)
    Rmr   =     (-0.153)*(r-i) - (0.117)
    Imi   =     (-0.386)*(i-z) - (0.397)

    B = Bmg + g
    U = UmB + B
    V = Vmg + g
    R = Rmr + r
    I = Imi + i

    return U, B, V, R, I


def get_old_hosts():
    """
    This loads hosts with the old name convention of "DLG#"
    """
    hostsd = {}

    hostsd['h1'] = NSAHost(76316, 'DLG1')
    hostsd['h2'] = NSAHost(46892, 'DLG2')
    hostsd['h3'] = NSAHost(133120, 'DLG3')
    hostsd['h4'] = NSAHost(156881, 'DLG4')
    hostsd['h5'] = NSAHost(159789, 'DLG5')
    hostsd['h6'] = NSAHost(140594, 'DLG6')

    return hostsd


def get_saga_hosts():
    """
    Returns a dictionary with the "official" SAGA names
    """
    hostsd = {}

    #note that the first name here is for the *host*, not the system, system is the second
    hostsd['odyssey'] = NSAHost(147100, ['Odysseus', 'Odyssey', 'NGC6181'])
    hostsd['iliad'] = NSAHost(150238, ['Achilles', 'Iliad', 'NGC7393'])
    hostsd['lortr'] = NSAHost(155005, ['Frodo Baggins', 'Lord of the Rings', 'NGC895'])
    hostsd['sw'] = NSAHost(155005, ['Luke Skywalker', 'Star Wars', 'NGC895'])

    return hostsd


#this adds the hosts from the get_saga_hosts function to the module's namespace
locals().update(get_saga_hosts())
