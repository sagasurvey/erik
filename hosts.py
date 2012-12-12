"""
Defines the hosts for the distant local group project
"""
from __future__ import division

import sys

import numpy as np
from matplotlib import pyplot as plt

SDSS_IMAGE_LIST_URL = 'http://skyserver.sdss3.org/dr9/en/tools/chart/list.asp'


class NSAHost(object):
    """
    A host for targetting extracted from the Nasa-Sloan Atlas.

    Parameters
    ----------
    nsaid : int
        The NSA ID# of this host
    name : str or None
        The name of this host (or None to go with "NSA###")
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

        self.nsaid = nsaid  # The NSA ID #
        self.name = 'NSA{0}'.format(self.nsaid) if name is None else name

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
            self.fnsdss = path.join('catalogs',
                'NSA{0}_sdss.dat'.format(self.nsaid))
        else:
            self.fnsdss = fnsdss

        if fnusnob is None:
            self.fnusnob = path.join('catalogs',
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
        Returns the angular distance (in degrees) given a projected physical distance (in kpc)
        """
        return np.degrees(distkpc / (1000 * self.distmpc)) * 60

    def projected_to_physical(self, angle):
        """
        Returns the projected physical distance (in kpc) given an angular distance (in arcmin)
        """
        if hasattr(angle, 'degrees'):
            angle = angle.degrees

        return np.radians(angle) * 1000 * self.distmpc

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
        from targeting import construct_sdss_query, download_sdss_query

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
        from targeting import construct_usnob_query

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


def sampled_imagelist(ras, decs, n=25, url=SDSS_IMAGE_LIST_URL, copytoclipboard=True):
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
        webbrowser.open(url)

    return text


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


if __name__ == '__main__':
    h1 = NSAHost(76316, 'DLG1')
    h2 = NSAHost(158901, 'DLG2')
    h3 = NSAHost(129387, 'DLG3')
