"""
Defines the hosts for the distant local group project
"""
from __future__ import division, print_function

import os
import sys

import numpy as np
from matplotlib import pyplot as plt

try:
    import six
except ImportErrir:
    from astropy.extern import six
urlopen = six.moves.urllib.request.urlopen

from astropy import units as u
from astropy.coordinates import SkyCoord


NSA_VERSION = '0.1.2'  # used to find the download location/file name
NSAFILENAME = 'nsa_v{0}.fits'.format(NSA_VERSION.replace('.', '_'))

SDSS_SQL_URL = 'http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'

USNOB_URL = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi'

TELL_IF_USING_CACHED = False  # If True, show an informational message about where the NSA is coming from


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
    environsradius : Quantity
        The distance (radius) to consider as the edge of the "environs" - can be
        given as an angle, or as a physical distance.
    fnsdss : str or None
        The filename for the SDSS data table.  If None, the default will
        be used.
    fnusnob : str or None
        The filename for the USNO-B data table.  If None, the default
        will be used.
    shortname : str or None
        Sets the `shortname` attribute (described below) - if None, a shortname
        will be created from a shortened version of the `name` attribute.

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
    shortname : str
        A name for the field that's 6 characters or less (needed for
        spectrographs like IMACS that have silly character size requirements).

    Notes
    -----
    This assumes the NSA catalog is sorted by NSAID.  This is True as of
    this writing but could in theory change.  In that case the catalog should
    be pre-sorted or something
    """
    def __init__(self, nsaid, name=None, environsradius=300*u.kpc,
                 fnsdss=None, fnusnob=None, shortname=None):
        from os import path
        from astropy.coordinates import ICRS, Galactic

        self.nsaid = nsaid  # The NSA ID #

        nsaname = 'NSA{0}'.format(self.nsaid)
        if name is None:
            self.name = nsaname
            self.altnames = []
        elif isinstance(name, six.string_types):
            self.name = name
            self.altnames = [nsaname]
        else:
            self.name = name[0]
            self.altnames = name[1:]
            if nsaname not in self.altnames:
                self.altnames.append(nsaname)

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

        galcoord = ICRS(self.ra*u.deg, self.dec*u.deg).transform_to(Galactic)
        self.l = galcoord.l.degree
        self.b = galcoord.b.degree

        for i, band in enumerate('FNugriz'):
            setattr(self, band, obj['ABSMAG'][i])

        if hasattr(environsradius, 'unit') and environsradius.unit.is_equivalent(u.kpc):
            self.environskpc = environsradius.to(u.kpc).value
        elif hasattr(environsradius, 'unit') and environsradius.unit.is_equivalent(u.arcmin):
            self.environsarcmin = environsradius.to(u.arcmin).value
        elif isinstance(environsradius, float) or isinstance(environsradius, int): #pre-Quantity behavior
            if environsradius > 0:
                self.environskpc = environsradius
            else:
                self.environsarcmin = -environsradius
        else:
            raise ValueError('invalid environsradius')

        if fnsdss is None:
            self.fnsdss = path.join('catalogs',
                '{0}_sdss.dat'.format(self.name))
            self.altfnsdss = [path.join('catalogs',
                '{0}_sdss.dat'.format(nm)) for nm in self.altnames]
        else:
            self.fnsdss = fnsdss
            self.altfnsdss = []

        if fnusnob is None:
            self.fnusnob = path.join('catalogs',
                '{0}_usnob.dat'.format(self.name))
            self.altfnusnob = [path.join('catalogs',
                '{0}_usnob.dat'.format(nm)) for nm in self.altnames]
        else:
            self.fnusnob = fnusnob
            self.altfnusnob = []

        self.sdssquerymagcut = None
        self.shortname = shortname

    @property
    def dist(self):
        """
        Distance to host (given WMAP7 cosmology/H0)
        """
        from astropy.cosmology import WMAP7

        return WMAP7.luminosity_distance(self.zdist)

    @property
    def disterr(self):
        """
        Distance error (given WMAP7 cosmology/H0)
        """
        from astropy.cosmology import WMAP7

        dist = WMAP7.luminosity_distance(self.zdist)
        dp = abs(dist - WMAP7.luminosity_distance(self.zdist + self.zdisterr))
        dm = abs(dist - WMAP7.luminosity_distance(self.zdist - self.zdisterr))

        return (dp + dm) / 2

    @property
    def distmpc(self):
        """
        `dist` in Mpc (for backwards-compatibility)
        """
        return self.dist.to(u.Mpc).value

    @property
    def disterrmpc(self):
        """
        `disterr` in Mpc (for backwards-compatibility)
        """
        return self.disterr.to(u.Mpc).value

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

    @property
    def coords(self):
        """
        The host's coordinates as an ICRS `SkyCoord` object.
        """
        return SkyCoord(self.ra*u.deg, self.dec*u.deg, frame='icrs')

    @property
    def shortname(self):
        if self._shortname is None:
            return self.name.replace('NGC', 'N')[:6]
        else:
            return self._shortname
    @shortname.setter
    def shortname(self, value):
        if value is not None and len(value) > 6:
            raise ValueError('shortname must be less than 6 characters!')
        self._shortname = value



    def physical_to_projected(self, dist):
        """
        Returns the angular distance (in arcmin) given a projected physical
        distance. `distkpc` must be a Quantity.
        E.g., ``physical_to_projected(30*u.arcmin)``
        """
        if u.kpc.is_equivalent(dist):
            dist = dist.to(u.Mpc).value
        else:
            raise ValueError('need to give physical_to_projected a Quantity.')

        return np.degrees(dist / self.distmpc) * 60 * u.arcmin

    def projected_to_physical(self, angle):
        """
        Returns the projected physical distance (in kpc) given an angular
        distance. `angle` must be a Quantity.
        E.g., ``projected_to_physical(30*u.arcmin)``
        """
        if u.degree.is_equivalent(angle):
            angle = angle.to(u.degree).value
        else:
            raise ValueError('need to give projected_to_physical a Quantity.')

        return np.radians(angle) * 1000 * self.distmpc * u.kpc

    def sdss_environs_query(self, dl=False, usecas=False, magcut=None,
                            inclphotzs=True, applyphotflags=False,
                            xmatchwise=False, usepost=False):
        """
        Constructs an SDSS query to get the SDSS objects around this
        host and possibly downloads the catalog.

        .. note ::
            If you do `usecas`=True with this, be sure to put the result
            in whatever file `fnsdss` points to.

        Parameters
        ----------
        dl : bool or str
            If True, download the catalog to `fnsdss`.  Otherwise, just
            return the query.  If 'overwrite', will overwrite `fnsdss` even if
            it already exists.
        usecas : bool
            If True, includes an `INTO` in the SQL for use with casjobs.
            Ignored if `dl` is True
        magcut : float or None
            `magcut` as accepted by `construct_sdss_query` or
            None to use `self.sdssquerymagcut`
        inclphotzs : bool
            If True, includes a further join to add phot-zs where present
        applyphotflags : bool
            If True, adds photometry flags to the query (See
            `construct_sdss_query` for exact flags)
        xmatchwise : bool
            If True, the query will also cross-match against WISE and include w1
            mags.
        usepost : bool
            If False, uses GET, otherwise POST

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

        raddeg = self.environsarcmin / 60.
        usecas = False if dl else usecas
        magcut = self.sdssquerymagcut if magcut is None else magcut

        intoname = self.name.replace(' ', '_') + ('_wwise' if xmatchwise else '')

        query = construct_sdss_query(self.ra, self.dec, raddeg,
            into=('{0}_environs'.format(intoname)) if usecas else None,
            magcut=magcut, inclphotzs=inclphotzs, applyphotflags=applyphotflags,
            xmatchwise=xmatchwise)

        if dl:
            altfns = [self.fnsdss]
            altfns.extend(self.altfnsdss)
            for fn in altfns:
                if exists(fn):
                    if dl == 'overwrite':
                        print('File', fn, 'exists, overwriting with new download.')
                    else:
                        print('File', fn, 'exists - not downloading anything.')
                        break
            else:
                msg = 'Downloading NSA ID{0} to {1}'.format(self.nsaid, self.fnsdss)
                download_sdss_query(query, fn=self.fnsdss, dlmsg=msg, usepost=usepost,
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

        raddeg = self.environsarcmin / 60.

        usnourl = construct_usnob_query(self.ra, self.dec, raddeg, verbosity=1)

        if dl:
            altfns = [self.fnusnob]
            altfns.extend(self.altfnusnob)
            for fn in altfns:
                if exists(fn):
                    print('File', fn, 'exists - not downloading anything.')
                    break
            else:
                u = urlopen(usnourl)
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
            from os.path import exists
            from astropy.io import ascii

            if exists(self.fnusnob):
                fn = self.fnusnob
            else:
                for fn in self.altfnusnob:
                    if exists(fn):
                        print('Could not find file "{0}" but did find "{1}" so using that.'.format(self.fnusnob, fn))
                        break
                else:
                    #didn't find one
                    raise IOError('Could not find file {0} nor any of {1}'.format(self.altfnusnob, self.altfnusnob))

            with open(fn) as f:
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

            self._cached_usnob = ascii.read(fn, names=colnames, guess=False, Reader=ascii.NoHeader)

        return self._cached_usnob

    def get_sdss_catalog(self):
        """
        Loads and retrieves the data for the SDSS environs catalog associated
        with this host.

        Note that this automatically converts all-upper tables to "mixed-case"
        for backwards-compatibility.

        Returns
        -------
        cat : astropy.table.Table
            The SDSS catalog
        """
        from os.path import exists

        if getattr(self, '_cached_sdss', None) is None:

            if exists(self.fnsdss):
                fn = self.fnsdss
            else:
                for fn in self.altfnsdss:
                    if exists(fn):
                        print('Could not find file "{0}" but did find "{1}" so using that.'.format(self.fnsdss, fn))
                        break
                else:
                    #didn't find one
                    raise IOError('Could not find file {0} nor any of {1}'.format(self.fnsdss, self.altfnsdss))
            self._cached_sdss = self.load_and_reprocess_sdss_catalog(fn)
        return self._cached_sdss


    catalog_cases_to_convert = ['ra', 'dec', 'rhost', 'type', 'phot_sg', 'rhost_kpc', 'objID']
    for band in 'ugriz':
        catalog_cases_to_convert.append(band)
        catalog_cases_to_convert.append('fibermag_' + band)
    catalog_aliases = {nm.upper():nm for nm in catalog_cases_to_convert}
    catalog_aliases['RHOST_ARCM'] = 'rhost'
    catalog_aliases['PHOTPTYPE'] = 'type'
    for band in 'ugriz':
        catalog_aliases['EXTINCTION_' + band.upper()] = 'A' + band
        catalog_aliases['PSFMAG_' + band.upper()] = 'psf_' + band
    del band, catalog_cases_to_convert # clean up the namespace


    def load_and_reprocess_sdss_catalog(self, fn):
        from astropy.io import ascii, fits
        from astropy.table import Table, Column, MaskedColumn
        from astropy.coordinates import SkyCoord

        if '.fits' in fn:
            tab = Table(fits.getdata(fn))
        else:
            tab = ascii.read(fn, delimiter=',')

        # if any of the upper-case versions are present, create mixed-case aliases
        for real, alias in self.catalog_aliases.items():
            if real in tab.colnames and alias not in tab.colnames:
                tab[alias] = tab[real]

        # add UBVRI converted from SDSS mags
        U, B, V, R, I = sdss_to_UBVRI(*[tab[b] for b in 'ugriz'])
        for b in 'UBVRI':
            dat = locals()[b]
            colcls = MaskedColumn if hasattr(dat, 'mask') else Column
            tab.add_column(colcls(name=b, data=dat))

        if 'psf_u' in tab.colnames:
            pU, pB, pV, pR, pI = sdss_to_UBVRI(*[tab['psf_' + b] for b in 'ugriz'])
            for b in 'UBVRI':
                dat = locals()['p' + b]
                colcls = MaskedColumn if hasattr(dat, 'mask') else Column
                tab.add_column(colcls(name='psf_' + b, data=dat))

        if 'rhost' not in tab.colnames:
            tabsc = SkyCoord(u.Quantity(tab['ra'], u.deg), u.Quantity(tab['dec'], u.deg))
            rhost = self.coords.separation(tabsc).degree
            colcls = MaskedColumn if hasattr(rhost, 'mask') else Column
            tab.add_column(colcls(name='rhost', data=rhost))
            tab.add_column(colcls(name='rhost_kpc', data=np.radians(rhost)*self.distmpc*1000))

        if 'type' not in tab.colnames or 'phot_sg' not in tab.colnames:
            # some catalogs had 'phot_sg' as an *integer* 3/6
            if 'phot_sg' in tab.colnames and tab['phot_sg'].dtype.kind == 'i':
                typeint = tab['phot_sg']
                tab.remove_column('phot_sg')
            elif 'type' in tab.colnames and tab['type'].dtype.kind == 'i':
                typeint = tab['type']
                tab.remove_column('type')
            else:
                raise ValueError("couldn't figure out what's up with type and phot_sg")

            sgstr = np.zeros(len(tab), dtype='S6')
            sgstr[typeint == 3] = 'GALAXY'
            sgstr[typeint == 6] = 'STAR'
            tab.add_column(colcls(name='type', data=typeint))
            tab.add_column(colcls(name='phot_sg', data=sgstr))

        tab['coord'] = SkyCoord(tab['ra']*u.deg, tab['dec']*u.deg)

        return tab


    def open_on_nsasite(self):
        """
        Uses `webbrowser` to open the page for this host on the NSA web site.
        """
        import webbrowser

        urltempl = 'http://www.nsatlas.org/getAtlas.html?search=nsaid&nsaID={0}&submit_form=Submit'

        webbrowser.open(urltempl.format(self.nsaid))

    def __repr__(self):
        clsname = super(NSAHost, self).__repr__().split()[0][1:]  # strips the "<" at the start
        altnamepart = '' if len(self.altnames) == 0 else (' AKA: ' + str(self.altnames))
        return "<{clsname} object w/ name '{name}'{altnamepart}>".format(clsname=clsname, name=self.name, altnamepart=altnamepart)

    def sdss_image_cutout(self, savefn=None, drorurl=10, imagesize=(512, 512),
                          opts='', scale = '0.396127', targets=None,
                          raoffset=0*u.deg, decoffset=0*u.deg):
        """
        Saves an SDSS image of the host, optionally with a list of targets
        marked.

        This requires the `requests` package (http://requests.readthedocs.org/).

        Parameters
        ----------
        savefn : None or str
            The file name to save to, or None to not save
        drorurl : str or int
            If an int, this is the data release to use, otherwise a string URL
            pointing to the image cutout URL.
        imagesize : 2-tuple of its
            Size of the image in (width, height)
        opts : str
            Options the SDSS server accepts like showing photometric objects or
            scale bars or the like.  Go to the finding chart tool to see
            exmaples.
        scale : number or angle Quantity
            If a raw number, it will be taken as the image scale in arcsec per
            pixel. If an angle Quantity, it's the size of the *whole* image.
        targets : None or table
            A table of targets to mark.  Should be the same sort of table output
            by `targeting.select_targets`.
        raoffset : angle Quantity
            How much to offset the image center in RA
        decoffset : angle Quantity
            How much to offset the image center in Dec

        Returns
        -------
        imagedata
            Either the image suitable for IPython if you're in IPython, or the
            raw image data if not.
        """
        import requests

        if isinstance(drorurl, six.string_types):
            url = drorurl
        else:
            url = 'http://skyservice.pha.jhu.edu/DR' + str(drorurl) + '/ImgCutout/getjpeg.aspx'

        if u.arcsec.is_equivalent(scale):
            avgimagesize = sum(imagesize) / len(imagesize)
            scale = scale.to(u.arcsec).value / avgimagesize

        if targets is None:
            qry = ''
        else:
            qry = ['id ra dec']
            for t in targets:
                qry.append('{0} {1} {2}'.format(t['objID'], t['ra'], t['dec']))
            qry = '\n'.join(qry)

        res = requests.post(url, data={'ra': (self.coords.ra + raoffset).to(u.deg).value,
                                       'dec': (self.coords.dec + decoffset).to(u.deg).value,
                                       'width': imagesize[0],
                                       'height': imagesize[1],
                                       'opt': opts,
                                       'scale': scale,
                                       'query': qry})

        if savefn:
            with open(savefn, 'w') as f:
                f.write(res.content)

        try:
            __IPYTHON__
            from IPython.display import Image
            return Image(data=res.content, format='jpeg')
        except NameError:
            return res.content

    def within_environs(self, scorcat):
        if isinstance(scorcat, SkyCoord):
            sc = scorcat
        else:
            try:
                sc = SkyCoord.guess_from_table(scorcat)
            except u.UnitsError:
                sc = SkyCoord.guess_from_table(scorcat, unit=u.deg)

        return sc.separation(self.coords) < self.environsarcmin * u.arcmin


def download_with_progress_updates(u, fw, nreports=100, msg=None, outstream=sys.stdout):
    """
    Download a file and give progress updates on the download.

    Parameters
    ----------
    u : result of `urlopen`
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

    from astropy.io import fits
    from hosts import download_with_progress_updates

    if fn is None:
        fn = NSAFILENAME

    if fn in _cachednsa:
        if TELL_IF_USING_CACHED:
            print('Using cached NSA for file', fn)
        return _cachednsa[fn]

    if os.path.exists(fn):
        if TELL_IF_USING_CACHED:
            print('Loading NSA from local file', fn)
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


def construct_sdss_query(ra, dec, radius=1*u.deg, into=None, magcut=None,
                         inclphotzs=True, applyphotflags=False, xmatchwise=False):
    """
    Generates the query to send to the SDSS to get the full SDSS catalog around
    a target.

    Parameters
    ----------
    ra : `Quantity` or float
        The center/host RA (in degrees if float)
    dec : `Quantity` or float
        The center/host Dec (in degrees if float)
    radius : `Quantity` or float
        The radius to search out to (in degrees if float)
    into : str or None
        The name of the table to construct in your `mydb` if you want to use
        this with CasJobs, or None to have no "into" in the SQL. This also
        adjust other parts of the query a little to work with CasJobs instead
        of the direct query.
    magcut : 2-tuple or None
        if not None, adds a magnitude cutoff.  Should be a 2-tuple
        ('magname', faintlimit). Ignored if None.
    inclphotzs : bool
        If True, includes a further join to add phot-zs where present
    applyphotflags: bool
        If True, the query will some basic photometric flags (see the code for
        the exact flags)
    xmatchwise : bool
        If True, the query will also cross-match against WISE and include w1
        mags.

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
    ISNULL(s.z, -1) as spec_z, ISNULL(s.zErr, -1) as spec_z_err, ISNULL(s.zWarning, -1) as spec_z_warn, s.class as spec_class,
    s.subclass as spec_subclass{photzdata}{wisedata}


    {intostr}
    FROM {funcprefix}fGetNearbyObjEq({ra}, {dec}, {radarcmin}) n, PhotoPrimary p
    {wisejoins}
    LEFT JOIN SpecObj s ON p.specObjID = s.specObjID
    {photzjoins}WHERE n.objID = p.objID{magcutwhere}
    {photflags}
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

    if inclphotzs:
        photzdata = ', ISNULL(pz.z,-1) as photz,ISNULL(pz.zerr,-1) as photz_err'
        photzjoins = 'LEFT JOIN PhotoZ pz ON p.ObjID = pz.ObjID\n'
    else:
        photzdata = photzjoins = ''

    if xmatchwise:
        #wisedata = ', ISNULL(W.w1mpro,-99) as w1, ISNULL(W.w1sigmpro, -99) as w1_err'
        wisedata = ', ' + dedent("""
        ISNULL(w.j_m_2mass,9999) as J, ISNULL(w.j_msig_2mass,9999) as Jerr,
        ISNULL(w.H_m_2mass,9999) as H, ISNULL(w.h_msig_2mass,9999) as Herr,
        ISNULL(w.k_m_2mass,9999) as K, ISNULL(w.k_msig_2mass,9999) as Kerr,
        ISNULL(w.w1mpro,9999) as w1, ISNULL(w.w1sigmpro,9999) as w1err,
        ISNULL(w.w2mpro,9999) as w2, ISNULL(w.w2sigmpro,9999) as w2err,
        ISNULL(w.w3mpro,9999) as w3, ISNULL(w.w3sigmpro,9999) as w3err,
        ISNULL(w.w4mpro,9999) as w4, ISNULL(w.w4sigmpro,9999) as w4err""")

        wisejoins = ('LEFT JOIN WISE_xmatch X ON p.objid = X.sdss_objid\n'
                     'LEFT JOIN wise_allsky W ON x.wise_cntr = W.cntr')
    else:
        wisejoins = wisedata = ''

    if applyphotflags:
        photflags = dedent("""
        AND (flags & dbo.fPhotoFlags('BINNED1')) != 0
        AND (flags & dbo.fPhotoFlags('SATURATED')) = 0
        AND (flags & dbo.fPhotoFlags('BAD_COUNTS_ERROR')) = 0 -- "interpolation affected many pixels; PSF flux error is inaccurate and likely underestimated."
        """[1:-1])
        # these were considered, but found to sometimes remove *maybe* galaxies
        #AND (flags & 0x8100000c00a0) == 0  -- not NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP, SATURATED, or BAD_COUNTS_ERROR
        #AND ((flags & 0x400000000000) == 0)  | (psfmagerr_g <= 0.2)  --  DEBLEND_NOPEAK on faint ones
    else:
        photflags = ''

    if isinstance(ra, float):
        ra = ra*u.deg
    ra = ra.to(u.deg).value

    if isinstance(dec, float):
        dec = dec*u.deg
    dec = dec.to(u.deg).value

    if isinstance(radius, float):
        radius = radius*u.deg
    radarcmin = radius.to(u.arcmin).value

    # ``**locals()`` means "use the local variable names to fill the template"
    return query_template.format(**locals())


def construct_usnob_query(ra, dec, radius=1*u.deg, verbosity=1, votable=False, baseurl=USNOB_URL):
    """
    Generate a USNO-B query for the area around a target.

    Parameters
    ----------
    ra : `Quantity` or float
        The center/host RA (in degrees if float)
    dec : `Quantity` or float
        The center/host Dec (in degrees if float)
    radius : `Quantity` or float
        The radius to search out to (in degrees if float)
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

    if isinstance(ra, float):
        ra = ra*u.deg
    if isinstance(dec, float):
        dec = dec*u.deg
    if isinstance(radius, float):
        radius = radius*u.deg

    parameters = urlencode([('CAT', 'USNO-B1'),
                            ('RA', ra.to(u.deg).value),
                            ('DEC', dec.to(u.deg).value),
                            ('SR', radius.to(u.deg).value),
                            ('VERB', verbosity),
                            ('cftype', 'XML/VO' if votable else 'ASCII'),
                            ('slf', 'ddd.ddd/dd.ddd'),
                            ('skey', 'Mag')])

    return baseurl + '?' + parameters


def download_sdss_query(query, fn=None, sdssurl=SDSS_SQL_URL, format='csv',
                   dlmsg='Downloading...', inclheader=True, usepost=False):
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
    usepost : bool
        If False, uses GET, otherwise POST

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
    import os
    import datetime
    from urllib import urlencode
    StringIO = six.moves.StringIO

    from hosts import download_with_progress_updates

    parameterstr = urlencode([('cmd', query.strip()), ('format', format)])

    #either open the requested file or a buffer to later return the values
    if fn is None:
        fw = StringIO()
    else:
        #make directories if needed
        fndir = os.path.split(fn)[0]
        if fndir and not os.path.isdir(fndir):
            os.makedirs(fndir)
        fw = open(fn, 'w')

    try:
        if usepost:
            import requests
            q = requests.post(sdssurl, data=parameterstr, stream=True).raw
        else:
            q = urlopen(sdssurl + '?' + parameterstr)
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
                if isinstance(inclheader, six.string_types):
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


_HIDE_GET_SAGA_HOSTS_WARN = False
def get_saga_hosts():
    """
    Returns a dictionary with the "official" SAGA names
    """
    from warnings import warn

    if not _HIDE_GET_SAGA_HOSTS_WARN:
        warn('This is out-of-date.  get_saga_hosts_from_google is a better choice')
    hostsd = {}


    #note that the first name here is for the *host*, not the system, system is the second
    hostsd['odyssey'] = NSAHost(147100, ['Odyssey', 'Odysseus', 'NGC6181']) #30.3 Mpc
    hostsd['iliad'] = NSAHost(150238, ['Iliad', 'Achilles', 'NGC7393'])
    hostsd['lotr'] = NSAHost(155005, ['LordoftheRings', 'FrodoBaggins',  'NGC895'])
    hostsd['starwars'] = NSAHost(53145, ['StarWars', 'LukeSkywalker', 'NGC5485'], shortname='SW')
    hostsd['aiw'] = NSAHost(140594, ['AliceInWonderland', 'Alice', 'NGC4030'], shortname='AIW')
    hostsd['beowulf'] = NSAHost(135667, ['Beowulf', 'Beowulf', 'NGC2750'], shortname='beo')
    hostsd['gilgamesh'] = NSAHost(166313, ['Gilgamesh', 'Gilgamesh', 'NGC5962'], shortname='gil')
    hostsd['hamlet'] = NSAHost(166035, ['Hamlet', 'Hamlet', 'NGC5899'], shortname='ham')

    return hostsd


def get_saga_hosts_from_google(googleusername=None, googlepasswd=None,
                               useobservingsummary=False, cachefile='hosts_dl.pkl',
                               clientsecretjsonorfn='client_secrets.json'):
    """
    Returns a lost of hosts obtained from querying the google spreadsheet
    "SAGA_Hosts+Satellites".  Currently only works for NSA hosts

    To get this to access google after they sunset the name/password login,
    you'll need to go to https://console.developers.google.com and create a
    project.   Click "credentials" under "APIs&Auth", click "Create new Client
    ID", and choose a "native" application.  After that's created, hit the
    "download JSON" button, and then point the `clientsecretjsonorfn` argument
    to the location of the downloaded file.

    Parameters
    ----------
    googleusername : str
        Your google login name, typically your gmail address.
    passwd : str or None, or None to use python's getpass to input it.
        Your google passwd, or None to use python's getpass to input it.
    useobservingsummary : bool
        If True, use the "SAGA Observing Summary" spreadsheet, instead.
    cachefile : str or bool
        If a string, specifies a filename to load the host list from, or if it
        is absent, the file the loaded data will get saved to.  To clear the
        cache, just delete the relevant file.  Note that if it ends in "pkl", a
        2 or 3 will be added for the python version.
    """
    import sys
    import getpass
    import gspread
    from utils import get_google_oauth2_credentials
    pickle = six.moves.cPickle

    if cachefile.endswith('.pkl'):
        cachefile = cachefile + str(sys.version_info[0])

    if os.path.isfile(cachefile):
        print('Using cached version of google hosts list from file "{0}"'.format(cachefile))
        with open(cachefile, 'rb') as f:
            return pickle.load(f)

    if useobservingsummary:
        ssname = 'SAGA Observing Summary'
    else:
        ssname = 'SAGA_Hosts+Satellites'

    if clientsecretjsonorfn:
        credentials = get_google_oauth2_credentials(clientsecretjsonorfn)
        try:
            c = gspread.authorize(credentials)
        except:
            print("It seems you had a problem authorizing.  You might want to "
                  "clearing the OAuthcredentials.txt file...")
            raise
    else:
        if googleusername is None:
            googleusername = getpass.getpass('Google username:')
        if googlepasswd is None:
            googlepasswd = getpass.getpass('Password for "{0}":'.format(googleusername))
        c = gspread.Client(auth=(googleusername, googlepasswd))
        c.login()

    ss = c.open(ssname)
    s = ss.get_worksheet(0)  # first worksheet

    hosts = []
    if useobservingsummary:
        col1 = s.col_values(1)
        startrow = col1.index('Summary of Observed Systems') + 3
        endrow = [i for i, v in enumerate(col1) if v is not None and v.startswith('these systems are currently deprecated')][0] + 1

        rowvals = [s.row_values(row) for row in range(startrow, endrow)]
        for r in rowvals:
            sysname = r[0]
            othernames = [nm.replace(' ', '') for nm in r[1].split(',')]

            nsanum = nsaidx = None
            for i, nm in enumerate(othernames):
                if nm.startswith('NSA'):
                    nsanum = int(nm[3:])
                    nsaidx = i
                    break
            if nsanum is None:
                continue  # skip this one, it's not an NSA object

            del othernames[i]
            if sysname is not None:
                othernames.insert(0, sysname)
            if not othernames:
                othernames = None

            hosts.append(NSAHost(nsanum, othernames))
    else:
        col1 = s.col_values(1)
        startrow = col1.index('SAGA') + 2
        endrow = startrow
        for v in col1[(startrow-1):]:
            endrow += 1
            if v.strip() == '':
                #we've hit the end
                break

        rowvals = [s.row_values(row) for row in range(startrow, endrow)]

        for r in rowvals:
            if ''.join(r).strip() == '':
                #blank line
                continue
            if not r[2].strip():
                print('Entry', r, 'does not have an NSA number, cannot use')
            nsanum = int(r[2])

            names = [r[0].strip()]
            if r[1] is not None:
                names.append('NGC' + r[1])
            names.append('NSA' + str(nsanum))

            hosts.append(NSAHost(nsanum, names))

    if cachefile:
        with open(cachefile, 'wb') as f:
            pickle.dump(hosts, f)

    return hosts

def use_base_catalogs(hostlst):
    """
    Update this list of hosts to use the "base catalog" version of the SDSS
    catalogs
    """
    for h in hostlst:
        h.fnsdss = 'base_sql_nsa{}.fits.gz'.format(h.nsaid)
        h._cached_sdss = None


#this adds the hosts from the get_saga_hosts function to the module's namespace
_HIDE_GET_SAGA_HOSTS_WARN = True
locals().update(get_saga_hosts())
_HIDE_GET_SAGA_HOSTS_WARN = False
