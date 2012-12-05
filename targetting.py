from __future__ import division

"""
These functions are for the "Distant local groups" project target selection.
"""

import sys

import numpy as np

NSA_VERSION = '0.1.2'  # used to find the download location/file name
NSAFILENAME = 'nsa_v{0}.fits'.format(NSA_VERSION.replace('.', '_'))

SDSS_SQL_URL = 'http://skyserver.sdss3.org/dr9/en/tools/search/x_sql.asp'


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


def construct_query(ra, dec, radius=1, into=None):
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
       p.ra, p.dec, p.type, p.flags,
       p.modelMag_u as u, p.modelMag_g as g, p.modelMag_r as r,p.modelMag_i as i,p.modelMag_z as z,
       p.modelMagErr_u as u_err, p.modelMagErr_g as g_err, p.modelMagErr_r as r_err,p.modelMagErr_i as i_err,p.modelMagErr_z as z_err,
       p.psfMag_u as psf_u, p.psfMag_g as psf_g, p.psfMag_r as psf_r,p.psfMag_i as psf_i,p.psfMag_z as psf_z,
       extinction_u as Au, extinction_g as Ag, extinction_r as Ar, extinction_i as Ai, extinction_z as Az

    {into}
    FROM {funcprefix}fGetNearbyObjEq({ra}, {dec}, {radarcmin}) n, PhotoPrimary p
    WHERE n.objID=p.objID
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


def download_query(query, fn=None, sdssurl=SDSS_SQL_URL, format='csv',
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


def



