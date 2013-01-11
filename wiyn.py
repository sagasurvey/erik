from __future__ import division

"""
These functions are for the "Distant local groups" project WIYN-related work.

How to design WIYN configurations:

1. Choose your hosts, and make sure you have the SDSS and USNO catalogs
downloaded for those hosts. (See NSAHost.sdss_environs_query and
NSAHost.usnob_environs_query.)

For each host:
2. Generate the master catalog with `construct_master_catalog`. If you
   want to adjust any selection parameters, do so in `select_targets` and
   pass the result into `construct_master_catalog`.
3. Create the first whydra input (".ast") file by calling
   `generate_ast_file` with the master catalog from #2.
4. Copy the .ast file from `wiyn_targets` to the `whydra` directory (on dept Linux machines)
5. Run whydra on this ast file (on dept Linux machines), producing a .hydra
   file with the same name.  Not that you may want to use the `do_whydra` script.
6. copy the .hydra file (on dept Linux machines) to the `wiyn_targets` directory
7. (optional) Use `imagelist_selected_fops` on the .hydra file to check if some
   of the FOP star are galaxies or double stars or something. if so, comment them
   out in the master list and they won't appear in susequent AST files. You might
   also want to re-generate the field you're working on if you're down to 2 or 3
   FOPs.
8. Repeat 3-7 until you have all the configurations you need
9. Copy them to the WIYN observing computer (oatmeal), and observe! (Note that
   you may need to go by way of cork to get to oatmeal)
"""

import numpy as np


def select_fops(host, faintlimit=13.5, brightlimit=12., randomize=True):
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


def select_sky_positions(host, nsky=250, sdsscat=None, usnocat=None, nearnesslimitarcsec=15):
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


def construct_master_catalog(host, fnout=None, targetcat=None, fopcat=None,
    skyradec=None, faintlimit=None):
    """
    This function produces the "master" catalog for each host for WIYN/hydra
    observations. The master catalog contains lines for all the objects/sky/fops
    and potential targets.  It is laid out such that the rows can simply be
    extracted as-is for .ast files for input to WIYN's `whydra` program.

    Parameters
    ----------
    host : NSAHost
        the host object to use to construct this catalog.  Note that SDSS and
        USNO data must be already downloaded for this host.
    fnout : str or None
        The filename to save as the master catalog or None to use a
        default based on the host's name
    targetcat : str or None
        The catalog of targets or None to use the default options to
        generate one.
    fopcat : str or None
        The catalog of FOP stars or None to use the default options to
        generate one.
    """
    import os
    from targeting import usno_vs_sdss_offset, select_targets

    if fnout is None:
        fnout = os.path.join('wiyn_targets', host.name + '.cat')

    if targetcat is None:
        if faintlimit is None:
            targetcat = select_targets(host)
        else:
            targetcat = select_targets(host, faintlimit=faintlimit)
        print len(targetcat), 'objects'
    if fopcat is None:
        fopcat = select_fops(host)
        print len(fopcat), 'FOPS'
    if skyradec is None:
        skyradec = select_sky_positions(host)

    if len(targetcat) > 1999:  # 1999 because the central object also gets one
        print('whydra cannot handle > 2000 objects, truncating')
        targetcat = targetcat[:1999]
    if len(fopcat) > 2000:
        print('whydra cannot handle > 2000 FOPS, truncating')
        fopcat = fopcat[:2000]
    if len(skyradec[0]) > 2000:
        print('whydra cannot handle > 2000 sky points, truncating')
        skyradec = skyradec[0][:2000], skyradec[1][:2000]

    #determine the SDSS vs. USNO offset
    dra, ddec = usno_vs_sdss_offset(host.get_sdss_catalog(), host.get_usnob_catalog())
    print 'USNO/SDSS offsets:', dra * 3600, ddec * 3600

    if os.path.exists(fnout):
        raise ValueError('File for generating master catalog {0} already exists!'.format(fnout))

    print 'Constucting catalog in', fnout
    with open(fnout, 'w') as fw:
        #first add host as center, and as object
        fw.write(_whydra_file_line(0001, 'Center'.format(host.nsaid), host.ra, host.dec, 'C'))
        fw.write('\n')

        fw.write(_whydra_file_line(1000, 'NSA{0}'.format(host.nsaid), host.ra, host.dec, 'O'))
        fw.write('\n')

        i = 2000
        for obj in targetcat:
            ln = _whydra_file_line(i, 'SDSS', obj['ra'], obj['dec'], 'O')
            i += 1
            fw.write(ln)
            fw.write(' # mag_r={0:.2f}'.format(obj['r']))
            fw.write('\n')

        i = 5000
        for fop in fopcat:
            ln = _whydra_file_line(i, 'USNO' + fop['id'], fop['RA'] - dra, fop['DEC'] - ddec, 'F')
            fw.write(ln)
            fw.write(' # mag_R2={0:.2f}'.format(fop['R2']))
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

    if decstr[0] == '-' and decstr[2] == ':':  # missing 0 pad
        decstr = '-0' + decstr[1:]

    if name == 'SDSS':
        name = 'J' + rastr[:-1].replace(':', '') + decstr[:-1].replace(':', '')

    return '{i:04} {name:20} {ra} {dec} {code}'.format(i=i, name=name, ra=rastr,
        dec=decstr, code=code)


def reprocess_master_catalog(mastercatfn, whydraoutputs=None):
    """
    Takes the requested master catalog, removes all *target* objects that have been
    assigned a fiber previously, and returns the resulting catalog.

    Parameters
    ----------
    mastercatfn : str
        Location of the master catalog file.
    whydraoutputs : list of str or None
        The names of ".hydra" files output by `whydra` or None to automatically
        locate them based on the name of `mastercatfn`

    Returns
    -------
    catlines : list of str
        The lines from the master catalog with already selected targets removed
    usedhydrafns : list of str
        The ".hydra" files actually used (useful when `whydraoutputs` is None).
    """
    from os import path
    from warnings import warn
    from glob import glob

    if whydraoutputs is None:
        basename = path.split(mastercatfn)[-1].split('.')[0]
        whydraoutputs = glob('wiyn_targets/{basenm}*.hydra'.format(basenm=basename))
    print 'Using existing catalogs', whydraoutputs, 'for removing from master'

    ids = []
    names = []
    catlines = []
    objcode = []
    with open(mastercatfn) as f:
        for l in f:
            lprecomment = l.split('#')[0]
            if lprecomment.strip() != '':
                ls = lprecomment.split()
                ids.append(int(ls[0]))
                names.append(ls[1])
                objcode.append(ls[4])
                catlines.append(l)

    #ids/names to *remove* because they are in a hydra output already
    ids2 = []
    names2 = []
    for fn in whydraoutputs:
        with open(fn) as f:
            for l in f:
                if 'STATUS=' in l:
                    status = l.split('STATUS=')[1].strip()
                    if status == 'OK' or status == 'EDGE':
                        continue  # skip

                    try:
                        int(status)  # check to make sure it's assigned a fiber
                        ls = l.split()
                        code = ls[8]
                        if code == 'O':  # only remove program science targets
                            ids2.append(int(ls[0]))
                            names2.append(ls[1])
                    except ValueError:
                        warn('Unrecognized status ' + status)

    # now use all lines expect those added above to be removed
    resultlines = []
    if len(set(names)) == len(catlines):
        #find which to remove based on *names*
        for l, nm in zip(catlines, names):
            if nm in names2:
                #this means *don't* put it in the output catalog
                names2.remove(nm)
            else:
                resultlines.append(l)
    else:
        warn('Names are not all unique in master catalog {0}!  Matching on IDs'.format(mastercatfn))

        #find which to remove based on hydra i-numbers
        for l, hid in zip(catlines, ids):
            if hid in ids2:
                #this means *don't* put it in the output catalog
                ids2.remove(hid)
            else:
                resultlines.append(l)

    return resultlines, whydraoutputs


def generate_ast_file(mastercatfn, lst, obsepoch=None, whydrafiles=None,
    texp=1.5, wl=7000, finame=None, outdir='wiyn_targets', faintmagcut=None):
    """
    Create the `.ast` file for input into the `whydra` program.

    This will normally automatically figure out what you want aside from the
    local sidereal time at observation, so usually you just need to give it
    that.

    Parameters
    ----------
    mastercatfn : str
        The path to the master catalog for this .ast to be generated from.
    lst : float
        Local sidereal time that at which this should be observed.  This
        is used during the observations to give a default for the LST at
        observation and can be updated, so just give a decent guess here.
    obsepoch : float (epoch) or `astropy.Time` object or None
        The epoch of observation (for astrometry purposes) - getting
        within a month or so should be more than good enough. If None,
        the `current` epoch when this function is called will be used.
    whydrafiles : list of str or None
        The list of previously-observed files. If None, this will be
        automatically guessed at based on files with the same base name
        as the master catalog.
    texp : float
        Estimated exposure time
    wl : float
        Center wavelength
    finame : str or None
        Name of the field. If None, the name is based on the master
        catalog and those with similar names in `whydrafiles` (e.g.
        if `DLG1_1.ast` exists, this will be `DLG1_2.ast`)
    outdir : str
        The directory to save `finame` in (and look for default file names in)
    faintmagcut : float or None
        A lower magnitude to cut off at if only bright objects are desired in
        this ast file.
    """

    import time
    import os
    import re

    from astropy.time import Time

    if not os.path.exists(mastercatfn):
        newpath = os.path.join(outdir, mastercatfn)
        if os.path.exists(newpath):
            mastercatfn = newpath
        else:
            raise ValueError('Master catalog ' + mastercatfn + ' does not exist')

    cataloglines, whydrafiles = reprocess_master_catalog(mastercatfn, whydrafiles)

    if finame is None:
        finame = os.path.split(mastercatfn)[-1].replace('.cat', '') + '_' + str(len(whydrafiles) + 1)

    fnout = os.path.join(outdir, finame + '.ast')

    if obsepoch is None:
        obsjepoch = Time(time.time(), format='unix', scale='utc').jyear
    elif isinstance(obsepoch, Time):
        obsjepoch = obsepoch.jyear
    else:
        obsjepoch = obsepoch

    if os.path.exists(fnout):
        raise ValueError(fnout + ' exists!')

    #used for faint mag cutting
    magre = re.compile(r'.*mag_(.+?)=([0-9]*?\.[0-9]*).*')

    print 'Writing to', fnout
    with open(fnout, 'w') as fw:
        #first do the header
        fw.write('FIELD NAME: {0}\n'.format(finame))
        fw.write('INPUT EPOCH: 2000.00\n')
        fw.write('CURRENT EPOCH: {0:.2f}\n'.format(obsjepoch))
        fw.write('SIDEREAL TIME: {0:.2f}\n'.format(lst))
        fw.write('EXPOSURE LENGTH: {0:.2f}\n'.format(texp))
        fw.write('WAVELENGTH: {0:f}\n'.format(int(wl)))
        fw.write('CABLE: RED\n')
        fw.write('WEIGHTING: WEAK\n')

        for l in cataloglines:
            if faintmagcut:
                m = magre.match(l)
                if m:
                    #magbang = m.group(1)
                    mag = float(m.group(2))
                    if mag > faintmagcut:
                        continue

            lprecomment = l.split('#')[0]
            if not lprecomment.strip() == '':
                fw.write(lprecomment)
                if not lprecomment.endswith('\n'):
                    fw.write('\n')


def imagelist_selected_fops(hydrafile, copytoclipboard=True):
    """
    Views the FOPs from `hydrafile` in the SDSS image list utility site.  See
    `targeting.sampled_imagelist` for more details.
    """
    from astropy.coordinates import Angle
    from targeting import sampled_imagelist

    ras = []
    decs = []
    names = []
    fibs = []
    with open(hydrafile) as f:
        for l in f:
            if 'STATUS=' in l and l[52] == 'F':
                statstr = l.split('STATUS=')[1].strip()
                if statstr == 'OK' or statstr == 'EDGE':
                    continue

                fibs.append(int(statstr))
                names.append(l[5:26])
                ras.append(Angle(l[26:38], unit='hour').degrees)
                decs.append(Angle(l[39:51], unit='deg').degrees)

    print 'FOP names', names
    print 'FOP fibers', fibs

    return sampled_imagelist(ras, decs, len(ras), names=names, copytoclipboard=copytoclipboard)
