from __future__ import division, print_function

from astropy import units as u

"""
These functions are for the SAGA project's WIYN/HYDRA-related work.

How to design HYDRA configurations:

1. Choose your hosts, and make sure you have the SDSS and USNO catalogs
downloaded for those hosts. (See NSAHost.sdss_environs_query and
NSAHost.usnob_environs_query.)

For each host:
2. Generate the master catalog with `construct_master_catalog`. If you
   want to adjust any selection parameters, do so in `select_targets` and
   pass the result into `construct_master_catalog`.
3. Create the first whydra input (".ast") file by calling
   `generate_ast_file` with the master catalog from #2.
4. Copy the .ast file from `hydra_targets` to the `hydra_simulator/whydra` directory
   (on Yale astro dept Linux machines)
5. Run whydra on this ast file (on dept Linux machines), producing a .hydra
   file with the same name.  Note that you may want to use the `do_whydra`
   script - copy it to the `whydra` directory and use it to automate running whydra.
6. copy the .hydra file (on dept Linux machines) to the `hydra_targets` directory
7. (optional) Use `imagelist_selected_fops` on the .hydra file to check if some
   of the FOP star are galaxies or double stars or something. if so, comment them
   out in the master list and they won't appear in susequent AST files. You might
   also want to re-generate the field you're working on if you're down to 2 or 3
   FOPs.
8. Repeat 3-7 until you have all the configurations you need
9. Run the `generate_wiyn_cache` function to generate the cache, and upload it via
   the WIYN web interface.
10. Copy .hydra files to the WIYN observing computer (oatmeal) (Note that
    you may need to go by way of cork to get to oatmeal)
11. Observe!
"""

import numpy as np

try:
    import six
except ImportErrir:
    from astropy.extern import six

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


def construct_master_catalog(host, fnout=None, targetcat={}, fopcat=None,
    skyradec=None, faintlimit=None, fibermaglimit=None, orderby=None,
    usnosdssoffsettol=0.5):
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
    targetcat : str or dict
        The catalog of targets or a dict to use `targeting.select_targets`
        to generate one (the dict is used as the kwargs).
    fopcat : str or None
        The catalog of FOP stars or None to use the default options to
        generate one.
    skyradec : 2xN array or None
        (ra, dec) for locations of sky fibers or if None, will use
        `select_sky_positions` to find good sky fiber locations.
    faintlimit : float or None
        The faint-end cutoff for r-band magnitudes in the target catalog,
        or None to have no cutoff
    fibermaglimit : float or None
        The faint-end cutoff for fiber2mag_r in the target catalog,
        or None to have no cutoff
    orderby: str or list of str
        Order the target list in order by field(2) in the target catalog.  If
        prefixed with '-', it's in *decreasing* order, otherwise increasing
        (i.e., without '-', smallest number first).  The special 'lowphotz'
        means to change order to blocks of [photz<.1, nophotz, photz>.1]
    usnosdssoffsettol : float
        The number of arcseconds median offset acceptable between the SDSS and
        USNO-B frames

    Returns
    -------
    targetcat
        The catalog of targets used to generate the object list
    """
    import os
    from collections import Mapping
    from targeting import usno_vs_sdss_offset, select_targets

    if fnout is None:
        fnout = os.path.join('hydra_targets', host.name + '.cat')

    if isinstance(targetcat, Mapping):
        if faintlimit is not None:
            if 'faintlimit' in targetcat:
                raise ValueError('cannot give both faintlimit kwarg and in targetcat dict')
            targetcat['faintlimit'] = faintlimit
        targetcat = select_targets(host, **targetcat)
    elif faintlimit is not None:
        #do the faintlimit cut manually if not using `select_targets`
        targetcat = targetcat[targetcat['r'] < faintlimit]

    if fibermaglimit is not None:
        fmsk = targetcat['fiber2mag_r'] < fibermaglimit
        print('Fiber mag cut removes', np.sum(~fmsk), 'of', len(fmsk), 'objects')
        targetcat = targetcat[fmsk]

    print(len(targetcat), 'objects')
    if fopcat is None:
        fopcat = select_fops(host)
    print(len(fopcat), 'FOPS')
    if skyradec is None:
        skyradec = select_sky_positions(host)

    if orderby:
        if isinstance(orderby, six.string_types):
            orderby = [orderby]

        for field in orderby:
            if field=='lowphotz':
                pz = targetcat['photz']
                highz = np.where(pz > .1)[0]
                lowz = np.where((-1 < pz)&(pz < .1))[0]
                noz = np.where(pz==-1)[0]
                print('Found {0} objects at low photz, {1} at high, and {2} '
                      'without photz'.format(len(lowz), len(highz), len(noz)))
                sorti = np.concatenate([lowz, noz, highz])
            elif field.startswith('-'):
                sorti = np.argsort(targetcat[field[1:]])[::-1]
            else:
                sorti = np.argsort(targetcat[field])
            targetcat = targetcat[sorti]

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
    dra, ddec = usno_vs_sdss_offset(host.get_sdss_catalog(),
                                    host.get_usnob_catalog(),
                                    raiseerror=usnosdssoffsettol)
    print('USNO/SDSS offsets:', dra * 3600, ddec * 3600)

    if os.path.exists(fnout):
        raise ValueError('File for generating master catalog {0} already exists!'.format(fnout))

    print('Constucting catalog in', fnout)
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

    return targetcat


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


def update_master_catalog(oldmastercatfn, newmastercatfn,
                          consolidategrps=1*u.arcsec,
                          removecoords=None,
                          removecoordtol=1*u.arcsec):
    """
    This takes a master catalog `oldmastercatfn` and writes it out to
    `newmastercatfn` with modifications:

    * if `consolidategrps` is not None, it will take clusterings of objects
      within `consolidategrps` distance and remove all except the brightest
    * if `removecoords` is not None, it should be a coordinate object that will
      be matched to the input catalog as objects to *remove* from the output
      catalog.  Typically for targets already successfully observed.

    Returns the catalog rows that were *removed*
    """
    from astropy.coordinates import ICRS

    with open(oldmastercatfn) as f:
        lines = [l for l in f]

    toremove = np.zeros(len(lines), dtype='S10')  # string gives the removal reason

    ras = []
    decs = []
    isobj = []
    mags = []
    for l in lines:
        ls = l.split()
        ras.append(ls[2])
        decs.append(ls[3])
        isobj.append(ls[4]=='O')
        if 'mag_' in ls[-1]:
            mags.append(float(ls[-1].split('=')[-1]))
        else:
            mags.append(999.)
    isobj = np.array(isobj, dtype=bool)
    mags = np.array(mags)

    coos = ICRS(ras, decs, unit=(u.hourangle, u.degree))
    objcoos = coos[isobj]

    if removecoords is not None:
        idx, dd, d3d = removecoords.match_to_catalog_sky(objcoos)

        toremove[np.arange(len(toremove))[isobj][idx[dd < removecoordtol]]] = 'matched'

        print(sum(dd < removecoordtol), 'objects in master catalog matched to be removed of', len(removecoords.ra), 'possible')
    if consolidategrps is not None:
        pairs = []
        n = 2
        while n > 0:
            idx, dd, d3d = objcoos.match_to_catalog_sky(objcoos, n)

            for i1, i2 in zip(np.where(dd <= consolidategrps)[0],
                              idx[dd <= consolidategrps]):
                pairs.append((i1, i2))

            if dd.arcsec.min() > consolidategrps.to(u.arcsec).value:
                #none are smaller than the minimum separation
                n = -1
            else:
                n += 1

        pairs = np.array(pairs)
        #make first in pair minimum
        toswap = np.argmin(pairs, 1).astype(bool)
        pairs[toswap] = pairs[toswap, ::-1]
        #sort in ascending order on min/0th
        pairs = pairs[np.argsort(pairs[:, 0])]

        grpnums = np.arange(len(objcoos.ra))

        i = 0
        while np.any(grpnums[pairs[:, 1]] != grpnums[pairs[:, 0]]):
            for i1, i2 in pairs:
                #go *in order* (which is why we can't do this with numpy directly)
                grpnums[i1] = grpnums[i2]
            i += 1
            if i > 100:
                raise ValueError('did>100 iterations, probably stuck?')

        #now gather the groups, and set all those *except* the brightest one to
        #False
        nremoved = ngrps = 0
        #note that the indexes in the groups are *object-only* indecies, not
        #everything including skys and FOPS
        objidxs = np.arange(len(toremove))[isobj]
        objmags = mags[isobj]
        for grp in np.unique(grpnums):
            msk = grpnums == grp
            if np.sum(msk) > 1:
                ngrps += 1
                magsort = np.argsort(objmags[msk])
                toremove[objidxs[msk][magsort[1:]]] = 'grpw' + str(objidxs[msk][magsort[0]])
                nremoved += len(magsort) - 1

        print('Removed', nremoved, 'from', ngrps, 'clustered groups')

    removedlines = []
    with open(newmastercatfn, 'w') as f:
        for i, l in enumerate(lines):
            if toremove[i]:
                removedlines.append(l)
            else:
                f.write(l)

    return removedlines, toremove[toremove != '']


def reprocess_master_catalog(mastercatfn, whydraoutputs=None):
    """
    Takes the requested master catalog, removes all *object* targets that have
    been assigned a fiber previously, and returns the resulting catalog.

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
        whydraoutputs = glob('hydra_targets/{basenm}*.hydra'.format(basenm=basename))
    print('Using existing catalogs', whydraoutputs, 'for removing from master')

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
    texp=1.5, wl=7000, finame=None, outdir='hydra_targets', faintmagcut=None,
    scpname='turtle', scpusername='ejt26'):
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
        if `DLG1_1.ast` exists, this will be `DLG1_2.ast`).  If there is
        an '%i' in this string, it will be replaced with the next in the
        sequence based on `whydrafiles`.
    outdir : str
        The directory to save `finame` in (and look for default file names in)
    faintmagcut : float or None
        A lower magnitude to cut off at if only bright objects are desired in
        this ast file.
    scpname : str
        The name of the computer to use in the scp commands displayed at the
        end.
    scpusername : str
        The username for the directories in the scp commands displayed at the end.
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
    elif '%i' in finame:
        finame = finame % (len(whydrafiles) + 1)

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

    print('Writing to', fnout)
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

    if scpname:
        print('SCP commands:')
        print('scp {0} {scpname}:/home/{scpusername}/hydra_simulator/whydra'.format(fnout, scpname=scpname, scpusername=scpusername))
        print('scp "{scpname}:/home/{scpusername}/hydra_simulator/whydra/{0}.hydra" hydra_targets'.format(finame, scpname=scpname, scpusername=scpusername))

def imagelist_selected_fops(hydrafile, copytoclipboard=True, openurl=True):
    """
    Short for `imagelist_fibers` with `objtype` set to 'fop' - this function
    mainly exists for backwards compatibility.
    """
    return imagelist_fibers(hydrafile, 'fop', copytoclipboard=copytoclipboard,
                            openurl=openurl)

def imagelist_fibers(hydrafile, objtype, copytoclipboard=True, openurl=True):
    """
    Views objects from `hydrafile` in the SDSS image list utility site.
    `objtype` can be 'target', 'fop', or 'sky'.

    See `targeting.sampled_imagelist` for more details.
    """
    from astropy.coordinates import Angle
    from targeting import sampled_imagelist


    if objtype == 'target':
        objcode = 'O'
    elif objtype == 'fop':
        objcode = 'F'
    elif objtype == 'sky':
        objcode = 'S'
    else:
        raise ValueError('Invalid objtype ' + str(objtype))

    ras = []
    decs = []
    names = []
    fibs = []
    with open(hydrafile) as f:
        for l in f:
            if 'STATUS=' in l and l[52] == objcode:
                statstr = l.split('STATUS=')[1].strip()
                if statstr == 'OK' or statstr == 'EDGE':
                    continue

                fibs.append(int(statstr))
                names.append(l[5:26])
                ras.append(Angle(l[26:38], unit='hour').degree)
                decs.append(Angle(l[39:51], unit='deg').degree)
    if openurl:
        return sampled_imagelist(ras, decs, len(ras), names=names, copytoclipboard=copytoclipboard)
    else:
        return sampled_imagelist(ras, decs, len(ras), names=names, copytoclipboard=copytoclipboard, url=None)


def parse_master(mastercatfn, objtype):
    """
    Returns a table of entries from a WIYN master catalog.

    Parameters
    ----------
    mastercatfn : str
        The file name of the master catalog
    objtype : str
        The type of object to show.  Valid options are 'target', 'fop', 'sky', or 'all'.

    Returns
    -------
    table : astropy.table.Table
        A table with the 'objID', 'ra', and 'dec' columns
    """
    from astropy.coordinates import Angle
    from astropy.table import Table, Column

    ras = []
    decs = []
    names = []

    if objtype == 'target':
        objcode = 'O'
    elif objtype == 'fop':
        objcode = 'F'
    elif objtype == 'sky':
        objcode = 'S'
    elif objtype == 'all':
        objcode = ''
    else:
        raise ValueError('Invalid objtype ' + str(objtype))

    with open(mastercatfn) as f:
        for l in f:
            ls = l.split()

            if l[0].startswith('#') or len(ls) < 5:
                continue

            if (not objcode) or ls[4] == objcode:
                names.append(l[5:26])
                ras.append(Angle(l[26:38], unit='hour').degree)
                decs.append(Angle(l[39:51], unit='deg').degree)

    tab = Table()
    tab.add_column(Column(data=names, name='objID'))
    tab.add_column(Column(data=ras, name='ra'))
    tab.add_column(Column(data=decs, name='dec'))
    return tab


def imagelist_from_master(mastercatfn, objtype, nobjs=None,
                          copytoclipboard=True, openurl=True,
                          nearloc=None):
    """
    Shows objects from a WIYN master catalog in the imagelist tool.

    Parameters
    ----------
    mastercatfn : str
        The file name of the master catalog
    objtype : str
        The type of object to show.  Valid options are 'target', 'fop', or 'sky'.
    nobjs : int or None
        The number of objects to show (randomly sampled) or None to show all
    copytoclipboard : bool
        If True, the imagelist table is copied to the clipboard
    openurl : bool
        If True, the imagelist url is opened in a web browser window
    nearloc : 2-tuple (coord, angle) or None
        Only show objects within `angle` of `coord` if not None
    """
    from astropy.coordinates import ICRS
    from targeting import sampled_imagelist

    tab = parse_master(mastercatfn, objtype)

    ras = tab['ra']
    decs = tab['dec']
    names = tab['objID']

    if nearloc:
        coord, angle = nearloc
        dsep = ICRS(ras*u.deg, decs*u.deg).separation(coord)
        msk = dsep < angle
        ras = ras[msk]
        decs = decs[msk]
        names = names[msk]

    if nobjs is None:
        nobjs = len(ras)

    if openurl:
        return sampled_imagelist(ras, decs, nobjs, names=names, copytoclipboard=copytoclipboard)
    else:
        return sampled_imagelist(ras, decs, nobjs, names=names, copytoclipboard=copytoclipboard, url=None)


def generate_wiyn_cache(outfn, infns='hydra_targets/*.hydra'):
    """
    Generate the file of target positions in the format WIYN wants.

    Submit it at the WIYN web page at http://www-kpno.kpno.noao.edu/Info/submitcache.html

    Parameters
    ----------
    outfn : str or None
        The name of a file to save the cache into, or None to return the file
        contents.
    infns : str or list of strings
        The files or file patterns of the hydra files to generate from

    Returns
    -------
    The file contents if `outfn` is None.
    """
    from glob import glob

    if isinstance(infns, six.string_types):
        infns = glob(infns)
    else:
        newinfns = []
        for infn in infns:
            newinfns.extend(glob(infn))
        infns = newinfns

    cenlines = {}
    for infn in infns:
        currnm = None
        with open(infn) as fr:
            for l in fr:
                if l.startswith('FIELD NAME'):
                    currnm = l[11:].strip()
                ls = l.split()
                if len(ls) > 8 and ls[8] == 'C':
                    if currnm is None:
                        raise ValueError('Field {0} does not have a name!'.format(infn))
                    cenlines[currnm] = l
                    break
            else:
                raise ValueError('File {0} did not have a center!'.format(infn))

    maxchars = max([len(nm) for nm in cenlines])
    paddednms = []
    for nm in sorted(cenlines):
        blankchars = maxchars + 1 - len(nm)
        paddednms.append(nm + ' ' * blankchars)

    outlines = []
    for nm, k in zip(paddednms, sorted(cenlines)):
        rastr = cenlines[k][26:38]
        decstr = cenlines[k][39:52]

        outlines.append('{nm} {rastr} {decstr}\n'.format(nm=nm, rastr=rastr, decstr=decstr))

    if outfn is None:
        return ''.join(outlines)
    else:
        with open(outfn, 'w') as fw:
            for l in outlines:
                fw.write(l)


def clean_master_catalog(mastercatfn, invalidnms, outfn=None):
    """
    This takes a WIYN master catalog and names of a bunch of invalid objects
    (e.g. initially selected FOPS that are actually galaxies), and removes
    them, writing a new catalog

    Parameters
    ----------
    mastercatfn : str
        The input master catalog
    invalidnms : list of str
        The names of the invalid objects (must match the second column
        of the master catalog, case-sensitive)
    outfn : str or None
        The name of the output file or None to overwrite the original `mastercatfn`
    """
    invalidnms = tuple(invalidnms)

    newlns = []
    skipped = []
    with open(mastercatfn) as fr:
        for l in fr:
            ls = l.split()
            if ls[1] in invalidnms:
                skipped.append(ls[1])
            else:
                newlns.append(l)

    if outfn is None:
        outfn = mastercatfn

    with open(outfn, 'w') as fw:
        for l in newlns:
            fw.write(l)
