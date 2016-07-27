"""
This file contains functions for targeting objects as part of the SAGA Survey's
AAT observations.

Note that priority levels 5 and 7 are only used if wise data is present
"""
from __future__ import print_function

import numpy as np

try:
    import six
except ImportErrir:
    from astropy.extern import six

from astropy import units as u
from astropy.coordinates import Angle


def prioritize_targets(targets, rvir=300*u.kpc, scheme='jun2014'):
    rvkpc = rvir.to(u.kpc).value

    pris = np.zeros(len(targets), dtype=int)
    if scheme == 'jul2014':
        g = targets['g']
        r = targets['r']
        i = targets['i']
        w1 = targets['w1'] if 'w1' in targets.colnames else None
        rkpc = targets['rhost_kpc']


        pris[(g-r < 1.2) & (r-i < 0.7) & (r<21) & (rkpc < rvkpc)] = 1
        pris[(g-r < 1.2) & (r-i < 0.7) & (r<20.5) & (rkpc < rvkpc)] = 2

        pris[(g-r < 1.0) & (r-i < 0.5) & (r<20.5) & (rkpc > rvkpc)] = 3

        #add priority level 5 and 7 *only* if WISE data are present
        pris[(g-r < 1.0) & (r-i < 0.5) & (r<21) & (rkpc < rvkpc)] = 4
        if w1 is not None:
            pris[(r-w1<2.5) & (g-r < 1.0) & (r-i < 0.5) & (r<21) & (rkpc < rvkpc)] = 5

        pris[(g-r < 1.0) & (r-i < 0.5) & (r<20.5) & (rkpc < rvkpc)] = 6
        if w1 is not None:
            pris[(r-w1<2.5) & (g-r < 1.0) & (r-i < 0.5) & (r<20.5) & (rkpc < rvkpc)] = 7
    elif scheme == 'jun2015baseline':
        sborder = np.argsort(targets['sb_petro_r'])
        # now split into 2 based on SB-ile
        pris[:] = 2
        pris[sborder[:len(pris)//2]] = 1

        #middle one for 3
        #pris[sborder[:2*len(pris)//3]] = 2

        pris[targets['rhost_kpc'] < rvkpc] += 2  # prefer close-in to outer
    else:
        raise ValueError('unrecognized scheme')

    return pris


def produce_master_fld(host, utcobsdate, catalog, pris, guidestars,
                       fluxstars, skyradec, outfn=None, randomizeorder=True,
                       fluxpri=8, inclhost=True, manualtargetlines=[]):
    """
    Priority of 1 to 9 (9 highest) means use, any other `pris` values skipped

    `skyradec` can either be (ra, dec) or a string to be a filename that will
    simply be pulled in wholesale.

    `inclhost` can give the priority, if desired - defaults to 9

    `manualtargetlines` is a list of lines (in string form) to add by hand
    """
    lines = []
    lines.append('LABEL ' + host.name + ' base catalog')
    lines.append('UTDATE  {yr} {mo:02} {day:02}'.format(yr=utcobsdate.year,
                                                        mo=utcobsdate.month,
                                                        day=utcobsdate.day))
    censtr = host.coords.to_string('hmsdms', sep=' ', precision=2, alwayssign=True)[1:]  # strips the leading '+'
    lines.append('CENTRE  ' + censtr)
    lines.append('EQUINOX J2000.0')
    #lines.append('WLEN1 3700')
    #lines.append('WLEN2 8800')
    lines.append('# End of Header')
    lines.append('')
    lines.append('# TargetName(unique for header) RA(h m s) Dec(d m s) TargetType(Program,Fiducial,Sky) Priority(9 is highest) Magnitude 0 TargetName')

    if randomizeorder:
        idxs = np.random.permutation(len(catalog))
    else:
        idxs = np.arange(len(catalog))

    if inclhost:
        if inclhost is True:
            inclhost = 9

        entry = []

        entry.append(host.name.replace(' ', ''))
        entry.append(host.coords.ra.to(u.hourangle).to_string(sep=' ', precision=2))
        entry.append(host.coords.dec.to_string(sep=' ', alwayssign=True, precision=2))
        entry.append('P')
        entry.append(str(int(inclhost)))
        entry.append('{0:0.2f}'.format(host.r + host.distmod))
        entry.append('0')
        entry.append('host galaxy')

        lines.append(' '.join(entry))

    if pris is None:
        pris = prioritize_targets(catalog)
    elif isinstance(pris, six.string_types):
        pris = prioritize_targets(catalog, scheme=pris)

    skippedbadpri = 0
    extranotes = 'extra_aat_notes' in catalog.colnames

    f2magnm = 'FIBER2MAG_R' if 'FIBER2MAG_R' in catalog.colnames else 'fiber2mag_r'

    for ci, pri in zip(catalog[idxs], pris[idxs]):
        if not 1 <= pri <= 9:
            skippedbadpri += 1
            continue

        entry = []

        entry.append(str(ci['objID']))
        entry.append(Angle(ci['ra'], u.deg).to(u.hourangle).to_string(sep=' ', precision=2))
        entry.append(Angle(ci['dec'], u.deg).to_string(sep=' ', alwayssign=True, precision=2))
        entry.append('P')
        entry.append(str(pri))
        entry.append('{0:0.2f}'.format(ci[f2magnm]))
        entry.append('0')
        entry.append('magcol=fiber2mag_r, model_r={0:.2f}'.format(ci['r']))
        if extranotes:
            entry[-1] = entry[-1] + ci['extra_aat_notes']

        lines.append(' '.join(entry))

    if skippedbadpri > 0:
        print('skipped', skippedbadpri, 'objects for priorities not in 1-9')

    #add any manually-added targets
    entry.extend(manualtargetlines)

    lines.append('\n#Flux stars')
    if randomizeorder:
        idxs = np.random.permutation(len(fluxstars))
    else:
        idxs = np.arange(len(fluxstars))
    for idx, fxs in enumerate(fluxstars[idxs]):
        entry = []

        entry.append('Flux' + str(idx))
        entry.append(Angle(fxs['ra'], u.deg).to(u.hourangle).to_string(sep=' ', precision=2))
        entry.append(Angle(fxs['dec'], u.deg).to_string(sep=' ', alwayssign=True, precision=2))
        entry.append('P')
        entry.append(str(int(fluxpri)))
        entry.append('{0:0.2f}'.format(fxs['r']))
        entry.append('0')
        entry.append('id=' + str(fxs['objID']))

        lines.append(' '.join(entry))

    lines.append('\n#Guide stars')
    if randomizeorder:
        idxs = np.random.permutation(len(guidestars))
    else:
        idxs = np.arange(len(guidestars))
    for idx, g in enumerate(guidestars[idxs]):
        entry = []

        entry.append('Guide' + str(idx))
        entry.append(Angle(g['ra'], u.deg).to(u.hourangle).to_string(sep=' ', precision=2))
        entry.append(Angle(g['dec'], u.deg).to_string(sep=' ', alwayssign=True, precision=2))
        entry.append('F')
        entry.append('9')
        entry.append('{0:0.2f}'.format(g['r']))
        entry.append('0')
        entry.append('id=' + str(g['objID']))

        lines.append(' '.join(entry))

    lines.append('\n#Sky positions')
    if isinstance(skyradec, six.string_types):
        with open(skyradec, 'r') as f:
            lines.extend(f.read().split('\n'))
    else:
        if randomizeorder:
            idxs = np.random.permutation(len(skyradec[0]))
        else:
            idxs = np.arange(len(skyradec[0]))
        skyradeczip = np.array(zip(*skyradec))
        for idx, (ra, dec) in enumerate(skyradeczip[idxs]):
            entry = []

            entry.append('Sky' + str(idx))
            entry.append(Angle(ra, u.deg).to(u.hourangle).to_string(sep=' ', precision=2))
            entry.append(Angle(dec, u.deg).to_string(sep=' ', alwayssign=True, precision=2))
            entry.append('S')
            entry.append('9')
            entry.append('20.00')
            entry.append('0')
            entry.append('sky')

            lines.append(' '.join(entry))

    if outfn is not None:
        with open(outfn, 'w') as f:
            for l in lines:
                f.write(l)
                f.write('\n')
    return lines


zlogcolnames = 'name,ra,dec,mag,z,sn,zqual,idx,sat?,star,unknown'.split(',')


def subsample_from_master_fld(masterfn, outfn, nperpri, nguides='all',
                              nflux='all', nsky='all', utcobsdate=None,
                              fieldname=None, listorem=None, dontrempri=None,
                              zlogfns=None, zltabkeepfunc=lambda entry: entry['zqual'] < 3,
                              guidemags='all'):
    """
    Selects from the master .fld and creates a smaller .fld file for consumption
    by configure.

    nperpri maps pri numbers to the number to do (any missing are
    treated as 0). 'all'

    `listorem` should be a list of ".lis" files of allocations as output by
    configure (or None)

    `dontrempri` is the priority to include even if it's in one of the remove
    lists.

    `zlogfns` should be a list of zlog files *which must match the corresponding
    .lis file* or None to not use the zlog file.

    `zltabkeepfunc` is a function that takes an entry in the zltab and returns
    True if the object should be *kept* even if its in listorem. (ignored if
    `zlogfns` is None).

    `guidemags` can be 'all' or a 2-tuple
    """
    from astropy.table import Table

    inhdr = True

    if fieldname is None:
        fieldname = 'subsampled'

    if nperpri == 'all':
        nperpri = {}
        for i in range(10):
            nperpri[i] = np.inf

    skydone = fluxdone = guidesdone = 0
    pridone = dict([(i, 0) for i in range(10) if i != 0])
    pritotal = dict([(i, 0) for i in range(10) if i != 0])

    namestoskip = []
    if listorem:
        if zlogfns and len(zlogfns) != len(listorem):
            raise ValueError('zlogfns and listorem do not match!')
        for i, lis in enumerate(listorem):
            prentslen = len(namestoskip)
            listab = load_lis_file(lis)[0]

            if zlogfns:
                zlogtab = Table.read(zlogfns[0], format='ascii', names=zlogcolnames)

                zlogfibnums = []
                for nm in zlogtab['name']:
                    if nm == 'noid':
                        zlogfibnums.append(-1)
                    else:
                        zlogfibnums.append(int(nm.split('_')[-1]))
            else:
                zlogtab = None

            keptfibers = []
            for entry in listab:
                nm = entry['ids']
                fibnum = entry['fibnums']

                if not (nm.startswith('Flux') or nm.startswith('Guide') or nm.startswith('Sky')):
                    if zlogtab is not None:
                        fibmatch = zlogfibnums == fibnum
                        if np.sum(fibmatch) == 0:
                            msg = 'Could not find a match in the zlog for fiber #{0}!'
                            print(msg.format(fibnum))
                            zlentry = None
                        elif np.sum(fibmatch) > 1:
                            msg = 'Found {0} zlog matches for fiber #{1}.  Using first one'
                            print(msg.format(np.sum(fibmatch), fibnum))
                            zlentry = zlogtab[fibmatch][0]
                        else:
                            zlentry = zlogtab[fibmatch][0]
                        if zlentry is not None and zltabkeepfunc(zlentry):
                            keptfibers.append(fibnum)
                            continue
                    namestoskip.append(nm)
            print('Kept the following fibers in due to zltabkeepfunc:', keptfibers)
            print('Found', len(namestoskip) - prentslen, 'objects to remove in', lis)

    with open(masterfn) as fr:
        with open(outfn, 'w') as fw:
            for l in fr:
                if inhdr:
                    if l.startswith('LABEL'):
                        fw.write(l.replace('base catalog', fieldname))
                    elif l.startswith('UTDATE'):
                        if utcobsdate is None:
                            fw.write(l)
                        else:
                            s = 'UTDATE  {yr} {mo:02} {day:02}\n'
                            fw.write(s.format(yr=utcobsdate.year,
                                              mo=utcobsdate.month,
                                              day=utcobsdate.day))
                    else:
                        fw.write(l)

                    if l.startswith('# End of Header'):
                        inhdr = False

                else:  # not in header
                    lst = l.strip()
                    if lst == '' or lst.startswith('#'):
                        fw.write(l)
                    elif l.startswith('Guide'):
                        if nguides == 'all' or (guidesdone < nguides):
                            if guidemags != 'all':
                                mag = float(l.split()[9])
                                if not (guidemags[0] < mag < guidemags[1]):
                                    continue  # outside of the valid mag range
                            fw.write(l)
                            guidesdone += 1
                    elif l.startswith('Flux'):
                        if nflux == 'all' or (fluxdone < nflux):
                            fw.write(l)
                            fluxdone += 1
                    elif l.startswith('Sky'):
                        if nsky == 'all' or (skydone < nsky):
                            fw.write(l)
                            skydone += 1
                    else:  # program target
                        ls = l.split()
                        pri = int(ls[8])

                        #skip *unless* dontrrempri == pri
                        if ls[0] in namestoskip:
                            del namestoskip[namestoskip.index(ls[0])]
                            if dontrempri != pri:
                                continue
                        ntodo = nperpri.get(pri, 0)
                        if pridone[pri] < nperpri.get(pri, 0):
                            fw.write(l)
                            pridone[pri] += 1
                        pritotal[pri] += 1
    if len(namestoskip) > 0:
        print('Had', len(namestoskip), 'unmatched list file objects:\n', namestoskip)

    msg = 'Total remaining in each priority ({0} fluxes, {1} guides, and {2} skies not included): {3}'
    print(msg.format(fluxdone, guidesdone, skydone, pritotal))


def imagelist_fld_targets(fldlinesorfn, ttype='all', **kwargs):
    from targeting import sampled_imagelist

    if isinstance(fldlinesorfn, six.string_types):
        with open(fldlinesorfn) as f:
            fldlines = f.read().split('\n')
    else:
        fldlines = fldlinesorfn

    prognum = -1
    if ttype.startswith('prog'):
        if ttype[4:] == '':
            prognum = int(ttype[4:])
        else:
            prognum = 0

    ras = []
    decs = []
    names = []
    for l in fldlines:
        if l.startswith('#') or l.startswith('*') or l.strip() == '':
            continue

        ls = l.split()
        if len(ls) > 9:
            if (ttype == 'all' or
               (ttype == 'sky' and ls[7] == 'S') or
               (ttype == 'guide' and ls[7] == 'F') or
               (ttype == 'flux' and ls[7] == 'P' and ls[0].startswith('Flux')) or
               (prognum>-1 and ls[7] == 'P' and not ls[0].startswith('Flux') and (prognum==0 or prognum==int(ls[8])))
               ):
                ras.append(Angle(ls[1]+'h'+ls[2]+'m'+ls[3]+'s').deg)
                decs.append(Angle(ls[4]+'d'+ls[5]+'m'+ls[6]+'s').deg)
                names.append(ls[0])

    kwargs['names'] = names
    return sampled_imagelist(ras, decs, **kwargs)


def select_guide_stars_usnob(host, faintlimit=13.5, brightlimit=12., randomize=True):
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

def select_guide_stars_sdss(cat, magrng=(12.5, 14)):
    msk = cat['type'] == 6 #type==6 means star

    magrng = min(*magrng), max(*magrng)

    msk = msk & (magrng[0] < cat['r']) & (cat['r'] < magrng[1])

    return cat[msk]

def select_sky_positions(host, nsky=250, sdsscat=None, usnocat=None,
                         nearnesslimitarcsec=15, outfn=None, rad=None):
    """
    Produces sky positions uniformly covering a circle centered at the host,
    with radius given by `environsarcmin`.

    Parameters
    ----------
    host : NSAHost
        The host to make sky
    nsky : int
        Number of sky positions to generate
    sdsscat : Table or None
        The SDSS catalog for this host or None to use `get_sdss_catalog`
    usnocat : Table or None
        The SDSS catalog for this host or None to use `get_usnob_catalog`
    nearnesslimitarcsec : float or Quantity
        How close a position has to be to a catalog entry to get eliminated
    outfn : None or str
        If given, a file to save the sky positions out suitable for use in an
        AAT fld file
    rad : angle quantity or None
        radius from host out to make sky positions.  If None, will use host
        `environsarcmin`.

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

    if not isinstance(nearnesslimitarcsec, u.Quantity):
        nearnesslimitarcsec = nearnesslimitarcsec * u.arcsec
    neardeg = nearnesslimitarcsec.to(u.deg).value

    skdt = cKDTree(np.array([sdsscat['ra'], sdsscat['dec']]).T)
    ukdt = cKDTree(np.array([usnocat['RA'], usnocat['DEC']]).T)

    if rad is None:
        raddeg = host.environsarcmin / 60.
    else:
        raddeg = rad.to(u.degree).value

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

    if outfn:
        with open(outfn, 'w') as f:
            for idx, (ra, dec) in enumerate(zip(ras[:nsky], decs[:nsky])):
                entry = []

                entry.append('Sky' + str(idx))
                entry.append(Angle(ra, u.deg).to(u.hourangle).to_string(sep=' ', precision=2))
                entry.append(Angle(dec, u.deg).to_string(sep=' ', alwayssign=True, precision=2))
                entry.append('S')
                entry.append('9')
                entry.append('20.00')
                entry.append('0')
                entry.append('sky')

                f.write(' '.join(entry))
                f.write('\n')

    return ras[:nsky], decs[:nsky]


def select_flux_stars(cat, magrng=(17, 17.7), extcorr=False, fluxfnout=None, onlyoutside=None):
    """
    Identifies flux calibration stars by choosing ~F stars based on color cuts.

    Uses the SDSS criteria for specphot standards:
    0.1 < (g-r) < 0.4

    their specphot stars are 16 < g < 17.1 and "reddening standards" are
    17.1 < g < 18.5

    `onlyoutside`, if not None, specifies that they should only be further from
    the given distance from the host. Should be an astropy quantity.
    """
    from astropy.units import degree, kpc
    from astropy.table import Table

    starmsk = cat['type'] == 6 #type==6 means star
    cat = cat[starmsk]


    u = cat['u']
    g = cat['g']
    r = cat['r']
    if extcorr:
        u = u - cat['Au']
        g = g - cat['Ag']
        r = r - cat['Ar']
    umg = u - g
    gmr = g - r

    #msk = (0.1 < gmr) & (gmr < 0.4) #old version

    std_color = (0.6 < umg) & (umg < 1.2)
    std_color = std_color & (0.0 < gmr) & (gmr < 0.6)
    std_color = std_color & (gmr > 0.75 * umg - 0.45)

    sp_std = std_color & (15.5 < g) & (g < 17.0)
    red_std = std_color & (17 < g) & (g < 18.5)

    minmag = min(*magrng)
    maxmag = max(*magrng)

    msk = (sp_std|red_std) & (minmag < r)& (r < maxmag)

    if onlyoutside:
        if onlyoutside.unit.is_equivalent(degree):
            msk = msk & ((degree * cat['rhost']) > onlyoutside)
        elif onlyoutside.unit.is_equivalent(kpc):
            msk = msk & ((kpc * cat['rhost_kpc']) > onlyoutside)
        else:
            raise ValueError('onlyoutside is not an angle or length')

    if fluxfnout:
        names = 'RA DEC u_psf g_psf r_psf i_psf z_psf extinction_r'.split()
        dat = [cat[msk]['ra'],
               cat[msk]['dec'],
               cat[msk]['psf_u'],
               cat[msk]['psf_g'],
               cat[msk]['psf_r'],
               cat[msk]['psf_i'],
               cat[msk]['psf_z'],
               cat[msk]['Ar']
              ]
        tab = Table(data=dat, names=names)
        with open(fluxfnout, 'w') as f:
            #f.write('#RA DEC u_psf g_psf r_psf i_psf z_psf extinction_r')
            tab.write(f, format='ascii.commented_header')

    return cat[msk]

def load_lis_file(fn):
    from astropy.coordinates import SkyCoord
    from astropy import table

    info = []

    comments = []

    fibnums = []
    ids = []
    ras = []
    decs = []
    codes = []
    pris = []
    mags = []

    with open(fn) as f:
        for l in f:
            if l.startswith('*'):
                l = l[1:]  # strip the *
                ls = l.split()
                if len(ls) < 5:
                    continue  # initial lines

                if ls[1] == 'Parked':
                    continue

                data = ls[:14]
                comments.append(' '.join(ls[14:]))

                fibnums.append(int(data[0]))
                ids.append(data[1])
                ras.append(':'.join(data[2:5]))
                decs.append(':'.join(data[5:8]))
                codes.append(data[8])
                pris.append(int(data[9]))
                mags.append(float(data[10]))
            else:
                info.append(l.strip())

    names = 'fibnums,ids,ras,decs,codes,pris,mags,comments'.split(',')
    tab = table.Table(names=names, data=[locals()[nm] for nm in names])

    sc = SkyCoord(ras, decs, unit=(u.hourangle, u.degree))

    return tab, sc, '\n'.join(info)


def load_fld(fn):
    from astropy.coordinates import SkyCoord
    from astropy import table

    header = []

    comment = []

    name = []
    ra = []
    dec = []
    code = []
    pri = []
    mag = []

    inheader = True

    with open(fn) as f:
        for l in f:
            if inheader:
                if l.startswith('# End of Header'):
                    inheader = False
                else:
                    header.append(l.strip())
            elif not l.startswith('#'):
                if l.strip() == '':
                    continue

                ls = l.split()
                if len(ls) < 5:
                    continue  # initial lines

                if ls[1] == 'Parked':
                    continue

                data = ls[:11]
                comment.append(' '.join(ls[11:]))

                name.append(data[0])
                ra.append(':'.join(data[1:4]))
                dec.append(':'.join(data[4:7]))
                code.append(data[7])
                pri.append(int(data[8]))
                mag.append(float(data[9]))

    names = 'name,ra,dec,code,pri,mag,comment'.split(',')
    tab = table.Table(names=names, data=[locals()[nm] for nm in names])
    sc = SkyCoord(ra, dec, unit=(u.hourangle, u.degree))

    return tab, sc, '\n'.join(header)
