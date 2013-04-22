"""
This file contains functions for targeting objects as part of the "distant
local group" project using MMT/Hectospec.
"""
import numpy as np

import targeting
from astropy import units as u


def make_catalog(host, fnout=None, targetfaintlim=(21.,22.), targetoutercutrad=300,
                 fluxrng=(17., 17.7), repeatflux=1, guidestarmagrng= (14, 15),
                 removegalsathighz=False, inclspecqsos=True, removespecstars=False,
                 fibermag_faintlimit=(22., 23.), fibermag_brightlimit=17,
                 fibermag_type='fibermag_r', useoutertargets=False):
    """
    Generates a catalog for the Hectospec "XFITFIBS" program.

    `targetfaintlim` can be just a mag to make everything "primary", or
    a 2-tuple that is (prisecboundary, faintestsec)

    `useoutertargets` means include targets beyond `targetoutercutrad`, but with
    a rank below where they would usually be

    Writes catalog to `fnout` if `fnout` is not None

    Returns the catalog as a `Table`

    `fibermag_faintlimit` is a lower limit *or* a two-tuple of limits for
    primary and secondary targets (secondary ignored if `targetfainlim` is just
    a single number) can also be None (`fibermag_brightlimit` can also be None)
    """
    from astropy.coordinates import Angle
    from astropy.table import Table

    fluxrank = 1
    hostrank = 2
    prisatrank = 3
    secsatrank = 5 if useoutertargets else 4 #the outer primaries get 4 if `useoutertargets` is True


    tabentries = []
    catlines = ['ra\tdec\tobject\trank\ttype\tmag',
                '--\t---\t------\t----\t----\t---']


    #add the host as the first entry
    rastr = Angle(host.ra, 'deg').format('hr', sep=':', precision=3)
    decstr = Angle(host.dec, 'deg').format('deg', sep=':', precision=3)
    mag = '{0:.2f}'.format(host.r + host.distmod)

    tabentries.append([rastr, decstr, host.name, str(hostrank), 'TARGET', mag])
    catlines.append('\t'.join(tabentries[-1]))

    try:
        prisecbounday, seclimit = targetfaintlim
    except TypeError:
        prisecbounday = seclimit = targetfaintlim


    if useoutertargets:
        targs = targeting.select_targets(host, faintlimit=seclimit, outercutrad=None, innercutrad=20, removegalsathighz=removegalsathighz, removespecstars=removespecstars, inclspecqsos=inclspecqsos)

        if targetoutercutrad < 0:  # arcmin
            outercutraddeg = -targetoutercutrad / 60.
        else:  # kpc
            outercutraddeg = np.degrees(targetoutercutrad / (1000 * host.distmpc))
        outermsk = targs['rhost'] > outercutraddeg
    else:
        targs = targeting.select_targets(host, faintlimit=seclimit, outercutrad=targetoutercutrad, removegalsathighz=removegalsathighz, removespecstars=removespecstars, inclspecqsos=inclspecqsos)
        outermsk = None

    print 'Found', len(targs), 'targets'
    if outermsk is not None:
        print (~outermsk).sum(), 'are in the inner zone'
    if prisecbounday != seclimit:
        print np.sum(targs['r'] < prisecbounday), 'Primaries and',np.sum((targs['r'] > prisecbounday) & (targs['r'] < seclimit)),'Secondaries'
        if outermsk is not None:
            inrtargs = targs[~outermsk]
            print np.sum(inrtargs['r'] < prisecbounday), 'Primaries and',np.sum((inrtargs['r'] > prisecbounday) & (inrtargs['r'] < seclimit)),'Secondaries are in the inner zone'



    if fibermag_faintlimit is not None:
        try:
            prifiberlimit, secfiberlimit = fibermag_faintlimit
        except TypeError:
            prifiberlimit = secfiberlimit = fibermag_faintlimit

        fibmag = targs[fibermag_type]

        #filter things that have faint fiber mags in the primary targets
        primsk = (targs['r'] < prisecbounday) & (fibmag < prifiberlimit)
        #filter things that have faint fiber mags in the secondart targets
        secmsk = (targs['r'] > prisecbounday) & (targs['r'] < seclimit) & (fibmag < secfiberlimit)

        print 'Accepting', primsk.sum(), 'primary targets and', secmsk.sum(), 'secondaries according to fiber faint limit'
        targs = targs[primsk | secmsk]

    if fibermag_brightlimit is not None:
        msk = targs[fibermag_type] > fibermag_brightlimit
        targs = targs[msk]
        print 'Filtering', msk.sum(),' targets due to fiber bright limit'

    for i, t in enumerate(targs):
        rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
        decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
        objnm = str(t['objID'])
        mag = '{0:.2f}'.format(t['r'])
        rank = prisatrank if t['r'] < prisecbounday else secsatrank
        if outermsk is not None and outermsk[i]:
            rank += 1


        tabentries.append([rastr, decstr, objnm, str(rank), 'TARGET', mag])
        catlines.append('\t'.join(tabentries[-1]))

    #calibration stars - don't obey outercutrad
    fluxtargs = select_flux_stars(host.get_sdss_catalog(), fluxrng)
    print 'Found', len(fluxtargs), 'Flux stars'
    for t in fluxtargs:

        rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
        decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
        objnm = str(t['objID'])
        mag = '{0:.2f}'.format(t['r'])

        for i in range(repeatflux):
            tabentries.append([rastr, decstr, objnm + '_' + str(i + 1), str(fluxrank), 'TARGET', mag])
            catlines.append('\t'.join(tabentries[-1]))

    #add guide stars - note these do *not* obey the outercutrad, intentionally
    guidestars = select_guide_stars(host.get_sdss_catalog(), guidestarmagrng)
    print 'Found', len(guidestars), 'guide stars'
    for t in guidestars:
        rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
        decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
        objnm = str(t['objID'])
        mag = '{0:.2f}'.format(t['r'])

        tabentries.append([rastr, decstr, objnm, '', 'guide', mag])
        catlines.append('\t'.join(tabentries[-1]))

    catstr = '\n'.join(catlines)
    if fnout:
        with open(fnout, 'w') as f:
            f.write(catstr)

    tab = Table(names=catlines[0].split('\t'), dtypes=['f', 'f', 'S50', 'S2', 'S20', 'f'])
    for e in tabentries:
        e[0] = Angle(e[0], unit=u.hour).degrees
        e[1] = Angle(e[1], unit=u.degree).degrees
        tab.add_row(e)
    return tab


def select_flux_stars(cat, magrng=(17, 17.7)):
    """
    Identifies flux calibration stars suitable for Hectospec.

    Primarily selects ~F stars based on color cut.

    Uses the SDSS criteria for specphot standards:
    0.1 < (g-r) < 0.4

    their specphot stars are 16 < g < 17.1 and "reddening standards" are
    17.1 < g < 18.5
    """
    starmsk = cat['type'] == 6 #type==6 means star
    cat = cat[starmsk]


    g = cat['g']
    r = cat['r']
    gmr = g - r

    msk = (0.1 < gmr) & (gmr < 0.4)

    minmag = min(*magrng)
    maxmag = max(*magrng)

    msk = msk & (minmag < r)& (r < maxmag)

    return cat[msk]

def select_guide_stars(cat, magrng=(14, 15)):
    msk = cat['type'] == 6 #type==6 means star

    magrng = min(*magrng), max(*magrng)

    msk = msk & (magrng[0] < cat['r']) & (cat['r'] < magrng[1])

    return cat[msk]



    return cat[msk]

def parse_cfg_file(fn):
    from astropy.coordinates import ICRSCoordinates
    from astropy.io import ascii

    intab = False
    inhdr = False

    coords = []
    targets = []
    ranks = []
    fields = []

    fi = 0
    with open(fn) as f:
        for l in f:
            if intab:
                if l.strip()=='':
                    intab = False
                else:
                    fiber, ra, dec, platex, platey, target, rank = [e.strip() for e in l.split('\t')]
                    coords.append(ICRSCoordinates(ra, dec, unit=('hr', 'deg')))
                    targets.append(int(float(target)))
                    ranks.append(int(float(rank)))
                    fields.append(fi)



            elif l.startswith('fiber'):
                hdrline = l
                inhdr = True
            elif inhdr and l.startswith('-'):
                hdrline2 = l
                inhdr = False
                intab = True
                fi += 1

    return np.array(coords), np.array(targets), np.array(ranks), np.array(fields)



