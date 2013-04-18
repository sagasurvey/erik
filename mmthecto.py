"""
This file contains functions for targeting objects as part of the "distant
local group" project using MMT/Hectospec.
"""
import targeting
from astropy import units as u


def make_catalog(host, fnout=None, targetfaintlim=21., targetoutercutrad=300,
                 fluxrng=(17., 17.7), repeatflux=1, guidestarmagrng= (14, 15),
                 removegalsathighz=False, inclspecqsos=True, removespecstars=False):
    """
    Generates a catalog for the Hectospec "XFITFIBS" program.

    `targetfaintlim` can be just a mag to make everything "primary", or
    a 2-tuple that is (prisecboundary, faintestsec)

    `usefaintflux` means go down to 18.5 instead of just 17.1 .  `repeatflux` is
    the number of times to repeat each flux star (e.g., the number of expected
    configs).

    Writes catalog to `fnout` if `fnout` is not None

    Returns the catalog as a `Table`
    """
    from astropy.coordinates import Angle
    from astropy.table import Table

    fluxrank = 1
    hostrank = 2
    prisatrank = 3
    secsatrank = 4


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
    targs = targeting.select_targets(host, faintlimit=seclimit, outercutrad=targetoutercutrad, removegalsathighz=removegalsathighz, removespecstars=removespecstars, inclspecqsos=inclspecqsos)

    print 'Found', len(targs), 'targets'
    for t in targs:
        rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
        decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
        objnm = str(t['objID'])
        mag = '{0:.2f}'.format(t['r'])
        rank = str(prisatrank if t['r'] < prisecbounday else secsatrank)

        tabentries.append([rastr, decstr, objnm, rank, 'TARGET', mag])
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
