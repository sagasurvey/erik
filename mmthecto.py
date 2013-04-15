"""
This file contains functions for targeting objects as part of the "distant
local group" project using MMT/Hectospec.
"""
import targeting


def make_catalog(host, fnout=None, targetfaintlim=21, targetoutercutrad=300,
                 fluxrank=10, usefaintflux=False, guidestarmagrng= (14, 15)):
    """
    Generates a catalog for the Hectospec "XFITFIBS" program.

    targets are all rank 1
    `stdrank` = None means don't add any flux standards

    Returns the catalog string and writes it to `fnout` unless `fnout` is None
    """
    from astropy.coordinates import Angle

    catlines = ['ra\tdec\tobject\trank\ttype\tmag',
                '--\t---\t------\t----\t----\t---']

    #add the host as the first entry in rank 1
    rastr = Angle(host.ra, 'deg').format('hr', sep=':', precision=3)
    decstr = Angle(host.dec, 'deg').format('deg', sep=':', precision=3)
    mag = '{0:.2f}'.format(host.r + host.distmod)
    catlines.append('\t'.join([rastr, decstr, host.name, '1', 'TARGET', mag]))


    targs = targeting.select_targets(host, faintlimit=targetfaintlim, outercutrad=targetoutercutrad)

    for t in targs:
        rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
        decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
        objnm = str(t['objID'])
        mag = '{0:.2f}'.format(t['r'])

        catlines.append('\t'.join([rastr, decstr, objnm, '1', 'TARGET', mag]))

    if fluxrank is not None:
        fluxtargs = select_flux_stars(host.get_sdss_catalog(), usefaintflux)
        for t in fluxtargs:

            rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
            decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
            objnm = str(t['objID'])
            mag = '{0:.2f}'.format(t['r'])

            catlines.append('\t'.join([rastr, decstr, objnm, str(fluxrank), 'TARGET', mag]))


    for t in select_guide_stars(host.get_sdss_catalog(), guidestarmagrng):
        rastr = Angle(t['ra'], 'deg').format('hr', sep=':', precision=3)
        decstr = Angle(t['dec'], 'deg').format('deg', sep=':', precision=3)
        objnm = str(t['objID'])
        mag = '{0:.2f}'.format(t['r'])

        catlines.append('\t'.join([rastr, decstr, objnm, '', 'guide', mag]))

    catstr = '\n'.join(catlines)
    if fnout:
        with open(fnout, 'w') as f:
            f.write(catstr)
    return catstr


def select_flux_stars(cat, use_fainter=True):
    """
    Identifies flux calibration stars suitable for Hectospec.

    Primarily selects ~F stars based on color cut.

    Uses the SDSS criteria for specphot standards:
    0.1 < (g-r) < 0.4
    16 < g < 17.1

    Can optionally include the stars that are "reddening standards":
    0.1 < (g-r) < 0.4
    17.1 < g < 18.5
    """
    starmsk = cat['type'] == 6 #type==6 means star
    cat = cat[starmsk]


    g = cat['g']
    r = cat['r']
    gmr = g - r

    msk = (0.1 < gmr) & (gmr < 0.4)
    msk = msk & (16 < r) & (r < 18.5)

    if not use_fainter:
        msk = msk & (r < 17.1)

    return cat[msk]

def select_guide_stars(cat, magrng=(14, 15)):
    msk = cat['type'] == 6 #type==6 means star

    magrng = min(*magrng), max(*magrng)

    msk = msk & (magrng[0] < cat['r']) & (cat['r'] < magrng[1])

    return cat[msk]



    return cat[msk]
