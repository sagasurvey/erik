from __future__ import division

"""
These functions are for BOSSANOVA (BOss Survey of Satellites Around Nearby Optically obserVable milky way Analogs)
"""

import numpy as np
from matplotlib import pyplot as plt

import targeting


def count_targets(hsts, verbose=True, remove_cached=True, rvir=300, targetingkwargs={}):
    """
    Generates a count of targets for each field.

    Parameters
    ----------
    hsts
        A list of `NSAHost` objects
    verbose : bool
        Whether or not to print a message when each host is examined
    remove_cached : bool
        Whether or not to remove the cached sdss catalog for each host
        after counting. This may be necessary to prevent running out of
        memory, depending on the number of hosts involved.
    rvir : float
        "virial radius" in kpc for the arcmin transform
    targetingkwargs : dict or list of dicts
        passed into  ` targeting.select_targets` if a single dictionary, otherwise
        the targeting will

    Returns
    -------
    ntargs : astropy.Table
        a table object with the names of the hosts and the target counts.

    """
    import sys
    import collections

    from astropy import table

    if isinstance(targetingkwargs, collections.Mapping):
        colnames = ['ntarg']
        targetingkwargs = [targetingkwargs.copy()]
    else:
        colnames = [('ntarg_' + t.get('colname', str(i))) for i, t in enumerate(targetingkwargs)]
        targetingkwargs = [t.copy() for t in targetingkwargs]

    for t in targetingkwargs:
        t.setdefault('outercutrad', 300)
        t.setdefault('removegama', False)
        if 'colname' in t:
            del t['colname']

    nms = []
    dists = []
    rvs = []
    cnts = [[] for t in targetingkwargs]

    for i, h in enumerate(hsts):
        if verbose:
            print 'Generating target count for', h.name, '#', i + 1, 'of', len(hsts)
            sys.stdout.flush()
        nms.append(h.name)
        dists.append(h.distmpc)
        rvs.append(h.physical_to_projected(300))

        for j, t in enumerate(targetingkwargs):
            if verbose:
                print 'Targeting parameters:', t
                sys.stdout.flush()

            tcat = targeting.select_targets(h, **t)
            cnts[j].append(len(tcat))

        if remove_cached:
            h._cached_sdss = None

    t = table.Table()
    t.add_column(table.Column(name='name', data=nms))
    t.add_column(table.Column(name='distmpc', data=dists, units='Mpc'))
    t.add_column(table.Column(name='rvirarcmin', data=rvs, units='arcmin'))

    for cnm, cnt in zip(colnames, cnts):
        t.add_column(table.Column(name=cnm, data=cnt))

    return t

_Vabs_mw_sats = {'Bootes I': -6.3099999999999987,
 'Bootes II': -2.7000000000000011,
 'Bootes III': -5.7500000000000018,
 'Canes Venatici I': -8.5900000000000016,
 'Canes Venatici II': -4.9199999999999982,
 'Canis Major': -14.389999999999999,
 'Carina': -9.1099999999999994,
 'Coma Berenices': -4.0999999999999996,
 'Draco': -8.7999999999999989,
 'Fornax': -13.44,
 'Hercules': -6.6000000000000014,
 'LMC': -18.120000000000001,
 'Leo I': -12.02,
 'Leo II': -9.8399999999999999,
 'Leo IV': -5.8400000000000016,
 'Leo V': -5.25,
 'Pisces II': -5.0,
 'SMC': -16.830000000000002,
 'Sagittarius dSph': -13.500000000000002,
 'Sculptor': -11.070000000000002,
 'Segue I': -1.5,
 'Segue II': -2.5,
 'Sextans I': -9.2700000000000014,
 'Ursa Major I': -5.5299999999999994,
 'Ursa Major II': -4.1999999999999993,
 'Ursa Minor': -8.7999999999999989,
 'Willman 1': -2.6999999999999993}

#now just assume they are all g-r=0.5, ~right for Draco ... Apply Jester+ transforms

_rabs_mw_sats = dict([(k, v + (-0.41 * (0.5) + 0.01)) for k, v in _Vabs_mw_sats.iteritems()])
_sorted_mw_rabs = np.sort(_rabs_mw_sats.values())


def count_mw_sats(h, maglim, mwsatsrmags=_sorted_mw_rabs):
    appmags = mwsatsrmags + h.distmod
    return np.sum(appmags < maglim)


def generate_count_table(hsts, fnout=None, maglims=[21, 20.5, 20], outercutrad=-90,remove_cached=True):
    from astropy.io import ascii
    from astropy import table

    targetingkwargs = []
    for m in maglims:
        targetingkwargs.append({'faintlimit': m, 'outercutrad': outercutrad, 'colname': str(m)})

    tab = count_targets(hsts, targetingkwargs=targetingkwargs, remove_cached=remove_cached)
    for m in maglims:
        satcnt = []
        for hs in hsts:
            satcnt.append(count_mw_sats(hs, m))
        tab.add_column(table.Column(name='nsat_' + str(m), data=satcnt))

    for m in maglims:
        nsatstr = 'nsat_' + str(m)
        ntargstr = 'ntarg_' + str(m)

        tab.add_column(table.Column(name='ntargpersat_' + str(m), data=tab[ntargstr] / tab[nsatstr]))

    if fnout:
        ascii.write(tab, fnout)

    return tab
