#!/usr/bin/env python
"""
Script/functions for generating the SAGA master host list.

Requires these input data files:
* 2MRS.csv (from the EDD: http://edd.ifa.hawaii.edu/)
* EDD.csv (from the EDD: http://edd.ifa.hawaii.edu/)
* LEDA.csv (from the EDD: http://edd.ifa.hawaii.edu/)
* upnearKK.csv (from the EDD: http://edd.ifa.hawaii.edu/)
* nsa_v0_1_2.fits (from http://www.nsatlas.org/)


Note that this script can be rather memory-intensive - it works fine for me with
16GB of memory, but much less than that might get pretty slow...
"""
from __future__ import division, print_function

import numpy as np

from astropy import units as u


def load_edd_csv(fn):
    from astropy.io import ascii

    return ascii.read(fn, delimiter=',', Reader=ascii.Basic, guess=False)


def load_nsa(fn='nsa_v0_1_2.fits', verbose=False):
    """
    This loads the NSA *and* breaks the FNugriz fields into distinct columns

    Not that it drops columns that have more complex dimensionality, like the
    radial profiles or stokes parameters
    """
    from astropy.io import fits
    from astropy import table

    tab = table.Table(fits.getdata(fn))

    newcols = []
    for nm in tab.colnames:
        col = tab[nm]
        if len(col.shape) == 2:
            if col.shape[1] == 7:
                for i, band in enumerate('FNugriz'):
                    newcols.append(table.Column(name=nm + '_' + band, data=col[:, i]))
            elif col.shape[1] == 5:
                for i, band in enumerate('ugriz'):
                    newcols.append(table.Column(name=nm + '_' + band, data=col[:, i]))
            else:
                if verbose:
                    print('Column', nm, "has weird dimension", col.shape, 'skipping')
        elif len(col.shape) > 2:
            if verbose:
                print('Column', nm, "has weird dimension", col.shape, 'skipping')
        else:
            newcols.append(col)

    newtab = table.Table()
    newtab.add_columns(newcols)
    return newtab


def generate_catalog(leda, twomass, edd, kknearby, nsa, matchtolerance=1*u.arcmin):
    from astropy import table
    from astropy.coordinates import ICRS

    #first join the small ones, because they both have "dist" columns
    small = table.join(edd, kknearby, keys=['pgc'], table_names=['edd', 'kk'], join_type='outer')

    #now join them with LEDA
    #we call the second one "kk" because the only thing shared with LEDA is 'm21' on the KK catalog
    ledaj = table.join(leda, small, keys=['pgc'], table_names=['leda', 'kk'], join_type='outer')

    #add in the 2mass stuff
    #call the first one "eddkk" because the shared columns are all either in the EDD or KK
    ledaj2 = table.join(ledaj, twomass, keys=['pgc'], table_names=['eddkk', '2mass'], join_type='outer')

    #now cross-match with NSA - need to match on RA/Dec because no PGC #s in NSA
    ral, decl = ledaj2['al2000'], ledaj2['de2000']
    lmsk = (~ral.mask) & (~decl.mask)
    lcoo = ICRS(u.hour * ral[lmsk], u.degree * decl[lmsk])
    nsacoo = ICRS(u.degree * nsa['RA'], u.degree * nsa['DEC'])

    idx, dd, dd3d = nsacoo.match_to_catalog_sky(lcoo)
    matchpgc = -np.ones(len(idx), dtype=int)  # non-matches get -1
    dmsk = dd < matchtolerance  # only match those with a closest neighbor w/i tol
    matchpgc[dmsk] = ledaj2['pgc'][lmsk][idx[dmsk]]

    if 'pgc' in nsa.colnames:
        nsa['pgc'] = matchpgc
    else:
        nsa.add_column(table.Column(name='pgc', data=matchpgc))

    mastercat = table.join(ledaj2, nsa, keys=['pgc'], table_names=['leda', 'nsa'], join_type='outer')

    return mastercat


def filter_master_catalog(mastercat, vcut=3000*u.km/u.s, musthavenirphot=False):
    """
    Removes entries in the simplified catalog to meet the  master catalog selection
    criteria
    """
    from astropy.constants import c

    ckps = c.to(u.km/u.s).value
    vcutkps = vcut.to(u.km/u.s).value

    vmsk = (mastercat['vhelio'] < vcutkps)

    msk = vmsk

    if musthavenirphot:
        msk = msk & ((~mastercat['i'].mask) | (~mastercat['z'].mask) | (~mastercat['I'].mask) | (~mastercat['K'].mask))

    return mastercat[msk]


def simplify_master_catalog(mastercat):
    """
    Removes most of the unnecessary columns from the master catalog and joins
    fields where relevant
    """
    from astropy import table

    from astropy.constants import c

    ckps = c.to(u.km/u.s).value

    tab = table.Table()

    #RADEC:
    ras = mastercat['al2000']*15
    ras[~mastercat['RA'].mask] = mastercat['RA'][~mastercat['RA'].mask]
    decs = mastercat['de2000']
    decs[~mastercat['DEC'].mask] = mastercat['DEC'][~mastercat['DEC'].mask]

    tab.add_column(table.MaskedColumn(name='RA', data=ras, units=u.deg))
    tab.add_column(table.MaskedColumn(name='Dec', data=decs, units=u.deg))

    #NAMES/IDs:
    #do these in order of how 'preferred' the object name is.
    nameorder = ('Objname', 'Name_eddkk', 'objname', 'Name_2mass')  # this is: EDD, KK, LEDA, 2MASS
    #need to figure out which has the *largest* name strings, because we have a fixed number of characters
    largestdt = np.dtype('S1')
    for nm in nameorder:
        if mastercat.dtype[nm] > largestdt:
            largestdt = mastercat.dtype[nm]
            largestdtnm = nm
    names = mastercat[largestdtnm].copy()  # these will all be overwritten - just use it for shape
    for nm in nameorder:
        msk = ~mastercat[nm].mask
        names[msk] = mastercat[nm][msk]
    tab.add_column(table.MaskedColumn(name='Name', data=names))

    pgc = mastercat['pgc'].copy()
    mastercat.mask = mastercat['pgc'] < 0
    tab.add_column(table.MaskedColumn(name='PGCNUM', data=pgc))

    tab.add_column(table.MaskedColumn(name='NSAID', data=mastercat['NSAID']))

    #VELOCITIES/redshifts and distances
    #start with LEDA
    vs = mastercat['v']
    v_errs = mastercat['e_v']
    #prefer NSA if available
    vs[~mastercat['ZDIST'].mask] = mastercat['ZDIST'][~mastercat['ZDIST'].mask] * ckps
    v_errs[~mastercat['ZDIST_ERR'].mask] = mastercat['ZDIST_ERR'][~mastercat['ZDIST_ERR'].mask] * ckps

    dist = mastercat['Dist_edd']
    dist[~mastercat['Dist_kk'].mask] = mastercat['Dist_kk'][~mastercat['Dist_kk'].mask]
    #for those without EDD or KK, use the redshift distance
    #msk = dist.mask.copy()
    #dist[msk] = luminosity_distance(vs/ckps)
    #dist.mask[msk] = False

    distmod = 5 * np.log10(dist) + 25  # used in phot section

    tab.add_column(table.MaskedColumn(name='vhelio', data=vs))
    tab.add_column(table.MaskedColumn(name='vhelio_err', data=v_errs))
    tab.add_column(table.MaskedColumn(name='distance', data=dist, units=u.Mpc))


    #NIR PHOTOMETRY
    tab.add_column(table.MaskedColumn(name='i', data=mastercat['ABSMAG_i'] + distmod))
    tab.add_column(table.MaskedColumn(name='z', data=mastercat['ABSMAG_z'] + distmod))
    tab.add_column(table.MaskedColumn(name='I', data=mastercat['it']))
    tab.add_column(table.MaskedColumn(name='K', data=mastercat['K_tc']))

    return tab

if __name__ == '__main__':
    print("Loading LEDA catalog...")
    leda = load_edd_csv('LEDA.csv')
    print("Loading 2MASS catalog...")
    twomass = load_edd_csv('2MRS.csv')
    print("Loading EDD catalog...")
    edd = load_edd_csv('EDD.csv')
    print("Loading KK nearby catalog...")
    kknearby = load_edd_csv('upnearKK.csv')
    nsa = load_nsa()

    #these variables are just for convinience in interactive work
    cats = [leda, twomass, edd, kknearby, nsa]
    eddcats = [leda, twomass, edd, kknearby]

    print('Generating master catalog...')
    mastercat = generate_catalog(*cats)
    print('Simplifying master catalog...')
    mastercat = simplify_master_catalog(mastercat)
    print('Filtering master catalog...')
    mastercat = filter_master_catalog(mastercat)

    outfn = 'mastercat.dat'
    #print('Writing master catalog to {outfn}...'.format(**locals()))
    #mastercat.write(outfn, format='ascii')
