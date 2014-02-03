#!/usr/bin/env python
"""
Script/functions for generating the SAGA master host list.

Requires these input data files:
* LEDA.csv (from the EDD: http://edd.ifa.hawaii.edu/ - the "LEDA" table with all columns, comma-separated)
* 2MRS.csv (from the EDD: http://edd.ifa.hawaii.edu/ - the "2MRS K<11.75 " table with all columns, comma-separated)
* EDD.csv (from the EDD: http://edd.ifa.hawaii.edu/ - the "EDD Distances" table with all columns, comma-separated)
* KKnearbygal.csv (from the EDD: http://edd.ifa.hawaii.edu/ - the "Updated Nearby Galaxy Catalog" table with all columns, comma-separated)
* nsa_v0_1_2.fits (from http://www.nsatlas.org/)


Note that this script can be rather memory-intensive - it works fine for me with
16GB of memory, but much less than that might get pretty slow...
"""
from __future__ import division, print_function

import numpy as np

from astropy import units as u
from astropy.cosmology import WMAP9  # for Z->distances


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


def filter_catalog(mastercat, vcut=3000*u.km/u.s, musthavenirphot=False):
    """
    Removes entries in the simplified catalog to meet the  master catalog selection
    criteria
    """

    #possible point of confusion:`msk` is True where we want to *keep* so
    #something, because it is used at the end as a bool index into the catalog.
    #The MaskedColumn `mask` attribute is the opposite - True is *masked*

    #filter anything without an RA or Dec
    msk = ~(mastercat['RA'].mask | mastercat['Dec'].mask)

    #also remove everything without a distance - this includes all w/velocities,
    #because for those the distance comes from assuming hubble flow
    msk = msk & (~mastercat['distance'].mask)


    # remove everything that has a velocity > `vcut`
    if vcut is not None:
        msk = msk & (mastercat['vhelio'] < vcut.to(u.km/u.s).value)

    if musthavenirphot:
        msk = msk & ((~mastercat['i'].mask) | (~mastercat['z'].mask) | (~mastercat['I'].mask) | (~mastercat['K'].mask))

    return mastercat[msk]


def simplify_catalog(mastercat, quickld=True):
    """
    Removes most of the unnecessary columns from the master catalog and joins
    fields where relevant

    Parameters
    ----------
    mastercat : astropy.table.Table
        The table from generate_catalog
    quickld : bool
        If True, means do the "quick" version of the luminosity distance
        calculation (takes <1 sec as opposed to a min or so, but is only good
        to a few kpc)
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

    #Names/IDs:
    pgc = mastercat['pgc'].copy()
    pgc.mask = mastercat['pgc'] < 0
    tab.add_column(table.MaskedColumn(name='PGC#', data=pgc))
    tab.add_column(table.MaskedColumn(name='NSAID', data=mastercat['NSAID']))
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
    tab.add_column(table.MaskedColumn(name='othername', data=names))

    #After this, everything should have either an NSAID, a PGC#, or a name (or more than one)

    #VELOCITIES/redshifts
    #start with LEDA
    vs = mastercat['v'].astype(float)
    v_errs = mastercat['e_v'].astype(float)

    #Now add vhelio from the the EDD
    eddvhel = mastercat['Vhel_eddkk']
    vs[~eddvhel.mask] = eddvhel[~eddvhel.mask]
    #EDD has no v-errors, so mask them
    v_errs[~eddvhel.mask] = 0
    v_errs.mask[~eddvhel.mask] = True

    #then the NSA, if available
    vs[~mastercat['ZDIST'].mask] = mastercat['ZDIST'][~mastercat['ZDIST'].mask] * ckps
    v_errs[~mastercat['ZDIST_ERR'].mask] = mastercat['ZDIST_ERR'][~mastercat['ZDIST_ERR'].mask] * ckps

    #finally, KK when present
    kkvh = mastercat['Vh']
    vs[~kkvh.mask] = kkvh[~kkvh.mask]
    #KK has no v-errors, so mask them
    v_errs[~kkvh.mask] = 0
    v_errs.mask[~kkvh.mask] = True

    #DISTANCES
    dist = mastercat['Dist_edd'].copy()
    dist[~mastercat['Dist_kk'].mask] = mastercat['Dist_kk'][~mastercat['Dist_kk'].mask]

    #for those *without* EDD or KK, use the redshift's luminosity distance
    premsk = dist.mask.copy()
    zs = vs[premsk]/ckps
    if quickld:
        ldx = np.linspace(zs.min(), zs.max(), 1000)
        ldy = WMAP9.luminosity_distance(ldx).to(u.Mpc).value
        ld = np.interp(zs, ldx, ldy)
    else:
        ld = WMAP9.luminosity_distance(zs).to(u.Mpc).value

    dist[premsk] = ld
    dist.mask[premsk] = vs.mask[premsk]

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


def load_master_catalog(fn='mastercat.dat'):
    from astropy.io import ascii
    return ascii.read(fn, delimiter=',')


def x_match_tests(cattomatch, tol=1*u.arcsec):
    """
    Does a bunch of cross-matches with other catalogs. `cattomatch` must be an
    `ICRS` object or a table.

    This depends on having a bunch of other catalogs in the current directory
    that aren't in the git repo, so you may want to just ask Erik if you want
    to run this.
    """
    import RC3
    from astropy.io import ascii
    from astropy.coordinates import ICRS

    if cattomatch.__class__.__name__.lower() == 'table':
        ra, dec = cattomatch['RA'], cattomatch['Dec']
        cattomatch = ICRS(u.Unit(ra.unit)*ra.view(np.ndarray), u.Unit(dec.unit)*dec.view(np.ndarray))

    rc3, rc3_coo = RC3.load_rc3()
    rc3wv = rc3[~rc3['cz'].mask]
    rc3wv_coo = rc3_coo[~rc3['cz'].mask]

    a3de = ascii.read('atlas3d_e.dat', data_start=3, format='fixed_width')
    a3dsp = ascii.read('atlas3d_sp.dat', data_start=3, format='fixed_width')
    a3de_coo = ICRS(u.deg*a3de['RA'], u.deg*a3de['DEC'])
    a3dsp_coo = ICRS(u.deg*a3dsp['RA'], u.deg*a3dsp['DEC'])

    nsah = ascii.read('hosts.dat')
    nsah_coo = ICRS(u.hourangle*nsah['RA'], u.deg*nsah['DEC'])

    rc3_idx, rc3_d = rc3_coo.match_to_catalog_sky(cattomatch)[:2]
    rc3wv_idx, rc3wv_d = rc3wv_coo.match_to_catalog_sky(cattomatch)[:2]
    a3de_idx, a3de_d = a3de_coo.match_to_catalog_sky(cattomatch)[:2]
    a3dsp_idx, a3dsp_d = a3dsp_coo.match_to_catalog_sky(cattomatch)[:2]
    nsah_idx, nsah_d = nsah_coo.match_to_catalog_sky(cattomatch)[:2]

    rc3_nomatch = rc3_d > tol
    rc3wv_nomatch = rc3wv_d > tol
    a3de_nomatch = a3de_d > tol
    a3dsp_nomatch = a3dsp_d > tol
    nsah_nomatch = nsah_d > tol

    print('RC3 non-matches: {0} of {1}'.format(np.sum(rc3_nomatch), len(rc3_nomatch)))
    print('RC3 with v non-matches: {0} of {1}'.format(np.sum(rc3wv_nomatch), len(rc3wv_nomatch)))
    print('ATLAS3D E non-matches: {0} of {1}'.format(np.sum(a3de_nomatch), len(a3de_nomatch)))
    print('ATLAS3D Spiral non-matches: {0} of {1}'.format(np.sum(a3dsp_nomatch), len(a3dsp_nomatch)))
    print('NSA Hosts non-matches: {0} of {1}'.format(np.sum(nsah_nomatch), len(nsah_nomatch)))

    dct = {}
    for nm in ('rc3', 'rc3wv', 'a3de', 'a3dsp', 'nsah'):
        dct[nm+'_match'] = (locals()[nm])[~locals()[nm + '_nomatch']]
        dct[nm+'_nomatch'] = (locals()[nm])[locals()[nm + '_nomatch']]
        dct[nm+'_catidx'] = (locals()[nm+'_idx'])
    return dct


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument('-q', '--quiet', action='store_true')
    p.add_argument('outfn', nargs='?', help="If this is not given, the catalog won't be saved", default=None)
    args = p.parse_args()

    if args.quiet:
        #silences the print function
        print = lambda s: None

    print("Loading LEDA catalog...")
    leda = load_edd_csv('LEDA.csv')
    print("Loading 2MASS catalog...")
    twomass = load_edd_csv('2MRS.csv')
    print("Loading EDD catalog...")
    edd = load_edd_csv('EDD.csv')
    print("Loading KK nearby catalog...")
    kknearby = load_edd_csv('KKnearbygal.csv')
    nsa = load_nsa()

    #these variables are just for convinience in interactive work
    cats = [leda, twomass, edd, kknearby, nsa]
    eddcats = [leda, twomass, edd, kknearby]

    print('Generating master catalog...')
    mastercat0 = generate_catalog(*cats)
    print('Simplifying master catalog...')
    mastercat1 = simplify_catalog(mastercat0)
    print('Filtering master catalog...')
    mastercat = filter_catalog(mastercat1)

    if args.outfn is not None:
        print('Writing master catalog to {outfn}...'.format(**locals()))

        oldmpo = str(np.ma.masked_print_option)
        try:
            np.ma.masked_print_option.set_display('')
            mastercat.write(args.outfn, format='ascii', delimiter=',')
        finally:
            np.ma.masked_print_option.set_display(oldmpo)
