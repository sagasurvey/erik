import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import table


def band_to_idx(band):
    return list('ugrizy').index(band)


def subselect_aperture(apmags, band, apsize=1*u.arcsec):
    if apsize is None:
        return apmags[:, band_to_idx(band), :]
    else:
        return apmags[:, band_to_idx(band), DECALS_AP_SIZES == apsize][:, 0]

def make_cutout_comparison_table(dcat, dhtml=True, doprint=True, inclres=False, inclmod=False, inclsdss=True, add_annotation=[]):
    """
    Produces a table comparing DECaLS and SDSS objects side-by-side

    `dcat` should be a *DECaLS* catalog, not SDSS
    """
    de_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer=decals-dr3&pixscale=0.1&bands=grz'
    demod_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer=decals-dr3-model&pixscale=0.1&bands=grz'
    deres_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer=decals-dr3-resid&pixscale=0.1&bands=grz'
    sd_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer=sdssco&pixscale=0.1&bands=gri'

    if doprint:
        print('put this into http://skyserver.sdss.org/dr13/en/tools/chart/listinfo.aspx')
        print('name ra dec')

    tabrows = []
    for row in dcat:
        dviewurl = 'http://legacysurvey.org/viewer?ra={}&dec={}'.format(row['ra'], row['dec'])
        sviewurl = 'http://skyserver.sdss.org/dr12/en/tools/chart/navi.aspx?ra={}&dec={}'.format(row['ra'], row['dec'])
        sc = SkyCoord(row['ra'], row['dec'], unit=u.deg)

        inforows = ['{}_{}'.format(row['brickid'], row['objid'])]
        for infocolnm in ['ra', 'dec', 'r', 'sb_r_0.5', 'decam_anymask', 'decam_allmask'] + list(add_annotation):
            if infocolnm in row.colnames:
                inforows.append('{}={}'.format(infocolnm, row[infocolnm]))
        objstr = '<br>'.join(inforows)

        deimg = '<a href="{}"><img src="{}"></a>'.format(dviewurl, de_cutout_url.format(sc))
        sdimg = '<a href="{}"><img src="{}"></a>'.format(sviewurl, sd_cutout_url.format(sc))

        headerelems = ['obj', 'DECALS']
        imgs = [objstr, deimg]

        if inclmod:
            imgs.append('<a href="{}"><img src="{}"></a>'.format(sviewurl, demod_cutout_url.format(sc)))
            headerelems.append('DECALS models')
        if inclres:
            imgs.append('<a href="{}"><img src="{}"></a>'.format(sviewurl, deres_cutout_url.format(sc)))
            headerelems.append('DECALS residuals')
        if inclsdss:
            imgs.append('<a href="{}"><img src="{}"></a>'.format(sviewurl, sd_cutout_url.format(sc)))
            headerelems.append('SDSS')


        tabrows.append('<tr><td>' + '</td><td>'.join(imgs) + '</td></tr>')
        if doprint:
            print(row['brickname']+'_'+str(row['objid']), row['ra'], row['dec'])

    headerelems

    htmlstr = """
    <table>

    <tr>
    <th>{}</th>
    </tr>

    {}
    </table>
    """.format('</th><th>'.join(headerelems), '\n'.join(tabrows))

    if dhtml:
        from IPython import display
        return display.HTML(htmlstr)
    else:
        return htmlstr

def fluxivar_to_mag_magerr(flux, ivar, deunit=True):
    """
    returns mag, mag_err as Quantities
    """
    flux = u.Quantity(flux, copy=False)
    ivar = u.Quantity(ivar, copy=False)
    if deunit:
        flux = u.Quantity(flux.value, u.dimensionless_unscaled)
        ivar = u.Quantity(ivar.value, u.dimensionless_unscaled)

    mag = np.array(22.5 - 2.5*np.log10(flux))
    flux_err = ivar**-0.5
    mag_err = 2.5/np.log(10) * flux_err/flux
    return mag*u.mag, mag_err.value*u.mag


DECALS_AP_SIZES = [0.5,0.75,1.0,1.5,2.0,3.5,5.0,7.0] * u.arcsec

def compute_sb(rad, apflux):
    """
    `rad` is the radius of the aperture
    `apflux` is the flux within the aperture

    returns SB as Quantity
    """
    if len(apflux.shape)==2:
        idxs = np.where(DECALS_AP_SIZES==rad)[0]
        assert len(idxs)==1, 'No aperture with size {}'.format(rad)

        apflux = apflux[:, idxs[0]]
    A = 2.5*np.log10(np.pi*(rad.to(u.arcsec).value)**2)
    return np.array(22.5 - 2.5*np.log10(apflux) + A) * u.mag * u.arcsec**-2


_NERSC_BASE = 'http://portal.nersc.gov/project/cosmo/data/legacysurvey/'
def brickname_to_catalog_url(brickname, drnum, baseurl=_NERSC_BASE):
    catalogurl = baseurl + '/dr' + str(drnum) + '/tractor'

    dirnm = brickname[:3]

    return catalogurl + '/' + dirnm + '/tractor-' + brickname + '.fits'


def show_in_decals(ra, dec=None,
                   urltempl='http://legacysurvey.org/viewer?ra={}&dec={}',
                   show=True):
    """
    Get the url for the supplied object to view in decals, and possibly actually
    display it in the browser.

    ra/dec can be in degrees, or if just ra is given it will be taken to be a
    host object, a singular SkyCoord,
    """
    import webbrowser

    if dec is None:
        obj = ra
        if hasattr(obj, 'coords'):
            coo = obj.coords
        elif hasattr(obj, 'ra') and hasattr(obj, 'dec'):
            coo = obj
        else:
            coo = SkyCoord.guess_from_table(obj, unit=u.deg)
    else:
        coo = SkyCoord(ra, dec, unit=u.deg)

    if coo.shape != tuple():
        raise ValueError('cant use show_in_decals with non-scalar input')

    url = urltempl.format(coo.ra.deg, coo.dec.deg)
    if show:
        webbrowser.open(url)
    return url


def find_host_bricks(hostlst, bricksdr, brickstab, environfactor=1.2, brick_check_func=np.any):
    """
    Find all the bricks around a set of hosts.  ``environfactor`` can be a
    multiple of the environs in the host objects, or a fixed value with degree
    units.
    """
    subbricks = brickstab['BRICKNAME', 'RA1', 'RA2', 'DEC1', 'DEC2']
    subbricks.rename_column('BRICKNAME', 'brickname')
    joined = table.join(subbricks, bricksdr)

    schosts = SkyCoord([h.coords for h in hostlst])
    hostnames = np.array([h.name for h in hostlst])
    if hasattr(environfactor, 'unit'):
        environrad = [environfactor.value for h in hostlst]*environfactor.unit
    else:
        environrad = [environfactor*h.environsarcmin for h in hostlst]*u.arcmin

    # we do the ravel because match_to_catalog_sky works best with 1d
    brickras = np.array([joined['RA1'], joined['RA1'], joined['RA2'], joined['RA2']]).ravel()
    brickdecs = np.array([joined['DEC1'], joined['DEC2'], joined['DEC1'], joined['DEC2']]).ravel()
    brickedge_scs = SkyCoord(brickras, brickdecs, unit=u.deg)

    idx, d2d, _ = brickedge_scs.match_to_catalog_sky(schosts)
    # reshape so that each *brick* is represented
    bricksidx = idx.reshape(4, idx.size//4)
    bricksd2d = d2d.reshape(4, idx.size//4)

    bricksin = brick_check_func(bricksd2d < environrad[bricksidx], axis=0)

    res = joined[bricksin]

    # do this matching b/c it's *possible* multiple matches exist if a brick is near-equidistant
    closest_host_idx = SkyCoord.guess_from_table(res, unit=u.deg).match_to_catalog_sky(schosts)[0]
    res['closest_host_idx'] = closest_host_idx
    res['closest_host_name'] = hostnames[closest_host_idx]

    return res
