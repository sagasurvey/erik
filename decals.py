import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

def make_cutout_comparison_table(dcat, dhtml=True, doprint=True):
    """
    Produces a table comparing DECaLS and SDSS objects side-by-side

    `dcat` should be a *DECaLS* catalog, not SDSS
    """
    de_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer=decals-dr3&pixscale=0.1&bands=grz'
    sd_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer=sdssco&pixscale=0.1&bands=gri'

    if doprint:
        print('put this into http://skyserver.sdss.org/dr13/en/tools/chart/listinfo.aspx')
        print('name ra dec')

    tabrows = []
    for row in dcat:
        dviewurl = 'http://legacysurvey.org/viewer?ra={}&dec={}'.format(row['ra'], row['dec'])
        sviewurl = 'http://skyserver.sdss.org/dr12/en/tools/chart/navi.aspx?ra={}&dec={}'.format(row['ra'], row['dec'])
        sc = SkyCoord(row['ra'], row['dec'], unit=u.deg)
        objstr = '{}_{}<br>RA={:.4f}<br>Dec={:.4f}<br>r={:.2f}<br>sb={:.2f}<br>anymsk={}<br>allmsk={}'.format(row['brickid'],
            row['objid'], row['ra'], row['dec'], row['r'], row['sb_r_0.5'], row['decam_anymask'], row['decam_allmask'])
        deimg = '<a href="{}"><img src="{}"></a>'.format(dviewurl, de_cutout_url.format(sc))
        sdimg = '<a href="{}"><img src="{}"></a>'.format(sviewurl, sd_cutout_url.format(sc))
        tabrows.append('<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format(objstr, deimg, sdimg))
        if doprint:
            print(row['brickname']+'_'+str(row['objid']), row['ra'], row['dec'])

    htmlstr = """
    <table>

    <tr>
    <th>obj</th>
    <th>DECALS</th>
    <th>SDSS</th>
    </tr>

    {}
    </table>
    """.format('\n'.join(tabrows))

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
