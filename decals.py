import os

import six

import numpy as np

from scipy import interpolate

from astropy.utils.data import download_file
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

def make_cutout_comparison_table(dcat, dhtml=True, doprint=True, inclres=False,
                                 inclmod=False, inclsdss=True,
                                 add_annotation=[], subsample=None):
    """
    Produces a table comparing DECaLS and SDSS objects side-by-side

    `dcat` should be a *DECaLS* catalog, not SDSS
    """
    if subsample:
        dcat = dcat[np.random.permutation(len(dcat))[:subsample]]
    else:
        dcat = dact.copy()

    for fromnm, tonm in [('RA', 'ra'), ('DEC', 'dec'), ('OBJID', 'objid')]:
        if tonm not in dcat.colnames and fromnm in dcat.colnames:
            dcat[tonm] = dcat[fromnm]
    if 'brickid' not in dcat.colnames:
        dcat['brickid'] = ['notDE']*len(dcat)
        dcat['brickname'] = ['notDE']*len(dcat)

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

def mags_catalog(cat, extcorr=False):
    if 'decam_flux' in cat.colnames:
        m, e = fluxivar_to_mag_magerr(cat['decam_flux'], cat['decam_flux_ivar'])
        cat['decam_mag'] = m
        cat['decam_mag_err'] = e
        cat['mag_r'] = m[:, 2]
        cat['mag_err_r'] = e[:, 2]
        if extcorr:
            raise NotImplementedError()
    else:
        for band in 'ugrizY':
            m, e = fluxivar_to_mag_magerr(cat['flux_'+band], cat['flux_ivar_'+band])
            if extcorr:
                m + 2.5 * np.log10(cat['mw_transmission_'+band])*u.mag
            cat['mag_'+band] = m
            cat['mag_err_'+band] = e



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

decam_band_name_to_idx = {nm:i for i, nm in enumerate('ugrizy')}

def aperture_sbs_catalog(cat, bandname='r'):
    """
    Takes a DECaLS tractor catalog and adds surface brightness profiles for the
    requested band to it for each of the decals flux apertures.
    """
    bandidx = decam_band_name_to_idx[bandname]

    if 'decam_apflux' in cat.colnames:
        ap_fluxes = cat['decam_apflux'][:, bandidx, :]
    elif 'apflux_' + bandname in cat.colnames:
        ap_fluxes = cat['apflux_' + bandname]

    for ap in DECALS_AP_SIZES:
        strnm = str(ap.value).replace('.0','').replace('.','p')
        cat['sb_{}_{}'.format(bandname, strnm)] = compute_sb(ap, ap_fluxes)

def interpolate_catalog_sb(cat, bandname='r', radtype='eff',
                           sbname='sbeff_r', radname='rad_sb',
                           loopfunc=lambda x:x):
    """
    Takes a DECaLS tractor catalog and adds r-band half-light surface brightness
    to it.

    ``radtype`` can be "eff" for model-determined reff, or a angle-unit quantity
    for a fixed aperture SB

    For details/tests that this function works, see the
    "DECALS low-SB_completeness figures" notebook.
    """
    bandidx = decam_band_name_to_idx[bandname]

    if 'decam_apflux' in cat.colnames:
        r_ap_fluxes = cat['decam_apflux'][:, bandidx, :]
    elif 'apflux_' + bandname in cat.colnames:
        r_ap_fluxes = cat['apflux_' + bandname]
    else:
        raise ValueError('found no valid {}-band apflux column!'.format(bandname))
    assert r_ap_fluxes.shape[-1] == 8, 'Column does not have 8 apertures'

    expflux_r = np.empty_like(r_ap_fluxes[:, 0])
    rad = np.empty(len(r_ap_fluxes[:, 0]))
    ap_sizesv = DECALS_AP_SIZES.to(u.arcsec).value

    intr = interpolate.BarycentricInterpolator(ap_sizesv, [0]*len(ap_sizesv))

    if loopfunc == 'ProgressBar':
        from astropy.utils.console import ProgressBar
        loopfunc = lambda x: ProgressBar(x)
    elif loopfunc == 'NBProgressBar':
        from astropy.utils.console import ProgressBar
        loopfunc = lambda x: ProgressBar(x, ipython_widget=True)
    elif loopfunc == 'tqdm':
        import tqdm
        loopfunc = lambda x: tqdm.tqdm(x)
    elif loopfunc == 'tqdm_notebook':
        import tqdm
        loopfunc = lambda x: tqdm.tqdm_notebook(x)

    for i in loopfunc(range(len(r_ap_fluxes))):
        f = r_ap_fluxes[i]

        if radtype != 'eff':
            r = radtype
        elif cat['type'][i] == 'PSF ':
            if 'decam_psfsize' in cat.colnames:
                r = cat['decam_psfsize'][i, bandidx]
            else:
                r = cat['psfsize_' + bandname][i]
        elif cat['type'][i] == 'DEV ':
            if 'shapeDev_r' in cat.colnames:
                r = cat['shapeDev_r' ][i]
            else:
                # DR4 changed to all lower-case... WWHHHHYYY!!?!?!??!?!?!?!?
                r = cat['shapedev_r'][i]
        else:
            if 'shapeExp_r' in cat.colnames:
                r = cat['shapeExp_r'][i]
            else:
                # DR4 changed to all lower-case... WWHHHHYYY!!?!?!??!?!?!?!?
                r = cat['shapeexp_r'][i]

        intr.set_yi(f)
        expflux_r[i] = intr(r)
        rad[i] = r

    cat[sbname] = compute_sb(rad*u.arcsec, np.array(expflux_r))
    cat[radname] = rad*u.arcsec


_NERSC_BASE = 'http://portal.nersc.gov/project/cosmo/data/legacysurvey/'
def brickname_to_catalog_url(brickname, drnum, baseurl=_NERSC_BASE):
    catalogurl = baseurl + '/dr' + str(drnum) + '/tractor'

    dirnm = brickname[:3]

    return catalogurl + '/' + dirnm + '/tractor-' + brickname + '.fits'

def download_bricks(bricktab, drnum, baseurl=_NERSC_BASE, update_tab=True):
    dled = []
    local_fns = []
    for row in bricktab:
        url = brickname_to_catalog_url(row['brickname'], 5)
        local_fn = os.path.join('decals_dr5', 'catalogs', 'tractor-{}.fits'.format(row['brickname']))
        if os.path.isfile(local_fn):
            print(local_fn, 'already exists, skipping...')
        else:
            dlfn = download_file(url)
            os.rename(dlfn, local_fn)
            dled.append(local_fn)
        local_fns.append(local_fn)
    if update_tab:
        bricktab['local_fn'] = local_fns
    return dled



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
    Find all the bricks around a set of hosts.

    Parameters
    ----------
    hostlst
        A list of host objects or a table with host info
    bricksdr
        The table from the "survey-bricks.fits.gz" DECaLS file
    brickstab
        The table from the "survey-bricks-dr#.fits.gz" DECaLS file
    environfactor
        The area to check.  If a raw number, will use that multiple of the
        environs, assuming `hostlst` is a Host object.  Otherwise, must be
        a quantity either distance or angular
    brick_check_func
        the function that maps a list of brick corner booleans to a single
        boolean

    Returns
    -------
    joinedtab
        A joined table with the contents of bricksdr and brickstab, with only
        bricks in the `environfactor` area
    """

    subbricks = brickstab['BRICKNAME', 'RA1', 'RA2', 'DEC1', 'DEC2']
    subbricks.rename_column('BRICKNAME', 'brickname')
    joined = table.join(subbricks, bricksdr)

    if isinstance(hostlst, table.Table):
        schosts = hostlst['coord']

        hostnames = np.choose(hostlst['SAGA_name']=='', [hostlst['SAGA_name'],
                              np.char.add('NSA', hostlst['NSAID'].astype('U'))])
        pgcs = hostnames == 'NSA-1'
        hostnames[pgcs] = np.char.add('PGC', hostlst['PGC#'].astype('U'))[pgcs]

        if hasattr(environfactor, 'unit'):
            if environfactor.unit.is_equivalent(u.kpc):
                environfactor = environfactor/(hostlst['distance']*u.Mpc)
                environrad = environfactor.to(u.arcmin, u.dimensionless_angles())
            elif environfactor.unit.is_equivalent(u.deg):
                environrad = environfactor.to(u.arcmin)
            else:
                raise u.UnitsError('environfactor not a distance or angle')
        else:
            raise TypeError('environfactor must be a Quantity with table hosts')

    else:
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


def show_decals_objects_in_nb(cat, nrows=3, dr=3, subsample=None, info_cols=[],
                              sdss_link=False, show_reticle=True):
    """
    Produces a table showing DECaLS cutouts from a decals catalog (must have
    'ra', 'dec', and 'objname' columns).

    If `dr` is 'fromcatalog' and catalog has a 'dr' the row will be used to choose which
    DR to show.

    ``subsample`` can be a number to randomly subsample to that many
    """
    from IPython import display

    if subsample:
        cat = cat[np.random.permutation(len(cat))[:subsample]]

    de_cutout_url = 'http://legacysurvey.org/viewer/jpeg-cutout/?ra={0.ra.deg}&dec={0.dec.deg}&layer={dr}&pixscale=0.1&bands=grz'
    sdss_cutout_url = 'http://skyserver.sdss.org/dr{dr}/SkyserverWS/ImgCutout/getjpeg?ra={0.ra.deg}&dec={0.dec.deg}&width=256&height=256&scale=0.15'

    entry_info = []
    rows = [[]]
    for row in cat:
        if 'RA' in row.colnames:
            sc = SkyCoord(row['RA'], row['DEC'], unit=u.deg)
        else:
            sc = SkyCoord(row['ra'], row['dec'], unit=u.deg)

        templ = ('<td style="text-align:center">{objstr}<a href="{dviewurl}"><{tagtype} '
                 'id="cnv{objnm}" width="256" height="256" src="{imgurl}"></a></td>')
        objnm = '{}'.format(row['objname' if 'objname' in row.colnames else 'name'])
        objstr = objnm

        needbr = True
        for colnm in info_cols:
            objstr += '<br>{} = {}'.format(colnm, row[colnm])
            needbr = False
        if needbr:
            objstr += '<br>'

        if sdss_link:
            sdss_navi_url = 'http://skyserver.sdss.org/dr13/en/tools/chart/navi.aspx?ra={0.ra.deg}&dec={0.dec.deg}'.format(sc)
            objstr += '<br><a href="{0}">SDSS DR13 Navigate</a>'.format(sdss_navi_url)

        dviewurl = 'http://legacysurvey.org/viewer?ra={0.ra.deg}&dec={0.dec.deg}&zoom=16'.format(sc)

        if dr == 'fromcatalog':
            thisdr = row['dr']
        else:
            thisdr = dr

        use_sdss = False
        if hasattr(thisdr, 'startswith') and thisdr.startswith('sdss'):
            use_sdss = True
        else:
            if thisdr == 4:
                thisdr = 'mzls+bass-dr4'
            elif thisdr > 2:
                thisdr = 'decals-dr' + str(thisdr)
        if use_sdss:
            imgurl = sdss_cutout_url.format(sc, dr=thisdr[4:])
        else:
            imgurl = de_cutout_url.format(sc, dr=thisdr)
        tagtype = 'canvas' if show_reticle else 'img'
        entry = templ.format(**locals())
        entry_info.append((objnm, imgurl))

        rows[-1].append(entry)
        if len(rows[-1]) >= nrows:
            rows[-1] = '<tr>' + ''.join(rows[-1]) + '</tr>'
            rows.append([])
    if not isinstance(rows[-1], six.string_types):
        rows[-1] = '<tr>' + ''.join(rows[-1]) + '</tr>'


    script = ''
    if show_reticle:
        script = """
        function draw(nm, imgurl){
          var img = new Image(256, 256);
          img.src = imgurl;
          img.onload = function(){
            var cnvs = document.getElementById("cnv"+nm);
            cnvs.style.z_index = 2;
            console.log("aha" + nm);

            var ctx = cnvs.getContext("2d");
            ctx.drawImage(img, 0, 0);
            ctx.beginPath();
            ctx.arc(128, 128, 15, 0, 2 * Math.PI, false);
            ctx.lineWidth = 1.5;
            ctx.strokeStyle = '#ddaa00';
            ctx.stroke();
          }


        }

        """
        for nm, url in entry_info:
            script += 'draw("{}","{}");\n'.format(nm, url)

    htmlstr = """
    <table>
    {}
    </table>

    <script type="text/javascript">
    {}
    </script>
    """.format('\n'.join(rows), script)

    return display.HTML(htmlstr)
