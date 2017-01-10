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
        objstr = '{}_{}<br>RA={:.4f}<br>Dec={:.4f}<br>r={:.2f}<br>sb={:.2f}'.format(row['brickid'], row['objid'], row['ra'], row['dec'], row['r'], row['sb_r_0.5'])
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
