from __future__ import print_function

from warnings import warn

import numpy as np
import SAGA as saga  # PEP8 name

from astropy import table


_bino_typemap ={'target': 1, 'sky': 2, 'standard': 3}


def write_bino_input(catalog, fn, overwrite=False):
    """
    Writes out a binospec input file. See
    http://mingus.mmto.arizona.edu/~bjw/mmt/binomask_info.html for more details.

    Note that if "mmtbino_type" is in the catalog, it will be used as "type".

    Parameters
    ----------
    catalog: astropy Table or list of Tables
        Must have 'coord' at a minimum (a SkyCoord). Recommended:
        'name','magnitude','priority','type'. 'priority' is 1 highest. 'type'
        column entries must be either “target”, “sky” or “standard”.

    Returns
    -------
    written_table
        The table that was actually written (with invalid cols pruned)
    """
    if not isinstance(catalog, table.Table):
        # assume list
        subtabs = []
        for elem in catalog:
            subtabs.append(write_bino_input(elem, None))
        to_write_catalog = table.vstack(subtabs)
    else:
        # actually parse it
        valid_col_names = 'name,magnitude,priority,pm-ra,pm-dec,epoch,type'.split(',')

        catalog_colnames = list(catalog.colnames)
        if 'mmtbino_type' in catalog_colnames:
            catalog_colnames.remove('mmtbino_type')
            catalog_colnames.insert(0, 'mmtbino_type')

        ras = catalog['coord'].ra.deg
        decs = catalog['coord'].dec.deg

        cols_to_write = []
        names_to_write = []
        for write_name in valid_col_names:
            found = False
            for colnm in catalog_colnames:
                if colnm.lower() == write_name or (colnm == 'mmtbino_type' and write_name == 'type'):
                    if found:
                        warn(f'Found column name {colnm}, which is a duplicate for '
                             f'{write_name}.  Only taking first matching column ("{found}")')
                    else:
                        found = colnm
                        if write_name == 'type':
                            col = [_bino_typemap[s] for s in catalog[colnm]]
                        else:
                            col = catalog[colnm]
                        cols_to_write.append(col)
                        names_to_write.append(write_name)

        idx = 1 if names_to_write[0] == 'name' else 0
        names_to_write.insert(idx, 'dec')
        cols_to_write.insert(idx, decs)
        names_to_write.insert(idx, 'ra')
        cols_to_write.insert(idx, ras)

        to_write_catalog = table.Table(cols_to_write, names=names_to_write)

    if fn is not None:
        to_write_catalog.write(fn, format='ascii', delimiter=',', overwrite=overwrite)
    return to_write_catalog

def select_flux_stars(cat, magrng=(17, 17.7), extcorr=False, cattype='sdss'):
    """
    Identifies flux calibration stars based on the SDSS approach

    Primarily selects ~F stars based on color cut.

    Uses the SDSS criteria for specphot standards:
    0.1 < (g-r) < 0.4

    their specphot stars are 16 < g < 17.1 and "reddening standards" are
    17.1 < g < 18.5

    Parameters
    ----------
    cat : astropy Table
        Catalog for extraction of flux stars
    magrng : 2 tuple
        Range of magnitudes to allow
    extcorr : bool
        Whether or not to apply extinction corrections to the observed values
    cattype : str
        Type of catalog.  Must be 'sdss', 'base' or 'decals5'
    """


    if cattype == 'sdss':
        starmsk = cat['type'] == 6 #type==6 means star
        u = cat['u']
        g = cat['g']
        r = cat['r']
        if extcorr:
            u = u - cat['Au']
            g = g - cat['Ag']
            r = r - cat['Ar']
    elif cattype == 'base':
        starmsk = cat['PHOT_SG'] == 'STAR'
        u = cat['u']
        g = cat['g']
        r = cat['r']
        if extcorr:
            u = u - cat['Au']
            g = g - cat['Ag']
            r = r - cat['Ar']
    elif cattype == 'decals5':
        starmsk = cat['type'] == 'PSF '
        u = np.array(22.5 - 2.5*np.log10(cat['flux_u']))
        g = np.array(22.5 - 2.5*np.log10(cat['flux_g']))
        r = np.array(22.5 - 2.5*np.log10(cat['flux_r']))
        if extcorr:
            u += 2.5 * np.log10(cat['mw_transmission_u'])
            g += 2.5 * np.log10(cat['mw_transmission_g'])
            r += 2.5 * np.log10(cat['mw_transmission_r'])

    umg = u - g
    gmr = g - r

    #msk = (0.1 < gmr) & (gmr < 0.4) #old version

    std_color = (0.6 < umg) & (umg < 1.2)
    std_color = std_color & (0.0 < gmr) & (gmr < 0.6)
    std_color = std_color & (gmr > 0.75 * umg - 0.45)

    sp_std = std_color & (15.5 < g) & (g < 17.0)
    red_std = std_color & (17 < g) & (g < 18.5)

    minmag = min(*magrng)
    maxmag = max(*magrng)

    msk = (sp_std|red_std) & (minmag < r)& (r < maxmag) & starmsk

    newcat = cat[msk]
    newcat['name'] = np.char.add('flux', (np.arange(len(newcat))+1).astype(str))
    newcat['magnitude'] = r[msk]
    newcat['type'] = ['standard']*len(newcat)

    return newcat
