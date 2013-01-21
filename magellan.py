from __future__ import division

import numpy as np
from matplotlib import pyplot as plt

"""
These functions are for the "Distant local groups" project Magellan-related work.

How to design IMACS masks:
"""


def plot_targets_and_imacs_fov(host, camera='short', offset=(0, 0), clf=True, **kwargs):
    """
    Plots the targets for the requested host and the IMACS FOV on top

    Parameters
    ----------
    host : NSAHost
        The host to
    camera : str
        The camera to make the FOV box for: 'short' (f/2) or 'long' (f/4)
    offset : 2-tuple
        RA/Dec offset from the host in arcmin

    kwargs are passed into `select_targets`
    """
    import targeting

    if camera == 'short':
        fov = 27.20
        fovrad = 36.4
    elif camera == 'long':
        fov = 15.46
        fovrad = None
    else:
        raise ValueError('unrecognized camera "{0}"'.format(camera))

    targets = targeting.select_targets(host, **kwargs)
    dra = (targets['ra'] - host.ra) * 60 * np.cos(np.radians(targets['dec']))
    ddec = (targets['dec'] - host.dec) * 60

    if clf:
        plt.clf()
    plt.plot(dra, ddec, '.', ms=1)

    xcorners = [offset[0] - fov / 2, offset[0] + fov / 2]
    xcorners = [xcorners[0], xcorners[1], xcorners[1], xcorners[0], xcorners[0]]
    ycorners = [offset[1] - fov / 2, offset[1] + fov / 2]
    ycorners = [ycorners[1], ycorners[1], ycorners[0], ycorners[0], ycorners[1]]

    plt.plot(xcorners, ycorners, 'k-')


def build_imacs_targetlists(host, onlygals=True):
    import os
    import targeting

    cat = host.get_sdss_catalog() # needed for ref stars

    targs = targeting.select_targets(host)

    if onlygals:
        galmsk = targs['type'] == 3
        targs = targs[galmsk]

    refmsk = (17 < cat['psf_r']) & (cat['psf_r'] < 19)
    for band in 'gri':
        refmsk & (np.abs(cat[band] - cat['psf_' + band]) < 0.25)

    fnobj = 'imacs_targets/{0}.dat'.format(host.name)

    if os.path.exists(fnobj):
        raise ValueError('Object catalog ({0}) already exists!'.format(fnobj))

    with open(fnobj, 'w') as f:
        f.write('&RADEGREE\n')
        f.write('@{0} {1} {2} -10\n'.format(host.name, host.ra, host.dec))
        for t in targs:
            f.write('@{0} {1} {2} {3}\n'.format(*[t[n] for n in 'objID,ra,dec,r'.split(',')]))

        #now the reference stars
        for t in cat[refmsk]:
            f.write('*{0} {1} {2} #r={3}\n'.format(*[t[n] for n in 'objID,ra,dec,r'.split(',')]))
