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



