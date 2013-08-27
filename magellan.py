from __future__ import division

import numpy as np
from matplotlib import pyplot as plt

"""
These functions are for the "Distant local groups" project Magellan-related work.

How to design IMACS masks:
1. Generate the target list using `build_imacs_targetlists`.  This will put a
   catalog and initial .obs file in the ``imacs_targets`` directory.
2. Open 2 terminals, and in both, cd into ``imacs_targets`` and do ``source ../imacs.envars``
3. Do "intgui -k <field>_ini.obs" in one of the terminals
4. Change the filename and title boxes to <field>_1
5. Use the GUI to pick the field location/rotation to optimize guide stars
   (the top one is the most important, I think), and perhaps try the maskgen
   button.
6. In the other terminal, do "maskgen <field>_1"
7. In the intgui, update the file and title to "<field>_2", and change the
   catalog to "<field>_1.obw".
8. Repeat 4-7 until done.
9. Have a look at what you selected with `plot_imacs_masks` and
   `imagelist_imacs_targets`.

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


#CENTER   24:00:00.000 -00:02:00.00
_obsfile_template = """OBSERVER  {observer}
TITLE   {title}
CENTER   {cenra} {cendec}
EQUINOX  2000.00000
POSITION 0.00000
WLIMIT  3500.0 9700.0
Wavelength  6749.28
TELESCOPE  Magellan
INSTRUMENT IMACS_sc
DISPERSER  IMACS_grism_300
SLITSIZE  1.000 6.000 6.000 0.000
REFHOLE   5.000 1 2.500 2.500 0.000
EXORDER  0
REPOBJ  0
REFLIMIT 12
DATE {obsmjd}
#  Object file list
OBJFILE  {catfile}
"""


def build_imacs_targeting_files(host, observername, date='2013-02-15',
                                onlygals=True, refmagrange={'r': (17, 19)},
                                overwrite=False):
    """
    Generates the target catalog and initial observation file for IMACS

    Parameters
    ----------
    host : NSAHost
        The host to target
    observername : str
        The name of the observer (for the .obs file)
    date : str
        The data of the observation - should be parsable by `astropy.time.Time`
    onlygals : bool
        If True, only generate targets that SDSS calls galaxies
    refmagrange: dict of 2-tuples.
        Magnitude range for reference stars - dict keys are band.
        http://www.lco.cl/telescopes-information/magellan/instruments/imacs/user-manual/the-imacs-user-manual#prepobs
        suggests 17 < R< 19.
    overwrite : bool
        If True, existing files will be overwritten.  Otherwise, an error will
        be raised if the files are present.

    """
    import os
    import targeting
    from astropy.coordinates import Angle
    from astropy.time import Time
    from astropy import units as u

    cat = host.get_sdss_catalog()  # needed for ref stars

    targs = targeting.select_targets(host)

    if onlygals:
        galmsk = targs['type'] == 3
        targs = targs[galmsk]

    refmsk = np.ones(len(cat), dtype=bool)
    for band in 'gri':
        refmsk & (np.abs(cat[band] - cat['psf_' + band]) < 0.25)
    for band, rng in refmagrange.iteritems():
        refmsk = refmsk & (min(rng) < cat['psf_' + band]) & (cat['psf_' + band] < max(rng))

    fncat = 'imacs_targets/{0}.cat'.format(host.name)
    fnobs = 'imacs_targets/{0}_ini.obs'.format(host.name)

    if not overwrite:
        if os.path.exists(fncat):
            if os.path.exists(fnobs):
                raise IOError('Object catalog and obs file ("{0}" and "{1}") already exist!'.format(fncat, fnobs))
            else:
                raise IOError('Object catalog ({0}) already exists!'.format(fncat))
        elif os.path.exists(fnobs):
            raise IOError('Observation file ({0}) already exists!'.format(fnobs))

    with open(fncat, 'w') as f:
        f.write('&RADEGREE\n')
        f.write('@{0} {1} {2} -10\n'.format(host.name, host.ra, host.dec))
        for t in targs:
            f.write('@{0} {1} {2} {3}\n'.format(*[t[n] for n in 'objID,ra,dec,r'.split(',')]))

        #now the reference stars
        for t in cat[refmsk]:
            f.write('*{0} {1} {2} #r={3}\n'.format(*[t[n] for n in 'objID,ra,dec,r'.split(',')]))
    print 'Wrote catalog to', fncat

    hra = Angle(host.ra, u.degree)
    hdec = Angle(host.dec, u.degree)

    obsmjd = int(Time(date, scale='utc').mjd)

    obsfile = _obsfile_template.format(title=host.name+'_ini',
        cenra=hra.format(u.hour, precision=3, sep=':', pad=True),
        cendec=hdec.format(u.degree, precision=2, sep=':', pad=True),
        observer=observername,
        obsmjd=obsmjd,
        catfile=host.name + '.cat')

    with open(fnobs, 'w') as f:
        f.write(obsfile)
    print 'Wrote obs file to', fnobs

def reprocess_catalog_for_prev_mmt_obs(fncat, hectocfg, fncatnew, rankcutoff=2, tolarcsec=1, magrng=None, hectofields='all'):
    """
    This removes everything from the mmt catalog in a given set of fields and magrngs
    """
    from scipy.spatial import cKDTree
    from mmthecto import parse_cfg_file

    hcoords, htargets, hranks, hfields = parse_cfg_file(hectocfg)

    cmsk = (hranks > rankcutoff)
    if hectofields != 'all':
        fimsk = np.zeros_like(cmsk)
        for fi in hectofields:
            fimsk = fimsk | (hfields == fi)
        print 'Removing', np.sum(~fimsk),'of', len(fimsk), 'objects due to not being in the requested fields'
        cmsk = cmsk & fimsk

    hras = np.array([c.ra.degrees for c in hcoords[cmsk]])
    hdecs = np.array([c.dec.degrees for c in hcoords[cmsk]])

    toldeg = tolarcsec / 3600.

    kdt = cKDTree(np.array([hras, hdecs], copy=False).T)

    newcatlines = []

    if magrng is None:
        check_mag = lambda mag: True
    else:
        magrng = min(magrng), max(magrng)
        check_mag = lambda mag: magrng[0] < mag < magrng[1]


    remras = []
    remdecs = []
    nremoved = 0
    nremaining = 0
    with open(fncat) as fr:
        for l in fr:
            if l.startswith('@'):
                nm, ra, dec, mag = l[1:-1].split()
                d, i = kdt.query([float(ra), float(dec)])
                if (d < toldeg) and check_mag(float(mag)):
                    remras.append(float(ra))
                    remdecs.append(float(dec))
                    nremoved += 1
                    continue  # means don't add to newcatlines
                else:
                    nremaining += 1

            newcatlines.append(l)

    print 'Removed', nremoved, 'from magellan catalog, leaving', nremaining,'\nwriting new catalog to', fncatnew
    with open(fncatnew, 'w') as f:
        for l in newcatlines:
            f.write(l)

    return np.array(remras), np.array(remdecs)


def get_smf_entries(fn, inclholes=False):
    from astropy.coordinates import Angle
    from astropy.units import hour, degree

    names = []
    radegs = []
    decdegs = []
    with open(fn) as f:
        for l in f:
            if l[:4] == 'SLIT' or (inclholes and (l[:4] == 'HOLE')):
                ls = l.split()
                names.append(ls[1])
                radegs.append(Angle(ls[2], hour).degrees)
                decdegs.append(Angle(ls[3], degree).degrees)
    return names, radegs, decdegs


def plot_imacs_masks(host, clf=True, save=False):
    from glob import glob

    smfs = glob('imacs_targets/{0}_?.SMF'.format(host.name))

    if clf:
        plt.clf()

    #catalog of everything
    ras = []
    decs = []
    with open('imacs_targets/{0}.cat'.format(host.name)) as f:
        for l in f:
            if l[0] == '@':
                ra, dec = l.split()[1:3]
                ras.append(float(ra))
                decs.append(float(dec))
    plt.plot(ras, decs, '.k', label='All', ms=1, alpha=.6)

    for fn in smfs:
        msknum = int(fn.split('_')[-1].split('.')[0])
        nmi, rai, deci = get_smf_entries(fn)

        plt.plot(rai, deci, '.', ms=5, alpha=.8, label=str(msknum))

    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.legend(loc=0)

    if save:
        plt.savefig('imacs_targets/targets_{0}.pdf'.format(host.name))


def imagelist_imacs_targets(smffn, n=100):
    import targeting

    nms, ras, decs = get_smf_entries(smffn)
    return targeting.sampled_imagelist(ras, decs, n, nms)


def load_ricardo_rvs_and_match(fn='rv_ricardo_jul32013.dat', fields=None, vtouse='VREL'):
    from os import path
    from astropy.io import ascii

    drv = ascii.read(fn, delimiter='\t', guess=False, Reader=ascii.CommentedHeader)
    ufis = np.unique(['_'.join(fn.split('_')[:-2]) for fn in drv['FILE']])

    fis = {}
    for fi in ufis:
        fis[fi] = get_smf_entries(path.join('imacs_targets', fi + '.SMF'), inclholes=True)

    rvs = []
    rverrs = []
    ras = []
    decs = []
    nms = []
    fns = []
    specclass = []

    for elem in drv:
        fi = '_'.join(elem['FILE'].split('_')[:-2])
        if fields is not None and fi not in fields:
            continue

        idx = int(elem['FILE'].split('_')[-2]) - 1

        rvs.append(elem[vtouse])
        rverrs.append(elem['VERR'])
        fns.append(elem['FILE'])
        specclass.append(elem['CLASS'])

        nms.append(fis[fi][0][idx])
        ras.append(fis[fi][1][idx])
        decs.append(fis[fi][2][idx])

    return dict([(enm, np.array(locals()[enm])) for enm in 'nms,ras,decs,rvs,rverrs,fns,specclass'.split(',')])





