from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u

"""
These functions are for the SAGA project Magellan-related spectroscopic mask
design.

Assumes you have the IMACS mask design software
(http://users.obs.carnegiescience.edu/clardy/imacs/maskmaking) installed.

How to design IMACS masks:
1. Generate the target list using `build_imacs_targeting_files`.  This will put
   a catalog and initial .obs file in the ``imacs_targets`` directory.
2. Open 2 terminals, and in both, cd into ``imacs_targets`` and do
   ``source ../imacs/imacs.envars``
   (Note that you will probably need to update that file, as it's currently
    based on paths to Magellan-related code for Erik's computer)
3. Do "intgui -k <field>_ini.obs" in one of the terminals
4. Change the filename and title boxes to <field>_1
5. Use the GUI to pick the field location/rotation to optimize guide stars
   you can use the middle mouse button to re-center the mask.  You can also use
   the "maskgen" button to do a test run of generating the mask (although you'll
   have to do it "for real" in the next step to get the mask files).
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
WLIMIT  4500.0 9000.0{dlimit}
Wavelength  6749.28{pdecide}
TELESCOPE  Magellan
INSTRUMENT IMACS_sc
DISPERSER  IMACS_grism_300
SLITSIZE  1.000 3.000 3.000 0.000
REFHOLE   5.000 1 2.500 2.500 0.000
EXORDER  3
REPOBJ  0
REFLIMIT 12
DATE {obsmjd}
#  Object file list
OBJFILE  {catfile}
"""


def build_imacs_targeting_files(host, observername, date=None,
                                onlygals=True, refmagrange={'r': (17, 19)},
                                overwrite=False, selectkws={}, targs=None,
                                pdecide=None, weakdlimit=False, inclhost=True):
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
    selectkws : dict
        Keywords to be passed into `targeting.select_targets`
    targs : astropy.table.Table
        The target list to use.  If given, will override `selectkws` and
        `targeting.select_targets` will be ignored.  Format should be the same
        as `targeting.select_targets`'s output - needs to have 'objID, 'type',
        'ra', 'dec', and 'r'.  If the column 'imacs_pri' is present, that
        will be included as a priority for the object
    pdecide : float or None
        The priority to use as 'PDECIDE', or None to go with the default
    weakdlimit : bool
        If True, set the 'DLIMIT' such that only the center wl (6750) is forced
        to be on the detector.
    inclhost : bool
        If True, always add a (priority 0) target for the host.
    """
    import os
    import targeting
    from astropy.coordinates import Angle
    from astropy.time import Time
    from astropy import units as u

    if targs is None:
        targs = targeting.select_targets(host, **selectkws)
    elif selectkws:
            raise ValueError('Cannot give both `targs` and `selectkws` arguments')

    cat = host.get_sdss_catalog()  # needed for ref stars

    if onlygals:
        if 'type' in targs.colnames:
            galmsk = targs['type'] == 3
            targs = targs[galmsk]
        else:
            print('onlygals was set, but no "type" in the catalog.  Not cutting.')

    refmsk = np.ones(len(cat), dtype=bool)
    for band in 'gri':
        refmsk & (np.abs(cat[band] - cat['psf_' + band]) < 0.25)
    for band, rng in refmagrange.items():
        refmsk = refmsk & (min(rng) < cat['psf_' + band]) & (cat['psf_' + band] < max(rng))

    fncat = 'imacs_targets/{0}.cat'.format(host.shortname)
    fnobs = 'imacs_targets/{0}_ini.obs'.format(host.shortname)

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

        if inclhost:
            f.write('@{0} {1} {2} -10 Pri=0\n'.format(host.name, host.ra, host.dec))

        if 'imacs_pri' in targs.colnames:
            for t in targs:
                f.write('@{0} {1} {2} {3} Pri={4}\n'.format(*[t[n] for n in
                        'objID,ra,dec,r,imacs_pri'.split(',')]))
        else:
            for t in targs:
                f.write('@{0} {1} {2} {3}\n'.format(*[t[n] for n in
                        'objID,ra,dec,r'.split(',')]))

        #now the reference stars
        for t in cat[refmsk]:
            f.write('*{0} {1} {2} #r={3}\n'.format(*[t[n] for n in 'objID,ra,dec,r'.split(',')]))
    print('Wrote catalog to', fncat)

    hra = Angle(host.ra, u.degree)
    hdec = Angle(host.dec, u.degree)

    obsmjd = int(Time(date, scale='utc').mjd)

    obsfile = _obsfile_template.format(title=host.shortname+'_ini',
        cenra=hra.to_string(u.hour, precision=3, sep=':', pad=True),
        cendec=hdec.to_string(u.degree, precision=2, sep=':', pad=True),
        observer=observername,
        obsmjd=obsmjd,
        catfile=host.shortname + '.cat',
        pdecide='' if pdecide is None else ('\nPDECIDE {0:.2f}'.format(pdecide)),
        dlimit='\nDLIMIT 5000 7500' if weakdlimit else '')

    with open(fnobs, 'w') as f:
        f.write(obsfile)
    print('Wrote obs file to', fnobs)

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
        print('Removing', np.sum(~fimsk),'of', len(fimsk), 'objects due to not being in the requested fields')
        cmsk = cmsk & fimsk

    hras = np.array([c.ra.degree for c in hcoords[cmsk]])
    hdecs = np.array([c.dec.degree for c in hcoords[cmsk]])

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

    print('Removed', nremoved, 'from magellan catalog, leaving', nremaining,'\nwriting new catalog to', fncatnew)
    with open(fncatnew, 'w') as f:
        for l in newcatlines:
            f.write(l)

    return np.array(remras), np.array(remdecs)


def add_use(oldfn, newfn, obwtoapply):
    """
    Adds "Use=N" columns to a catalog as output by `build_imacs_targeting_files`
    based on a supplied set of previous observations.

    Parameters
    ----------
    oldfn : str
        The ".cat" file to start from.
    newfn : str
        The new file to output based on `oldfn` but with "Use=N" at the end of
        some lines.
    obwtoapply : list of str
        A list of ".obw" files to parse for "Use=N" entries.  Must have the same
        object names as in `oldfn`.
    Returns
    -------
    nnotfound : dict
        A dictionary mapping from object name to the number of uses.  Only
        includes those that were in one of `obwtoapply` but not found in
        `oldfn`.
    """
    from shutil import move

    objnameton = {}
    for obwfn in obwtoapply:
        with open(obwfn) as f:
            for l in f:
                if l.startswith('@'):  # indicates it's a target object
                    ls = l.strip().split()
                    if ls[-1].startswith('Use='):
                        usen = int(ls[-1].replace('Use=', ''))
                        #if it's already present and a larger usen, don't add it
                        if not (ls[0] in objnameton and usen > objnameton[ls[0]]):
                            objnameton[ls[0]] = usen

    print('Trying to mark', len(objnameton), 'objects as used in', oldfn, 'to make', newfn)

    if newfn == oldfn:
        realnewfn = newfn
        newfn += '.tmp'
    else:
        realnewfn = None
    with open(newfn, 'w') as fw:
        with open(oldfn) as fr:
            for l in fr:
                if l.startswith('@'):  # indicates it's a target object
                    ls = l.strip().split()
                    if ls[0] in objnameton and ' Use=' not in l:
                        fw.write(l[:-1])  # strips the newline
                        fw.write('  Use=' + str(objnameton.pop(ls[0])))
                        fw.write('\n')
                        continue

                fw.write(l)

    if realnewfn:
        move(newfn, realnewfn)

    if len(objnameton) > 0:
        print(len(objnameton), 'objects in obw files', obwtoapply, 'were not matched in', oldfn)

    return objnameton


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
                radegs.append(Angle(ls[2], hour).degree)
                decdegs.append(Angle(ls[3], degree).degree)
    return names, radegs, decdegs


def plot_imacs_masks(host, clf=True, save=False, eastleft=False, altname=None,
                     skipnums=[], plotpris=False, showrvir=False,
                     showfootprint=False):
    import itertools
    import matplotlib
    from glob import glob

    smfs = glob('imacs_targets/{0}_*.SMF'.format(host.shortname))
    if altname:
        smfs.extend(glob('imacs_targets/{0}_*.SMF'.format(altname)))
    smfs = dict([(int(smf.split('_')[-1].split('.')[0]), smf) for smf in smfs])
    smfs = [smfs[i] for i in sorted(smfs)]
    if len(smfs)>0:
        print(smfs)

    if clf:
        plt.clf()

    #catalog of everything
    nms = []
    ras = []
    decs = []
    pris = []
    with open('imacs_targets/{0}.cat'.format(host.shortname)) as f:
        for l in f:
            if l[0] == '@':
                nm, ra, dec = l.split()[:3]
                nms.append(nm[1:])
                ras.append(float(ra))
                decs.append(float(dec))
                if 'Pri=' in l:
                    pris.append(float(l[l.index('Pri=')+4:].split()[0]))
                else:
                    pris.append(None)
    nms = np.array(nms)
    ras = np.array(ras)
    decs = np.array(decs)
    pris = np.array(pris)

    if np.any(plotpris):
        if plotpris is True:
            plotpris = np.unique(pris)
        ccycle = itertools.cycle(matplotlib.rcParams['axes.color_cycle'])
        for pri in np.sort(plotpris):
            msk = pris == pri
            plt.scatter(ras[msk], decs[msk], s=6, c=ccycle.next(), edgecolor='none', alpha=.7, label='Pri={0}, n={1}'.format(pri, msk.sum()))
    else:
        plt.scatter(ras, decs, label='All', s=4, c='k', edgecolor='none',alpha=.6)

    n = 0
    label = 'Have slits ($n_{{\\rm mask}}={0}$)'.format(len(smfs))
    allnms = []
    for fn in smfs:
        msknum = int(fn.split('_')[-1].split('.')[0])
        if msknum in skipnums:
            continue
        nmi, rai, deci = get_smf_entries(fn)
        allnms.extend(nmi)
        if np.any(plotpris):
            plt.scatter(rai, deci, color='r', s=10, alpha=.8, label=label, edgecolor='k', lw=1)
            label = ''
        else:
            label = str(msknum)
            plt.scatter(rai, deci, color='r', s=3, alpha=.7, label=label, edgecolor='none')

    print('Total targets already observed=', len(allnms))
    if np.any(plotpris):
        obspris = pris[np.in1d(nms, allnms)]
        print('Priority distribution of slitted targets:')
        ipri = 1
        for pri in np.unique(pris):
            ntarg = np.sum(obspris==pri)
            nrem = np.sum(pris==pri) - np.sum(obspris==pri)
            print('Priority', pri, ':', ntarg, 'targeted,', nrem, 'remaining:', 100*ntarg/(ntarg+nrem),'%')
            # if pri<=2:
            #     msg = 'Targeted {0} of {1} Pri {3} objects: {2:.0f}%'
            #     txt = msg.format(ntarg, ntarg+nrem, 100*ntarg/(ntarg+nrem), pri)
            #     plt.text(0.05, 0.025*ipri, txt, transform=plt.gca().transAxes)
            #     ipri+=1
        ntarg = len(obspris)
        nrem = len(pris) - len(obspris)
        print('All priorities:', ntarg, 'targeted,', nrem, 'remaining:', 100*ntarg/(ntarg+nrem),'%')

    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.legend(loc='upper left')

    if showrvir:
        circle = plt.Circle((host.ra, host.dec), host.physical_to_projected(300*u.kpc).to(u.deg).value, color='k', lw=2)
        plt.gca().add_artist(circle)

    if showfootprint:
        fov = (27.2*u.arcmin).to(u.deg).value
        footprint = plt.Rectangle((host.ra-fov/2, host.dec-fov/2), fov, fov, color='r', lw=2, fc='none')
        plt.gca().add_artist(footprint)

    if eastleft:
        plt.xlim(max(plt.xlim()), min(plt.xlim()))
    else:
        plt.xlim(min(plt.xlim()), max(plt.xlim()))

    if save:
        plt.savefig('imacs_targets/targets_{0}.pdf'.format(host.shortname))


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
