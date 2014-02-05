# coding: utf-8
from __future__ import print_function, division

"""
A python function to grab the RC3 catalog from vizier and parse it.

This does not include the actual rc3 data file.  The python `load_rc3` function
will get the data file automatically from
ftp://cdsarc.u-strasbg.fr/pub/cats/VII/155/rc3.gz ... but if that server is down
or something and you can find an RC3 data file from somewhere else, that should
work too - just be sure to rename it "rc3.dat" or pass in the filename to the
`load_rc3` function
"""


import os
import numpy as np
from astropy import units as u

_DEFAULT_RC3_URL = "ftp://cdsarc.u-strasbg.fr/pub/cats/VII/155/rc3.gz"

RC3_COLUMN_DEFINITIONS = """
RAh 1 2 i2
RAm 3 4 i2
RAs 5 8 f4
Decsign 10 10 a1
Decd 11 12 i2
Decm 13 14 i2
Decs 15 16 i2
RA1950h 18 19 i2
RA1950m 20 21 i2
RA1950s 22 25 f4
Dec1950d 27 29 i2
Dec1950m 30 31 i2
Dec1950s 32 33 i2
Gallong 35 40 f
Gallat 42 47 f
Sgallong 49 54 f
Sgallat 56 61 f
name 63 74 a12
altname 75 89 a15
desig 91 104 a14
PGCnum 106 116 a11
type 118 124 a7
typesrc 126 130 a5
Ttype 132 135 f
Ttypeerr 137 139 f
lumcl 141 144 f
lumclerr 146 148 f
lumclnest 150 150 i
logd25 152 155 f
logd25unc 156 156
logd25err 158 160 f
logr25 162 165 f
logr25unc 166 166
logr25err 168 170 f
logDo 172 175 f
logAe 177 180 f
logAeerr 182 184 f
PA 186 188  i
BT 190 194 f
BTcode 195 195
BTerr 197 199 f
Bmag 201 205 f
Bmagerr 207 209 f
BoT 211 215 g
mp25 217 221 f
mp25err 223 226 f
mpe 228 232 f
mpeerr 234 236 f
mFIR 238 242 f
m21 244 248 f
m21err 250 251 f
B-VT 253 256 f
B-VTerr  258 260 f
B-Ve 262 265 f
B-Veerr 267 269 f
B-VoT 271 274 f
U-BT 276 280 f
U-BTerr 282 284 f
U-Be 286 289 f
U-Beerr 291 293 f
U-BoT 295 298 f
HI 300 304 f
Ai 306 309 f
A21 311 313 f
Ag 315 318 f
W20 320 322 i
W20err 324 325 i
W50 327 329 i
W50err 331 332 i
V21 334 338 i
V21err 340 341 i
cz 343 347 i
czerr 349 351 i
Vgsr 353 357 i
V3K 359 363 i
"""[1:-1]


def load_rc3(datfn='rc3.dat'):
    """
    Loads the RC3 catalog

    Parameters
    ----------
    datfn : str or None
        If a string, will look for that file and try to load it.  If the file is
        not present (or `datfn` is None), this will look for the RC3 data file
        online on vizier and try to download it. Or, if the string starts with
        "http://" or "ftp://", it will use that as the download url.

    Returns
    -------
    tab : astropy.table.table
        The RC3 table
    coo : astropy.coordinates.ICRS
        An object with coordinates aligned to the table
    """
    from astropy.coordinates import ICRS, Latitude, Longitude
    from astropy.io import ascii
    from astropy.utils.data import download_file

    if not os.path.exists(datfn):
        if datfn.startswith('http://') or datfn.startswith('ftp://'):
            dlurl = datfn
        else:
            dlurl = _DEFAULT_RC3_URL

        datfn = download_file(dlurl, cache=True)

    names = []
    colstarts = []
    colends = []
    for l in RC3_COLUMN_DEFINITIONS.split('\n'):
        ls = l.strip().split()
        names.append(ls[0])
        colstarts.append(int(ls[1]) - 1)
        colends.append(int(ls[2]) - 1)
    tab = ascii.read(datfn, format='fixed_width_no_header',
                     names=names, col_starts=colstarts, col_ends=colends)
    nosec = tab['Decs'].mask.copy()
    tab['Decs'].mask[nosec] = False
    tab['Decs'][nosec] = 0.0
    dsign = np.ones(len(tab))
    dsign[tab['Decsign'] == '-'] = -1.0

    ra = Longitude(tab['RAh'] + tab['RAm']/60. + tab['RAs']/3600., u.hourangle)
    dec = Latitude(dsign*(tab['Decd'] + tab['Decm']/60. + tab['Decs']/3600.), u.degree)
    coo = ICRS(ra, dec)

    return tab, coo

if __name__ == '__main__':
    rc3 = load_rc3()
    print('Loaded', len(rc3), 'columns of the RC3')
