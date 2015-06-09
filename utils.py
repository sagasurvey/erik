"""
Various utility functions that have no other obvious place to live.
"""
from __future__ import division, print_function

import numpy as np


def int_to_roman(input):
    """
    Convert an integer to Roman numerals.

    From http://code.activestate.com/recipes/81611-roman-numerals/

    Examples:
    >>> int_to_roman(0)
    Traceback (most recent call last):
    ValueError: Argument must be between 1 and 3999

    >>> int_to_roman(-1)
    Traceback (most recent call last):
    ValueError: Argument must be between 1 and 3999

    >>> int_to_roman(1.5)
    Traceback (most recent call last):
    TypeError: expected integer, got <type 'float'>

    >>> for i in range(1, 21): print int_to_roman(i)
    ...
    I
    II
    III
    IV
    V
    VI
    VII
    VIII
    IX
    X
    XI
    XII
    XIII
    XIV
    XV
    XVI
    XVII
    XVIII
    XIX
    XX
    >>> print int_to_roman(2000)
    MM
    >>> print int_to_roman(1999)
    MCMXCIX
    """
    if type(input) != type(1):
        raise TypeError("expected integer, got %s" % type(input))
    if not 0 < input < 4000:
        raise ValueError("Argument must be between 1 and 3999")
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
    result = ""
    for i in range(len(ints)):
        count = int(input / ints[i])
        result += nums[i] * count
        input -= ints[i] * count
    return result


def roman_to_int(input):
    """
    From http://code.activestate.com/recipes/81611-roman-numerals/

    Convert a roman numeral to an integer.

    >>> r = range(1, 4000)
    >>> nums = [int_to_roman(i) for i in r]
    >>> ints = [roman_to_int(n) for n in nums]
    >>> print r == ints
    1

    >>> roman_to_int('VVVIV')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: VVVIV
    >>> roman_to_int(1)
    Traceback (most recent call last):
     ...
    TypeError: expected string, got <type 'int'>
    >>> roman_to_int('a')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: A
    >>> roman_to_int('IL')
    Traceback (most recent call last):
     ...
    ValueError: input is not a valid roman numeral: IL
    """
    if type(input) != type(""):
        raise TypeError("expected string, got %s" % type(input))
    input = input.upper()
    nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
    ints = [1000, 500, 100, 50,  10,  5,   1]
    places = []
    for c in input:
        if not c in nums:
            raise ValueError("input is not a valid roman numeral: %s" % input)
    for i in range(len(input)):
        c = input[i]
        value = ints[nums.index(c)]
        # If the next place holds a larger number, this value is negative.
        try:
            nextvalue = ints[nums.index(input[i + 1])]
            if nextvalue > value:
                value *= -1
        except IndexError:
            # there is no next place.
            pass
        places.append(value)
    sum = 0
    for n in places:
        sum += n
    # Easiest test for validity...
    if int_to_roman(sum) == input:
        return sum
    else:
        raise ValueError('input is not a valid roman numeral: %s' % input)


def get_mcconn_table(fnorurl='http://www.astro.uvic.ca/~alan/Nearby_Dwarfs_Database_files/NearbyGalaxies.dat', dropmw=True, sanitizenames=False):
    """
    Gets the Mcconnachie 12 table from the file.

    Parameters
    ----------
    fnorurl : str, optional
        Local file or URL.  If URL, will be cached.
    dropmw : bool, optional
        If true, the entry for the MW will be skipped.
    sanitizenames : bool, optional
        If True, names like "gal (I)" will be changed to "gal I"

    Returns
    -------
    mctab : astropy.table.Table
        The Mcconnachie 12 "database" as an astropy table
    """
    import os

    from astropy.coordinates import Angle, ICRS, SkyCoord
    from astropy.table import Table, Column
    from astropy.utils import data
    from astropy import units as u

    if fnorurl is None:
        if os.path.isfile('NearbyGalaxies.dat'):
            fnorurl = 'NearbyGalaxies.dat'
        else:
            fnorurl = 'https://www.astrosci.ca/users/alan/Nearby_Dwarfs_Database_files/NearbyGalaxies.dat'

    datstr = data.get_file_contents(fnorurl, cache=True)
    datlines = datstr.split('\n')

    # pull from the header.
    for hdri, hdrstr in enumerate(datlines):
        if hdrstr.startswith('GalaxyName'):
            break
    else:
        raise ValueError('No field name line found')
    # hdrstr now is the line with the field names

    # type is 'pm' for float w/ err+/- or 'coord'
    fieldinfo = [('name', 'S19', None),
                 ('center', 'coord', None),
                 ('EBmV', float, u.mag),
                 ('distmod', 'pm', u.mag),
                 ('vh', 'pm', u.km / u.s),
                 ('Vmag', 'pm', u.mag),
                 ('PA', 'pm', u.degree),
                 ('e', 'pm', u.dimensionless_unscaled),
                 ('muV0', 'pm', u.mag * u.arcsec ** -2),
                 ('rh', 'pm', u.arcmin),
                 ('sigma_s', 'pm', u.km / u.s),
                 ('vrot_s', 'pm', u.km / u.s),
                 ('MHI', float, u.Unit('1e6 solMass')),
                 ('sigma_g', 'pm', u.km / u.s),
                 ('vrot_g', 'pm', u.km / u.s),
                 ('[Fe/H]', 'pm', u.dimensionless_unscaled),
                 ('F', int, None),
                 ('References', 'S40', None)]

    fieldnames = []
    fielddtypes = []
    fieldunits = []
    for nm, tp, un in fieldinfo:
        if tp == 'coord':
            fieldnames.append(nm)
            fielddtypes.append(object)
            fieldunits.append(None)
        elif tp == 'pm':
            fieldnames.append(nm)
            fielddtypes.append(float)
            fieldunits.append(un)
            fieldnames.append(nm + '+')
            fielddtypes.append(float)
            fieldunits.append(un)
            fieldnames.append(nm + '-')
            fielddtypes.append(float)
            fieldunits.append(un)
        else:
            fieldnames.append(nm)
            fielddtypes.append(tp)
            fieldunits.append(un)

    t = Table(names=fieldnames, dtype=fielddtypes)
    for nm, un in zip(fieldnames, fieldunits):
        t[nm].units = un

    for l in datlines[(hdri + 2 + int(dropmw)):]:
        if l.strip() == '':
            continue

        vals = l[19:].split()
        ra = Angle(tuple([float(v) for v in vals[0:3]]), unit=u.hour)
        dec = Angle(tuple([float(v) for v in vals[3:6]]), unit=u.degree)
        vals = vals[6:]
        vals.insert(0, ICRS(ra, dec))
        vals.insert(0, l[:19].strip())
        if '(' not in vals[-1]:
            vals.append('')
        t.add_row(vals)

    # now add derived columns
    t.add_column(Column(name='Vabs', data=t['Vmag'] - t['distmod'], unit=u.mag))  # Vmag apparently already includes dereddening?
    t.add_column(Column(name='logLV', data=(t['Vabs'] - 4.83) / -2.5, unit=u.solLum))
    t.add_column(Column(name='distance', data=10 ** (t['distmod'] / 5. - 2), unit=u.kpc))
    t.add_column(Column(name='distance+', data=t['distmod+'] * t['distance'] * np.log(10) / 5, unit=u.kpc))
    t.add_column(Column(name='distance-', data=t['distmod-'] * t['distance'] * np.log(10) / 5, unit=u.kpc))
    t.add_column(Column(name='rh_phys', data=t['distance'] * np.radians(t['rh'] / 60.), unit=u.kpc))
    t.add_column(Column(name='radeg', data=[((ti.ra.degree + 180) % 360 - 180) for ti in t['center']], unit=u.degree))
    t.add_column(Column(name='decdeg', data=[ti.dec.degree for ti in t['center']], unit=u.degree))

    # Now replace the 'center' with a version that includes the distance.
    newscs = []
    xs = []
    ys = []
    zs = []
    for ic in t['center']:
        newscs.append(SkyCoord(ic.ra, ic.dec, distance=t['distance']))
        xs.append(newscs[-1].cartesian.x)
        ys.append(newscs[-1].cartesian.y)
        zs.append(newscs[-1].cartesian.z)
    t.add_column(Column(name='x', data=xs, unit=u.kpc))
    t.add_column(Column(name='y', data=ys, unit=u.kpc))
    t.add_column(Column(name='z', data=zs, unit=u.kpc))

    # andromeda dSph numbers
    andnum = []
    for nm in t['name']:
        if nm.startswith('Andromeda '):
            andnum.append(roman_to_int(nm.replace('Andromeda ', '')))
        elif nm == 'Andromeda':
            andnum.append(0)
        elif nm in ('M32', 'NGC 205', 'NGC 185', 'NGC 147'):
            andnum.append(-1)
        elif nm in ('IC 10', 'LGS 3'):
            andnum.append(-2)
        elif nm in ('IC 1613', 'Pegasus dIrr'):
            andnum.append(-3)
        else:
            andnum.append(-99)
    t.add_column(Column(name='and_number', data=andnum, unit=None))

    if sanitizenames:
        t['name'] = [nm.replace('(I)', 'I') for nm in t['name']]

    return t
