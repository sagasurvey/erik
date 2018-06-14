"""
This script uses 2dfdr to reduce 2dF/AAOMega data as part of the SAGA Survey. It
assumes:

* That the raw data from any given night is all in a directory that also contains
  "ccd_1" and "ccd_2" directories.
* That any given set of exposures contains a red-appropriate flat, a
  blue-appropriate flat, an arc, and some number of science exposures. (The
  default is to assume that order, although that's customizable)
* That the 2dfdr binaries are somewhere on your current path.

Usage is e.g.:
python ../reduce_saga.py pgc64427_2 25

Run from the "180614" which contains the raw data for that night
(where the pgc64427_2 data begins with flats on exposure 25)

"""

import os
import re
import time
import argparse
import subprocess

CCDNUM_TO_NAME = {1: 'blue', 2: 'red'}
POLL_PERIOD = 1

def check_make_biases(basepath):
    bias1 = os.path.join(basepath, 'bias/ccd1')
    bias2 = os.path.join(basepath, 'bias/ccd2')

    if not (os.path.isdir(bias1) and os.path.isdir(bias2)):
        os.makedirs(bias1)
        os.makedirs(bias2)

    biases_present = True
    biasfn1 = os.path.join(bias1, 'BIAScombined.fits')
    biasfn2 = os.path.join(bias2, 'BIAScombined.fits')
    if not os.path.exists(biasfn1):
        print('ccd1 bias missing.')
        biases_present = False
    if not os.path.exists(biasfn2):
        print('ccd2 bias missing.')
        biases_present = False

    if not biases_present:
        raise NotImplementedError('cannot make biases yet... the directories '
            'are made ({} & {}) but you will have to link the files and combine'
            ' them by hand'.format(bias1, bias2))

    return biasfn1, biasfn2

def check_or_make_symlink(src, dest):
    if os.path.exists(dest):
        if os.path.islink(dest):
            if os.readlink(dest) != src:
                print(dest, 'is a link, but it points to an unexpected place:', src, 'might be an issue?')
            else:
                print(dest, 'is the correct link, doing nothing')
        else:
            raise OSError(dest, 'already exists and not a symlink.  Aborting.')
    else:
        print('symlinking',dest ,'to', src)
        os.symlink(src, dest)

def determine_rawfns(rawdir, ccdnum):
    rex = re.compile(r'(.*?'+str(ccdnum)+r')\d{4}\.fits')

    raw_bases =[]
    for fn in os.listdir(rawdir):
        match = rex.match(fn)
        if match:
            raw_bases.append(match.group(1))
    if raw_bases:
        for b in raw_bases[1:]:
            if b != raw_bases[0]:
                raise ValueError('at least two different base file patterns '
                                 'found in one directory:{} and {}.  Not sure '
                                 'what to do, so aborting.'.format(b, raw_bases[0]))
        return raw_bases[0]
    else:
        return None


def reduce_field(fieldname, idx_template, basepath, redflat_expnum,
                 blueflat_expnum, arc_expnum, sci_expnums):
    raw1 = os.path.join(basepath, 'ccd_1')
    raw2 = os.path.join(basepath, 'ccd_2')
    if not (os.path.exists(raw1) and os.path.exists(raw2)):
        raise ValueError('raw data directories are missing')
    raw_base1 = determine_rawfns(raw1, 1)
    raw_base2 = determine_rawfns(raw2, 2)

    biasfn1, biasfn2 = check_make_biases(basepath)

    fieldpath = os.path.join(basepath, fieldname)
    ccd1 = os.path.join(fieldpath, 'ccd1')
    ccd2 = os.path.join(fieldpath, 'ccd2')
    if not os.path.isdir(ccd1):
        os.makedirs(ccd1)
    if not os.path.isdir(ccd2):
        os.makedirs(ccd2)

    blueflat_fn = '{:04}.fits'.format(blueflat_expnum)
    check_or_make_symlink(os.path.relpath(os.path.join(raw1, raw_base1 + blueflat_fn), ccd1),
               os.path.join(ccd1, raw_base1 + blueflat_fn))

    redflat_fn = '{:04}.fits'.format(redflat_expnum)
    check_or_make_symlink(os.path.relpath(os.path.join(raw2, raw_base2 + redflat_fn), ccd2),
               os.path.join(ccd2, raw_base2 + redflat_fn))

    arc_fn = '{:04}.fits'.format(arc_expnum)
    check_or_make_symlink(os.path.relpath(os.path.join(raw1, raw_base1 + arc_fn), ccd1),
               os.path.join(ccd1, raw_base1 + arc_fn))
    check_or_make_symlink(os.path.relpath(os.path.join(raw2, raw_base2 + arc_fn), ccd2),
               os.path.join(ccd2, raw_base2 + arc_fn))

    for enum in sci_expnums:
        sci_fn = '{:04}.fits'.format(enum)
        check_or_make_symlink(os.path.relpath(os.path.join(raw1, raw_base1 + sci_fn), ccd1),
                   os.path.join(ccd1, raw_base1 + sci_fn))
        check_or_make_symlink(os.path.relpath(os.path.join(raw2, raw_base2 + sci_fn), ccd2),
                   os.path.join(ccd2, raw_base2 + sci_fn))

    check_or_make_symlink(os.path.relpath(biasfn1, ccd1), os.path.join(ccd1, 'BIAScombined.fits'))
    check_or_make_symlink(os.path.relpath(biasfn2, ccd2), os.path.join(ccd2, 'BIAScombined.fits'))

    cmdline_templ = 'aaorun reduce_run {0} -idxfile {1}'
    cmdline1 = cmdline_templ.format(raw_base1, idx_template.format(CCDNUM_TO_NAME[1]))
    print('ccd1 reducing as: "{}" in {}'.format(cmdline1, ccd1))
    cmdline2 = cmdline_templ.format(raw_base2, idx_template.format(CCDNUM_TO_NAME[2]))
    print('ccd2 reducing as: "{}" in {}'.format(cmdline2, ccd2))

    log1fn = os.path.join(fieldpath, 'ccd1.log')
    log2fn = os.path.join(fieldpath, 'ccd2.log')
    print('Spawning two aaoruns. Logs will be in {} and {}.'.format(log1fn, log2fn))

    p1_passed = p2_passed = None
    bluecombinedfn = os.path.join(ccd1, raw_base1 + '_combined.fits')
    if os.path.exists(bluecombinedfn):
        print('ccd1 results already done:', bluecombinedfn)
        p1_passed = True
    redcombinedfn = os.path.join(ccd2, raw_base2 + '_combined.fits')
    if os.path.exists(redcombinedfn):
        print('ccd2 results already done:', redcombinedfn)
        p2_passed = True

    with open(log1fn, 'r' if p1_passed else 'w') as log1f:
        with open(log2fn, 'r' if p2_passed else 'w') as log2f:
            if not p1_passed:
                p1 = subprocess.Popen(cmdline1, shell=True, cwd=ccd1,
                                      stdout=log1f, stderr=subprocess.STDOUT)
            else:
                p1 = None
            if not p2_passed:
                p2 = subprocess.Popen(cmdline2, shell=True, cwd=ccd2,
                                      stdout=log2f, stderr=subprocess.STDOUT)
            else:
                p2 = None

            while not p1_passed and not p2_passed:
                time.sleep(POLL_PERIOD)
                if not p1_passed and p1.poll() is not None:
                    if p1.returncode == 0:
                        print('ccd1 completed sucessfully!')
                        p1_passed = True
                    else:
                        print('ccd1 failed')
                        p1_passed = False
                if not p2_passed and p2.poll() is not None:
                    if p2.returncode == 0:
                        print('ccd2 completed sucessfully!')
                        p2_passed = True
                    else:
                        print('ccd2 failed')
                        p2_passed = False

    if p1_passed and p2_passed:
        splicelogfn = os.path.join(fieldpath, 'splice.log')
        splice_outfn = os.path.join(basepath, fieldname+'_spliced.fits')

        print('Splicing together the results.  Log file:', splicelogfn)

        splice_cmdline_templ = 'aaorun splice "{blue} {red}" -idxfile {idx} -output_file {out}'
        splice_cmdline = splice_cmdline_templ.format(blue=bluecombinedfn,
                                                     red=redcombinedfn,
                                                     idx=idx_template.format(CCDNUM_TO_NAME[1]),
                                                     out=splice_outfn)
        with open(splicelogfn, 'w') as splicef:
            ps = subprocess.Popen(cmdline1, shell=True, cwd=fieldpath,
                                  stdout=splicef, stderr=subprocess.STDOUT)
            retcode = ps.wait()
            if retcode == 0:
                print('Splice succeeded: {}!'.format(splice_outfn))
            else:
                print('Splice failed')




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('fieldname')
    parser.add_argument('firstexposurenum')
    parser.add_argument('-n', '--n-science', default=4)
    parser.add_argument('-r', '--redflat-index', default=0)
    parser.add_argument('-b', '--blueflat-index', default=1)
    parser.add_argument('-a', '--arc-index', default=2)
    parser.add_argument('-o', '--sci-offset', default=3)
    parser.add_argument('--basepath', default = '.')
    parser.add_argument('--idx-template', default=os.path.join(os.environ['HOME'], 'observing/aat/ozdes_{}_mod.idx'))

    args = parser.parse_args()

    enum1 = int(args.firstexposurenum)
    redexp = enum1 + int(args.redflat_index)
    blueexp = enum1 + int(args.blueflat_index)
    arcexp = enum1 + int(args.arc_index)
    sciexp_start = enum1 + int(args.sci_offset)
    sciexps = [i+sciexp_start for i in range(args.n_science)]

    reduce_field(args.fieldname, args.idx_template, args.basepath, redexp,
                 blueexp, arcexp, sciexps)
