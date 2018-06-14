"""
This script uses 2dfdr to reduce 2dF/AAOMega data as part of the SAGA Survey. It
assumes:
* That the raw data from any given night is all in a directory that also contains
  "ccd_1" and "ccd_2" directories.
* That any given set of exposures contains a red-appropriate flat, a
  blue-appropriate flat, an arc, and some number of science exposures.
* That the 2dfdr binaries are somewhere on your current path.

Additional *default* assumptions that can be relaxed via appropriate
command-line arguments are:
* the first 10 exposures are biases
* the order of each exposure set is "redflat, blueflat, arc, science x n"
* there are 4 science exposures
* the current working directory is the one with the data and that should get
  all the reduce files
* idx files are "$HOME/observing/aat/ozdes_blue_mod.idx" and
  "$HOME/observing/aat/ozdes_red_mod.idx"

So then basic usage is e.g.:
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
REDUCE_TEMPL = 'aaorun reduce_run {0} -idxfile {1}'

def check_make_biases(basepath, idx_template, raw_bases, biasexpnums):
    biasbase = os.path.join(basepath, 'bias')
    bias1 = os.path.join(biasbase, 'ccd1')
    bias2 = os.path.join(biasbase, 'ccd2')
    biasi = (bias1, bias2)

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
        print('No combined biases present.  Need to build them.')

        for i, biasdir in enumerate((bias1, bias2)):
            ccdnum = i + 1
            rawdir = 'ccd_{}'.format(ccdnum)

            #need to link biases
            for expnum in biasexpnums:
                fn = raw_bases[i] + '{:04}.fits'.format(expnum)
                check_or_make_symlink(os.path.relpath(os.path.join(rawdir, fn), biasi[i]),
                                      os.path.join(biasi[i], fn))


            idxfile = idx_template.format(CCDNUM_TO_NAME[ccdnum])
            bias_cmdline = REDUCE_TEMPL.format(raw_bases[i], idxfile)
            biaslogfn = os.path.join(biasbase, 'bias{}.log'.format(ccdnum))
            print('Making bias in', biasdir, 'as',bias_cmdline,  'with log file:', biaslogfn)
            with open(biaslogfn, 'w') as biaslogf:
                ps = subprocess.Popen(bias_cmdline, shell=True, cwd=biasdir,
                                      stdout=biaslogf, stderr=subprocess.STDOUT)
                retcode = ps.wait()
                if retcode != 0:
                    raise ValueError('Bias {} failed'.format(biasdir))

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
                 blueflat_expnum, arc_expnum, sci_expnums, bias_expnums):
    raw1 = os.path.join(basepath, 'ccd_1')
    raw2 = os.path.join(basepath, 'ccd_2')
    if not (os.path.exists(raw1) and os.path.exists(raw2)):
        raise ValueError('raw data directories are missing')
    raw_base1 = determine_rawfns(raw1, 1)
    raw_base2 = determine_rawfns(raw2, 2)

    biasfn1, biasfn2 = check_make_biases(basepath, idx_template,
                                         (raw_base1, raw_base2),
                                         bias_expnums)

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

    cmdline1 = REDUCE_TEMPL.format(raw_base1, idx_template.format(CCDNUM_TO_NAME[1]))
    cmdline2 = REDUCE_TEMPL.format(raw_base2, idx_template.format(CCDNUM_TO_NAME[2]))
    log1fn = os.path.join(fieldpath, 'ccd1.log')
    log2fn = os.path.join(fieldpath, 'ccd2.log')

    p1_passed = p2_passed = None
    bluecombinedfn = os.path.join(ccd1, raw_base1 + '_combined.fits')
    if os.path.exists(bluecombinedfn):
        print('ccd1 results already done:', bluecombinedfn)
        p1_passed = True
    else:
        print('ccd1 reducing as: "{}" in {}. Log file:{}'.format(cmdline1, ccd1, log1fn))
    redcombinedfn = os.path.join(ccd2, raw_base2 + '_combined.fits')
    if os.path.exists(redcombinedfn):
        print('ccd2 results already done:', redcombinedfn)
        p2_passed = True
    else:
        print('ccd2 reducing as: "{}" in {}. Log file:{}'.format(cmdline2, ccd2, log2fn))

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

            while p1_passed is None or p2_passed is None:
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

        splice_cmdline_templ = 'aaorun splice "{blue} {red}" -idxfile {idx} -output_file {out}'
        splice_cmdline = splice_cmdline_templ.format(blue=bluecombinedfn,
                                                     red=redcombinedfn,
                                                     idx=idx_template.format(CCDNUM_TO_NAME[1]),
                                                     out=splice_outfn)

        print('Splicing together the results using ', splice_cmdline,'.  Log file:', splicelogfn)

        with open(splicelogfn, 'w') as splicef:
            ps = subprocess.Popen(splice_cmdline, shell=True, cwd=basepath,
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
    parser.add_argument('--bias-start', default=1)
    parser.add_argument('--n-bias', default=10)
    parser.add_argument('--basepath', default = '.')
    parser.add_argument('--idx-template', default=os.path.join(os.environ['HOME'], 'observing/aat/ozdes_{}_mod.idx'))

    args = parser.parse_args()

    enum1 = int(args.firstexposurenum)
    redexp = enum1 + int(args.redflat_index)
    blueexp = enum1 + int(args.blueflat_index)
    arcexp = enum1 + int(args.arc_index)
    sciexp_start = enum1 + int(args.sci_offset)
    sciexps = [i+sciexp_start for i in range(args.n_science)]
    biasexps = [i+args.bias_start for i in range(args.n_bias)]

    reduce_field(args.fieldname, args.idx_template, args.basepath, redexp,
                 blueexp, arcexp, sciexps, biasexps)
