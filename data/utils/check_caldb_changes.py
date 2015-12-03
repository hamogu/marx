'''This script helps to spot CALDB changes that effect MARX calibration files.

It is not foolproof, it just **helps** spotting changes!
Here is what it does:

- Find CALDB files of the same name as a MARX cal file.
- Compare those two files with the astropy FITSDiff object.

Files that are generated from CALDB files by reordering are
generated using ISIS and the newly generated versions are compared to the
current MARX version.

Sometimes, CALDB gets updated in a dramatically new way, e.g. the file format
changes. This script will miss those changes, they have to be found by hand.
'''

import os
from glob import glob
from subprocess import call

from astropy.io.fits import FITSDiff

marxcaldb = glob('../caldb/*.fits')
CALDB = os.environ['CALDB']
ignore_key = ['HISTORY', 'HISTNUM', 'CHECKSUM']

if CALDB is None:
    raise Exception('Environment variable CALDB must be set before calling this script.')


def find(name, path):
    for root, dirs, files in os.walk(path):
        if os.path.basename(name) in files:
            # Check if this is the last version number of a CALDB file.
            marxver = int(name[-9:-5])
            for f in files:
                if (int(f[-9:-5]) > marxver) and (f[:-9] == os.path.basename(name)[:-9]):
                    print('Check if {0} superceedes {1}'.format(f, name))
            return os.path.join(root, os.path.basename(name))


for f in marxcaldb:
    diff = None
    caldbfile = find(f, CALDB)
    if caldbfile is not None:
        diff = FITSDiff(f, caldbfile, ignore_keywords=ignore_key)
    # special cases
    else:
        fbase = os.path.basename(f)
        if os.path.basename(f) == 'acisfef.fits':
            caldbfile = find('acisD2000-01-29fef_phaN0005.fits', CALDB)
            call('../caldb/fixfef.sl {0} {1}_new'.format(caldbfile, f), shell=True)
            diff = FITSDiff(f, f + '_new', ignore_keywords=ignore_key)
        elif fbase.startswith('acisD1999-07-22subpixN00'):
            caldbfile = find(fbase.replace('_marx', ''), CALDB)
            call('../utils/mksubpix.sl {0} {1}_new'.format(caldbfile, f), shell = True)
            diff = FITSDiff(f, f + '_new', ignore_keywords=ignore_key)
        elif fbase.startswith('acisD1999-08-13contamN00'):
            caldbfile = find(fbase.replace('_marx', ''), CALDB)
            call('python ../utils/mkcontam.py {0} {1}_new'.format(caldbfile, f), shell = True)
            diff = FITSDiff(f, f + '_new', ignore_keywords=ignore_key)
        else:
            print("Don't know how to compare {0} - check manually".format(fbase))

    if caldbfile is None:
        print('No comparison for {0} found'.format(f))
        print('Check CALDB manually!')
    else:
        if not diff.identical:
            print diff.report()
        else:
            print('{0} is up to date.' .format(f))

print('''MARX uses some CALDB files directly, others are derived
from CALDB (e.g. by reordering).
This script has made new versions of those generated files, called *.fits_new
If no differences were found, those file can be deleted.

If any of these calibration files have changed, that probably also has an
effect on /hrma/corr*.dat.
See utils/Makefile for instructions on how to redo those.

If any file name changes, edit marxcaldb.dat, Markefile.in, and marx.lis accordingly; then do git rm, git add
''')
