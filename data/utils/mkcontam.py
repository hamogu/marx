import argparse

import numpy as np
from astropy.io import fits

'''
On 20 Feb 2010, JD asked Alexey about the uncertaintly on the fxy
image values.  Here is his response:

    "I'd say, a 10% relative uncertainty is a good estimate, probably
    slightly on the concervative side. Or, 10% relative and ~ +- 0.05
    absolute, whichever is greater. (The absolute uncertainty is in
    units of = tau for 0.67 keV at the given date and location).

On 23 Jan 2015, HMG discussed this with Herman Marshall looking at how
the change of the ACIS contamination has changed this picture. While the
functional form is still the same, the contamination is now much higher.
The final uncertainty does not depend linearly on fxy, since the absorption is
:math: `abs \propto exp(-\tau_1 f(x,y))`.
If :math:`\tau_1` is large, even a small fractional change in :math:`f(x,y)`
can make a big difference in the number of detected photons.
Thus, we decided to significantly reduce the relative tolerance.
'''


def rebin_factor(a, n):
    '''Rebin an array to a new shape.
    newshape must be a factor of a.shape.
    '''
    assert len(a.shape) == 2
    x, y = a.shape
    assert x == y

    return a.reshape((x, y / n, n)).mean(axis=2).reshape((x / n, n, y / n)).mean(axis=1)


def binup_factor(a, n):
    '''Increase size of array by n in each dimension.

    Each pixel will become and n*n island of pixels with identical values.
    '''
    return a.repeat(n, axis=1).repeat(n, axis=0)


def max_block_image(img, rel_tolerance, abs_tolerance):
    '''Find the maximal acceptable blocking factor'''
    for n in np.arange(1, 11):
        binned_img = rebin_factor(img, 2**n)
        expanded_binned = binup_factor(binned_img, 2**n)
        if not np.allclose(img, expanded_binned, rel_tolerance, abs_tolerance):
            return 2**(n - 1)
    return 2**n


def block_imgs(imgstack, rel_tolerance=0.01, abs_tolerance=0.01):
    '''Downsample images as much as possible

    Parameters
    ----------
    imgstack : np.ndarray
        Array of n 2d images (of m*m pixels) with shape (n,m,m)
    rel_tolerance : float
        Maximum relative tolerance
    abs_tolerance : float
        Maximum relative tolerance

    Returns
    -------
    newstack : np.ndarray
        downsampled images
    maxb : int
        downsampling factor
    '''
    n_comp = imgstack.shape[0]
    maxb = np.zeros(n_comp)

    for i in range(n_comp):
        maxb[i] = max_block_image(imgstack[i,:,:], rel_tolerance, abs_tolerance)
    maxb = int(np.min(maxb))

    newstack = np.zeros((n_comp, imgstack.shape[1] / maxb, imgstack.shape[2] / maxb))
    for i in range(n_comp):
        newstack[i, :, :] = rebin_factor(imgstack[i, :, :], maxb)
    return newstack, maxb


def marx_contam_file(filein, fileout, rel_tolerance=0.01, abs_tolerance=0.01):
    '''Change ACIS Contam file from CALDB to MARX format

    Parameters
    ----------
    filein : string
        filename and path to CALDB file
    fileout : sring
        filename and path where MARX file should be written
    rel_tolerance : float
        Maximum relative tolerance
    abs_tolerance : float
        Maximum relative tolerance
    '''
    hdus = fits.open(filein)

    prihdr = fits.Header()
    prihdr['HISTORY'] = 'This file was created using the marx script mkcontam.py {0} {1}'.format(filein, fileout)

    listhdus = [fits.PrimaryHDU(header=prihdr)]
    for i in range(10):
        h = hdus[i + 1]
        data = h.data

        if np.any(data['component'] != 0) or np.any(data['weight'] != 1):
            raise NotImplementedError("This script can deal with one contamination component only.")

        newstack, maxb = block_imgs(data['fxy'])

        hdr = fits.Header()
        hdr['mission'] = "AXAF"
        hdr['instrum'] = "ACIS"
        hdr['detnam'] = "ACIS-{0}".format(i)
        hdr['ccd_id'] = i
        hdr['marxvers'] = 5.3
        hdr['fxyblk'] = maxb

        cols = []
        for n in ["component", "n_energy", "energy", "mu", "n_time", "time",
                  "tau0", "tau1"]:

            if n not in data.dtype.names:
                raise ValueError('Input file is missing column {0}'.format(n))
            colindex = np.where(np.array(data.dtype.names) == n)[0][0]
            cols.append(h.columns[colindex])
        dim1 = newstack.shape[1]
        dim2 = newstack.shape[2]
        cols.append(fits.Column(array=newstack, name='fxy',
                                format='{0}E'.format(dim1 * dim2),
                                dim='({0},{1})'.format(dim1, dim2)))

        tab = fits.BinTableHDU.from_columns(fits.ColDefs(cols),
                                            header=hdr,
                                            name='ACIS{0}_CONTAM'.format(i))
        listhdus.append(tab)
    newhdulist = fits.HDUList(listhdus)
    newhdulist.writeto(fileout, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Change ACIS Contam file from CALDB to MARX format')
    parser.add_argument('infile', help='filename and path to CALDB file')
    parser.add_argument('outfile', help='filename and path where MARX file should be written')
    parser.add_argument('-r', '--relative-tolerance', default=0.01, type=float,
                        help='relative tolerance allowed in rebinning (default=0.01)')
    parser.add_argument('-a', '--absolute-tolerance', default=0.01, type=float,
                        help='absolute tolerance allowed in rebinning (default=0.01)')

    args = parser.parse_args()
    marx_contam_file(args.infile, args.outfile,
                     rel_tolerance=args.relative_tolerance,
                     abs_tolerance=args.absolute_tolerance)
