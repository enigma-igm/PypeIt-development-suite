"""
Combine boxcar spectra from several spec1d files that share one wavelength grid (same slit, same calibrations)
with an inverse-variance weighted mean, and convert from summed MJy/sr to flux density.

The reduced images carry the calwebb units (BUNIT in the _cal.fits SCI header): MJy/sr for source_type
EXTENDED, where the boxcar sum times the pixel solid angle (PIXAR_SR) is the flux density in MJy, or MJy for
source_type POINT, where the boxcar sum already is. The output is in microjansky (1e-29 erg/s/cm^2/Hz).

Usage: python3 -m pypeitdev.jwst.combine_spec1d spec1d_a.fits spec1d_b.fits --cal_file x_cal.fits --outfile combined.fits
"""
import argparse
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table
from pypeit import specobjs
from pypeit.core import coadd
from pypeit.utils import inverse


def slit_conversion(cal_file, source_name):
    """Return the factor from boxcar-summed calwebb units to microjansky, the units, and PIXAR_SR for one slit."""
    with fits.open(cal_file) as hdul:
        sci_hdrs = [h.header for h in hdul if h.name == 'SCI' and str(h.header.get('SRCNAME')).strip() == source_name]
    if len(sci_hdrs) == 0:
        raise ValueError('No SCI extension with SRCNAME {0} in {1}'.format(source_name, cal_file))
    sci_hdr = sci_hdrs[0]
    bunit = sci_hdr['BUNIT'].strip()
    pixar_sr = sci_hdr['PIXAR_SR']
    if bunit == 'MJy/sr':
        to_ujy = pixar_sr * 1e12  # MJy/sr summed over pixels -> MJy -> microjansky
    elif bunit == 'MJy':
        to_ujy = 1e12  # MJy -> microjansky
    else:
        raise ValueError('Unexpected BUNIT {0} in {1}'.format(bunit, cal_file))
    print('{0} slit {1}: calwebb units {2}, PIXAR_SR {3:.4g}; conversion to microjansky {4:.4g}'.format(
        os.path.basename(cal_file), sci_hdr.get('SLTNAME'), bunit, pixar_sr, to_ujy))
    return to_ujy, bunit, pixar_sr


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('spec1d_files', type=str, nargs='+', help='spec1d files, one object each')
    parser.add_argument('--cal_file', type=str, nargs='+', required=True,
                        help='calwebb _cal.fits file(s) to read BUNIT and PIXAR_SR from: one for all spectra, '
                             'or one per spec1d file (same order) when they come from different slits')
    parser.add_argument('--source_name', type=str, required=True,
                        help='Source name of the slit in the _cal.fits file, e.g. 4713_6733; units and pixel area are per slit')
    parser.add_argument('--outfile', type=str, required=True, help='Output FITS table')
    parser.add_argument('--extraction', type=str, default='BOX', choices=['BOX', 'OPT'], help='Extraction to combine')
    parser.add_argument('--resample', action='store_true',
                        help='Linearly interpolate every spectrum onto the wavelength grid of the first one '
                             '(pypeit.core.coadd.interp_oned), e.g. for dithers along a slit')
    args = parser.parse_args()

    if len(args.cal_file) == 1:
        cal_files = args.cal_file * len(args.spec1d_files)
    elif len(args.cal_file) == len(args.spec1d_files):
        cal_files = args.cal_file
    else:
        raise ValueError('Give one --cal_file, or one per spec1d file')
    conversions = {f: slit_conversion(f, args.source_name) for f in set(cal_files)}

    waves, fluxes, ivars, masks = [], [], [], []
    for f, cal_file in zip(args.spec1d_files, cal_files):
        to_ujy, bunit, pixar_sr = conversions[cal_file]
        sobjs = specobjs.SpecObjs.from_fitsfile(f, chk_version=False)
        if sobjs.nobj != 1:
            raise ValueError('{0} has {1} objects; expected 1'.format(f, sobjs.nobj))
        sobj = sobjs[0]
        waves.append(getattr(sobj, args.extraction + '_WAVE'))
        fluxes.append(getattr(sobj, args.extraction + '_COUNTS') * to_ujy)
        ivars.append(getattr(sobj, args.extraction + '_COUNTS_IVAR') / to_ujy**2)
        masks.append(getattr(sobj, args.extraction + '_MASK').astype(bool))

    if args.resample:
        for i in range(1, len(waves)):
            out = coadd.interp_oned(waves[0], waves[i], fluxes[i], ivars[i], masks[i], kind='linear')
            fluxes[i], ivars[i], masks[i] = out[0], out[1], out[2]
            waves[i] = waves[0]
    waves, fluxes, ivars, masks = map(np.array, (waves, fluxes, ivars, masks))
    good_wave = np.all(waves > 1.0, axis=0)
    max_wave_diff = np.max(np.abs(waves[:, good_wave] - waves[0, good_wave]))
    dwave = np.median(np.diff(waves[0, good_wave]))
    print('largest wavelength difference between spectra: {0:.3g} Angstrom ({1:.2%} of a pixel)'.format(
        max_wave_diff, max_wave_diff / dwave))
    # Traces at slightly different spatial positions give slightly different boxcar wavelengths
    if max_wave_diff > 0.05 * dwave:
        raise ValueError('Spectra differ by more than 5% of a pixel in wavelength; resample them first.')

    weights = ivars * masks
    wsum = weights.sum(axis=0)
    flux = (weights * fluxes).sum(axis=0) * inverse(wsum)
    nused = masks.sum(axis=0)
    mask = (wsum > 0) & good_wave

    out = Table({'wave': waves[0], 'flux': flux, 'ivar': wsum, 'mask': mask, 'nexp': nused})
    out.meta['FLUXUNIT'] = 'microjansky'
    out.meta['NSPEC1D'] = len(args.spec1d_files)
    for i, (f, cal_file) in enumerate(zip(args.spec1d_files, cal_files)):
        out.meta['SPEC1D{0}'.format(i)] = os.path.basename(f)
        out.meta['CALBUN{0}'.format(i)] = conversions[cal_file][1]
        out.meta['PIXAR{0}'.format(i)] = conversions[cal_file][2]
    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    out.write(args.outfile, overwrite=True)
    print('Wrote', args.outfile)


if __name__ == '__main__':
    main()
