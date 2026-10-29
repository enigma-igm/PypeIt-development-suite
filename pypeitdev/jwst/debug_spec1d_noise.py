"""
Compare the quoted errors of a combined spectrum (from combine_spec1d) with its actual pixel-to-pixel scatter,
and report the boxcar extraction radius of the input spec1d files.

Usage: python3 -m pypeitdev.jwst.debug_spec1d_noise combined.fits spec1d_a.fits --wave_min 10000 --wave_max 12000
"""
import argparse

import numpy as np
from astropy.table import Table
from pypeit import specobjs
from pypeit.utils import fast_running_median, nan_mad_std


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('combined_file', type=str, help='Output of combine_spec1d')
    parser.add_argument('spec1d_files', type=str, nargs='*', help='spec1d files to read the boxcar radius from')
    parser.add_argument('--wave_min', type=float, default=10000.0, help='Minimum wavelength [Angstrom]')
    parser.add_argument('--wave_max', type=float, default=12000.0, help='Maximum wavelength [Angstrom]')
    parser.add_argument('--smooth', type=int, default=25, help='Running-median width for the continuum [pixels]')
    args = parser.parse_args()

    t = Table.read(args.combined_file)
    use = np.asarray(t['mask'], dtype=bool) & (t['wave'] > args.wave_min) & (t['wave'] < args.wave_max)
    flux = np.asarray(t['flux'])[use]
    sigma = 1.0 / np.sqrt(np.asarray(t['ivar'])[use])
    resid = flux - fast_running_median(flux, args.smooth)
    scatter = nan_mad_std(resid)
    print('{0:.0f}-{1:.0f} A: median flux {2:.3f}, median quoted error {3:.3f}, actual scatter {4:.3f} '
          '(scatter / quoted error = {5:.2f})'.format(args.wave_min, args.wave_max, np.median(flux),
                                                      np.median(sigma), scatter, scatter / np.median(sigma)))
    print('signal-to-noise per pixel: {0:.1f} from quoted errors, {1:.1f} from actual scatter'.format(
        np.median(flux) / np.median(sigma), np.median(flux) / scatter))

    # Correlated noise makes the pixel-to-pixel scatter underestimate the true error per pixel.
    # Check the correlation of neighbouring residuals and how the scatter of binned residuals scales.
    r = resid - np.median(resid)
    for lag in [1, 2, 3]:
        print('  correlation of residuals {0} pixel(s) apart: {1:.3f}'.format(lag, np.corrcoef(r[:-lag], r[lag:])[0, 1]))
    for nbin in [1, 2, 4, 8]:
        n = (r.size // nbin) * nbin
        binned = r[:n].reshape(-1, nbin).mean(axis=1)
        print('  bins of {0}: scatter {1:.4f}, expected for uncorrelated noise {2:.4f}, ratio {3:.2f}'.format(
            nbin, nan_mad_std(binned), scatter / np.sqrt(nbin), nan_mad_std(binned) / (scatter / np.sqrt(nbin))))

    for f in args.spec1d_files:
        sobj = specobjs.SpecObjs.from_fitsfile(f, chk_version=False)[0]
        print('{0}: boxcar width (median pixels summed per spectral pixel) {1:.1f}, object FWHM {2:.2f} pixels'.format(
            f.split('spec1d_')[-1][:25], np.median(sobj.BOX_NPIX[sobj.BOX_NPIX > 0]), np.median(sobj.FWHMFIT)
            if sobj.FWHMFIT is not None else sobj.FWHM))


if __name__ == '__main__':
    main()
