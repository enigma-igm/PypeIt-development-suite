"""
Summarize the extracted objects in a folder of spec1d files: position, signal-to-noise, wavelength range,
and median signal-to-noise per pixel blueward and redward of a wavelength (e.g. the Lyman-alpha break).

Usage: python3 -m pypeitdev.jwst.debug_spec1d_summary /path/to/Science --break_wave 9440
"""
import argparse
import glob
import os

import numpy as np
from pypeit import specobjs


def window_snr(wave, flux, ivar, mask, wmin, wmax):
    use = mask & (wave > wmin) & (wave < wmax) & (ivar > 0)
    if use.sum() == 0:
        return np.nan, np.nan
    snr = flux[use] * np.sqrt(ivar[use])
    return np.median(snr), np.median(flux[use])


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('science_dir', type=str, help='Directory holding the spec1d files')
    parser.add_argument('--break_wave', type=float, default=9440.0, help='Wavelength to compare either side of [Angstrom]')
    parser.add_argument('--width', type=float, default=400.0, help='Width of each comparison window [Angstrom]')
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.science_dir, 'spec1d_*.fits')))
    print('{0} spec1d files'.format(len(files)))
    for f in files:
        sobjs = specobjs.SpecObjs.from_fitsfile(f, chk_version=False)
        for sobj in sobjs:
            print(os.path.basename(f)[7:32], 'objects: {0}, spatial pixel {1:.1f}, find signal-to-noise {2:.1f}'.format(
                sobjs.nobj, np.median(sobj.TRACE_SPAT), sobj.smash_snr))
            for kind in ['OPT', 'BOX']:
                wave = getattr(sobj, kind + '_WAVE')
                if wave is None:
                    print('    {0}: not extracted'.format(kind))
                    continue
                flux = getattr(sobj, kind + '_COUNTS')
                ivar = getattr(sobj, kind + '_COUNTS_IVAR')
                mask = getattr(sobj, kind + '_MASK')
                good = mask & (wave > 1.0)
                blue_snr, blue_flux = window_snr(wave, flux, ivar, mask, args.break_wave - args.width, args.break_wave - 20)
                red_snr, red_flux = window_snr(wave, flux, ivar, mask, args.break_wave + 50, args.break_wave + args.width)
                print('    {0}: wavelength {1:.0f}-{2:.0f} A, good pixels {3}/{4}; median signal-to-noise per pixel '
                      'blueward {5:.2f}, redward {6:.2f}; median counts blueward {7:.3f}, redward {8:.3f}'.format(
                          kind, wave[good].min(), wave[good].max(), good.sum(), wave.size,
                          blue_snr, red_snr, blue_flux, red_flux))


if __name__ == '__main__':
    main()
