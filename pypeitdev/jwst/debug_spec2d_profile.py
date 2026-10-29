"""
Print the slit width and the wavelength-collapsed spatial profile of sky-subtracted spec2d files.

Usage: python3 -m pypeitdev.jwst.debug_spec2d_profile /path/to/Science --wave_min 9600 --wave_max 12800
"""
import argparse
import glob
import os

import numpy as np
from pypeit import spec2dobj


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('science_dir', type=str, help='Directory holding the spec2d files')
    parser.add_argument('--wave_min', type=float, default=9600.0, help='Minimum wavelength to collapse over [Angstrom]')
    parser.add_argument('--wave_max', type=float, default=12800.0, help='Maximum wavelength to collapse over [Angstrom]')
    parser.add_argument('--raw', action='store_true', help='Use the science image without subtracting the sky model')
    args = parser.parse_args()

    for f in sorted(glob.glob(os.path.join(args.science_dir, 'spec2d_*.fits'))):
        s = spec2dobj.Spec2DObj.from_file(f, 'DET01')
        left, right, _ = s.slits.select_edges()
        spec_mid = left.shape[0] // 2
        print(os.path.basename(f))
        print('  slit edges at spectral middle: left {0:.1f}, right {1:.1f}, width {2:.1f} pixels'.format(
            left[spec_mid, 0], right[spec_mid, 0], right[spec_mid, 0] - left[spec_mid, 0]))

        resid = s.sciimg if args.raw else s.sciimg - s.skymodel
        good = s.bpmmask.mask == 0
        inwave = (s.waveimg > args.wave_min) & (s.waveimg < args.wave_max)
        use = good & inwave & (s.ivarraw > 0)
        nspat = resid.shape[1]
        spat_lo = int(np.floor(np.nanmin(left[:, 0])))
        spat_hi = int(np.ceil(np.nanmax(right[:, 0])))
        profile, profile_sig = [], []
        for ispat in range(max(spat_lo, 0), min(spat_hi + 1, nspat)):
            col = use[:, ispat]
            if col.sum() < 10:
                continue
            flux = np.sum(resid[col, ispat])
            sig = np.sqrt(np.sum(1.0 / s.ivarraw[col, ispat]))
            profile.append((ispat, flux, sig))
        print('  spatial pixel, summed sky-subtracted counts, signal-to-noise')
        for ispat, flux, sig in profile:
            print('    {0:4d} {1:12.2f} {2:8.1f}'.format(ispat, flux, flux / sig))


if __name__ == '__main__':
    main()
