"""
Print the metadata of the slits in calwebb _cal.fits files that carry a given source name.

Usage: python3 -m pypeitdev.jwst.debug_cal_slit_info 4713_6733 /path/to/level2/*_cal.fits
"""
import argparse
import os

import numpy as np
from stdatamodels.jwst import datamodels


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('source_name', type=str, help='Source name, e.g. 4713_6733')
    parser.add_argument('cal_files', type=str, nargs='+', help='calwebb _cal.fits files')
    args = parser.parse_args()

    for f in args.cal_files:
        with datamodels.open(f) as model:
            print(os.path.basename(f))
            for slit in model.slits:
                if str(slit.source_name).strip() != args.source_name:
                    continue
                wave = np.asarray(slit.wavelength)
                print('  slit name {0}, source id {1}, source ra {2:.6f}, dec {3:.6f}, shutter state {4}'.format(
                    slit.name, slit.source_id, slit.source_ra, slit.source_dec, slit.shutter_state))
                print('  source position in shutter x {0:.3f}, y {1:.3f}; quadrant {2}, row {3}, column {4}'.format(
                    slit.source_xpos, slit.source_ypos, slit.quadrant, slit.xcen, slit.ycen))
                print('  subimage origin x {0}, y {1}, size {2}x{3}; wavelength {4:.0f}-{5:.0f} A'.format(
                    slit.xstart, slit.ystart, slit.xsize, slit.ysize,
                    1e4 * np.nanmin(wave), 1e4 * np.nanmax(wave)))


if __name__ == '__main__':
    main()
