"""
Report the ratio of the uniform (extended) to point-source pathloss for one MSA source in calwebb _cal.fits files.
For the same count rate, a POINT reduction would be brighter than an EXTENDED one by about this ratio.

Usage: python3 -m pypeitdev.jwst.debug_pathloss_ratio 4713_6733 /path/to/level2/*_cal.fits
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

    windows = [(9600, 10000), (10000, 12000), (12000, 12600)]
    for f in args.cal_files:
        with datamodels.open(f) as model:
            for slit in model.slits:
                if str(slit.source_name).strip() != args.source_name:
                    continue
                wave = 1e4 * np.asarray(slit.wavelength, dtype=float)
                pl_un = np.asarray(slit.pathloss_uniform, dtype=float) if slit.pathloss_uniform is not None else None
                pl_ps = np.asarray(slit.pathloss_point, dtype=float) if slit.pathloss_point is not None else None
                if pl_un is None or pl_ps is None or pl_un.size == 0 or pl_ps.size == 0:
                    print(os.path.basename(f), 'pathloss arrays missing')
                    continue
                ratio = pl_un / pl_ps
                parts = []
                for lo, hi in windows:
                    use = np.isfinite(ratio) & (wave > lo) & (wave < hi)
                    parts.append('{0}-{1} A: {2:.3f}'.format(lo, hi, np.median(ratio[use]) if use.any() else np.nan))
                print('{0}  uniform/point pathloss  {1}'.format(os.path.basename(f), ';  '.join(parts)))


if __name__ == '__main__':
    main()
